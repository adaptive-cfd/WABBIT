!> \brief wrapper for RHS call in time step function, computes RHS in work array
!! (inplace)
!!
!! calls RHS depending on physics
!!
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!**********************************************************************************************

subroutine RHS_wrapper(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID, stage)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_rhs  (:, :, :, :, :)       !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)      !> hvy_mask: the mask function, usx,usy,usz,color,sponge
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik), intent(in), optional        :: stage

    real(kind=rk), dimension(3)         :: volume_int                   !> global integral
    real(kind=rk), dimension(3)         :: dx, x0                       !> spacing and origin of a block
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, hvy_id  ! loop variables
    integer(kind=ik)                    :: g, hvy_id_mask, Bs(1:3), grhs, Neqn_RHS, Jmin_active, Jmax_active, stage_set
    integer(kind=2)                     :: n_domain(1:3)
    real(kind=rk)                       :: t0, t1
    logical                             :: is_equidistant

    Neqn_RHS = params%n_eqn_rhs
    Bs   = params%Bs
    g    = params%g
    grhs = params%g_rhs
    t0       = MPI_wtime()
    n_domain = 0
    x0       = 0.0_rk
    dx       = 0.0_rk
    hvy_ID   = 1

    Jmin_active = minActiveLevel_tree(tree_ID_flow)
    Jmax_active = maxActiveLevel_tree(tree_ID_flow)
    is_equidistant = (Jmin_active == Jmax_active)

    ! assume this is stage 1 if we don't know. This is just important for optimizations for the Poisson solver
    stage_set = 1
    if (present(stage)) stage_set = stage
 
    !-------------------------------------------------------------------------
    ! CoarseExtension update of input data
    !-------------------------------------------------------------------------
    ! we assume data are sync'ed on call
    ! call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
    ! hvy_n(tree_ID), inputDataSynced=.true. )

    !-------------------------------------------------------------------------
    ! create mask function at current time
    !-------------------------------------------------------------------------
    t1 = MPI_wtime()
    call createMask_tree(params, time, hvy_mask, hvy_tmp)
    call toc( "RHS_wrapper::createMask_tree", 31, MPI_wtime()-t1 )


    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    t1 = MPI_wtime()
    ! performs initializations in the RHS module, such as resetting integrals
    ! for this stage, just pass any block (result does not depend on block), hvy_id=1 and set x0=dx=0
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,1), g, x0, dx, &
         hvy_rhs(:,:,:,:,1), hvy_mask(:,:,:,:,1), "init_stage" )

    !-------------------------------------------------------------------------
    ! 2nd stage: integral_stage. (called for all blocks)
    !-------------------------------------------------------------------------
    ! For some RHS, the eqn depend not only on local, block based qtys, such as
    ! the state vector, but also on the entire grid, for example to compute a
    ! global forcing term (e.g. in FSI the forces on bodies). As the physics
    ! modules cannot see the grid, (they only see blocks), in order to encapsulate
    ! them nicer, two RHS stages have to be defined: integral / local stage.
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
        ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
        ! Note: hvy_mask is not used if in this case by the RHS routines...
        hvy_id_mask = hvy_id
        if (size(hvy_mask,5) == 1) hvy_id_mask = 1

        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, x0, dx,&
        hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id_mask), "integral_stage", n_domain )
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used ror example for MPI_REDUCES
    ! for this stage, just pass any block (result does not depend on block), hvy_id=1 and set x0=dx=0
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, 1), g, x0, dx, &
    hvy_rhs(:,:,:,:,1), hvy_mask(:,:,:,:,1), "post_stage" )
    call toc( "RHS_wrapper::integral-stage", 32, MPI_wtime()-t1 )


    !-------------------------------------------------------------------------
    ! 4th stage: local evaluation of RHS on all blocks (called for all blocks)
    !-------------------------------------------------------------------------
    ! the second stage then is what you would usually do: evaluate local differential
    ! operators etc.

    t1 = MPI_wtime()
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
        ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
        ! Note: hvy_mask is not used if in this case by the RHS routines...
        hvy_id_mask = hvy_id
        if (size(hvy_mask,5) == 1) hvy_id_mask = 1

        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, &
        x0, dx, hvy_rhs(:,:,:,:, hvy_id), hvy_mask(:,:,:,:,hvy_id_mask), "local_stage", n_domain)
    enddo
    call toc( "RHS_wrapper::local-stage", 33, MPI_wtime()-t1 )

    ! For purely hyperbolic equations, at this stage the RHS is completely ready.
    ! However, for NSPP we have to actually solve the pressure-poisson equation to get the correct pressure gradient term.
    ! This has to be done at exactly this position. However, we will need to sync some variables, which you usually would love to avoid in the RHS_wrapper
    if (params%PDE_type == "mixed_hyperbolic_parabolic") then

        ! TBH, this makes less sense to me somehow, I am skipping it for now as for several stages we always alternate between computing P and dP
        ! ! incremental version, add old pressure gradient so that new pressure to be solved is only Dp - this reduces the amount of iterations for projection needed to converge, however it only really helps for adaptive grids. We assume that hvy_block contains the pressure of the last MultiGrid solve
        ! if (stage > 1 .and. .not. is_equidistant) then
        !     do k = 1, hvy_n(tree_ID)
        !         hvy_id = hvy_active(k, tree_ID)
        !         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        !         call get_block_spacing_origin( params, lgt_id, x0, dx )

        !         ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        !         if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        !         ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
        !         ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
        !         ! Note: hvy_mask is not used if in this case by the RHS routines...
        !         hvy_id_mask = hvy_id
        !         if (size(hvy_mask,5) == 1) hvy_id_mask = 1

        !         ! RHS_meta will assume the pressure to be in component dim+1, it then will update the hvy_rhs with the pressure gradient term
        !         call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
        !             hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "pressure_gradient_stage", n_domain )
        !     enddo
        ! endif

        ! sync RHS
        t1 = MPI_wtime()
        call sync_ghosts_RHS_tree( params, hvy_rhs(:,:,:,1:Neqn_RHS, :), tree_ID, g_minus=grhs, g_plus=grhs )
        call toc( "RHS_wrapper::sync-RHS", 10035, MPI_wtime()-t1 )

        ! ! debug save
        ! call saveHDF5_tree("RHSx_0.h5", 0.0_rk, 0, 1, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)
        ! call saveHDF5_tree("RHSy_0.h5", 0.0_rk, 0, 2, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)

        ! divergence of RHS
        t1 = MPI_wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            if ( .not. All(params%periodic_BC) ) then
                ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
            endif
            ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
            ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
            ! Note: hvy_mask is not used if in this case by the RHS routines...
            hvy_id_mask = hvy_id
            if (size(hvy_mask,5) == 1) hvy_id_mask = 1
            ! RHS_meta will compute the divergence and store it in the first component of hvy_block

            call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
             hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "divergence_stage", n_domain )
        enddo
        call toc( "RHS_wrapper::divergence-stage", 10036, MPI_wtime()-t1 )

        ! ! debug save
        ! call saveHDF5_tree("divRHS_0.h5", 0.0_rk, 0, 1, params, hvy_block(:,:,:,:,:), tree_ID, no_sync=.false.)

        ! MultiGrid solve, we sync before and after in multigrid_solve
        ! Logic here: We can reuse the latest solution as initialization if the grid is equidistant and we are sure it is the same as before
        t1 = MPI_wtime()
        call multigrid_solve(params, hvy_block(:,:,:,params%dim+1:params%dim+1,:), hvy_block(:,:,:,1:1,:), hvy_block(:,:,:,2:2,:), tree_ID, init_0=.not.(is_equidistant .and. (stage > 1 .or. (stage == 1 .and. .not. params%adapt_tree))), verbose=.false., hvy_full=hvy_block, time=time)

        ! ! debug save
        ! call saveHDF5_tree("phi_0.h5", 0.0_rk, 0, params%dim+1, params, hvy_block(:,:,:,:,:), tree_ID, no_sync=.false.)

        ! update RHS with pressure gradient
        t1 = MPI_wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            if ( .not. All(params%periodic_BC) ) then
                ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
            endif
            ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
            ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
            ! Note: hvy_mask is not used if in this case by the RHS routines...
            hvy_id_mask = hvy_id
            if (size(hvy_mask,5) == 1) hvy_id_mask = 1

            ! RHS_meta will assume the pressure to be in component dim+1, it then will update the hvy_rhs with the pressure gradient term
            call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
             hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "pressure_gradient_stage", n_domain )
        enddo
        call toc( "RHS_wrapper::NSPP-pressure-poisson-solver", 10037, MPI_wtime()-t1 )

        ! ! debug save
        ! call saveHDF5_tree("RHS2x_0.h5", 0.0_rk, 0, 1, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)
        ! call saveHDF5_tree("RHS2y_0.h5", 0.0_rk, 0, 2, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)

    endif

    ! This subroutine is actually a hack. The RHS of the ACM eqn can be equipped with a pressure gradient term that is supposed 
    ! to keep the mean flow constant - but that seems to work only approximately. I do not know where the 
    ! difference comes from. As the mean pressure forcing does not perfectly set the mean flow constant,
    ! we remove in this function the part that the forcing did not manage to remove completely. As it was now useless, we removed
    ! the mean pressure gradient from the RHS_ACM to simplify the code.
    ! Because it is a tree-level routine, it cannot be merged easily in the ACM module.
    ! -TE
    if (params_acm%use_channel_forcing) then
        call ACM_remove_channel_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
    endif
    if (params_acm%HIT_linear_forcing) then
        call HIT_remove_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
    endif

    call toc( "RHS_wrapper (TOTAL)", 30, MPI_wtime()-t0 )
end subroutine RHS_wrapper



subroutine ACM_remove_channel_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_rhs  (:, :, :, :, :)     !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)      !> hvy_mask: the mask function, usx,usy,usz,color,sponge
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk)                :: rhs_mean, dV, V_channel, y
    integer(kind=ik)             :: ix, iy, iz, k, hvy_id, mpierr, lgt_id, g, Bs(1:3)
    real(kind=rk), dimension(3)  :: dx, x0

    ! This subroutine is actually a hack. The RHS of the ACM eqn can be equipped with a pressure gradient term that is supposed 
    ! to keep the mean flow constant - but that seems to work only approximately. I do not know where the 
    ! difference comes from. As the mean pressure forcing does not perfectly set the mean flow constant,
    ! we remove in this function the part that the forcing did not manage to remove completely. As it was now useless, we removed
    ! the mean pressure gradient from the RHS_ACM to simplify the code.
    ! Because it is a tree-level routine, it cannot be merged easily in the ACM module.
    ! -TE
    g = params%g 
    Bs = params%Bs

    ! --- Step1: compute the mean rhs
    rhs_mean = 0.0_rk
    
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )
        ! volume element
        dV = product(dx)

        do iy = g+1, Bs(2)+g
            y = x0(2) + dble(iy-(g+1)) * dx(2)

            ! mean of rhs inside the fluid (excluding the channel walls) 
            if (( y>params_acm%h_channel).and.(y<params_acm%domain_size(2)-params_acm%h_channel)) then
                rhs_mean = rhs_mean + sum(hvy_rhs(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 1, hvy_id))*dV
            endif
        enddo
    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE, rhs_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
    
    ! volume of fluid in the channel
    V_channel = params_acm%domain_size(1)*params_acm%domain_size(3)*(params_acm%domain_size(2)-2.0_rk*params_acm%h_channel)

    ! --- Step 2: remove it from RHS, so that the zero mode is not drifting
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        do iy = g+1, Bs(2)+g
            y = x0(2) + dble(iy-(g+1)) * dx(2)

            ! exclude channel walls
            if (( y>params_acm%h_channel).and.(y<params_acm%domain_size(2)-params_acm%h_channel)) then
                hvy_rhs(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 1, hvy_id) = hvy_rhs(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 1, hvy_id) - rhs_mean / V_channel
            end if
        end do
    end do

end subroutine


subroutine HIT_remove_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_rhs  (:, :, :, :, :)     !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)      !> hvy_mask: the mask function, usx,usy,usz,color,sponge
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk)                :: rhs_mean, dV, V_domain, u_RMS
    integer(kind=ik)             :: ix, iy, iz, k, hvy_id, mpierr, lgt_id, g, Bs(1:3)
    real(kind=rk), dimension(3)  :: dx, x0

    ! This subroutine is actually a hack. It enforces zero mean flow for HIT simulations with linear forcing.
    ! Because it is a tree-level routine, it cannot be merged easily in the ACM module.
    ! -TE
    g = params%g 
    Bs = params%Bs

    ! --- mean flow was already computed in RHS_ACM, we just need to remove it here
    
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! remove mean flow
        hvy_block(:,:,:,1,hvy_id) = hvy_block(:,:,:,1,hvy_id) - params_acm%mean_flow(1)
        hvy_block(:,:,:,2,hvy_id) = hvy_block(:,:,:,2,hvy_id) - params_acm%mean_flow(2)
        if (params%dim == 3) hvy_block(:,:,:,3,hvy_id) = hvy_block(:,:,:,3,hvy_id) - params_acm%mean_flow(3)
    enddo

end subroutine HIT_remove_meanrhs