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

subroutine RHS_wrapper(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_rhs  (:, :, :, :, :)       !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)      !> hvy_mask are qtys that depend on grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk), dimension(3)         :: volume_int                   !> global integral
    real(kind=rk), dimension(3)         :: dx, x0                       !> spacing and origin of a block
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, hvy_id  ! loop variables
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=2)                     :: n_domain(1:3)
    real(kind=rk)                       :: t0, t1

    ! grid parameter
    Bs = params%Bs
    g  = params%g
    t0 = MPI_wtime()
    n_domain = 0

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
    hvy_id = hvy_active(1, tree_ID) ! for this stage, just pass any block (result does not depend on block)
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
         hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "init_stage" )

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

        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
        endif

        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, x0, dx,&
        hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "integral_stage", n_domain )
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used ror example for MPI_REDUCES
    hvy_id = hvy_active(1, tree_ID) ! for this stage, just pass any block (result does not depend on block)
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, x0, dx, &
    hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "post_stage" )
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

        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
        endif

        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, &
        x0, dx, hvy_rhs(:,:,:,:, hvy_id), hvy_mask(:,:,:,:, hvy_id), "local_stage", n_domain)
    enddo
    call toc( "RHS_wrapper::local-stage", 33, MPI_wtime()-t1 )

    call toc( "RHS_wrapper (TOTAL)", 30, MPI_wtime()-t0 )
end subroutine RHS_wrapper
