!> \file
!> \callgraph
!********************************************************************************************
! WABBIT
! ============================================================================================
!> \name     adapt_mesh.f90
!! \version  0.4
!! \author   msr, engels
!
!> \brief This routine performs the coarsing of the mesh, where possible. For the given mesh
!! we compute the details-coefficients on all blocks. If four sister blocks have maximum
!! details below the specified tolerance, (so they are insignificant), they are merged to
!! one coarser block one level below. This process is repeated until the grid does not change
!! anymore.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!
!> \note The block thresholding is done with the restriction/prediction operators acting on the
!! entire block, INCLUDING GHOST NODES. Ghost node syncing is performed in threshold_block.
!
!> \note It is well possible to start with a very fine mesh and end up with only one active
!! block after this routine. You do *NOT* have to call it several times.
!
!> \details
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
!!
!> = log ======================================================================================
!!
!! 10/11/16 - switch to v0.4
! ==========================================================================================
!********************************************************************************************
!> \image html adapt_mesh.svg width=400

subroutine adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, indicator, hvy_synch, hvy_work )

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    real(kind=rk), intent(in)           :: time
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)   :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> coarsening indicator
    character(len=*), intent(in)        :: indicator
    ! loop variables
    integer(kind=ik)                    :: lgt_n_old, iteration, k, max_neighbors
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1, t_misc
    ! MPI error variable
    integer(kind=ik)                    :: ierr
    logical::test

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    t0 = MPI_Wtime()
    t1 = t0
    t_misc = 0.0_rk
    lgt_n_old = 0
    iteration = 0

    if ( params%threeD_case ) then
        max_neighbors = 56
    else
        max_neighbors = 12
    end if

!---------------------------------------------------------------------------------------------
! main body

    !> we iterate until the number of blocks is constant (note: as only coarseing
    !! is done here, no new blocks arise that could compromise the number of blocks -
    !! if it's constant, its because no more blocks are refined)
    do while ( lgt_n_old /= lgt_n )

        lgt_n_old = lgt_n

        !> (a) check where coarsening is possible
        ! ------------------------------------------------------------------------------------
        ! first: synchronize ghost nodes - thresholding on block with ghost nodes
        ! synchronize ghostnodes, grid has changed, not in the first one, but in later loops
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, hvy_synch )

        !! calculate detail on the entire grid. Note this is a wrapper for block_coarsening_indicator, which
        !! acts on a single block only
        call grid_coarsening_indicator( time, params, lgt_block, hvy_block, hvy_work, lgt_active, lgt_n, &
        hvy_active, hvy_n, indicator, iteration, hvy_neighbor, hvy_synch)


        !> (b) check if block has reached maximal level, if so, remove refinement flags
        t0 = MPI_Wtime()
        call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (min/max)", MPI_Wtime()-t0 )


        !> (c) unmark blocks that cannot be coarsened due to gradedness
        t0 = MPI_Wtime()
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (gradedness)", MPI_Wtime()-t0 )


        !> (d) ensure completeness
        t0 = MPI_Wtime()
        call ensure_completeness( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (completeness)", MPI_Wtime()-t0 )


        !> (e) adapt the mesh, i.e. actually merge blocks
        t0 = MPI_Wtime()
        call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (coarse mesh)", MPI_Wtime()-t0 )


        ! the following calls are indeed required (threshold->ghosts->neighbors->active)
        ! update lists of active blocks (light and heavy data)
        ! update list of sorted nunmerical treecodes, used for finding blocks
        t0 = MPI_Wtime()
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        t_misc = t_misc + (MPI_Wtime() - t0)


        t0 = MPI_Wtime()
        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (update neighbors)", MPI_Wtime()-t0 )

        iteration = iteration + 1
    end do

    ! The grid adaptation is done now, the blocks that can be coarsened are coarser.
    ! If a block is on Jmax now, we assign it the status +11.
    ! NOTE: Consider two blocks, a coarse on Jmax-1 and a fine on Jmax. If you refine only
    ! the coarse one (Jmax-1 -> Jmax), because you cannot refine the other one anymore
    ! (by defintion of Jmax), then the redundant layer in both blocks is different.
    ! To corrent that, you need to know which of the blocks results from interpolation and
    ! which one has previously been at Jmax. This latter one gets the 11 status.
    do k = 1, lgt_n
        if ( lgt_block( lgt_active(k), params%max_treelevel+1) == params%max_treelevel ) then
            lgt_block( lgt_active(k), params%max_treelevel+2 ) = 11
        end if
    end do

    !> At this point the coarsening is done. All blocks that can be coarsened are coarsened
    !! they may have passed several level also. Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now
    call balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, hvy_work )


    !> load balancing destroys the lists again, so we have to create them one last time to
    !! end on a valid mesh
    !! update lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    t0 = MPI_wtime()
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
    t_misc = t_misc + (MPI_Wtime() - t0)


    ! update neighbor relations
    t0 = MPI_Wtime()
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
    ! CPU timing (only in debug mode)
    call toc( params, "adapt_mesh (update neighbors) ", MPI_Wtime()-t0 )


    ! time remaining parts of this routine.
    call toc( params, "adapt_mesh (...)", t_misc )
    call toc( params, "adapt_mesh (TOTAL)", MPI_wtime()-t1)
end subroutine adapt_mesh





! ============================================================================================
!> \name coarsening_indicator.f90
!> \version 0.5
!> \author engels
!> \brief Set coarsening status for all active blocks, different methods possible
!
!> \details This routine sets the coarsening flag for all blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very complex, but expected to grow
!! in the future.
!! \n
!! ------------------ \n
!! Refinement status: \n
!! ------------------ \n
!! +1 refine \n
!! 0 do nothing \n
!! -1 block wants to coarsen (ignoring other constraints, such as gradedness) \n
!! -2 block will coarsen and be merged with her sisters \n
!! ------------------ \n
!! \n
!! = log ======================================================================================
!! \n
!! 29/05/2018 create
! ********************************************************************************************
subroutine grid_coarsening_indicator( time, params, lgt_block, hvy_block, hvy_work, lgt_active, lgt_n, &
  hvy_active, hvy_n, indicator, iteration, hvy_neighbor, hvy_synch)
  !---------------------------------------------------------------------------------------------
  ! modules
    use module_indicators

    implicit none
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    !> coarsening iteration index. coarsening is done until the grid has reached
    !! the steady state; therefore, this routine is called several times during the
    !! mesh adaptation. Random coarsening (used for testing) is done only in the first call.
    integer(kind=ik), intent(in)        :: iteration
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)


    ! local variables
    integer(kind=ik) :: k, Jmax, neq, lgt_id, Bs, g, mpierr
    ! local block spacing and origin
    real(kind=rk) :: dx(1:3), x0(1:3), tmp(1:params%number_data_fields)
    real(kind=rk) :: norm(1:params%number_data_fields)


    Jmax = params%max_treelevel
    neq = params%number_data_fields
    Bs = params%number_block_nodes
    g = params%number_ghost_nodes

    !> reset refinement status to "stay" on all blocks
    do k = 1, lgt_n
        lgt_block( lgt_active(k), Jmax+2 ) = 0
    enddo


    !> Compute vorticity (if required)
    !! it is a little unfortunate that we have to compute the vorticity here and not in
    !! block_coarsening_indicator, where it should be. but after computing it, we have to synch
    !! its vorticity ghost nodes in order to apply the detail operator to the entire
    !! vorticity field (incl gost nodes)
    if (params%coarsening_indicator=="threshold-vorticity") then
        ! loop over my active hvy data
        do k = 1, hvy_n
            ! get lgt id of block
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

            ! some indicators may depend on the grid (e.g. to compute the vorticity), hence
            ! we pass the spacing and origin of the block
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            ! actual computation of vorticity on the block
            call compute_vorticity(hvy_block(:,:,:,1,hvy_active(k)), hvy_block(:,:,:,2,hvy_active(k)), hvy_block(:,:,:,3,hvy_active(k)), &
            dx, Bs, g, params%order_discretization, hvy_work(:,:,:,1:3,hvy_active(k)))
        enddo

        ! note here we synch hvy_work (=vorticity) and not hvy_block
        call sync_ghosts( params, lgt_block, hvy_work, hvy_neighbor, hvy_active, hvy_n, hvy_synch )
    endif

    !> Compute normalization for eps, if desired.
    !! versions <30.05.2018 used fixed eps for all qtys of the state vector, but that is not very smart
    !! as each qty can have different mangitudes. If the switch eps_normalized is on, we compute here
    !! the vector of normalization factors for each qty that adaptivity will bebased on (state vector
    !! or vorticity). We currently use the L_infty norm. I have made bad experience with L_2 norm
    !! (to be checked...)
    if (params%eps_normalized .and. params%physics_type=="ACM-new") then
        norm = 0.0_rk
        if (params%coarsening_indicator=="threshold-vorticity") then
            ! loop over my active hvy data
            do k = 1, hvy_n
                norm(1) = max(norm(1),  maxval(abs(hvy_work(:, :, :, 1, hvy_active(k)))) )
            enddo
        else
            do k = 1, hvy_n
                ! call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
                ! call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                ! norm(1) = norm(1) + sum(hvy_block(:,:,:,1,hvy_active(k))**2)*dx(1)*dx(2)

                ! max over ux, uy
                norm(1) = max(norm(1), maxval(abs(hvy_block(:,:,:,1:2,hvy_active(k)))) )
                ! pressure
                norm(3) = max(norm(3), maxval(abs(hvy_block(:,:,:,3,hvy_active(k)))) )
            enddo
            ! isotropy: uy=ux
            norm(2) = norm(1)
            ! norm(2) = sqrt( norm(1) )
        endif

        tmp = norm
        call MPI_ALLREDUCE(tmp, norm, neq, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
    else
        norm = 1.0_rk
    endif

    ! during dev is it useful to know what the normalization is, if that is actives
    if (params%rank == 0 .and. params%eps_normalized) then
        open(14,file='eps_norm.t',status='unknown',position='append')
        write (14,'(5(g15.8,1x))') time, norm, params%eps
        close(14)
    endif

    !> evaluate coarsing criterion on all blocks
    ! loop over all my blocks
    do k = 1, hvy_n
        ! get lgt id of block
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

        ! some indicators may depend on the grid, hence
        ! we pass the spacing and origin of the block (as we have to compute vorticity
        ! here, this can actually be omitted.)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! evaluate the criterion on this block.
        call block_coarsening_indicator( params, hvy_block(:,:,:,1:neq,hvy_active(k)), &
        hvy_work(:,:,:,1:neq,hvy_active(k)), dx, x0, indicator, iteration, &
        lgt_block(lgt_id, Jmax+2), norm )
    enddo


    !> after modifying all refinement statusses, we need to synchronize light data
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.true. )

end subroutine grid_coarsening_indicator
