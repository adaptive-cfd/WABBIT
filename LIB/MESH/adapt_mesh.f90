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

subroutine adapt_mesh( params, lgt_block, hvy_block, hvy_synch, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, indicator, com_lists, com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

!---------------------------------------------------------------------------------------------
! modules
    use module_indicators

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy synch array
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)
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
    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)
    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)
    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)         :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! loop variables
    integer(kind=ik)                    :: lgt_n_old, iteration, k, max_neighbors

    ! MPI error variable
    integer(kind=ik)                    :: ierr

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0

!---------------------------------------------------------------------------------------------
! variables initialization

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

        ! timing
        sub_t0 = MPI_Wtime()

        lgt_n_old = lgt_n

        !> (a) check where coarsening is possible
        ! ------------------------------------------------------------------------------------
        ! first: synchronize ghost nodes - thresholding on block with ghost nodes
        ! synchronize ghostnodes, grid has changed, not in the first one, but in later loops
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_synch, hvy_neighbor, hvy_active, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:), com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

        ! timing
        call toc( params, "-adapt mesh: synch ghost", MPI_wtime()-sub_t0, .true. )
        sub_t0 = MPI_Wtime()

        ! calculate detail
        call coarsening_indicator( params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_active, hvy_n, indicator, iteration)

        !> (b) check if block has reached maximal level, if so, remove refinement flags
        call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )

        !> (c) unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )

        !> (d) ensure completeness
        call ensure_completeness( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )

        ! timing
        call toc( params, "-adapt mesh: indicator, gradedness, completeness", MPI_wtime()-sub_t0, .true. )
        sub_t0 = MPI_Wtime()

        !> (e) adapt the mesh, i.e. actually merge blocks
        call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )

        ! timing
        call toc( params, "-adapt mesh: coarse mesh", MPI_wtime()-sub_t0, .true. )
        sub_t0 = MPI_Wtime()

        ! the following calls are indeed required (threshold->ghosts->neighbors->active)
        ! update lists of active blocks (light and heavy data)
        ! update list of sorted nunmerical treecodes, used for finding blocks
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )

        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

        ! timing
        call toc( params, "-adapt mesh: update neighbors", MPI_wtime()-sub_t0, .true. )
        sub_t0 = MPI_Wtime()

        iteration = iteration + 1
    end do

!    !> At this point the coarsening is done. All blocks that can be coarsened are coarsened
!    !! they may have passed several level also. Now, the distribution of blocks may no longer
!    !! be balanced, so we have to balance load now
!    call balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )
!
!
!    !> load balancing destroys the lists again, so we have to create them one last time to
!    !! end on a valid mesh
!    !! update lists of active blocks (light and heavy data)
!    ! update list of sorted nunmerical treecodes, used for finding blocks
!    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
!
!    ! update neighbor relations
!    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

end subroutine adapt_mesh
