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

subroutine adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, indicator )

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
    integer(kind=ik)                    :: lgt_n_old, iteration, k

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    sub_t0 = MPI_Wtime()

    time_sum = 0.0_rk

    lgt_n_old = 0
    iteration = 0

!---------------------------------------------------------------------------------------------
! main body

    !> we iterate until the number of blocks is constant (note: as only coarseing
    !! is done here, no new blocks arise that could compromise the number of blocks -
    !! if it's constant, its because no more blocks are refined)
    do while ( lgt_n_old /= lgt_n )

        lgt_n_old = lgt_n

        ! end time
        sub_t1 = MPI_Wtime()
        time_sum = time_sum + (sub_t1 - sub_t0)

        !> (a) check where coarsening is possible
        call coarsening_indicator( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, indicator, iteration )

        ! start time
        sub_t0 = MPI_Wtime()

        !> (b) check if block has reached maximal level, if so, remove refinement flags
        call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
        !> (c) unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )
        !> (d) ensure completeness
        call ensure_completeness( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )

        !> (e) adapt the mesh, i.e. actually merge blocks
        call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )

        ! the following calls are indeed required (threshold->ghosts->neighbors->active)
        ! update lists of active blocks (light and heavy data)
        call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
        call create_hvy_active_list( lgt_block, hvy_active, hvy_n )
        ! update list of sorted nunmerical treecodes, used for finding blocks
        call create_lgt_sortednumlist( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )
        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        iteration = iteration + 1

    end do

    ! end time
    sub_t1 = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)

    !> At this point the coarsening is done. All blocks that can be coarsened are coarsened
    !! they may have passed several level also. Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now
    if ( params%threeD_case ) then
        ! 3D:
        call balance_load_3D( params, lgt_block, hvy_block, lgt_active, lgt_n )
    else
        ! 2D:
        call balance_load_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )
    end if

    ! start time
    sub_t0 = MPI_Wtime()

    !> load balancing destroys the lists again, so we have to create them one last time to
    !! end on a valid mesh
    !! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_lgt_sortednumlist( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )
    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    ! end time
    sub_t1 = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "adapt_mesh (w/o load balanc./thresh.)" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "adapt_mesh (w/o load balanc./thresh.)"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + time_sum
    end if

end subroutine adapt_mesh
