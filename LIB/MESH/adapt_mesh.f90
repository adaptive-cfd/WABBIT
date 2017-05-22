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
!> \image html adapt_mesh.png width=400

subroutine adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, indicator )

!---------------------------------------------------------------------------------------------
! modules

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
    integer(kind=ik)                    :: j, ierr, Jmax, lgt_n_old, iteration
    ! random variable for coarsening
    real(kind=rk)                       :: r

!---------------------------------------------------------------------------------------------
! variables initialization
  Jmax = params%max_treelevel
  lgt_n_old = 0
  iteration = 0
!---------------------------------------------------------------------------------------------
! main body

    ! we iterate until the number of blocks is constant (note: as only coarseing
    ! is done here, no new blocks arise that could compromise the number of blocks -
    ! if it's constant, its because no more blocks are refined)
    do while ( lgt_n_old /= lgt_n )
        lgt_n_old = lgt_n

        ! check where to coarsen (refinement done with safety zone)
        if ( indicator == "threshold") then
          ! use wavelet indicator to check where to coarsen. threshold_block performs
          ! the required ghost node sync and loops over all active blocks.
          call threshold_block( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

        elseif (indicator == "random") then
          ! randomly coarse some blocks. used for testing. note we tag for coarsening
          ! only once in the first iteration.
          if (iteration == 0) then
            call init_random_seed()
            ! unset all refinement flags
            lgt_block( :,Jmax+2 ) = 0
            ! only root rank sets the flag, then we sync. It is messy if all procs set a
            ! random value which is not sync'ed
            if (params%rank == 0) then
              do j = 1, lgt_n
                ! random number
                call random_number(r)
                ! set refinement status to coarsen
                if ( r <= 0.25_rk ) then
                    lgt_block( lgt_active(j), Jmax+2 ) = -1
                end if
              end do
            endif
            ! sync light data, as only root sets random coarsening
            call MPI_BCAST( lgt_block(:,params%max_treelevel+2), size(lgt_block,1), MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr )
          endif
        else
            call error_msg("ERROR: unknown coarsening operator")

        endif

        ! update lists of active blocks (light and heavy data)
        call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
        ! hvy_active list is required for update_neighbors
        call create_hvy_active_list( lgt_block, hvy_active, hvy_n )
        ! update list of sorted nunmerical treecodes, used for finding blocks
        call create_lgt_sortednumlist( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )
        ! update neighbor relations, required for gradedness
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

        ! unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )

        ! ensure completeness
        call ensure_completeness( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )

        ! adapt the mesh
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

    ! At this point the coarsening is done. All blocks that can be coarsened are coarsened
    ! they may have passed several level also. Now, the distribution of blocks may no longer
    ! be balanced, so we have to balance load now
    if ( params%threeD_case ) then
        ! 3D:
        call balance_load_3D( params, lgt_block, hvy_block, lgt_active, lgt_n )
    else
        ! 2D:
        call balance_load_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )
    end if

    ! load balancing destroys the lists again, so we have to create them one last time to
    ! end on a valid mesh
    ! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_lgt_sortednumlist( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )
    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

end subroutine adapt_mesh
