!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name refine_mesh.f90
!> \version 0.5
!> \author msr, engels
!
!> \brief This routine first sets the refinement flag for all blocks to +1
!! and then executes the refinement directly. Blocks that cannot be refined because they
!! are already on the finest allowed level are unaltered.
!! 
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!
!> 
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
!! \n
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 03/02/17 - 3D heavy data \n
!! 04/04/17 - include active lists, neighbor relations and load balancing, symmetrical to adapt_mesh \n
!! 05/04/17 - Provide an interface to use different criteria for refinement, rename routines
! ********************************************************************************************

subroutine refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, indicator )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)        :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)           :: indicator

    ! loop variables
    integer(kind=ik)                    :: k, ierr
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1
    ! chance for block refinement, random number
    real(kind=rk)                       :: ref_chance, r
    ! shortcuts for levels
    integer(kind=ik)                    :: Jmin, Jmax, max_blocks, d
!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization
    Jmin = params%min_treelevel
    Jmax = params%max_treelevel
    ! set data dimension
    if ( params%threeD_case ) then
        d = 3
    else
        d = 2
    endif


!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! loop over the blocks and set their refinement status.
    ! NOTE: refinement is an absolute statement, that means once set, the block will be refined
    ! (which is not the case in block coarsening), it may even entrail other blocks in
    ! its vicinity to be refined as well.
    select case (indicator)
        case ("everywhere")
          ! set status "refine" for all active blocks, which is just setting the
          ! last index in the light data block list to +1. This indicator is used
          ! to refine the entire mesh at the beginning of a time step, if error
          ! control is desired.
          do k = 1, lgt_n
              lgt_block( lgt_active(k), params%max_treelevel+2 ) = +1
          end do

      case ("random")
          ! randomized refinement. This can be used to generate debug meshes for
          ! testing purposes. For example the unit tests use that
          ref_chance = 0.25_rk
          ! random refinement can set at most this many blocks to refine (avoid errors
          ! sue to insufficient memory) (since we already have lgt_n blocks we can set the status
          ! at most for Nmax-lgt_n blocks, whcih genrate ech 2**d new blocks)
          max_blocks = (size(lgt_block,1)-lgt_n-10) / 2**d
          ! set random seed
          call init_random_seed()
          ! unset all refinement flags
          lgt_block( :,Jmax+2 ) = 0
          ! only root rank sets the flag, then we sync. It is messy if all procs set a
          ! random value which is not sync'ed
          if (params%rank == 0) then
            do k = 1, lgt_n
              ! random number
              call random_number(r)
              ! set refinement status to refine
              if ( r <= ref_chance .and. sum(lgt_block( :,Jmax+2 )) <= max_blocks) then
                  lgt_block( lgt_active(k), Jmax+2 ) = 1
              else
                  lgt_block( lgt_active(k), Jmax+2 ) = 0
              end if
            end do
          endif
          ! sync light data, as only root sets random refinement
          call MPI_BCAST( lgt_block(:,params%max_treelevel+2), size(lgt_block,1), MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr )

      case default
          call error_msg("ERROR: refine_mesh: the refinement indicator is unkown")

    end select

    ! check if block has reached maximal level, if so, remove refinement flags
    call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )

    ! ensure gradedness of mesh. If the refinement is done everywhere, there is
    ! no way gradedness can be damaged, so we skip the call in this case. However,
    ! in all other indicators, this step is very important.
    if ( indicator /= "everywhere") then
      call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )
    endif

    ! execute refinement, interpolate the new mesh. All blocks go one level up
    ! except if they are already on the highest level.
    ! FIXME: For consistency, it would be better to always refine (allowing one level
    ! beyond maxlevel), but afterwards coarsen to fall back to maxlevel again
    if ( params%threeD_case ) then
        ! 3D:
        call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
    else
        ! 2D:
        call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_active, hvy_n )
    end if

    ! as the grid changed now with the refinement, we have to update the list of
    ! active blocks so other routines can loop just over these active blocks
    ! and do not have to ensure that the active list is up-to-date
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

    ! balance load
    if ( params%threeD_case ) then
        ! 3D:
        call balance_load_3D( params, lgt_block, hvy_block, lgt_active, lgt_n )
    else
        ! 2D:
        call balance_load_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )
    end if

    ! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

    ! end time
    sub_t1 = MPI_Wtime()

    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "refine_mesh" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "refine_mesh"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine refine_mesh
