!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name refine_mesh.f90
!> \version 0.5
!> \author msr, engels
!> \brief Refine the mash (tag / interpolate / update lists)
!
!> \details This routine first sets the refinement flag for all blocks to +1
!! and then executes the refinement directly. Blocks that cannot be refined because they
!! are already on the finest allowed level are unaltered.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!!
!! Note we assume, on call, that active lists / neighbor are up-to-date
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

subroutine refine_mesh( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, indicator  )

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
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)     :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)           :: indicator

    ! cpu time variables for running time calculation
    real(kind=rk)                          :: t0
    integer(kind=ik)                       :: k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    t0 = MPI_Wtime()

    !> (a) loop over the blocks and set their refinement status.
    call refinement_indicator( params, lgt_block, lgt_active, lgt_n, indicator )

    !> (b) check if block has reached maximal level, if so, remove refinement flags
    call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )

    !> (c) ensure gradedness of mesh. If the refinement is done everywhere, there is
    !! no way gradedness can be damaged, so we skip the call in this case. However,
    !! in all other indicators, this step is very important.
    if ( indicator /= "everywhere") then
      call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )
    endif

    !> (d) execute refinement, interpolate the new mesh. All blocks go one level up
    !! except if they are already on the highest level. Note that those blocks have
    !! the status +11
    if ( params%threeD_case ) then
        ! 3D:
        call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
    else
        ! 2D:
        call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_active, hvy_n )
    end if

    !> (e) as the grid changed now with the refinement, we have to update the list of
    !! active blocks so other routines can loop just over these active blocks
    !! and do not have to ensure that the active list is up-to-date
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
         hvy_active, hvy_n, lgt_sortednumlist, .true. )
    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
         lgt_sortednumlist, hvy_active, hvy_n )

    !> At this point the refinement is done. Since not all blocks are refined, namely only those
    !! that were not on Jmax, Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now
    call balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
    hvy_active, hvy_n, hvy_work )


    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
         hvy_active, hvy_n, lgt_sortednumlist, .true. )
    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
         lgt_sortednumlist, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! End of routine
    call toc( params, "refine_mesh", MPI_wtime()-t0 )

end subroutine refine_mesh
