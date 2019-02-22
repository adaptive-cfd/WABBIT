! ----------------------------------------------------------------------------------------
!> \file
!> \brief Reset hvy and ligth data arrays
!> \version 0.5
!> \author msr
!! \details
!! \date 04/11/16 - switch to v0.4, now run complete initialization
!!            within these subroutine and return initialized block data to main program
!! \date 07/12/16 - now uses heavy work data array
!! \date 25/01/17 - switch to 3D, v0.5
! ----------------------------------------------------------------------------------------

!> \brief reset grid, set all blocks to empty
subroutine reset_grid(params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active, lgt_n,&
     hvy_active, hvy_n, lgt_sortednumlist, verbosity )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)       :: params
    !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: 00, k1, k2 etc)
    real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik),  intent(inout)        :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)         :: lgt_n
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)         :: hvy_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)      :: lgt_sortednumlist(:,:)
    !> write output
    logical, intent(in)                     :: verbosity
! ----------------------------------------------------------------------------------------

    if ( (params%rank == 0) .and. verbosity ) then
      write(*,'(A)') "RESET-GRID: resetting grid to empty (deactivate all blocks)."
    endif

    call reset_lgt_data(lgt_block, lgt_active, params%max_treelevel, lgt_n, lgt_sortednumlist)


    ! reset data
    ! hvy_block = 9.99e99_rk ! Thomas (30-11-2018): do not reset anymore, use 'pseudo-dynamic' memory management
    ! hvy_work = 9.99e99_rk
    ! hvy_tmp =  9.99e99_rk
    
    hvy_neighbor = -1
    ! as the grid has changed (we deleted it here), we now update the heavy and light
    ! active lists
    ! update list of sorted nunmerical treecodes, used for finding blocks
    ! set number of active blocks to maximum number, to reset everything
    lgt_n = size(lgt_active,1)
    hvy_n = size(hvy_active,1)
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, &
         hvy_n, lgt_sortednumlist, .true. )

end subroutine reset_grid




!> Resets the light data. After calling this function
!> only one tree is left in the forest, and all blocks are inactive with
!> refinement status 0
subroutine reset_lgt_data(lgt_block, lgt_active,max_treelevel, lgt_n, lgt_sortednumlist)

    implicit none

    !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout),optional:: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout),optional :: lgt_sortednumlist(:,:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik)  , intent(in)         :: max_treelevel


    ! lgt_block = -1 ! Thomas (30-11-2018): do not reset anymore, use 'pseudo-dynamic' memory management
    ! lgt_active = -1

    ! reset data:
    ! all blocks are inactive, reset treecode
    lgt_block(:,1:max_treelevel)   = -1
    ! start with only one tree in the forest (i.e. all blocks are associated with only one tree)
    lgt_block(:, max_treelevel+idx_tree_id) = -1
    ! all blocks are inactive, reset mesh level
    lgt_block(:, max_treelevel+idx_mesh_lvl) = -1
    ! set refinement to 0
    lgt_block(:, max_treelevel+idx_refine_sts) = 0

    ! reset sorted list of numerical treecodes
    if ( present(lgt_sortednumlist) ) then
      lgt_sortednumlist = -1
    end if
    if ( present(lgt_n) ) then
      lgt_n = size(lgt_active,1)
    end if

end subroutine reset_lgt_data
