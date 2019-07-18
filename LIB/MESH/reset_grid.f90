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
subroutine reset_tree(params, lgt_block, lgt_active, lgt_n,&
     hvy_active, hvy_n, lgt_sortednumlist, verbosity, tree_ID )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)       :: params
    !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)         :: lgt_n
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)         :: hvy_n
    !> which tree to reset?
    integer(kind=ik), intent(in)            :: tree_ID
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)      :: lgt_sortednumlist(:,:)
    !> write output
    logical, intent(in)                     :: verbosity
    integer(kind=ik) :: k

    if ( (params%rank == 0) .and. verbosity ) then
      write(*,'(A)') "RESET-GRID: resetting grid to empty (deactivate all blocks)."
    endif

    do k = 1, size(lgt_block,1)
        if (lgt_block(k,params%max_treelevel+IDX_TREE_ID) == tree_ID) then
            ! this block is a part of this tree -> delete it
            lgt_block(k,:) = -1
        endif
    enddo


    ! as the grid has changed (we deleted it here), we now update the heavy and light
    ! active lists
    ! update list of sorted nunmerical treecodes, used for finding blocks

    ! set number of active blocks to maximum number, to reset everything
    lgt_n = size(lgt_active,1)
    hvy_n = size(hvy_active,1)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, &
         hvy_n, lgt_sortednumlist, tree_ID )

end subroutine reset_tree




!> Resets the light data. After calling this function
!> only one tree is left in the forest, and all blocks are inactive with
!> refinement status 0
subroutine reset_forest(lgt_block, lgt_active, max_treelevel, lgt_n, lgt_sortednumlist)

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
    lgt_block(:,:) = -1
    ! reset tree_id
    lgt_block(:, max_treelevel+IDX_TREE_ID) = -1
    ! all blocks are inactive, reset mesh level
    lgt_block(:, max_treelevel+IDX_MESH_LVL) = -1
    ! set refinement to 0
    lgt_block(:, max_treelevel+IDX_REFINE_STS) = 0

    ! reset sorted list of numerical treecodes
    if ( present(lgt_sortednumlist) ) then
      lgt_sortednumlist = -1
    end if
    if ( present(lgt_n) ) then
      lgt_n = size(lgt_active,1)
    end if

end subroutine reset_forest
