!> \brief reset grid, set all blocks to empty
subroutine reset_tree(params, verbosity, tree_ID )
    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas 
    implicit none

    type (type_params), intent(in)          :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)            :: tree_ID                  !> which tree to reset?
    logical, intent(in)                     :: verbosity
    integer(kind=ik)                        :: k

    if ( (params%rank == 0) .and. verbosity ) then
      write(*,'(A)') "RESET-GRID: resetting grid to empty (deactivate all blocks)."
    endif

    do k = 1, size(lgt_block,1)
        if (lgt_block(k,params%Jmax+IDX_TREE_ID) == tree_ID) then
            ! this block is a part of this tree -> delete it
            lgt_block(k,:) = -1
        endif
    enddo


    ! as the grid has changed (we deleted it here), we now update the heavy and light
    ! active lists
    ! update list of sorted nunmerical treecodes, used for finding blocks

    ! set number of active blocks to maximum number, to reset everything
    lgt_n(tree_ID) = size(lgt_active,1)
    hvy_n(tree_ID) = size(hvy_active,1)

    call createActiveSortedLists_tree( params, tree_ID )

end subroutine reset_tree


!> Resets the light data. After calling this function
!> only one tree is left in the forest, and all blocks are inactive with
!> refinement status 0
subroutine reset_forest(params)
    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas
    implicit none

    type (type_params), intent(in)          :: params                     !> user defined parameter structure

    integer(kind=ik)                        :: i, max_treelevel, n_trees

    ! lgt_block = -1 ! Thomas (30-11-2018): do not reset anymore, use 'pseudo-dynamic' memory management
    ! lgt_active = -1
    n_trees = params%forest_size
    max_treelevel = params%Jmax
    ! reset data:
    ! all blocks are inactive, reset treecode
    lgt_block(:,:) = -1
    ! reset tree_ID
    lgt_block(:, max_treelevel+IDX_TREE_ID) = -1
    ! all blocks are inactive, reset mesh level
    lgt_block(:, max_treelevel+IDX_MESH_LVL) = -1
    ! set refinement to 0
    lgt_block(:, max_treelevel+IDX_REFINE_STS) = 0

    ! set TC to -1
    lgt_block(:, max_treelevel+IDX_TC_1) = -1
    lgt_block(:, max_treelevel+IDX_TC_2) = -1

    ! reset sorted list of numerical treecodes
    do i = 1,size(lgt_active,2)
            lgt_sortednumlist(:,:,i) = -1
            lgt_n(i) = size(lgt_active(:,i),1)
            hvy_n(i) = size(hvy_active(:,i),1)
    end do

    call createActiveSortedLists_forest(params)
end subroutine reset_forest
