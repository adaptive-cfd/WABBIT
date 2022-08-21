!> \brief reset grid, set all blocks to empty
subroutine reset_tree(params, lgt_block, lgt_active, lgt_n, &
     hvy_active, hvy_n, lgt_sortednumlist, verbosity, tree_ID )

    implicit none

    type (type_params), intent(in)          :: params                   !> user defined parameter structure
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)          !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_active(:,:)            !> list of active blocks (light data)
    integer(kind=ik), intent(inout)         :: lgt_n(:)                    !> number of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:,:)            !> list of active blocks (light data)
    integer(kind=ik), intent(inout)         :: hvy_n(:)                    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)            :: tree_ID                  !> which tree to reset?
    integer(kind=tsize), intent(inout)      :: lgt_sortednumlist(:,:,:)   !> sorted list of numerical treecodes, used for block finding
    !> write output
    logical, intent(in)                     :: verbosity
    integer(kind=ik)                        :: k

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
    lgt_n(tree_ID) = size(lgt_active,1)
    hvy_n(tree_ID) = size(hvy_active,1)

    call createActiveSortedLists_tree( params, lgt_block, lgt_active, lgt_n, hvy_active, &
         hvy_n, lgt_sortednumlist, tree_ID )

end subroutine reset_tree


!> Resets the light data. After calling this function
!> only one tree is left in the forest, and all blocks are inactive with
!> refinement status 0
subroutine reset_forest(params, lgt_block, lgt_active, lgt_n,&
     hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    implicit none

    type (type_params), intent(in)          :: params                     !> user defined parameter structure
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)            !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_active(:,:)            !> list of active blocks (light data)
    integer(kind=ik), intent(inout),optional:: lgt_n(:)                   !> number of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:,:)            !> list of active blocks (light data)
    integer(kind=ik), intent(inout)         :: hvy_n(:)                   !> number of active blocks (heavy data)
    integer(kind=tsize), intent(inout)      :: lgt_sortednumlist(:,:,:)   !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik),  intent(inout),optional        :: tree_n
    integer(kind=ik)                        :: i, max_treelevel, n_trees

    ! lgt_block = -1 ! Thomas (30-11-2018): do not reset anymore, use 'pseudo-dynamic' memory management
    ! lgt_active = -1
    n_trees = params%forest_size
    max_treelevel = params%max_treelevel
    ! reset data:
    ! all blocks are inactive, reset treecode
    lgt_block(:,:) = -1
    ! reset tree_ID
    lgt_block(:, max_treelevel+IDX_TREE_ID) = -1
    ! all blocks are inactive, reset mesh level
    lgt_block(:, max_treelevel+IDX_MESH_LVL) = -1
    ! set refinement to 0
    lgt_block(:, max_treelevel+IDX_REFINE_STS) = 0

    ! reset sorted list of numerical treecodes
    do i =1,size(lgt_active,2)
            lgt_sortednumlist(:,:,i) = -1
            lgt_n(i) = size(lgt_active(:,i),1)
            hvy_n(i) = size(hvy_active(:,i),1)
    end do

    if (present(tree_n)) then
        call createActiveSortedLists_forest( params, lgt_block, lgt_active, &
               lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
    else
        call createActiveSortedLists_forest( params, lgt_block, lgt_active, &
               lgt_n, hvy_active, hvy_n, lgt_sortednumlist, n_trees)
    end if
end subroutine reset_forest
