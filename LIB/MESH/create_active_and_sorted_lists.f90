!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_active_and_sorted_lists.f90
!> \version 0.5
!> \author msr
!
!> \brief create all active (lgt/hvy) lists, create also sorted list of active
!! light blocks with numerical treecodes
!! \n
!! input:
!!           - light data
!!
!! output:
!!           - list of active blocks
!!           - number of active blocks
!!           - sorted light data with numerical treecodes
!!
!> = log ======================================================================================
!! \n
!! 14/06/17 - create subroutine
!!
! ********************************************************************************************




! ################################################################################
!> Updates active lgt/hvy lists from lgt_block data. 
!> \author PKrah
subroutine create_active_and_sorted_lists_tree( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, create_sorted_list)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)      :: params    !< user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)!< light data array
    integer(kind=ik), intent(inout)     :: lgt_active(:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n        !< number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_n        !< number of active blocks (light data)
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)!< sorted light data with numerical treecodes
    logical, intent(in)                 :: create_sorted_list!< switch for sorted list creation
    !-----------------------------------------------------------------
    
    ! loop variables
    integer(kind=ik)                    :: k, N, heavy_id, block_rank
    ! process rank
    integer(kind=ik)                    :: rank, tree_id, tree_id_idx

    ! reset old lists, use old numbers of active blocks. If the old numbers are
    ! invalid (too small too large), then we delete everything in the lists
    !> \todo Check if resetting the arrays is not a waste of time in any case!
    if (lgt_n>size(lgt_active)) lgt_n = size(lgt_active)
    if (hvy_n>size(hvy_active)) hvy_n = size(hvy_active)

    if (lgt_n<=0) lgt_n = size(lgt_active)
    if (hvy_n<=0) hvy_n = size(hvy_active)

    ! reset the active lists
    lgt_active(1:lgt_n)          = -1
    hvy_active(1:hvy_n)          = -1
    lgt_sortednumlist(1:lgt_n,:) = -1

    ! reset active block numbers
    lgt_n = 0
    hvy_n = 0

    rank = params%rank
    N    = params%number_blocks
    tree_id_idx = params%max_treelevel + idx_tree_id

    ! =======================================================
    ! loop over all light data
    ! =======================================================
    do k = 1, size(lgt_block, 1)
        ! block is active
        if ( lgt_block(k, 1) /= -1 ) then

            ! which tree id has the current block k?
            tree_id = lgt_block(k,tree_id_idx)
            ! ---------------------------
            ! update light active
            ! ---------------------------
            ! save lgt id as active block
            lgt_n      = lgt_n + 1
            lgt_active( lgt_n ) = k

            ! ---------------------------
            ! update hvy active
            ! ---------------------------
            ! save heavy id, only if proc responsable for block
            call lgt_id_to_proc_rank( block_rank, k, N )
            if ( rank == block_rank ) then
                ! convert light data id into heavy data id
                call lgt_id_to_hvy_id( heavy_id, k, rank, N)
                hvy_active( hvy_n + 1 ) = heavy_id
                hvy_n                   = hvy_n + 1
            end if

            if (create_sorted_list) then
                ! sorted list
                ! first index stores the light id of the block
                lgt_sortednumlist(lgt_n, 1) = k
                ! second index stores the numerical treecode
                ! + the tree index
                lgt_sortednumlist(lgt_n, 2) = treecode2int( &
                        lgt_block(k, 1:params%max_treelevel), tree_id )
            end if

        end if
    end do
    ! =======================================================
    if (create_sorted_list) then
        ! sort list
        if (lgt_n > 1) then
            call quicksort(lgt_sortednumlist, 1, lgt_n, 2)
        end if
    end if

end subroutine create_active_and_sorted_lists_tree



















! ################################################################################
!> Updates active lgt/hvy lists from lgt_block data. 
!> Returns active lists for each tree in the forest and lgt_n for each tree.
!
!> \details
!> -------------------------------------------------------------
!>     code                    | explanation
!> -------------------------------------------------------------
!> lgt_active(:,tree_id)       | active block list of tree
!> lgt_n(:,tree_id)            | number of active blocks in tree
!> lgt_active(:,fsize+1 )| list of all active blocks
!> lgt_n(:,fsize+1)      | total number of active blocks
!> -------------------------------------------------------------
!> \author PKrah
subroutine create_active_and_sorted_lists_forest( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, create_sorted_list, tree_n)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)      :: params    !< user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)!< light data array
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)!< number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_n!< number of active blocks (light data)
    integer(kind=ik), intent(out) :: tree_n!< highest tree id
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)!< sorted light data with numerical treecodes
    logical, intent(in)                 :: create_sorted_list!< switch for sorted list creation
    !-----------------------------------------------------------------
    
    ! loop variables
    integer(kind=ik)                    :: k, N, heavy_id, block_rank, fsize
    ! process rank
    integer(kind=ik)                    :: rank, tree_id, tree_id_idx, lgt_n_sum

    rank = params%rank
    N    = params%number_blocks
    fsize= params%forest_size
    tree_id_idx = params%max_treelevel + idx_tree_id

    ! Note: that the total number of active blocks in all trees is stored in 
    ! the last element of lgt_n and can also be calculated from the lgt_n sum
    ! over all active trees
    if (tree_n>1) then
      lgt_n_sum = max(sum( lgt_n(1:tree_n) ),lgt_n(fsize +1))
    else
      lgt_n_sum = lgt_n(fsize+1)
    endif

    ! reset old lists, use old numbers of active blocks. If the old numbers are
    ! invalid (too small too large), then we delete everything in the lists
    !> \todo Check if resetting the arrays is not a waste of time in any case!
    if (lgt_n_sum > size(lgt_active(:,fsize+1))) lgt_n_sum = size(lgt_active(:,fsize+1))
    if (hvy_n>size(hvy_active)) hvy_n = size(hvy_active)

    if (lgt_n_sum <= 0) lgt_n_sum = size(lgt_active(:,fsize+1))
    if (hvy_n<=0) hvy_n = size(hvy_active)

    ! reset the active lists 
    lgt_active = -1
    hvy_active(1:hvy_n)    = -1
    lgt_sortednumlist(1:lgt_n_sum,:) = -1

    ! reset active block numbers
    lgt_n_sum = 0 ! lgt_n_sum = sum(lgt_n(1:tree_n))
    lgt_n = 0
    hvy_n = 0
    tree_n= 0

    ! =======================================================
    ! loop over all light data
    ! =======================================================
    do k = 1, size(lgt_block, 1)
        ! block is active
        if ( lgt_block(k, 1) /= -1 ) then

            ! which tree id has the current block k?
            tree_id = lgt_block(k,tree_id_idx)
            ! find the highest tree number tree_n
            if (tree_id > tree_n) tree_n = tree_id

            ! ---------------------------
            ! update light active
            ! ---------------------------
            ! save lgt id as active block
            lgt_n(tree_id) = lgt_n(tree_id) + 1
            lgt_n_sum      = lgt_n_sum + 1
            lgt_active( lgt_n(tree_id), tree_id) = k
            lgt_active( lgt_n_sum, fsize+1) = k

            ! ---------------------------
            ! update hvy active
            ! ---------------------------
            ! save heavy id, only if proc responsable for block
            call lgt_id_to_proc_rank( block_rank, k, N )
            if ( rank == block_rank ) then
                ! convert light data id into heavy data id
                call lgt_id_to_hvy_id( heavy_id, k, rank, N)
                hvy_active( hvy_n + 1 ) = heavy_id
                hvy_n                   = hvy_n + 1
            end if

            if (create_sorted_list) then
                ! sorted list
                ! first index stores the light id of the block
                lgt_sortednumlist(lgt_n_sum, 1) = k
                ! second index stores the numerical treecode
                lgt_sortednumlist(lgt_n_sum, 2) = treecode2int( &
                        lgt_block(k, 1:params%max_treelevel), tree_id )
            end if

        end if
    end do
    ! =======================================================
    
    lgt_n(fsize + 1) = lgt_n_sum
    if (create_sorted_list) then
        ! sort list
        if (lgt_n_sum > 1) then
            call quicksort(lgt_sortednumlist, 1, lgt_n_sum, 2)
        end if
    end if

    ! check if the number of trees is not bigger then the size of the forest
    ! The forest size is defined as the maximum number of trees in the forest.
    if (tree_n > fsize) call abort(1402192, "To many trees in the forest!!" )

end subroutine create_active_and_sorted_lists_forest


