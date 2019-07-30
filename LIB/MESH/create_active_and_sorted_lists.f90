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
!> Updates active lgt/hvy lists from lgt_block data FOR ONE TREE
!> \author PKrah
subroutine create_active_and_sorted_lists_tree( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)      :: params    !< user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)!< light data array
    integer(kind=ik), intent(inout)     :: lgt_active(:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n        !< number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_n        !< number of active blocks (light data)
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)!< sorted light data with numerical treecodes
    !-----------------------------------------------------------------

    ! loop variables
    integer(kind=ik)                    :: k, N, heavy_id, block_rank
    ! process rank
    integer(kind=ik)                    :: rank, TREE_ID_IDX
    real(kind=rk) :: t0

    t0 = MPI_wtime()


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
    TREE_ID_IDX = params%max_treelevel + IDX_TREE_ID

    ! =======================================================
    ! loop over all light data
    ! =======================================================
    do k = 1, size(lgt_block, 1)
        ! block is active
        if ( lgt_block(k, 1) /= -1 .and. lgt_block(k,TREE_ID_IDX)==tree_ID) then
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

            ! sorted list
            ! first index stores the light id of the block
            lgt_sortednumlist(lgt_n, 1) = k
            ! second index stores the numerical treecode
            ! + the tree index
            lgt_sortednumlist(lgt_n, 2) = treecode2int(lgt_block(k, 1:params%max_treelevel), tree_id )

        end if
    end do
    ! =======================================================
    ! sort list
    if (lgt_n > 1) then
        call quicksort(lgt_sortednumlist, 1, lgt_n, 2)
    end if
    call toc("create_active_and_sorted_lists_tree", MPI_wtime()-t0)
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
!> -------------------------------------------------------------
!> \author PKrah
subroutine create_active_and_sorted_lists_forest( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)      :: params    !< user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)!< light data array
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)!< number of ACTIVE blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_active(:,:)!< list of active blocks for each tree (light data)
    integer(kind=ik), intent(inout)     :: hvy_n(:)!< number of ACTIVE blocks (light data) in each tree
    integer(kind=ik), intent(inout)     :: tree_n!< number of ACTIVE trees in forest
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:,:)!< sorted light data with numerical treecodes
    !-----------------------------------------------------------------

    ! loop variables
    integer(kind=ik)                    :: k, N, heavy_id, block_rank, fsize
    ! process rank
    integer(kind=ik)                    :: rank, tree_id, TREE_ID_IDX, lgt_n_sum, hvy_n_sum
    integer(kind=tsize) :: treecode_int
    real(kind=rk) :: t0

    t0 = MPI_wtime()


    rank = params%rank
    N    = params%number_blocks
    fsize = params%forest_size
    TREE_ID_IDX = params%max_treelevel + IDX_TREE_ID

    ! =======================================================
    ! Reset active lists of all trees
    ! =======================================================
    ! note: this seems to be a complicated way of reseting the
    !       active lists, but it is very crucial for performance!
    !       NEVER RESET the full array without reasons!!!
    do tree_id = 1, tree_n 
      ! check if lgt_n or hvy_n of tree is valid (not too small or too large)
      if (lgt_n(tree_id)>size(lgt_active(:,tree_id)) .or. lgt_n(tree_id) <=0) &
        lgt_n(tree_id) = size(lgt_active(:,tree_id))
      if (hvy_n(tree_id)>size(hvy_active(:,tree_id)) .or. hvy_n(tree_id) <=0) &
        hvy_n(tree_id) = size(hvy_active(:,tree_id))

      ! reset the active lists
      lgt_active(1:lgt_n(tree_id),tree_id) = -1
      hvy_active(1:hvy_n(tree_id),tree_id) = -1
      lgt_sortednumlist(1:lgt_n(tree_id),:,tree_id) = -1
    end do


    ! reset active block numbers
    lgt_n = 0
    hvy_n = 0
    tree_n = 0

    ! =======================================================
    ! loop over all light data
    ! =======================================================
    do k = 1, size(lgt_block, 1)
        ! block is active
        if ( lgt_block(k, 1) /= -1 ) then

            ! which tree id has the current block k?
            tree_id = lgt_block(k, TREE_ID_IDX)
            ! find the highest tree number. this is should be the same as
            ! the number of active trees: tree_n
            tree_n = max(tree_id, tree_n)

            ! ---------------------------
            ! update light active
            ! ---------------------------
            ! save lgt id as active block
            lgt_n(tree_id) = lgt_n(tree_id) + 1
            lgt_active( lgt_n(tree_id), tree_id) = k

            ! ---------------------------
            ! update hvy active
            ! ---------------------------
            ! save heavy id, only if proc responsable for block
            call lgt_id_to_proc_rank( block_rank, k, N )
            if ( rank == block_rank ) then
                ! convert light data id into heavy data id
                call lgt_id_to_hvy_id( heavy_id, k, rank, N)
                hvy_n(tree_id) = hvy_n(tree_id) + 1
                hvy_active( hvy_n(tree_id) , tree_id ) = heavy_id
            end if

            ! sorted list
            treecode_int = treecode2int(lgt_block(k, 1:params%max_treelevel), tree_id )
            ! first index stores the light id of the block
            lgt_sortednumlist(lgt_n(tree_id), 1, tree_id) = k
            ! second index stores the numerical treecode
            lgt_sortednumlist(lgt_n(tree_id), 2, tree_id) = treecode_int

        end if
    end do

    ! =======================================================
    ! sort list of every single tree
    do tree_id = 1, tree_n
        if (lgt_n(tree_id) > 1) then
            call quicksort(lgt_sortednumlist(:,:,tree_id), 1, lgt_n(tree_id), 2)
        end if
    end do

    ! check if the number of trees is not bigger then the size of the forest
    ! The forest size is defined as the maximum number of trees in the forest.
    if (tree_n > fsize) call abort(1402192, "Too many trees in the forest!!" )

    call toc("create_active_and_sorted_lists_forest", MPI_wtime()-t0)
end subroutine create_active_and_sorted_lists_forest
