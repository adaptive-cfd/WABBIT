!> \brief create all active (lgt/hvy) lists, create also sorted list of active
!! light blocks with numerical treecodes
!! input:    - light data
!! output:   - list of active blocks
!!           - number of active blocks
!!           - sorted light data with numerical treecodes
! ********************************************************************************************
!> Updates active lgt/hvy lists from lgt_block data FOR ONE TREE by looping
!> over the hvy ids and synchronizing the active lists with the other procs
subroutine create_active_and_sorted_lists_tree_comm( params, lgt_block, lgt_active, &
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
    integer(kind=ik)                    :: k, N, hvy_id, block_rank, ierr, lgt_id
    ! process rank
    integer(kind=ik)                    :: rank, TREE_ID_IDX, N_blocks_tot, mpisize
    real(kind=rk) :: t0,t(5)

    integer(kind=ik), allocatable, save     :: proc_lgt_n(:), proc_lgt_start(:)

    t0 = MPI_wtime()

    mpisize = params%number_procs
    rank = params%rank
    N    = params%number_blocks
    TREE_ID_IDX = params%max_treelevel + IDX_TREE_ID

    if (.not.allocated(proc_lgt_n)) allocate( proc_lgt_n(1:mpisize) )
    if (.not.allocated(proc_lgt_start)) allocate( proc_lgt_start(1:mpisize) )

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
    t(1) = MPI_wtime()

    ! =======================================================
    ! loop over all hvy_data
    ! =======================================================
    do hvy_id = 1, N
        call hvy2lgt( k, hvy_id, rank, N )
        ! block is active
        if (lgt_block(k,TREE_ID_IDX)==tree_ID) then
             ! assuming that there are more blocks per rank then blocks with the corresponding tree ID
            if ( lgt_block(k, 1) /= -1) then
            ! ---------------------------
            ! update hvy active
            ! ---------------------------
            ! save heavy id, only if proc responsable for block
            hvy_active( hvy_n + 1 ) = hvy_id
            hvy_n                   = hvy_n + 1
            end if
        end if
    end do
    t(2) = MPI_wtime()
    ! ==========================================================================
    ! MPI Synchronization
    ! ==========================================================================
    ! Note: lgt_n is hear clearly the same as hvy_n since we have looped
    ! over lgt ids of this proc
    !n sum up lgt_n over all procs
    call MPI_allgather(hvy_n, 1, MPI_INTEGER4, proc_lgt_n, 1, MPI_INTEGER4, WABBIT_COMM, ierr )
    ! this is the global buffer size
    lgt_n = sum(proc_lgt_n)
    do k = 1, mpisize
        proc_lgt_start(k) = sum(proc_lgt_n(1:k-1))! + 1
    enddo

    call MPI_allgatherv( rank*N + hvy_active(1:hvy_n), hvy_n, MPI_INTEGER4, &
    lgt_active( 1:lgt_n), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
    WABBIT_COMM, ierr)

    ! call MPI_allgatherv( lgt_sortednumlist(1:lgt_n, 2), lgt_n, MPI_INTEGER4, &
    ! my_lgt_recv_buffer( 1:lgt_n,2 ), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
    ! WABBIT_COMM, ierr)
    t(3) = MPI_wtime()
    ! UNPACK SYNCHRONIZED DATA
    do k = 1, lgt_n
        lgt_id = lgt_active(k)
        lgt_sortednumlist(k, 1) = lgt_id
        lgt_sortednumlist(k, 2) = treecode2int(lgt_block(lgt_id, 1:params%max_treelevel), tree_id )
        lgt_active(k)=lgt_id
    end do

    t(4) = MPI_wtime()
    ! =======================================================
    ! sort list
    if (lgt_n > 1) then
        call quicksort(lgt_sortednumlist, 1, lgt_n, 2)
    end if
    t(5) = MPI_wtime()
    call toc("create_active_and_sorted_lists_tree (reset lgt_n)", t(1)-t0)
    call toc("create_active_and_sorted_lists_tree (loop over hvy_n)", t(2)-t(1))
    call toc("create_active_and_sorted_lists_tree (comm)", t(3)-t(2))
    call toc("create_active_and_sorted_lists_tree (loop over lgt_n_sum)", t(4)-t(3))
    call toc("create_active_and_sorted_lists_tree (quicksort)", t(5)-t(4))
    call toc("create_active_and_sorted_lists_tree", MPI_wtime()-t0)
end subroutine create_active_and_sorted_lists_tree_comm


! ################################################################################
!> Updates active lgt/hvy lists from lgt_block data FOR ONE TREE by looping
!> over the hvy ids and synchronizing the active lists with the other procs
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
    integer(kind=ik)                    :: k, N, hvy_id, block_rank, ierr,lgt_n_sum, lgt_id
    ! process rank
    integer(kind=ik)                    :: rank, TREE_ID_IDX, N_blocks_tot, mpisize
    real(kind=rk) :: t0,t(5)

    !integer(kind=ik), allocatable, save     :: my_lgt_recv_buffer(:,:)
    integer(kind=tsize), allocatable, save     :: my_lgt_send_buffer(:,:)
    integer(kind=ik), allocatable, save     :: proc_lgt_n(:), proc_lgt_start(:)

    t0 = MPI_wtime()

    mpisize = params%number_procs
    rank = params%rank
    N    = params%number_blocks
    TREE_ID_IDX = params%max_treelevel + IDX_TREE_ID

    if (.not.allocated(proc_lgt_n)) allocate( proc_lgt_n(1:mpisize) )
    if (.not.allocated(proc_lgt_start)) allocate( proc_lgt_start(1:mpisize) )
    !if (.not.allocated(my_lgt_recv_buffer)) allocate( my_lgt_recv_buffer( size(lgt_active),1) ) ! 2 components: lgt_sortednumlist(:,1:2)
    if (.not.allocated(my_lgt_send_buffer)) then
        allocate( my_lgt_send_buffer( N,2) ) ! 2 components: lgt_sortednumlist(:,1:2)
        my_lgt_send_buffer=-1
    endif
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
    my_lgt_send_buffer(1:hvy_n,:) = -1

    ! reset active block numbers
    lgt_n = 0
    hvy_n = 0
    t(1) = MPI_wtime()

    ! =======================================================
    ! loop over all hvy_data
    ! =======================================================
    do hvy_id = 1, N
        call hvy2lgt( k, hvy_id, rank, N )
        ! block is active
        if (lgt_block(k,TREE_ID_IDX)==tree_ID) then
            if ( lgt_block(k, 1) /= -1 ) then
                ! ---------------------------
                ! update hvy active
                ! ---------------------------
                ! save heavy id, only if proc responsable for block
                hvy_active( hvy_n + 1 ) = hvy_id
                hvy_n                   = hvy_n + 1
                ! sorted list
                ! first index stores the light id of the block
                my_lgt_send_buffer(hvy_n, 1) = k
                ! second index stores the numerical treecode
                ! + the tree index
                my_lgt_send_buffer(hvy_n, 2) = treecode2int(lgt_block(k, 1:params%max_treelevel), tree_id )
            end if
        end if
    end do
    t(2) = MPI_wtime()
    ! ==========================================================================
    ! MPI Synchronization
    ! ==========================================================================
    ! Note: lgt_n is hear clearly the same as hvy_n since we have looped
    ! over lgt ids of this proc
    !n sum up lgt_n over all procs
    call MPI_allgather(hvy_n, 1, MPI_INTEGER4, proc_lgt_n, 1, MPI_INTEGER4, WABBIT_COMM, ierr )
    ! this is the global buffer size
    lgt_n = sum(proc_lgt_n)
    do k = 1, mpisize
        proc_lgt_start(k) = sum(proc_lgt_n(1:k-1))! + 1
    enddo

    call MPI_allgatherv( rank*N + hvy_active(1:hvy_n), hvy_n, MPI_INTEGER4, &
    lgt_active( 1:lgt_n), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
    WABBIT_COMM, ierr)
    call MPI_allgatherv( my_lgt_send_buffer(1:hvy_n, 1), hvy_n, MPI_INTEGER8, &
    lgt_sortednumlist( 1:lgt_n, 1 ), proc_lgt_n, proc_lgt_start, MPI_INTEGER8, &
    WABBIT_COMM, ierr)
    call MPI_allgatherv( my_lgt_send_buffer(1:hvy_n, 2), hvy_n, MPI_INTEGER8, &
    lgt_sortednumlist( 1:lgt_n, 2 ), proc_lgt_n, proc_lgt_start, MPI_INTEGER8, &
    WABBIT_COMM, ierr)

    t(3) = MPI_wtime()

    ! no loop anymore at this point
    t(4) = MPI_wtime()
    ! =======================================================
    ! sort list
    if (lgt_n > 1) then
        call quicksort(lgt_sortednumlist, 1, lgt_n, 2)
    end if
    t(5) = MPI_wtime()
    call toc("create_active_and_sorted_lists_tree (reset lgt_n)", t(1)-t0)
    call toc("create_active_and_sorted_lists_tree (loop over hvy_n)", t(2)-t(1))
    call toc("create_active_and_sorted_lists_tree (comm)", t(3)-t(2))
    call toc("create_active_and_sorted_lists_tree (loop over lgt_n)", t(4)-t(3))
    call toc("create_active_and_sorted_lists_tree (quicksort)", t(5)-t(4))
    call toc("create_active_and_sorted_lists_tree", MPI_wtime()-t0)
end subroutine create_active_and_sorted_lists_tree

! ################################################################################
!> Updates active lgt/hvy lists from lgt_block data FOR ONE TREE
!> \author PKrah
subroutine create_active_and_sorted_lists_tree_old( params, lgt_block, lgt_active, &
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
    real(kind=rk) :: t0,t(5)

    t0 = MPI_wtime()


    ! reset old lists, use old numbers of active blocks. If the old numbers are
    ! invalid (too small too large), then we delete everything in the lists
    !> \todo Check if resetting the arrays is not a waste of time in any case!
    if (lgt_n>size(lgt_active)) lgt_n = size(lgt_active)
    if (hvy_n>size(hvy_active)) hvy_n = size(hvy_active)

    if (lgt_n<=0) lgt_n = size(lgt_active)
    if (hvy_n<=0) hvy_n = size(hvy_active)
    t(1) = MPI_wtime()

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

    t(2) = MPI_WTIME()
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
            call lgt2proc( block_rank, k, N )
            if ( rank == block_rank ) then
                ! convert light data id into heavy data id
                call lgt2hvy( heavy_id, k, rank, N)
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
    t(3) = MPI_wtime()
    ! =======================================================
    ! sort list
    if (lgt_n > 1) then
        call quicksort(lgt_sortednumlist, 1, lgt_n, 2)
    end if
    t(4) = MPI_wtime()
    call toc("create_active_and_sorted_lists_tree (set lgt_n)", t(1)-t0)
    call toc("create_active_and_sorted_lists_tree (reset lgt,hvy list)", t(2)-t(1))
    call toc("create_active_and_sorted_lists_tree (loop over lgt_n)", t(3)-t(2))
    call toc("create_active_and_sorted_lists_tree (quicksort)", t(4)-t(3))
    call toc("create_active_and_sorted_lists_tree", MPI_wtime()-t0)
end subroutine create_active_and_sorted_lists_tree_old


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
subroutine create_active_and_sorted_lists_forest_comm( params, lgt_block, lgt_active, &
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
    integer(kind=ik)                    :: k, N, hvy_id, block_rank, fsize
    ! process rank
    integer(kind=ik)                    :: rank, tree_id, TREE_ID_IDX, lgt_n_sum, hvy_n_sum
    integer(kind=tsize) :: treecode_int
    real(kind=rk) :: t0

    integer(kind=ik) :: ierr, mpisize
    integer(kind=tsize), allocatable, save     :: my_lgt_send_buffer( :, :, :)
    integer(kind=ik), allocatable, save     :: proc_lgt_n(:), proc_lgt_start(:)


    t0 = MPI_wtime()

    mpisize = params%number_procs
    rank = params%rank
    N    = params%number_blocks
    fsize = params%forest_size
    TREE_ID_IDX = params%max_treelevel + IDX_TREE_ID

    if (.not.allocated(proc_lgt_n)) allocate( proc_lgt_n(1:mpisize) )
    if (.not.allocated(proc_lgt_start)) allocate( proc_lgt_start(1:mpisize) )
    !if (.not.allocated(my_lgt_recv_buffer)) allocate( my_lgt_recv_buffer( size(lgt_active),1) ) ! 2 components: lgt_sortednumlist(:,1:2)
    if (.not.allocated(my_lgt_send_buffer)) then
        allocate( my_lgt_send_buffer( N, 2, fsize) ) ! 2 components: lgt_sortednumlist(:,1:2)
        my_lgt_send_buffer = -1
    endif


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
      lgt_sortednumlist(1:lgt_n(tree_id), :, tree_id) = -1
      my_lgt_send_buffer(1:hvy_n(tree_id), :, tree_id) = -1
    end do


    ! reset active block numbers
    lgt_n = 0
    hvy_n = 0
    tree_n = 0

    ! =======================================================
    ! loop over all light data
    ! =======================================================
    do hvy_id = 1, N
        call hvy2lgt( k, hvy_id, rank, N )
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
            hvy_n(tree_id) = hvy_n(tree_id) + 1
            hvy_active( hvy_n(tree_id) , tree_id ) = hvy_id

            ! sorted list
            treecode_int = treecode2int(lgt_block(k, 1:params%max_treelevel), tree_id )
            ! first index stores the light id of the block
            my_lgt_send_buffer(lgt_n(tree_id), 1, tree_id) = k
            ! second index stores the numerical treecode
            my_lgt_send_buffer(lgt_n(tree_id), 2, tree_id) = treecode_int

        end if
    end do

    ! ==========================================================================
    ! MPI Synchronization
    ! ==========================================================================
    ! Note: lgt_n is hear clearly the same as hvy_n since we have looped
    ! over lgt ids of this proc
    call MPI_ALLREDUCE(MPI_IN_PLACE,tree_n,1,MPI_INTEGER4, MPI_MAX, WABBIT_COMM, ierr)

    do tree_id = 1, tree_n
    !n sum up lgt_n over all procs
        call MPI_allgather(hvy_n(tree_id), 1, MPI_INTEGER4, proc_lgt_n, 1, MPI_INTEGER4, WABBIT_COMM, ierr )
        ! this is the global buffer size
        lgt_n(tree_id) = sum(proc_lgt_n)
        do k = 1, mpisize
            proc_lgt_start(k) = sum(proc_lgt_n(1:k-1))! + 1
        enddo

        call MPI_allgatherv( rank*N + hvy_active(1:hvy_n(tree_id), tree_id), hvy_n(tree_id), MPI_INTEGER4, &
        lgt_active( 1:lgt_n(tree_id), tree_id), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
        WABBIT_COMM, ierr)

        call MPI_allgatherv( my_lgt_send_buffer(1:hvy_n(tree_id), 1, tree_id), hvy_n(tree_id), MPI_INTEGER8, &
        lgt_sortednumlist( 1:lgt_n(tree_id), 1, tree_id), proc_lgt_n, proc_lgt_start, MPI_INTEGER8, &
        WABBIT_COMM, ierr)
        call MPI_allgatherv( my_lgt_send_buffer(1:hvy_n(tree_id), 2, tree_id), hvy_n(tree_id), MPI_INTEGER8, &
        lgt_sortednumlist( 1:lgt_n(tree_id), 2, tree_id ), proc_lgt_n, proc_lgt_start, MPI_INTEGER8, &
        WABBIT_COMM, ierr)

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
end subroutine create_active_and_sorted_lists_forest_comm



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
            call lgt2proc( block_rank, k, N )
            if ( rank == block_rank ) then
                ! convert light data id into heavy data id
                call lgt2hvy( heavy_id, k, rank, N)
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
