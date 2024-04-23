! !> \brief create all active (lgt/hvy) lists, create also sorted list of active
! !! light blocks with numerical treecodes
! !! input:    - light data
! !! output:   - list of active blocks
! !!           - number of active blocks
! !!           - sorted light data with numerical treecodes
! ! ********************************************************************************************
! !> Updates active lgt/hvy lists from lgt_block data FOR ONE TREE by looping
! !> over the hvy ids and synchronizing the active lists with the other procs
! ! This routine updates active lists for A SINGLE TREE inside a forest
! subroutine createActiveSortedLists_tree_comm( params, tree_ID)

!     ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
!     use module_params
    
!     implicit none
!     !-----------------------------------------------------------------
!     type (type_params), intent(in)      :: params
!     integer(kind=ik), intent(in)        :: tree_ID
!     !-----------------------------------------------------------------

!     integer(kind=ik)                    :: k, N, hvy_id, block_rank, ierr, lgt_id
!     integer(kind=ik)                    :: rank, TREE_ID_IDX, N_blocks_tot, mpisize
!     real(kind=rk) :: t0,t(5)

!     integer(kind=ik), allocatable, save     :: proc_lgt_n(:), proc_lgt_start(:)

!     ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
!     ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
!     ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
!     ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
!     ! to the last subroutine.)  -Thomas

!     t0 = MPI_wtime()

!     mpisize = params%number_procs
!     rank = params%rank
!     N    = params%number_blocks
!     TREE_ID_IDX = IDX_TREE_ID

!     if (.not.allocated(proc_lgt_n)) allocate( proc_lgt_n(1:mpisize) )
!     if (.not.allocated(proc_lgt_start)) allocate( proc_lgt_start(1:mpisize) )

!     ! reset old lists, use old numbers of active blocks. If the old numbers are
!     ! invalid (too small too large), then we delete everything in the lists
!     !> \todo Check if resetting the arrays is not a waste of time in any case!
!     if (lgt_n(tree_ID) > size(lgt_active, 1)) lgt_n(tree_ID) = size(lgt_active, 1)
!     if (hvy_n(tree_ID) > size(hvy_active, 1)) hvy_n(tree_ID) = size(hvy_active, 1)

!     if (lgt_n(tree_ID)<=0) lgt_n(tree_ID) = size(lgt_active, 1)
!     if (hvy_n(tree_ID)<=0) hvy_n(tree_ID) = size(hvy_active, 1)

!     ! reset the active lists
!     lgt_active(1:lgt_n(tree_ID), tree_ID) = -1
!     hvy_active(1:hvy_n(tree_ID), tree_ID) = -1
!     lgt_sortednumlist(1:lgt_n(tree_ID),:, tree_ID) = -1

!     ! reset active block numbers
!     lgt_n(tree_ID) = 0
!     hvy_n(tree_ID) = 0
!     t(1) = MPI_wtime()

!     ! =======================================================
!     ! loop over all hvy_data
!     ! =======================================================
!     do hvy_id = 1, N
!         call hvy2lgt( k, hvy_id, rank, N )
!         ! block is part of the tree in question
!         if (lgt_block(k, TREE_ID_IDX) == tree_ID) then
!             ! assuming that there are more blocks per rank then blocks with the corresponding tree ID
!             if ( lgt_block(k, 1) /= -1) then ! active ?
!                 ! ---------------------------
!                 ! update hvy active
!                 ! ---------------------------
!                 ! save heavy id, only if proc responsable for block
!                 hvy_active( hvy_n(tree_ID) + 1, tree_ID ) = hvy_id
!                 hvy_n(tree_ID) = hvy_n(tree_ID) + 1
!             end if
!         end if
!     end do
!     t(2) = MPI_wtime()
!     ! ==========================================================================
!     ! MPI Synchronization
!     ! ==========================================================================
!     ! Note: lgt_n is hear clearly the same as hvy_n since we have looped
!     ! over lgt ids of this proc
!     !n sum up lgt_n over all procs
!     call MPI_allgather(hvy_n(tree_ID), 1, MPI_INTEGER4, proc_lgt_n, 1, MPI_INTEGER4, WABBIT_COMM, ierr )
!     ! this is the global buffer size
!     lgt_n(tree_ID) = sum(proc_lgt_n)
!     do k = 1, mpisize
!         proc_lgt_start(k) = sum(proc_lgt_n(1:k-1))! + 1
!     enddo

!     call MPI_allgatherv( rank*N + hvy_active(1:hvy_n(tree_ID), tree_ID), hvy_n(tree_ID), MPI_INTEGER4, &
!     lgt_active(1:lgt_n(tree_ID), tree_ID), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
!     WABBIT_COMM, ierr)

!     ! call MPI_allgatherv( lgt_sortednumlist(1:lgt_n, 2), lgt_n, MPI_INTEGER4, &
!     ! my_lgt_recv_buffer( 1:lgt_n,2 ), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
!     ! WABBIT_COMM, ierr)
!     t(3) = MPI_wtime()
!     ! UNPACK SYNCHRONIZED DATA
!     do k = 1, lgt_n(tree_ID)
!         lgt_id = lgt_active(k, tree_ID)
!         lgt_sortednumlist(k, 1, tree_ID) = lgt_id
!         lgt_sortednumlist(k, 2, tree_ID) = treecode2int(lgt_block(lgt_id, 1:params%Jmax), tree_ID )
!         ! lgt_sortednumlist(k, 2, tree_ID) = treearray2bid(lgt_block(lgt_id, 1:params%Jmax), &
!         !     dim=params%dim, tree_ID=tree_ID, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)
!     end do

!     t(4) = MPI_wtime()
!     ! =======================================================
!     ! sort list
!     if (lgt_n(tree_ID) > 1) then
!         call quicksort(lgt_sortednumlist(:,:,tree_ID), 1, lgt_n(tree_ID), 2)
!     end if
!     t(5) = MPI_wtime()
!     call toc("createActiveSortedLists_tree (reset lgt_n)", t(1)-t0)
!     call toc("createActiveSortedLists_tree (loop over hvy_n)", t(2)-t(1))
!     call toc("createActiveSortedLists_tree (comm)", t(3)-t(2))
!     call toc("createActiveSortedLists_tree (loop over lgt_n_sum)", t(4)-t(3))
!     call toc("createActiveSortedLists_tree (quicksort)", t(5)-t(4))
!     call toc("createActiveSortedLists_tree", MPI_wtime()-t0)
! end subroutine createActiveSortedLists_tree_comm


! ################################################################################
!> Updates active lgt/hvy lists from lgt_block data FOR ONE TREE by looping
!> over the hvy ids and synchronizing the active lists with the other procs
!> \author PKrah
! This routine updates active lists for A SINGLE TREE inside a forest
subroutine createActiveSortedLists_tree( params, tree_ID)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: tree_ID
    !-----------------------------------------------------------------

    integer(kind=ik)                    :: k, N, hvy_id, block_rank, ierr,lgt_n_sum, lgt_id
    integer(kind=ik)                    :: rank, N_blocks_tot, mpisize
    integer(kind=tsize)                 :: treecode_ID
    real(kind=rk) :: t0,t(5)

    !integer(kind=ik), allocatable, save     :: my_lgt_recv_buffer(:,:)
    integer(kind=ik), allocatable, save     :: my_lgt_send_buffer(:,:)
    integer(kind=ik), allocatable, save     :: proc_lgt_n(:), proc_lgt_start(:)

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    t0 = MPI_wtime()

    mpisize = params%number_procs
    rank = params%rank
    N    = params%number_blocks

    if (.not.allocated(proc_lgt_n)) allocate( proc_lgt_n(1:mpisize) )
    if (.not.allocated(proc_lgt_start)) allocate( proc_lgt_start(1:mpisize) )
    
    ! 4 components: lgt_sortednumlist(1:4, :)
    ! ID, TC1, TC2, Tree_ID * 100 + Level
    if (.not.allocated(my_lgt_send_buffer)) then
        allocate( my_lgt_send_buffer( 4,N) )
        my_lgt_send_buffer=-1
    endif

    ! reset old lists, use old numbers of active blocks. If the old numbers are
    ! invalid (too small too large), then we delete everything in the lists
    !> \todo Check if resetting the arrays is not a waste of time in any case!
    if (lgt_n(tree_ID) > size(lgt_active, 1)) lgt_n(tree_ID) = size(lgt_active, 1)
    if (hvy_n(tree_ID) > size(hvy_active, 1)) hvy_n(tree_ID) = size(hvy_active, 1)

    if (lgt_n(tree_ID)<=0) lgt_n(tree_ID) = size(lgt_active, 1)
    if (hvy_n(tree_ID)<=0) hvy_n(tree_ID) = size(hvy_active, 1)

    ! reset the active lists
    lgt_active(1:lgt_n(tree_ID), tree_ID) = -1
    hvy_active(1:hvy_n(tree_ID), tree_ID) = -1
    lgt_sortednumlist(:, 1:lgt_n(tree_ID), tree_ID) = -1
    my_lgt_send_buffer(:, 1:hvy_n(tree_ID)) = -1

    ! reset active block numbers
    lgt_n(tree_ID) = 0
    hvy_n(tree_ID) = 0
    t(1) = MPI_wtime()

    ! =======================================================
    ! loop over all hvy_data
    ! =======================================================
    do hvy_id = 1, N
        call hvy2lgt( k, hvy_id, rank, N )

        if (lgt_block(k, IDX_TREE_ID) == tree_ID) then
            ! check if block is active: TC > 0
            ! performance: don't construct tc and check only first int
            if (lgt_block(k, IDX_TC_1 ) >= 0) then
                ! ---------------------------
                ! update hvy active
                ! ---------------------------
                ! save heavy id, only if proc responsable for block
                hvy_active( hvy_n(tree_ID) + 1, tree_ID ) = hvy_id
                hvy_n(tree_ID) = hvy_n(tree_ID) + 1
                ! sorted list
                ! first index stores the light id of the block
                my_lgt_send_buffer(1, hvy_n(tree_ID)) = k

                ! second and third index stores the numerical treecode
                my_lgt_send_buffer(2:3, hvy_n(tree_ID)) = lgt_block(k, IDX_TC_1 : IDX_TC_2)

                ! tree_ID and level are combined into one number for performance purposes
                ! as max_level is 31 for 2D, we shift tree_ID by 100
                ! this might be confusing but can be entangled easily
                ! fourth index stores the level and tree_id
                my_lgt_send_buffer(4, hvy_n(tree_ID)) = lgt_block(k, IDX_MESH_LVL) + tree_ID * 100

                ! ! fourth index stores the level
                ! my_lgt_send_buffer(hvy_n(tree_ID), 4) = lgt_block(k, IDX_MESH_LVL)
                ! ! fifth index stores the tree_id
                ! my_lgt_send_buffer(hvy_n(tree_ID), 5) = tree_ID
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
    call MPI_allgather(hvy_n(tree_ID), 1, MPI_INTEGER4, proc_lgt_n, 1, MPI_INTEGER4, WABBIT_COMM, ierr )

    ! this is the global buffer size
    lgt_n(tree_ID) = sum(proc_lgt_n)
    do k = 1, mpisize
        proc_lgt_start(k) = sum(proc_lgt_n(1:k-1))! + 1
    enddo


    call MPI_allgatherv( rank*N + hvy_active(1:hvy_n(tree_ID), tree_ID), hvy_n(tree_ID), MPI_INTEGER4, &
    lgt_active( 1:lgt_n(tree_ID), tree_ID), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
    WABBIT_COMM, ierr)

    call MPI_allgatherv( my_lgt_send_buffer(1:4, 1:hvy_n(tree_ID)), 4*hvy_n(tree_ID), MPI_INTEGER4, &
    lgt_sortednumlist(1:4, 1:lgt_n(tree_ID), tree_ID), 4*proc_lgt_n, 4*proc_lgt_start, MPI_INTEGER4, &
    WABBIT_COMM, ierr)

    t(3) = MPI_wtime()

    ! no loop anymore at this point
    t(4) = MPI_wtime()
    ! =======================================================
    ! sort list
    if (lgt_n(tree_ID) > 1) then
        call quicksort(lgt_sortednumlist(:,:,tree_ID), 1, lgt_n(tree_ID), 4)
    end if
    t(5) = MPI_wtime()
    call toc("createActiveSortedLists_tree (reset lgt_n)", t(1)-t0)
    call toc("createActiveSortedLists_tree (loop over hvy_n)", t(2)-t(1))
    call toc("createActiveSortedLists_tree (comm)", t(3)-t(2))
    call toc("createActiveSortedLists_tree (loop over lgt_n)", t(4)-t(3))
    call toc("createActiveSortedLists_tree (quicksort)", t(5)-t(4))
    call toc("createActiveSortedLists_tree", MPI_wtime()-t0)
end subroutine createActiveSortedLists_tree

! ! ################################################################################
! !> Updates active lgt/hvy lists from lgt_block data FOR ONE TREE
! !> \author PKrah
! ! This routine updates active lists for A SINGLE TREE inside a forest
! subroutine createActiveSortedLists_tree_old( params, tree_ID)

!     implicit none
!     !-----------------------------------------------------------------
!     type (type_params), intent(in)      :: params
!     integer(kind=ik), intent(in)        :: tree_ID
!     !-----------------------------------------------------------------

!     integer(kind=ik)                    :: k, N, heavy_id, block_rank
!     integer(kind=ik)                    :: rank, TREE_ID_IDX
!     real(kind=rk) :: t0,t(5)

!     ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
!     ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
!     ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
!     ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
!     ! to the last subroutine.)  -Thomas

!     t0 = MPI_wtime()


!     ! reset old lists, use old numbers of active blocks. If the old numbers are
!     ! invalid (too small too large), then we delete everything in the lists
!     !> \todo Check if resetting the arrays is not a waste of time in any case!
!     if (lgt_n(tree_ID) > size(lgt_active, 1)) lgt_n(tree_ID) = size(lgt_active, 1)
!     if (hvy_n(tree_ID) > size(hvy_active, 1)) hvy_n(tree_ID) = size(hvy_active, 1)

!     if (lgt_n(tree_ID)<=0) lgt_n(tree_ID) = size(lgt_active, 1)
!     if (hvy_n(tree_ID)<=0) hvy_n(tree_ID) = size(hvy_active, 1)


!     ! reset the active lists
!     lgt_active(1:lgt_n(tree_ID), tree_ID) = -1
!     hvy_active(1:hvy_n(tree_ID), tree_ID) = -1
!     lgt_sortednumlist(1:lgt_n(tree_ID),:, tree_ID) = -1

!     ! reset active block numbers
!     lgt_n(tree_ID) = 0
!     hvy_n(tree_ID) = 0

!     rank = params%rank
!     N    = params%number_blocks
!     TREE_ID_IDX = IDX_TREE_ID

!     ! =======================================================
!     ! loop over all light data
!     ! =======================================================
!     t(2) = MPI_WTIME()
!     do k = 1, size(lgt_block, 1)
!         ! block is active
!         if ( lgt_block(k, 1) /= -1 .and. lgt_block(k,TREE_ID_IDX)==tree_ID) then
!             ! ---------------------------
!             ! update light active
!             ! ---------------------------
!             ! save lgt id as active block
!             lgt_n(tree_ID) = lgt_n(tree_ID) + 1
!             lgt_active(lgt_n(tree_ID), tree_ID) = k

!             ! ---------------------------
!             ! update hvy active
!             ! ---------------------------
!             ! save heavy id, only if proc responsable for block
!             call lgt2proc( block_rank, k, N )
!             if ( rank == block_rank ) then
!                 ! convert light data id into heavy data id
!                 call lgt2hvy( heavy_id, k, rank, N)

!                 hvy_active( hvy_n(tree_ID) + 1, tree_ID) = heavy_id
!                 hvy_n(tree_ID) = hvy_n(tree_ID) + 1
!             end if

!             ! sorted list
!             ! first index stores the light id of the block
!             lgt_sortednumlist(lgt_n(tree_ID), 1, tree_ID) = k
!             ! second index stores the numerical treecode
!             ! + the tree index
!             lgt_sortednumlist(lgt_n(tree_ID), 2, tree_ID) = treecode2int(lgt_block(k, 1:params%Jmax), tree_ID )
!             ! lgt_sortednumlist(lgt_n(tree_ID), 2, tree_ID) = treearray2bid(lgt_block(k, 1:params%Jmax), &
!             !     dim=params%dim, tree_ID=tree_ID, level=lgt_block(k, IDX_MESH_LVL), max_level=params%Jmax)
!         end if
!     end do
!     t(3) = MPI_wtime()
!     ! =======================================================
!     ! sort list
!     if (lgt_n(tree_ID) > 1) then
!         call quicksort(lgt_sortednumlist(:,:,tree_ID), 1, lgt_n(tree_ID), 2)
!     end if
!     t(4) = MPI_wtime()
!     call toc("createActiveSortedLists_tree (set lgt_n)", t(1)-t0)
!     call toc("createActiveSortedLists_tree (reset lgt,hvy list)", t(2)-t(1))
!     call toc("createActiveSortedLists_tree (loop over lgt_n)", t(3)-t(2))
!     call toc("createActiveSortedLists_tree (quicksort)", t(4)-t(3))
!     call toc("createActiveSortedLists_tree", MPI_wtime()-t0)
! end subroutine createActiveSortedLists_tree_old


! ! ################################################################################
! !> Updates active lgt/hvy lists from lgt_block data.
! !> Returns active lists for each tree in the forest and lgt_n for each tree.
! !
! !> \details
! !> -------------------------------------------------------------
! !>     code                    | explanation
! !> -------------------------------------------------------------
! !> lgt_active(:,tree_ID)       | active block list of tree
! !> lgt_n(:,tree_ID)            | number of active blocks in tree
! !> -------------------------------------------------------------
! !> \author PKrah
! ! This routine updates active lists for ALL TREES inside a forest
! subroutine createActiveSortedLists_forest_comm( params)

!     implicit none
!     !-----------------------------------------------------------------
!     type (type_params), intent(in) :: params
!     !-----------------------------------------------------------------

!     integer(kind=ik) :: k, N, hvy_id, block_rank, fsize
!     integer(kind=ik) :: rank, tree_ID, TREE_ID_IDX, lgt_n_sum, hvy_n_sum
!     integer(kind=tsize) :: treecode_int
!     real(kind=rk) :: t0

!     integer(kind=ik) :: ierr, mpisize
!     integer(kind=tsize), allocatable, save     :: my_lgt_send_buffer( :, :, :)
!     integer(kind=ik), allocatable, save     :: proc_lgt_n(:), proc_lgt_start(:)

!     ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
!     ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
!     ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
!     ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
!     ! to the last subroutine.)  -Thomas


!     t0 = MPI_wtime()

!     mpisize = params%number_procs
!     rank = params%rank
!     N    = params%number_blocks
!     fsize = params%forest_size
!     TREE_ID_IDX = IDX_TREE_ID

!     if (.not.allocated(proc_lgt_n)) allocate( proc_lgt_n(1:mpisize) )
!     if (.not.allocated(proc_lgt_start)) allocate( proc_lgt_start(1:mpisize) )
!     !if (.not.allocated(my_lgt_recv_buffer)) allocate( my_lgt_recv_buffer( size(lgt_active),1) ) ! 2 components: lgt_sortednumlist(:,1:2)
!     if (.not.allocated(my_lgt_send_buffer)) then
!         allocate( my_lgt_send_buffer( N, 2, fsize) ) ! 2 components: lgt_sortednumlist(:,1:2)
!         my_lgt_send_buffer = -1
!     endif


!     ! =======================================================
!     ! Reset active lists of all trees
!     ! =======================================================
!     ! note: this seems to be a complicated way of reseting the
!     !       active lists, but it is very crucial for performance!
!     !       NEVER RESET the full array without reasons!!!
!     do tree_ID = 1, tree_n
!         ! check if lgt_n or hvy_n of tree is valid (not too small or too large)
!         if (lgt_n(tree_ID)>size(lgt_active(:,tree_ID)) .or. lgt_n(tree_ID) <=0) &
!         lgt_n(tree_ID) = size(lgt_active(:,tree_ID))
!         if (hvy_n(tree_ID)>size(hvy_active(:,tree_ID)) .or. hvy_n(tree_ID) <=0) &
!         hvy_n(tree_ID) = size(hvy_active(:,tree_ID))

!         ! reset the active lists
!         lgt_active(1:lgt_n(tree_ID), tree_ID) = -1
!         hvy_active(1:hvy_n(tree_ID), tree_ID) = -1
!         lgt_sortednumlist(1:lgt_n(tree_ID), :, tree_ID) = -1
!         my_lgt_send_buffer(1:hvy_n(tree_ID), :, tree_ID) = -1
!     end do


!     ! reset active block numbers
!     lgt_n = 0
!     hvy_n = 0
!     tree_n = 0

!     ! =======================================================
!     ! loop over all light data
!     ! =======================================================
!     do hvy_id = 1, N
!         call hvy2lgt( k, hvy_id, rank, N )
!         ! block is active
!         if ( lgt_block(k, 1) /= -1 ) then
!             ! which tree id has the current block k?
!             tree_ID = lgt_block(k, TREE_ID_IDX)
!             ! find the highest tree number. this is should be the same as
!             ! the number of active trees: tree_n
!             tree_n = max(tree_ID, tree_n)

!             ! ---------------------------
!             ! update light active
!             ! ---------------------------
!             ! save lgt id as active block
!             lgt_n(tree_ID) = lgt_n(tree_ID) + 1
!             lgt_active( lgt_n(tree_ID), tree_ID) = k

!             ! ---------------------------
!             ! update hvy active
!             ! ---------------------------
!             hvy_n(tree_ID) = hvy_n(tree_ID) + 1
!             hvy_active( hvy_n(tree_ID) , tree_ID ) = hvy_id

!             ! sorted list
!             treecode_int = treecode2int(lgt_block(k, 1:params%Jmax), tree_ID )
!             ! treecode_int = treearray2bid(lgt_block(k, 1:params%Jmax), &
!             !     dim=params%dim, tree_ID=tree_ID, level=lgt_block(k, IDX_MESH_LVL), max_level=params%Jmax)
            
!             ! first index stores the light id of the block
!             my_lgt_send_buffer(lgt_n(tree_ID), 1, tree_ID) = k
!             ! second index stores the numerical treecode
!             my_lgt_send_buffer(lgt_n(tree_ID), 2, tree_ID) = treecode_int

!         end if
!     end do

!     ! ==========================================================================
!     ! MPI Synchronization
!     ! ==========================================================================
!     ! Note: lgt_n is hear clearly the same as hvy_n since we have looped
!     ! over lgt ids of this proc
!     call MPI_ALLREDUCE(MPI_IN_PLACE,tree_n,1,MPI_INTEGER4, MPI_MAX, WABBIT_COMM, ierr)

!     do tree_ID = 1, tree_n
!         !n sum up lgt_n over all procs
!         call MPI_allgather(hvy_n(tree_ID), 1, MPI_INTEGER4, proc_lgt_n, 1, MPI_INTEGER4, WABBIT_COMM, ierr )
!         ! this is the global buffer size
!         lgt_n(tree_ID) = sum(proc_lgt_n)
!         do k = 1, mpisize
!             proc_lgt_start(k) = sum(proc_lgt_n(1:k-1))! + 1
!         enddo

!         call MPI_allgatherv( rank*N + hvy_active(1:hvy_n(tree_ID), tree_ID), hvy_n(tree_ID), MPI_INTEGER4, &
!         lgt_active( 1:lgt_n(tree_ID), tree_ID), proc_lgt_n, proc_lgt_start, MPI_INTEGER4, &
!         WABBIT_COMM, ierr)

!         call MPI_allgatherv( my_lgt_send_buffer(1:hvy_n(tree_ID), 1, tree_ID), hvy_n(tree_ID), MPI_INTEGER8, &
!         lgt_sortednumlist( 1:lgt_n(tree_ID), 1, tree_ID), proc_lgt_n, proc_lgt_start, MPI_INTEGER8, &
!         WABBIT_COMM, ierr)
!         call MPI_allgatherv( my_lgt_send_buffer(1:hvy_n(tree_ID), 2, tree_ID), hvy_n(tree_ID), MPI_INTEGER8, &
!         lgt_sortednumlist( 1:lgt_n(tree_ID), 2, tree_ID ), proc_lgt_n, proc_lgt_start, MPI_INTEGER8, &
!         WABBIT_COMM, ierr)

!     end do


!     ! =======================================================
!     ! sort list of every single tree
!     do tree_ID = 1, tree_n
!         if (lgt_n(tree_ID) > 1) then
!             call quicksort(lgt_sortednumlist(:,:,tree_ID), 1, lgt_n(tree_ID), 2)
!         end if
!     end do

!     ! check if the number of trees is not bigger then the size of the forest
!     ! The forest size is defined as the maximum number of trees in the forest.
!     if (tree_n > fsize) call abort(1402192, "Too many trees in the forest!!" )

!     call toc("createActiveSortedLists_forest", MPI_wtime()-t0)
! end subroutine createActiveSortedLists_forest_comm



! ################################################################################
!> Updates active lgt/hvy lists from lgt_block data.
!> Returns active lists for each tree in the forest and lgt_n for each tree.
!
!> \details
!> -------------------------------------------------------------
!>     code                    | explanation
!> -------------------------------------------------------------
!> lgt_active(:,tree_ID)       | active block list of tree
!> lgt_n(:,tree_ID)            | number of active blocks in tree
!> -------------------------------------------------------------
!> \author PKrah
! This routine updates active lists for ALL TREES in a forest
subroutine createActiveSortedLists_forest(params)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in) :: params
    !-----------------------------------------------------------------

    integer(kind=ik) :: k, N, heavy_id, block_rank, fsize
    integer(kind=ik) :: rank, tree_ID, lgt_n_sum, hvy_n_sum
    integer(kind=tsize) :: treecode_int
    integer(kind=ik) :: tc_ik(2)
    real(kind=rk) :: t0

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    t0 = MPI_wtime()


    rank  = params%rank
    N     = params%number_blocks
    fsize = params%forest_size

    ! =======================================================
    ! Reset active lists of all trees
    ! =======================================================
    ! note: this seems to be a complicated way of reseting the
    !       active lists, but it is very important for performance!
    !       NEVER RESET the full array without reasons!!!
    do tree_ID = 1, tree_n
        ! check if lgt_n or hvy_n of tree is valid (not too small or too large)
        if (lgt_n(tree_ID)>size(lgt_active(:,tree_ID)) .or. lgt_n(tree_ID) <=0) &
        lgt_n(tree_ID) = size(lgt_active(:,tree_ID))
        if (hvy_n(tree_ID)>size(hvy_active(:,tree_ID)) .or. hvy_n(tree_ID) <=0) &
        hvy_n(tree_ID) = size(hvy_active(:,tree_ID))

        ! reset the active lists
        lgt_active(1:lgt_n(tree_ID), tree_ID) = -1
        hvy_active(1:hvy_n(tree_ID), tree_ID) = -1
        lgt_sortednumlist(:, 1:lgt_n(tree_ID), tree_ID) = -1
    end do


    ! reset active block numbers
    lgt_n  = 0
    hvy_n  = 0
    tree_n = 0

    ! =======================================================
    ! loop over all light data
    ! =======================================================
    do k = 1, size(lgt_block, 1)
        ! check if block is active: TC > 0
        ! performance: don't construct tc and check only first int
        if (lgt_block(k, IDX_TC_1 ) >= 0) then

            ! which tree id has the current block k?
            tree_ID = lgt_block(k, IDX_TREE_ID)
            ! find the highest tree number. this is should be the same as
            ! the number of active trees: tree_n
            tree_n = max(tree_ID, tree_n)

            ! ---------------------------
            ! update light active
            ! ---------------------------
            ! save lgt id as active block
            lgt_n(tree_ID) = lgt_n(tree_ID) + 1
            lgt_active( lgt_n(tree_ID), tree_ID) = k

            ! ---------------------------
            ! update hvy active
            ! ---------------------------
            ! save heavy id, only if proc responsable for block
            call lgt2proc( block_rank, k, N )
            if ( rank == block_rank ) then
                ! convert light data id into heavy data id
                call lgt2hvy( heavy_id, k, rank, N)
                hvy_n(tree_ID) = hvy_n(tree_ID) + 1
                hvy_active( hvy_n(tree_ID) , tree_ID ) = heavy_id
            end if

            ! sorted list

            ! first index stores the light id of the block
            lgt_sortednumlist(1, lgt_n(tree_ID), tree_ID) = k
            ! second index stores the numerical treecode
            lgt_sortednumlist(2:3, lgt_n(tree_ID), tree_ID) = lgt_block(k, IDX_TC_1 : IDX_TC_2)

            ! tree_ID and level are combined into one number for performance purposes
            ! as max_level is 31 for 2D, we shift tree_ID by 100
            ! this might be confusing but can be entangled easily
            ! fourth index stores the level and tree_id
            lgt_sortednumlist(4, lgt_n(tree_ID), tree_ID) = lgt_block(k, IDX_MESH_LVL) + tree_ID * 100

            ! ! fourth index stores the level
            ! lgt_sortednumlist(lgt_n(tree_ID), 4, tree_ID) = lgt_block(k, IDX_MESH_LVL)
            ! ! fifth index stores the tree_ID
            ! lgt_sortednumlist(lgt_n(tree_ID), 5, tree_ID) = tree_ID 
        end if
    end do

    ! =======================================================
    ! sort list of every single tree
    do tree_ID = 1, tree_n
        if (lgt_n(tree_ID) > 1) then
            call quicksort(lgt_sortednumlist(:,:,tree_ID), 1, lgt_n(tree_ID), 4)
        end if
    end do

    ! check if the number of trees is not bigger then the size of the forest
    ! The forest size is defined as the maximum number of trees in the forest.
    if (tree_n > fsize) call abort(1402192, "Too many trees in the forest!!" )

    call toc("createActiveSortedLists_forest", MPI_wtime()-t0)
end subroutine createActiveSortedLists_forest
