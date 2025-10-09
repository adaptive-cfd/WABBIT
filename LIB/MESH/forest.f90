! \brief Delete blocks of tree where all point values are 0, usefull for mask
subroutine prune_tree( params, hvy_block, tree_ID)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)    :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data

    integer(kind=ik) :: k, lgt_id, hvy_id, rank, N, g ,Bs(3)

    rank = params%rank
    N = params%number_blocks
    g = params%g
    Bs = params%Bs


    if (rank==0) write(*,'("Tree-pruning, before Nb=",i7)') lgt_n(tree_ID)

    if (params%dim == 3) then
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, rank, N )

            ! pruning condition: all entries of the block are zero (or below)
            if (.not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id) > 0.0_rk) ) then
                ! pruning: delete the block from the tree
                lgt_block(lgt_id, :) = -1_ik
            endif
        end do
    else
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, rank, N )

            ! pruning condition: all entries of the block are zero (or below)
            if ( .not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1, hvy_id) > 0.0_rk) ) then
                ! pruning: delete the block from the tree
                lgt_block(lgt_id, :) = -1_ik
            endif
        end do
    endif

    call synchronize_lgt_data( params, refinement_status_only=.false. )

    call createActiveSortedLists_forest(params)
    ! do not call updateNeighbors_tree on the pruned tree..

    if (rank==0) write(*,'("Tree-pruning, pruned to Nb=",i7)') lgt_n(tree_ID)
end subroutine

!---------------------------------------------------------------------------------

! subroutine add_pruned_to_full_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
!     hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_ID_pruned, tree_ID_full)
!     implicit none
!     !-----------------------------------------------------------------
!     type (type_params), intent(in) :: params   !< params structure
!     integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
!     integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
!     integer(kind=ik), intent(in)      :: tree_ID_pruned, tree_ID_full
!     integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
!     integer(kind=ik), intent(inout)   :: lgt_block(:, :)  !< light data array
!     real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
!     integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
!     integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
!     integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
!
!     integer(kind=ik) :: k, lgt_id, Jmax, hvy_id, rank, N, fsize, i
!     integer(kind=ik) :: lgt_id1, lgt_id2, hvy_id1, hvy_id2
!     integer(kind=ik) :: level1, level2, rank_pruned, rank_full, n_comm
!     logical :: exists
!     integer(kind=ik), allocatable, save :: comm_list(:,:)
!     integer(kind=tsize) :: treecode1, treecode2
!
!     fsize = params%forest_size
!     Jmax = params%Jmax ! max treelevel
!     rank = params%rank
!     N = params%number_blocks
!
!     ! NOTES:
!     !   We assume that the pruned tree is on JMAX (rhs level) and we assume that
!     !   the full one is not finer, only coarser.
!
!     if (.not.allocated(comm_list)) allocate( comm_list( 3, params%number_procs*N) )
!
!     call createActiveSortedLists_forest( params, lgt_block, lgt_active, lgt_n, &
!     hvy_active, hvy_n, lgt_sortednumlist, tree_n)
!
!     ! a pruned tree has fewer entries: loop over it instead of the other one?
!     ! if you find a block in the full tree -> well then that's good, copy.
!     ! else: the target grid is either refined or coarsened at this position.
!
!     ! Step 1: XFER. we look for blocks that exist in both pruned and full tree and
!     ! if they are on different mpiranks, we xfer the pruned trees block to the rank
!     ! holding the corresponding full trees block.
!
!     ! Step 1a: prepare for xfer, gather all required xfers
!     n_comm = 0
!     do k = 1, lgt_n(tree_ID_pruned)
!
!         lgt_id1 = lgt_active(k, tree_ID_pruned)
!         level1  = lgt_block(lgt_id1, IDX_MESH_LVL)
!
!         ! does the block from the pruned tree exist in the full one?
!         call doesBlockExist_tree( lgt_block(lgt_id1, 1:level1), exists, lgt_id2, &
!         lgt_sortednumlist(:,:,tree_ID_full), lgt_n(tree_ID_full), tree_ID_full)
!
!         if (exists) then
!             ! we found the pruned trees block in the full tree : we happily copy
!             ! its heavy data
!
!             ! now check on which CPU this block is currently
!             call lgt2proc( rank_pruned, lgt_id1, N)
!             call lgt2proc( rank_full, lgt_id2, N)
!
!             if (rank_full /= rank_pruned) then
!                 n_comm = n_comm + 1
!                 comm_list(1, n_comm) = rank_pruned   ! sender mpirank
!                 comm_list(2, n_comm) = rank_full   ! receiver mpirank
!                 comm_list(3, n_comm) = lgt_id1 ! block lgt_id to send
!                 write(*,*) "found on different rank", n_comm, rank_pruned, rank_full
!             endif
!         else
!             ! it does not exist, but maybe it is coarsened by one level?
!             ! NOTE: code will not detect if coarsened by more than one level
!             call doesBlockExist_tree( lgt_block(lgt_id1, 1:level1-1), exists, lgt_id2, &
!             lgt_sortednumlist(:,:,tree_ID_full), lgt_n(tree_ID_full), tree_ID_full)
!
!             if (exists) then
!
!             else
!                 ! nothing...the pruned tree block is not there, not even on a coarser level
!             endif
!
!         endif
!     enddo
!
!     ! Step 1b: actual xfer.
!     call block_xfer( params, comm_list, n_comm, lgt_block, hvy_block )
!
!     call createActiveSortedLists_forest( params, lgt_block, lgt_active, &
!     lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
!
!
!     ! Step 2: ADDITION. now we're sure that blocks existing in both trees are on the
!     ! same mpirank. therefore, the responsible rank can just add them together.
!     do k = 1, hvy_n(tree_ID_pruned)
!
!         hvy_id1  = hvy_active(k, tree_ID_pruned)
!         call hvy2lgt(lgt_id1, hvy_id1, params%rank, params%number_blocks)
!         level1  = lgt_block(lgt_id1, IDX_MESH_LVL)
!
!         call doesBlockExist_tree(lgt_block(lgt_id1, 1:level1), exists, lgt_id2, &
!         lgt_sortednumlist(:,:,tree_ID_full), lgt_n(tree_ID_full), tree_ID_full)
!
!         ! we found the pruned trees block in the full tree : we happily copy
!         ! its heavy data
!         if (exists) then
!             ! now check on which CPU this block is currently
!             call lgt2proc( rank_pruned, lgt_id1, N)
!             call lgt2proc( rank_full, lgt_id2, N)
!
!             if (rank_full /= rank_pruned) then
!                 ! This is an erorr one should never see: the block xfer is done, and coexisting blocks should
!                 ! be on the same mpirank now.
!                 call abort(030719,"pruned_to_full_tree: although we xferred, we found a coexisting block on a different rank.")
!             endif
!
!             ! only responsible mpirank can perform the addition
!             if (params%rank==rank_full) then
!                 ! get both heavy ids
!                 call lgt2hvy( hvy_id1, lgt_id1, rank_pruned, N )
!                 call lgt2hvy( hvy_id2, lgt_id2, rank_full  , N )
!                 ! actual addition
!                 hvy_block(:,:,:,1,hvy_id2) = hvy_block(:,:,:,1,hvy_id2) + hvy_block(:,:,:,1,hvy_id1)
!             endif
!         else
!             ! we did not find it. The grid has changed in the interior of the
!             ! obstacle, and we can set those interior blocks to constant 1. But
!             ! we need to figure out which blocks to set to 1.
!
!             ! we can further assume that all blocks in the pruned tree are the minimum:
!             ! only the interface is very important, interior blocks are dictated by
!             ! gradedness. Hence: assuming that the fluid/solid interface is always
!             ! on the finest level, only a refinement in the interior is possible (and
!             ! not a coarsening) => look for sister blocks on higher levels
!
!             ! find all sister blocks and add 1 to them. no xfer required.
!             treecode1 = treecode2int( lgt_block(lgt_id1, 1:level1) )
!
!             do i = 1, lgt_n(tree_ID_full)
!                 lgt_id2 = lgt_active(i, tree_ID_full)
!                 treecode2 = treecode2int( lgt_block(lgt_id2, 1:level1) )
!
!                 if (treecode1 == treecode2) then
!                     ! this is one of the sisters
!                     call lgt2proc( rank_full, lgt_id2, N)
!                     if (params%rank == rank_full) then
!                         call lgt2hvy( hvy_id2, lgt_id2, rank_full, N)
!                         hvy_block(:,:,:,1,hvy_id2) = hvy_block(:,:,:,1,hvy_id2) + 1.0_rk
!                     endif
!                 endif
!             enddo
!         endif
!     enddo
!
! end subroutine
subroutine add_pruned_to_full_tree( params, hvy_block, tree_ID_pruned, tree_ID_full)
    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)    :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID_pruned, tree_ID_full
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data

    integer(kind=ik) :: k, lgt_id, Jmax, hvy_id, rank, N, fsize, i, ierr
    integer(kind=ik) :: lgt_id1, lgt_id2, hvy_id1, hvy_id2
    integer(kind=ik) :: level1, level2, rank_pruned, rank_full, n_comm
    integer(kind=ik) :: num_blocks_count(3), i_var  ! debug counters: [send, recv, keep]
    logical :: exists
    integer(kind=ik), allocatable, save :: comm_list(:,:)
    integer(kind=tsize) :: treecode1, treecode2
    real(kind=rk) :: t_cycle, t_block
    character(len=4) :: string_kind(3)
    character(len=cshort) :: format_string

    fsize = params%forest_size
    Jmax = params%Jmax ! max treelevel
    rank = params%rank
    N = params%number_blocks
    t_cycle = MPI_Wtime()

    ! debug counters: [send, recv, keep]
    num_blocks_count(1:3) = 0

    ! init
    treecode1 = 0_tsize; treecode2 = 0_tsize

    ! NOTES:
    !   We assume that the pruned tree is on JMAX (rhs level) and we assume that
    !   the full one is finer, never coarser.

    if (.not.allocated(comm_list)) allocate( comm_list( 3, params%number_procs*N ) )

    t_block = MPI_Wtime()
    call createActiveSortedLists_forest(params)
    call toc("add_pruned_to_full_tree (createActiveSortedLists_forest)", 261, MPI_Wtime() - t_block)

    ! a pruned tree has fewer entries: loop over it instead of the fuller (fluid) one!
    !
    ! if you find a pruned tree block in the full (fluid) tree -> well then that's good, copy.
    ! else: We just assume that this block does not exist, as it maybe is on another level, but we don't care
    !       Why? Well because sometimes we have full grid formulation (mask on JMax and JMax-1)
    !       and sometimes only leaf grid (mask on JMax only). For the second point the JMax-1 blocks are ignored.
    ! The implication is: we do not use geometries where the fluid grid is COARSENED inside the obstacle. In such a
    ! case, the pruned tree block is not found on the full fluid tree, but not because the geometry is not relevant,
    ! rather because the fluid grid is coarsened. In such a case, the block should be set entirely to ones.
    ! With out STL file generation, a band around the surface (a shell) is created in the mask function - and this
    ! function here is used for STL files. In this case, coarsening inside the solid shell would require very thick 
    ! shells (at least four blocks on the finest level).

    ! Step 1: XFER. we look for blocks that exist in both pruned and full tree and
    ! if they are on different mpiranks, we xfer the pruned trees block to the rank
    ! holding the corresponding full trees block.

    ! Step 1a: prepare for xfer, gather all required xfers
    t_block = MPI_Wtime()
    n_comm = 0
    do k = 1, lgt_n(tree_ID_pruned)

        lgt_id1 = lgt_active(k, tree_ID_pruned)
        level1  = lgt_block(lgt_id1, IDX_MESH_LVL)
        treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

        ! do we find the treecode of the pruned tree-block in the full grid?
        call doesBlockExist_tree(treecode1, exists, lgt_id2, dim=params%dim, level=level1, tree_id=tree_ID_full, max_level=params%Jmax)

        if (exists) then
            ! we found the pruned trees block in the full tree : we happily copy
            ! its heavy data

            ! now check on which CPU this block is currently
            call lgt2proc( rank_pruned, lgt_id1, N)
            call lgt2proc( rank_full, lgt_id2, N)

            if (rank_full /= rank_pruned) then
                n_comm = n_comm + 1
                comm_list(1, n_comm) = rank_pruned   ! sender mpirank
                comm_list(2, n_comm) = rank_full   ! receiver mpirank
                comm_list(3, n_comm) = lgt_id1 ! block lgt_id to send
                
                ! count transferred blocks for this rank - used for development debugging output
                if (params%rank == rank_pruned) then
                    num_blocks_count(1) = num_blocks_count(1) + 1  ! send
                endif
                if (params%rank == rank_full) then
                    num_blocks_count(2) = num_blocks_count(2) + 1  ! recv
                endif
            else
                if (params%rank == rank_full) then
                    num_blocks_count(3) = num_blocks_count(3) + 1  ! keep
                endif
            endif
        endif
    enddo
    call toc("add_pruned_to_full_tree (gather comms)", 262, MPI_Wtime() - t_block)

    ! Step 1b: actual xfer.
    t_block = MPI_Wtime()
    call block_xfer( params, comm_list, n_comm, hvy_block )
    call toc("add_pruned_to_full_tree (block_xfer)", 263, MPI_Wtime() - t_block)

    ! As some blocks have been transferred, the active lists are outdated.
    t_block = MPI_Wtime()
    call createActiveSortedLists_tree(params, tree_ID_pruned)
    call toc("add_pruned_to_full_tree (createActiveSortedLists_tree)", 264, MPI_Wtime() - t_block)

    ! Step 2: ADDITION. now we're sure that blocks existing in both trees are on the
    ! same mpirank. therefore, the responsible mpirank can just add them together.
    t_block = MPI_Wtime()
    do k = 1, hvy_n(tree_ID_pruned)

        hvy_id1  = hvy_active(k, tree_ID_pruned)
        call hvy2lgt(lgt_id1, hvy_id1, params%rank, params%number_blocks)
        level1  = lgt_block(lgt_id1, IDX_MESH_LVL)
        treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

        call doesBlockExist_tree(treecode1, exists, lgt_id2, dim=params%dim, level=level1, tree_id=tree_ID_full, max_level=params%Jmax)

        ! we found the pruned trees block in the full tree : we happily copy
        ! its heavy data
        if (exists) then
            ! now check on which CPU this block is currently
            call lgt2proc( rank_pruned, lgt_id1, N)
            call lgt2proc( rank_full, lgt_id2, N)

            if (rank_full /= rank_pruned) then
                ! This is an error one should never see: the block xfer is done, and coexisting blocks should
                ! be on the same mpirank now.
                call abort(030719,"pruned_to_full_tree: although we xferred, we found a coexisting block on a different rank.")
            endif

            ! only responsible mpirank can perform the addition
            if (params%rank==rank_full) then
                ! get both heavy ids
                call lgt2hvy( hvy_id1, lgt_id1, rank_pruned, N )
                call lgt2hvy( hvy_id2, lgt_id2, rank_full  , N )

                ! actual addition
                ! hvy_block(:,:,:,:,hvy_id2) = hvy_block(:,:,:,:,hvy_id2) + hvy_block(:,:,:,:,hvy_id1)
                ! hvy_block(:,:,:,:,hvy_id2) = max(hvy_block(:,:,:,:,hvy_id2),hvy_block(:,:,:,:,hvy_id1))

                where ( hvy_block(:,:,:,1,hvy_id2)<=hvy_block(:,:,:,1,hvy_id1) )
                    hvy_block(:,:,:,1,hvy_id2) = hvy_block(:,:,:,1,hvy_id1)
                    hvy_block(:,:,:,2,hvy_id2) = hvy_block(:,:,:,2,hvy_id1)
                    hvy_block(:,:,:,3,hvy_id2) = hvy_block(:,:,:,3,hvy_id1)
                    hvy_block(:,:,:,4,hvy_id2) = hvy_block(:,:,:,4,hvy_id1)
                    hvy_block(:,:,:,5,hvy_id2) = hvy_block(:,:,:,5,hvy_id1)
                end where

            endif

        else
            ! We did not find it. However, we do not have blocks with flat 1 in our current implentation of a tube around the border
            ! So, this block might be a mother or daughter block of another border block, in leaf grids it might just not exist at the current moment
            ! But that is no problem - We just ignore it
        endif
    enddo
    call toc("add_pruned_to_full_tree (addition)", 265, MPI_Wtime() - t_block)

    if (params%debug_pruned2full) then
        ! development output, gather information on rank 0 and print to file
        string_kind = (/ "send", "recv", "keep" /)
        
        do i_var = 1, 3
            ! gather number of blocks transferred/added/ignored on rank 0
            call MPI_GATHER(num_blocks_count(i_var), 1, MPI_INTEGER, comm_list, 1, MPI_INTEGER, 0, WABBIT_COMM, ierr)
            if (params%rank == 0) then
                ! Single IO operation with dynamic format
                open(unit=99, file="debug_pruned2full.csv", status="unknown", position="append")
                if (params%number_procs == 1) then
                    write(99, '(A,i0)') string_kind(i_var)//",", comm_list(1,1)
                else
                    write(format_string, '("(A, i0,",i0,"("","",i0))")') params%number_procs - 1
                    write(99, format_string) string_kind(i_var)//",", comm_list(1,1), comm_list(1,2:params%number_procs)
                endif
                close(99)
            endif
        enddo
    endif

    call toc("add_pruned_to_full_tree (TOTAL)", 260, MPI_Wtime() - t_cycle)

end subroutine



!##############################################################
!> deletes the lgt_block data of tree with given tree_ID
!> CAUTION: active lists will be outdated!!!
subroutine delete_tree(params, tree_ID)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)   :: params    !< user defined parameter structure
    integer(kind=ik), intent(in)     :: tree_ID!< highest tree id
    !-----------------------------------------------------------------
    integer(kind=ik)                 :: k, lgt_id, hvy_id

    ! loop over active list of tree
    do k = 1, lgt_n(tree_ID)
        lgt_id = lgt_active(k, tree_ID)
        lgt_block(lgt_id, :) = -1_ik
    end do
    lgt_active(1:lgt_n(tree_ID),tree_ID)=-1_ik
    lgt_n(tree_ID)=0_ik
    hvy_active(1:hvy_n(tree_ID),tree_ID)=-1_ik
    hvy_n(tree_ID)=0_ik
end subroutine delete_tree
!##############################################################



!##############################################################
!> copy hvy and lgt data of tree_ID_source to tree_ID_dest
!> after copy the active and sorted lists are updated, as well
!> as neighbors and ghost nodes
subroutine copy_tree(params, hvy_block, tree_ID_dest, tree_ID_source, skip_sync_ghosts)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)    :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID_dest, tree_ID_source !< all data from tree_ID_source gets copied to destination_tree_ID
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    logical, intent(in), optional     :: skip_sync_ghosts  !< sometimes synching ghosts is not desired so we skip it
    !-----------------------------------------------------------------
    integer(kind=ik)    :: Jmax, lgt_id_dest, lgt_id_source, hvy_id_dest, hvy_id_source, fsize
    integer(kind=ik)    :: k, N, rank
    real(kind=rk) :: t_elapse
    logical skipSyncGhosts

    Jmax = params%Jmax ! max treelevel
    fsize= params%forest_size   ! maximal number of trees in forest
    rank = params%rank
    N = params%number_blocks

    skipSyncGhosts = .false.
    if (present(skip_sync_ghosts)) skipSyncGhosts = skip_sync_ghosts

    if (tree_ID_dest > fsize) call abort(0403191,"tree-copy: destination treeID ou of valid range")

    ! first we delete tree_ID_dest if it is already allocated.
    ! tree_ID_dest only exists if it is in the list of active trees, i.e. tree_ID_dest <= tree_n
    if (lgt_n(tree_ID_dest) > 0 ) then
        call delete_tree(params, tree_ID_dest)
    end if

    call createActiveSortedLists_forest(params)

    ! Loop over the active hvy_data
    t_elapse = MPI_WTIME()
    do k = 1, hvy_n(tree_ID_source)
        hvy_id_source = hvy_active(k, tree_ID_source)
        call hvy2lgt( lgt_id_source, hvy_id_source, rank, N )

        ! first we have to find out if the hvy data belongs to the tree we want to copy
        call get_free_local_light_id( params, rank, lgt_id_dest)
        call lgt2hvy( hvy_id_dest, lgt_id_dest, rank, N )
        !--------------------
        ! Light DATA
        !--------------------
        lgt_block(lgt_id_dest, :) = lgt_block(lgt_id_source, :)  ! copy light data
        lgt_block(lgt_id_dest, IDX_TREE_ID) = tree_ID_dest  ! asign new tree_ID
        !--------------------
        ! Heavy DATA
        !--------------------
        hvy_block( :, :, :, :, hvy_id_dest) = hvy_block( :, :, :, :, hvy_id_source)
    end do ! loop over source tree
    call toc( "copy_tree (copy heavy_data)", 250, MPI_Wtime()-t_elapse )


    ! always synchronize lgt_data when you have changed lgt_block locally
    ! (i.e. when looping over lgt_ids which originate from a hvy_id)
    t_elapse = MPI_WTIME()
    call synchronize_lgt_data( params, refinement_status_only=.false. )
    call updateMetadata_tree( params, tree_ID_dest )

    if (.not. skipSyncGhosts) call sync_ghosts_tree( params, hvy_block, tree_ID_dest)

    call toc( "copy_tree (copy synchronize hvy and lgt)", 251, MPI_Wtime()-t_elapse )

end subroutine
!##############################################################

!##############################################################
!> multiply every block of a tree with a given value alpha
subroutine multiply_tree_with_scalar(params, hvy_block, tree_ID, alpha, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in) :: params   !< params structure
    integer(kind=ik), intent(in)   :: tree_ID !< all data from tree_ID2 gets copied to tree_ID1
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(in)      :: alpha !<prefactor
    logical, intent(in),optional   :: verbosity !< if true aditional stdout is printed
    !-----------------------------------------------------------------
    integer(kind=ik)    :: hvy_id, k
    logical :: verbose = .false.

    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose ) write(*,'("scalar multiplication tree: ",i3)') tree_ID

    ! Loop over the active hvy_data, multiply block data by alpha
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        hvy_block( :, :, :, :, hvy_id) = alpha * hvy_block( :, :, :, :, hvy_id)
    end do


end subroutine
!##############################################################


!##############################################################
!> scalar product assoziated L2 norm ||u||_L2 = sqrt(<u,u>_L2)
!> note that this not implies a simple sum over the squares of u_ijk
!> it is rather a weighted sum depending on the wavelet basis.
function compute_tree_L2norm(params, hvy_block, hvy_tmp, tree_ID, verbosity ) &
    result(L2norm)

    implicit none
    !-----------------------------------------------------------------
    ! inputs:
    type (type_params), intent(inout)  :: params
    integer(kind=ik), intent(in)      :: tree_ID
    !> heavy data array - block data
    real(kind=rk),  intent(inout)   :: hvy_block(:, :, :, :, :)
    !> heavy temp data: needed in blockxfer which is called in add_two_trees
    real(kind=rk),   intent(inout)  :: hvy_tmp(:, :, :, :, :)
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    ! result
    real(kind=rk) :: L2norm
    !-----------------------------------------------------------------
    logical :: verbose=.false.
    if (present(verbosity)) verbose = verbosity

    L2norm = scalar_product_two_trees(params, hvy_block, hvy_tmp, tree_ID1=tree_ID, tree_ID2=tree_ID)
    L2norm = sqrt(L2norm)
    if (params%rank == 0 .and. verbose ) write(*,'("L2 norm: ",e13.6)') L2norm
end function
!##############################################################


!##############################################################
!> This function coarses a tree to its reference mesh.
!> The reference mesh is saved in lgt_..._ref and stores the treecode and
!> all necessary light data.
!> All blocks on the actual tree are coarsened to have the same
!> treecode structure as the reference.
!> This routine is usefull, when trying to go back to a old
!> treestructure, after calling for example refine_trees2sametreelevel.
subroutine coarse_tree_2_reference_mesh(params, lgt_block_ref, lgt_active_ref, lgt_n_ref, &
    hvy_block, hvy_tmp, tree_ID, ref_can_be_finer, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID !< number of the tree
    integer(kind=ik), intent(inout)   :: lgt_n_ref !< number of light active blocks
    integer(kind=ik), intent(inout)   :: lgt_block_ref(:, : )  !< light data array
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), intent(inout)   :: lgt_active_ref(:) !< active lists
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
    !> if true, the routine will also work if the reference mesh has some finer blocks, needs refine_tree_2_reference_mesh afterwards
    logical, intent(in),optional      :: ref_can_be_finer
    logical, intent(in),optional      :: verbosity
    !-----------------------------------------------------------------
    integer(kind=ik)    :: rank, level_ref, level, Jmax, lgt_id_ref, lgt_id, fsize
    integer(kind=ik)    :: k1, k2, Nblocks_2coarsen, level_min
    integer(kind=tsize) :: treecode_ref, treecode
    logical :: verbose, exists, refCanBeFiner

    Jmax = params%Jmax ! max treelevel
    fsize= params%forest_size   ! maximal number of trees in forest
    verbose = .false.
    refCanBeFiner = .false.
    if (present(ref_can_be_finer)) refCanBeFiner = ref_can_be_finer

    if (present(verbosity)) verbose=verbosity

    call updateMetadata_tree(params, tree_ID)

    ! Operations (addition, multiplication, etc) using two trees (=grids) require both trees (=grids)
    ! to be identical. Therefore, this routine first initializes a full grid structure with -1 on the actual tree,
    ! gives all blocks with identical TC/level to reference tree the status 0 and deletes all finer blocks, voila

    ! init full grid structure on the actual tree
    call init_full_tree(params, tree_ID, set_ref=-1, verbose_check=.false.)
    if (params%rank==0 .and. verbose) write(*,'("coarse_tree_2_reference_mesh : init_full_tree completed, size Nb=",i7)') lgt_n(tree_ID)

    ! init all ref stats to -1
    do k1 = 1, lgt_n(tree_ID)
        lgt_id = lgt_active(k1, tree_ID)
        lgt_block(lgt_id, IDX_REFINE_STS) = -1
    end do

    !----------------------------
    ! loop over all blocks of the reference tree and find the corresponding block with our sorted list to be kept
    do k1 = 1, lgt_n_ref
        lgt_id_ref = lgt_active_ref(k1)
        level_ref  = lgt_block_ref(lgt_id_ref, IDX_MESH_LVL)

        ! get the treecode of the reference block
        treecode_ref = get_tc(lgt_block_ref(lgt_id_ref, IDX_TC_1 : IDX_TC_2))

        call doesBlockExist_tree(treecode_ref, exists, lgt_id, dim=params%dim, level=level_ref, tree_id=tree_ID, max_level=params%Jmax)

        ! if the block does not exists the reference block is finer than the actual tree (or Jmin is not respected)
        ! we have two options:
        !    1) panic, scream and abort, because something went wrong (default)
        !    2) if refCanBeFiner is set to true, we coarsen the reference block until we find a block on the actual tree
        !       or we fall below Jmin, means we just keep the finest mesh possible at that position and afterwards need to refine as well
        if (.not. exists .and. .not. refCanBeFiner) then
            call abort(250930, "Something went wrong: Block on reference tree1 is too fine (or lower than Jmin)!!!")
        endif
        if (.not. exists .and. refCanBeFiner) then
            do while( .not. exists .and. level_ref >= params%Jmin )
                level_ref = level_ref - 1
                treecode_ref = tc_clear_until_level_b( treecode_ref, params%dim, level_ref, params%Jmax)
                call doesBlockExist_tree(treecode_ref, exists, lgt_id, dim=params%dim, level=level_ref, tree_id=tree_ID, max_level=params%Jmax)
            end do
            if (.not. exists) then
                call abort(250930, "Something went wrong: I have no Idea how you ended up here")
            endif
        endif

        ! set lgt_id to be kept
        lgt_block(lgt_id, IDX_REFINE_STS) = 0

    end do ! loop over reference blocks
    if (params%rank==0 .and. verbose) write(*,'("coarse_tree_2_reference_mesh : marking blocks completed")')
    !----------------------------

    !----------------------------
    ! Coarsen the tagged blocks
    !    ensure_gradedness will take care of keeping mother blocks
    !    This will also apply CE to the grid if set
    call adapt_tree(0.0_rk, params, hvy_block, tree_ID, "nothing (external)", hvy_tmp, ignore_coarsening=.false., init_full_tree_grid=.false.)
    !----------------------------
end subroutine



!##############################################################
!> This function refines a tree to its reference mesh.
!> The reference mesh is saved in lgt_..._ref and stores the treecode and
!> all necessary light data.
!> All blocks on the actual tree are refined to have the same
!> treecode structure as the reference.
subroutine refine_tree_2_reference_mesh(params, lgt_block_ref, lgt_active_ref, lgt_n_ref, &
    hvy_block, hvy_tmp, tree_ID, ref_can_be_coarser, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID !< number of the tree
    integer(kind=ik), intent(inout)   :: lgt_n_ref !< number of light active blocks
    integer(kind=ik), intent(inout)   :: lgt_block_ref(:, : )  !< light data array
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), intent(inout)   :: lgt_active_ref(:) !< active lists
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
    !> if true, the routine will also work if the reference mesh has some coarser blocks, needs coarse_tree_2_reference_mesh afterwards
    logical, intent(in),optional      :: ref_can_be_coarser
    logical, intent(in),optional      :: verbosity
    !-----------------------------------------------------------------
    integer(kind=ik)    :: level_ref, level, Jmax, lgt_id_ref, lgt_id, fsize
    integer(kind=ik)    :: k1, k2, Nblocks_2coarsen, level_min, mpierr, loop_i
    integer(kind=tsize) :: treecode_ref, treecode
    logical :: verbose, exists, refCanBeCoarser, continue_loop

    Jmax = params%Jmax ! max treelevel
    fsize= params%forest_size   ! maximal number of trees in forest
    verbose = .false.
    refCanBeCoarser = .false.
    if (present(ref_can_be_coarser)) refCanBeCoarser = ref_can_be_coarser

    if (present(verbosity)) verbose=verbosity

    call updateMetadata_tree(params, tree_ID)

    ! Operations (addition, multiplication, etc) using two trees (=grids) require both trees (=grids)
    ! to be identical. Therefore, this routine flags all blocks for refinement which are not on the reference tree, but a higher one is

    ! init full grid structure on the actual tree, so that we also check if reference mesh is coarser
    if (refCanBeCoarser) call init_full_tree(params, tree_ID, set_ref=0, verbose_check=.false.)

    continue_loop = .true.
    loop_i = 1
    do while (continue_loop)
        continue_loop = .false.
        if (params%rank==0 .and. verbose) write(*,'("refine_tree_2_reference_mesh : loop ",i3)') loop_i

        ! init all ref stats to 0
        do k1 = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k1, tree_ID)
            lgt_block(lgt_id, IDX_REFINE_STS) = 0
        end do

        !----------------------------
        ! loop over all blocks of the reference tree and find the corresponding block with our sorted list to be kept
        do k1 = 1, lgt_n_ref
            lgt_id_ref = lgt_active_ref(k1)
            level_ref  = lgt_block_ref(lgt_id_ref, IDX_MESH_LVL)

            ! get the treecode of the reference block
            treecode_ref = get_tc(lgt_block_ref(lgt_id_ref, IDX_TC_1 : IDX_TC_2))

            call doesBlockExist_tree(treecode_ref, exists, lgt_id, dim=params%dim, level=level_ref, tree_id=tree_ID, max_level=params%Jmax)

            ! if the block does not exists the reference block is finer than the actual tree (or Jmin is not respected)
            ! we now go downwards in levels until we find the mother block and set it to be refined.
            ! If this does not exist, then the reference mesh is coarser than the actual mesh or Jmin is not respected
            if (.not. exists) then
                do while( .not. exists .and. level_ref >= params%Jmin )
                    level_ref = level_ref - 1
                    treecode_ref = tc_clear_until_level_b( treecode_ref, params%dim, level_ref, params%Jmax)
                    call doesBlockExist_tree(treecode_ref, exists, lgt_id, dim=params%dim, level=level_ref, tree_id=tree_ID, max_level=params%Jmax)
                end do
                if (.not. exists) then
                    call abort(250930, "Something went wrong: Block on reference tree2 is too coarse (or lower than Jmin)!!!")
                endif
                ! set lgt_id to be refined
                lgt_block(lgt_id, IDX_REFINE_STS) = 1
            endif

        end do ! loop over reference blocks
        !----------------------------

        !----------------------------
        ! refine the tagged blocks
        call refine_tree( params, hvy_block, indicator='nothing (external)', tree_id=tree_ID, error_OOM=exists, check_full_tree=.false.)
        !----------------------------

        ! gather continue_loop to see if any mpirank wants to continue
        call MPI_Allreduce(MPI_IN_PLACE, continue_loop, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, mpierr)
        
    end do

    if (params%rank==0 .and. verbose) write(*,'("refine_tree_2_reference_mesh : refinement completed, size Nb=",i7)') lgt_n(tree_ID)

    ! remove full tree
    if (refCanBeCoarser) call prune_fulltree2leafs(params, tree_ID)

end subroutine

!> \brief Stores the lgt data of two trees in an extra array as backup
subroutine store_ref_meshes(lgt_block_ref, lgt_active_ref, lgt_n_ref, tree_ID1, tree_ID2)

    implicit none
    !-----------------------------------------------------------------
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2
    integer(kind=ik), allocatable, intent(inout) ::  lgt_block_ref(:, : ), lgt_active_ref(:,:)
    integer(kind=ik), intent(inout)   ::lgt_n_ref(2)
    !-----------------------------------------------------------------
    integer(kind=ik) :: k1
    if (maxval(lgt_n_ref) < max(lgt_n(tree_ID1),lgt_n(tree_ID2))) then
        lgt_n_ref(1) = lgt_n(tree_ID1)
        lgt_n_ref(2) = lgt_n(tree_ID2)
        if (allocated(lgt_active_ref)) deallocate(lgt_block_ref,lgt_active_ref)
        allocate(lgt_block_ref(2*maxval(lgt_n_ref), size(lgt_block,2)))
        allocate(lgt_active_ref(2*maxval(lgt_n_ref), 2))
    else
        lgt_n_ref(1) = lgt_n(tree_ID1)
        lgt_n_ref(2) = lgt_n(tree_ID2)
    endif

    ! store lgt data of tree 1
    do k1 = 1, lgt_n(tree_ID1)
        lgt_block_ref(k1,:) = lgt_block(lgt_active(k1,tree_ID1), :)
        lgt_active_ref(k1,1)    = k1
    end do
    ! store lgt data of tree 2
    do k1 = 1, lgt_n(tree_ID2)
        lgt_block_ref(k1+lgt_n_ref(1),:) = lgt_block(lgt_active(k1,tree_ID2), :)
        lgt_active_ref(k1,2)    = k1 + lgt_n_ref(1)
    end do

end subroutine

!##############################################################
!> This function takes two trees and refines their treestructures
!> until they are refined to the same max level.
!> The resulting trees are identical in treestructure. Hence
!> they have the same treecodes and number of blocks.
!> Remark: You may want to balance the load after using this routine.
subroutine refine_trees2same_lvl(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)    :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity
    !-----------------------------------------------------------------
    integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
    integer(kind=ik)    :: k1, k2, Nblocks_2refine, level_min
    integer(kind=tsize) :: treecode1, treecode2
    logical :: verbose = .false.

    Jmax = params%Jmax ! max treelevel
    fsize= params%forest_size   ! maximal number of trees in forest

    if (present(verbosity)) verbose=verbosity
    if ( params%rank == 0 .and. verbose ) write(*,'("Refining trees to same level: ",i9,",",i9)') tree_ID1, tree_ID2

    ! The Trees can only be added when their grids are identical. At present we
    ! try to keep the finest levels of all trees. This means we refine
    ! all blocks which are not on the same level.

    ! Operations (addition, multiplication, etc) using two trees (=grids) require both trees (=grids)
    ! to be identical. Therefore, this routine modifies both trees such that for any position x in space
    ! the resolution is the higher


    ! loop until both trees have the same grid structure, meaning that no block has to be refined anymore
    do while( .true. )
        Nblocks_2refine = 0

        ! Loop over the active light data of both trees and compare the meshlevels
        do k1 = 1, lgt_n(tree_ID1)
            !--------
            ! TREE 1
            !--------
            lgt_id1 = lgt_active(k1, tree_ID1)
            level1  = lgt_block(lgt_id1, IDX_MESH_LVL)

            do k2 = 1, lgt_n(tree_ID2)
                !--------
                ! TREE 2
                !--------
                lgt_id2  = lgt_active(k2, tree_ID2)
                ! Skip this block if it is already tagged for refinement.
                if ( lgt_block(lgt_id2, IDX_REFINE_STS) == +1 ) then
                    cycle
                endif
                level2 = lgt_block(lgt_id2, IDX_MESH_LVL)

                ! The treecodes can only be compared if they have the same size
                ! (i.e. treecode length). Therefore we have to find the minimum of both
                ! levels.
                level_min = min(level1, level2)

                treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                treecode1 = tc_clear_until_level_b(treecode1, params%dim, level_min, Jmax)
                treecode2 = tc_clear_until_level_b(treecode2, params%dim, level_min, Jmax)

                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                ! Comparison of the trees:
                ! - Compare the treecodes up to the smallest level both blocks share (level_min)
                ! - check if they are already on the same level. in this case, they are identical
                !   and nothing needs to be done (the block exists in both trees already)
                ! - If the treecodes agree and they are not on both the same level, mark the coarser
                !   block for refinement. As there may well be more than one level difference between
                !   both blocks, this process potentially has to be repeated afterwards.
                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                if ((treecode1 == treecode2) .and. (level1 == level2)) then
                    ! treecode comparison okay and both blocks on same level: they're the same block
                    exit  ! exit the inner loop (TREE2)
                endif

                if ((treecode1 == treecode2) .and. (level1 /= level2)) then
                    Nblocks_2refine = Nblocks_2refine + 1
                    ! now tag the coarser block to be refined
                    if (level1 > level2) then
                        ! mark the block on tree 2 for refinement
                        lgt_block(lgt_id2, IDX_REFINE_STS) = +1
                    else
                        ! mark the block on tree 1 for refinement
                        lgt_block(lgt_id1, IDX_REFINE_STS) = +1
                    end if

                    exit  ! exit the inner loop (TREE2)

                end if
                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            end do ! loop over tree2
        end do ! loop over tree1

        ! Decide if trees have the same treestructure or not (i.e. exit loop or refine)
        if (Nblocks_2refine == 0) then
            exit   ! EXIT the (while true) loop when nothing has to be refined anymore
        else
            !----------------------------
            ! refine the tagged blocks
            !----------------------------
            if (params%rank == 0 .and. verbose) write(*,'("Number of blocks marked for refinement: ",i9)') Nblocks_2refine

            ! 1) check gradedness of the grid (meshlevel of adjacent blocks should not differe by more than 1
            call ensureGradedness_tree( params, tree_ID1 )

            call ensureGradedness_tree( params, tree_ID2 )


            ! 2) refine blocks
            call refinement_execute_tree( params, hvy_block, tree_ID1 )
            call refinement_execute_tree( params, hvy_block, tree_ID2 )

            ! since lgt_block was synced we have to create the active lists again
            call createActiveSortedLists_tree( params, tree_ID1 )
            call createActiveSortedLists_tree( params, tree_ID2 )

            ! update neighbor relations and synchronice ghosts of 1st tree
            call updateNeighbors_tree( params, tree_ID1 , search_overlapping=.false.)

            call sync_ghosts_tree( params, hvy_block, tree_ID1)

            ! update neighbor relations and synchronice ghosts of 2nd tree
            call updateNeighbors_tree( params, tree_ID2, search_overlapping=.false.)

            call sync_ghosts_tree( params, hvy_block, tree_ID2)

        endif
    end do
end subroutine

!##############################################################
!> \brief Refines all active blocks lower than lvl 12??
subroutine refine_tree2(params, hvy_block, hvy_tmp, tree_ID)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)    :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
    !-----------------------------------------------------------------
    integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
    integer(kind=ik)    :: k
    integer(kind=tsize) :: treecode1, treecode2
    logical :: verbose = .true.

    Jmax = params%Jmax ! max treelevel
    fsize = params%forest_size   ! maximal number of trees in forest

    do k = 1, lgt_n(tree_ID)
        if (lgt_block( lgt_active(k,tree_ID), IDX_MESH_LVL)<12) then
            lgt_block( lgt_active(k,tree_ID), IDX_REFINE_STS) = +1
        endif
    enddo

    call respectJmaxJmin_tree( params,  tree_ID )

    ! 2) refine blocks
    call refinement_execute_tree( params, hvy_block, tree_ID )
    ! since lgt_block was synced we have to create the active lists again
    ! update neighbor relations
    call updateMetadata_tree(params, tree_ID)

    ! synchronize ghosts
    call sync_ghosts_tree( params, hvy_block, tree_ID)
end subroutine
!##############################################################



!#############################################################
!> Perform pointwise operations (+, -, /, *) for two given trees.
!> The trees will be refined to same treestructure and then the
!> operation will be executed on each block pointwise.
!> Its like the Matlabs .* or ./ for trees.
!> NOTE: if dest_tree_ID is present then the result will be saved
!>       to the given dest_tree_ID. Otherwise the operation will
!>       be saved on the first tree_ID1
!                                                      @@@
!                       @                              @@@                               @
!         @             @    @ @                                          @             @    @ @
!           @    @@      @ @                     @ @ @ @ @ @ @ @             @    @@      @ @
!           @  @          @        @@                                       @  @          @        @@
!      @     @@        @  @      @@@                   @@@              @     @@        @  @      @@@
!      @ @     @@@ @    @ @     @    @                 @@@             @ @     @@@ @    @ @     @    @
!       @@        @@     @     @@@@@@@                                  @@        @@     @     @@@@@@@
!         @@@      @@  @@    @@                         @                 @@@      @@  @@    @@
!           @@@    @@@@@ @@@@@@@         @@             @                   @@@    @@@@@ @@@@@@@         @@
!    @@@@@@          @@@@       @@@@@@@@                @            @@@@@@          @@@@       @@@@@@@@
!        @           @@@         @               @ @ @ @ @ @ @ @         @           @@@         @
!      @             @@@                                @              @             @@@
!                   @@@                                 @                           @@@
!                   @@@                                 @                           @@@
!                   @@@                                                             @@@
!                   @@@                             @       @                       @@@
!                   @@@                              @     @                        @@@
!                   @@@                               @   @                         @@@
!                   @@@@                                @                           @@@@
!                    @@@@                        @ @ @ @ @ @ @ @                    @@@@
!                  @@@@@@@@    @                        @                         @@@@@@@@    @
!             @@@     @    @@@                        @   @                   @@@     @    @@@
!          @@@  @      @       @                     @     @               @@@  @      @       @
!                                                   @       @
! (done with png to ascii)
subroutine tree_pointwise_arithmetic(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2, &
    operation, dest_tree_ID,a,b)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    character (len=*), intent(in)  :: operation !< which arithmetical operation (+,-,*,/) which is applied
    integer(kind=ik),optional, intent(in)::dest_tree_ID !< optional for saving results to destination tree id
    real(kind=rk), intent(in),optional      :: a,b ! scalar factors for addition: result = a * tree_ID1 + b *tree_ID2
    !-----------------------------------------------------------------
    integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
    integer(kind=ik)    :: k1, k2, N, iq, iz, iy, ix, &
    hvy_id1, hvy_id2, Bs(3), g, lgt_id_dest, hvy_id_dest
    integer(kind=tsize) :: treecode1, treecode2
    real (kind=rk) :: t_elapse, alpha(2)
    character (len=5) :: op !< which arithmetical operation (+,-,*,/) which is applied
    integer(kind=ik) , save, allocatable :: lgt_active_ref(:,:), lgt_block_ref(:,:)
    integer(kind=ik) , save :: lgt_n_ref(2)=0_ik
    character(len=cshort) :: MR
    Jmax = params%Jmax ! max treelevel
    fsize= params%forest_size   ! maximal number of trees in forest
    N    = params%number_blocks ! number of blocks per rank
    rank = params%rank       ! proc rank
    g = params%g         ! number of ghost nodes
    Bs= params%Bs               ! number of grid points per block

    ! decide if inplace or out of place operation
    if(present(dest_tree_ID))then
        ! operation will be out of place: result saved in dest_id
        op = trim(operation)//"out"
        ! delete tree_ID_dest if it is already allocated.
        ! tree_ID_dest only exists if it is in the list of active trees, i.e. tree_ID_dest <= tree_n
        if (dest_tree_ID <= tree_n) then
            ! Caution: active lists will be outdated
            call delete_tree(params, dest_tree_ID)

        end if
    else
        op = trim(operation)
    endif

    alpha = (/1 , 1/)
    if(present(a)) alpha(1) = a
    if(present(b)) alpha(2) = b
    !=============================================
    ! Prepare the trees for pointwise arithmentic
    !=============================================
    ! The Trees can only be added when their grids are identical. At present we
    ! try to keep the finest levels of all trees. This means we refine
    ! all blocks which are not on the same level.

    if (tree_ID1 .ne. tree_ID2 .and. params%Jmax .ne. params%Jmin) then
        t_elapse = MPI_WTIME()
        call store_ref_meshes(lgt_block_ref, lgt_active_ref, lgt_n_ref, tree_ID1, tree_ID2)

        call refine_trees2same_lvl(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2)
        call toc( "pointwise_tree_arithmetic (refine_trees2same_lvl)", 252, MPI_Wtime()-t_elapse )
    end if

    if (tree_ID1 .ne. tree_ID2) then
        ! because all trees have the same treestructrue their hilbertcurve is identical
        ! and therefore balance the load will try to distribute blocks with the same
        ! treecode (but on different trees) at the same rank.
        t_elapse = MPI_WTIME()
        call balanceLoad_tree( params, hvy_block, tree_ID1 )
        call balanceLoad_tree( params, hvy_block, tree_ID2 )
    end if

    call toc( "pointwise_tree_arithmetic (balancing after refine_trees2same_lvl)", 253, MPI_Wtime()-t_elapse )

    !=================================================
    ! Decide which pointwice arithmetic shell be used
    !=================================================
    !      _
    !     / \    from here we assume that trees have the same grid structure
    !    /   \   and the blocks of the same treecode are on the same rank !
    !   /     \
    !  /caution\
    ! +--------+
    t_elapse = MPI_WTIME()
    select case(op)
    case("+out")
        !         #
        !         #
        !    ###########           ADDITION out of place
        !         #
        !         #
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1,tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2,tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )
                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    ! copy light data from one of the added trees
                    ! first we have to find out if the hvy data belongs to the tree we want to copy
                    call get_free_local_light_id( params, rank, lgt_id_dest)
                    call lgt2hvy( hvy_id_dest, lgt_id_dest, rank, N )
                    !--------------------
                    ! Light DATA
                    !--------------------
                    lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                    lgt_block(lgt_id_dest, IDX_TREE_ID) = dest_tree_ID  ! asign new tree_ID

                    if (params%dim==2) then
                        !##########################################
                        ! 2d
                        do iq = 1, params%N_eqn
                            do iy = g+1, Bs(2) + g
                                do ix = g+1, Bs(1) + g
                                    hvy_block(ix,iy,1,iq,hvy_id_dest) =alpha(1) * hvy_block(ix,iy,1,iq,hvy_id1) + &
                                    alpha(2) * hvy_block(ix,iy,1,iq,hvy_id2)
                                end do
                            end do
                        end do
                        !##########################################
                    else
                        !##########################################
                        ! 3D
                        do iq = 1, params%N_eqn
                            do iz = g+1, Bs(3) + g
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,iz,iq,hvy_id_dest) = alpha(1) * hvy_block(ix,iy,iz,iq,hvy_id1) + &
                                        alpha(2) * hvy_block(ix,iy,iz,iq,hvy_id2)
                                    end do
                                end do
                            end do
                        end do
                        !##########################################
                    endif
                    exit
                end if
            end do
        end do
    case("+")
        !         #
        !         #
        !    ###########           ADDITION
        !         #
        !         #
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1,tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2,tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )
                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    if (params%dim==2) then
                        !##########################################
                        ! 2d
                        do iq = 1, params%N_eqn
                            do iy = g+1, Bs(2) + g
                                do ix = g+1, Bs(1) + g
                                    hvy_block(ix,iy,1,iq,hvy_id1) = alpha(1) * hvy_block(ix,iy,1,iq,hvy_id1) + &
                                    alpha(2) * hvy_block(ix,iy,1,iq,hvy_id2)
                                end do
                            end do
                        end do
                        !##########################################
                    else
                        !##########################################
                        ! 3D
                        do iq = 1, params%N_eqn
                            do iz = g+1, Bs(3) + g
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,iz,iq,hvy_id1) = alpha(1) * hvy_block(ix,iy,iz,iq,hvy_id1) + &
                                        alpha(2) * hvy_block(ix,iy,iz,iq,hvy_id2)
                                    end do
                                end do
                            end do
                        end do
                        !##########################################
                    endif
                    exit
                end if
            end do
        end do
    case("-")

        !
        !
        !    ###########           SUBSTRACTION
        !
        !
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1, tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            ! we want to add everything to tree1
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2,tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )
                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    if (params%dim==2) then
                        !##########################################
                        ! 2d
                        do iq = 1, params%N_eqn
                            do iy = g+1, Bs(2) + g
                                do ix = g+1, Bs(1) + g
                                    hvy_block(ix,iy,1,iq,hvy_id1) = hvy_block(ix,iy,1,iq,hvy_id1) - &
                                    hvy_block(ix,iy,1,iq,hvy_id2)
                                end do
                            end do
                        end do
                        !##########################################
                    else
                        !##########################################
                        ! 3D
                        do iq = 1, params%N_eqn
                            do iz = g+1, Bs(3) + g
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,iz,iq,hvy_id1) = hvy_block(ix,iy,iz,iq,hvy_id1) - &
                                        hvy_block(ix,iy,iz,iq,hvy_id2)
                                    end do
                                end do
                            end do
                        end do
                        !##########################################
                    endif
                    exit
                end if
            end do
        end do
    case("/")
        !         #
        !
        !    ###########          Division
        !
        !         #
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1, tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            ! we want to add everything to tree1
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2, tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )
                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    !> \todo make this one faster by looping only over inner grid points
                    hvy_block(:,:,:,:,hvy_id1) = hvy_block(:,:,:,:,hvy_id1) / &
                    hvy_block(:,:,:,:,hvy_id2)
                    exit
                end if
            end do
        end do
    case("/out")
        !         #
        !
        !    ###########           Division out of place
        !
        !         #
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1,tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2,tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )
                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    ! copy light data from one of the added trees
                    ! first we have to find out if the hvy data belongs to the tree we want to copy
                    call get_free_local_light_id( params, rank, lgt_id_dest)
                    call lgt2hvy( hvy_id_dest, lgt_id_dest, rank, N )
                    !--------------------
                    ! Light DATA
                    !--------------------
                    lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                    lgt_block(lgt_id_dest, IDX_TREE_ID) = dest_tree_ID  ! asign new tree_ID

                    if (params%dim==2) then
                        !##########################################
                        ! 2d
                        do iq = 1, params%N_eqn
                            do iy = g+1, Bs(2) + g
                                do ix = g+1, Bs(1) + g
                                    hvy_block(ix,iy,1,iq,hvy_id_dest) = hvy_block(ix,iy,1,iq,hvy_id1) / &
                                    hvy_block(ix,iy,1,iq,hvy_id2)
                                end do
                            end do
                        end do
                        !##########################################
                    else
                        !##########################################
                        ! 3D
                        do iq = 1, params%N_eqn
                            do iz = g+1, Bs(3) + g
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,iz,iq,hvy_id_dest) = hvy_block(ix,iy,iz,iq,hvy_id1) / &
                                        hvy_block(ix,iy,iz,iq,hvy_id2)
                                    end do
                                end do
                            end do
                        end do
                        !##########################################
                    endif
                    exit
                end if
            end do
        end do
    case("-out")
        !
        !
        !    ###########           Substraction out of place
        !
        !
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1,tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2,tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )
                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    ! copy light data from one of the added trees
                    ! first we have to find out if the hvy data belongs to the tree we want to copy
                    call get_free_local_light_id( params, rank, lgt_id_dest)
                    call lgt2hvy( hvy_id_dest, lgt_id_dest, rank, N )
                    !--------------------
                    ! Light DATA
                    !--------------------
                    lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                    lgt_block(lgt_id_dest, IDX_TREE_ID) = dest_tree_ID  ! asign new tree_ID

                    if (params%dim==2) then
                        !##########################################
                        ! 2d
                        do iq = 1, params%N_eqn
                            do iy = g+1, Bs(2) + g
                                do ix = g+1, Bs(1) + g
                                    hvy_block(ix,iy,1,iq,hvy_id_dest) = hvy_block(ix,iy,1,iq,hvy_id1) - &
                                    hvy_block(ix,iy,1,iq,hvy_id2)
                                end do
                            end do
                        end do
                        !##########################################
                    else
                        !##########################################
                        ! 3D
                        do iq = 1, params%N_eqn
                            do iz = g+1, Bs(3) + g
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,iz,iq,hvy_id_dest) = hvy_block(ix,iy,iz,iq,hvy_id1) - &
                                        hvy_block(ix,iy,iz,iq,hvy_id2)
                                    end do
                                end do
                            end do
                        end do
                        !##########################################
                    endif
                    exit
                end if
            end do
        end do
    case("*")
        !        # #
        !         #
        !     #########        Multiplication
        !         #
        !        # #
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1,tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            ! we want to add everything to tree1
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2,tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )
                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    if (params%dim==2) then
                        !##########################################
                        ! 2d
                        do iq = 1, params%N_eqn
                            do iy = g+1, Bs(2) + g
                                do ix = g+1, Bs(1) + g
                                    hvy_block(ix,iy,1,iq,hvy_id1) = hvy_block(ix,iy,1,iq,hvy_id1) * &
                                    hvy_block(ix,iy,1,iq,hvy_id2)
                                end do
                            end do
                        end do
                        !##########################################
                    else
                        !##########################################
                        ! 3D
                        do iq = 1, params%N_eqn
                            do iz = g+1, Bs(3) + g
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,iz,iq,hvy_id1) = hvy_block(ix,iy,iz,iq,hvy_id1) * &
                                        hvy_block(ix,iy,iz,iq,hvy_id2)
                                    end do
                                end do
                            end do
                        end do
                        !##########################################
                    endif
                    exit
                end if
            end do
        end do
    case("*out")
        !        # #
        !         #
        !     #########        Multiplication out of place
        !         #
        !        # #
        do k1 = 1, hvy_n(tree_ID1)
            hvy_id1 = hvy_active(k1,tree_ID1)
            call hvy2lgt(lgt_id1, hvy_id1, rank, N )
            level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))

            do k2 = 1, hvy_n(tree_ID2)
                hvy_id2 = hvy_active(k2,tree_ID2)
                call hvy2lgt(lgt_id2, hvy_id2, rank, N )

                level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
                treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))

                if (treecode1 .ne. treecode2 ) then
                    cycle
                else
                    ! copy light data from one of the added trees
                    ! first we have to find out if the hvy data belongs to the tree we want to copy
                    call get_free_local_light_id( params, rank, lgt_id_dest)
                    call lgt2hvy( hvy_id_dest, lgt_id_dest, rank, N )
                    !--------------------
                    ! Light DATA
                    !--------------------
                    lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                    lgt_block(lgt_id_dest, IDX_TREE_ID) = dest_tree_ID  ! asign new tree_ID

                    if (params%dim==2) then
                        !##########################################
                        ! 2d
                        do iq = 1, params%N_eqn
                            do iy = g+1, Bs(2) + g
                                do ix = g+1, Bs(1) + g
                                    hvy_block(ix,iy,1,iq,hvy_id_dest) = hvy_block(ix,iy,1,iq,hvy_id1) * &
                                    hvy_block(ix,iy,1,iq,hvy_id2)
                                end do
                            end do
                        end do
                        !##########################################
                    else
                        !##########################################
                        ! 3D
                        do iq = 1, params%N_eqn
                            do iz = g+1, Bs(3) + g
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,iz,iq,hvy_id_dest) = hvy_block(ix,iy,iz,iq,hvy_id1) * &
                                        hvy_block(ix,iy,iz,iq,hvy_id2)
                                    end do
                                end do
                            end do
                        end do
                        !##########################################
                    endif
                    exit
                endif
            end do
        end do
    case default
        call abort(135,"Operation:  "//op//"  unknown")
    end select
    call toc( "pointwise_tree_arithmetic (hvy_data operation)", 254, MPI_Wtime()-t_elapse )

    t_elapse = MPI_Wtime()

    !!! Attention we should not use prefilter of CDF44 wavelets here. Because
    !!! it will change the original data!!!
    MR = params%wavelet
    params%wavelet = "CDF40"

    if (present(dest_tree_ID)) then
        ! we have to synchronize lgt data since we were updating it locally on this procesor
        call synchronize_lgt_data( params, refinement_status_only=.false. )

        call createActiveSortedLists_forest(params)
        !
        if (tree_ID1 .ne. tree_ID2 .and. params%Jmax .ne. params%Jmin) then
            call coarse_tree_2_reference_mesh(params, lgt_block_ref, lgt_active_ref(:,1), lgt_n_ref(1), hvy_block, hvy_tmp, tree_ID1, verbosity=.False.)
            call coarse_tree_2_reference_mesh(params, lgt_block_ref, lgt_active_ref(:,2), lgt_n_ref(2), hvy_block, hvy_tmp, tree_ID2, verbosity=.False.)
        endif
    else
        if (tree_ID1 .ne. tree_ID2 .and. params%Jmax .ne. params%Jmin) then
            call coarse_tree_2_reference_mesh(params, lgt_block_ref, lgt_active_ref(:,2), lgt_n_ref(2), hvy_block, hvy_tmp, tree_ID2, verbosity=.False.)
        endif

        call updateMetadata_tree( params, tree_ID1 )

        call sync_ghosts_tree( params, hvy_block, tree_ID1)
    endif

    params%wavelet = MR



    !   if (present(dest_tree_ID)) then
    !     ! we have to synchronize lgt data since we were updating it locally on this procesor
    !     call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )
    !     call createActiveSortedLists_forest( params, lgt_block, lgt_active, &
    !     lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
    !     !
    !     ! update neighbor relations and synchronice ghosts of 1st tree
    !     call updateNeighbors_tree( params, lgt_block, hvy_neighbor, lgt_active(:, dest_tree_ID), &
    !     lgt_n(dest_tree_ID), lgt_sortednumlist(:,:,dest_tree_ID), hvy_active(:, dest_tree_ID), hvy_n(dest_tree_ID) )
    !     call sync_ghosts_tree( params, hvy_block, dest_tree_ID)
    !     if (tree_ID1 .ne. tree_ID2) then
    !         call coarse_tree_2_reference_mesh(params, tree_n, &
    !             lgt_block, lgt_active(:,tree_ID1),lgt_n, lgt_sortednumlist(:,:,tree_ID1), &
    !             lgt_block_ref, lgt_active_ref,lgt_n_ref, &
    !             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_ID, verbosity)
    !
    !         call adapt_tree( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_ID2), &
    !         lgt_n(tree_ID2), lgt_sortednumlist(:,:,tree_ID2), hvy_active(:,tree_ID2), &
    !         hvy_n(tree_ID2), tree_ID2, params%coarsening_indicator, hvy_tmp )
    !     endif
    ! else
    !
    !     ! update neighbor relations and synchronice ghosts of 1st tree
    !     call updateNeighbors_tree( params, lgt_block, hvy_neighbor, lgt_active(:, tree_ID1), &
    !     lgt_n(tree_ID1), lgt_sortednumlist(:,:,tree_ID1), hvy_active(:, tree_ID1), hvy_n(tree_ID1) )
    !     call sync_ghosts_tree( params, hvy_block, tree_ID1)
    !
    !     call adapt_tree( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_ID2), &
    !     lgt_n(tree_ID2), lgt_sortednumlist(:,:,tree_ID2), hvy_active(:,tree_ID2), &
    !     hvy_n(tree_ID2), tree_ID2, params%coarsening_indicator, hvy_tmp )
    ! end if
    call toc( "pointwise_tree_arithmetic (coarse to reference mesh)", 255, MPI_Wtime()-t_elapse )
end subroutine
!##############################################################

!##############################################################
! This routine sums up all trees in treeid_list and saves the sum
! in dest_tree_ID
subroutine sum_trees(params, hvy_block, hvy_tmp, treeid_list, dest_tree_ID, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: treeid_list(:), dest_tree_ID !< number of the tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    integer(kind=ik):: i
    logical :: verbose=.false.

    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose) write(*,'("Adding trees: ",i4)') treeid_list

    call copy_tree(params, hvy_block, dest_tree_ID, treeid_list(1))

    do i = 2, size(treeid_list)
        call add_two_trees(params, hvy_block, hvy_tmp, dest_tree_ID, treeid_list(i))
    enddo

end subroutine
!##############################################################

!##############################################################
!> This routine computes the average of all trees in treeid_list and
!> saves the result in dest_tree_ID
!> result = sum(trees) / number_trees
!> to construct a tree_ID list use for example: tree_ID_list= (/(id, id= 1, N)/)
subroutine average_trees(params, hvy_block, hvy_tmp, treeid_list, dest_tree_ID, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: treeid_list(:) !< List of tree ids you want to average
    integer(kind=ik), intent(in)      :: dest_tree_ID !< tree id of the averaged tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    integer(kind=ik):: N
    logical :: verbose=.false.

    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose) write(*,'("Averaging trees: ",i4)') treeid_list

    N = size(treeid_list)

    ! first sum up all trees from the given list
    call sum_trees(params, hvy_block, hvy_tmp, treeid_list, dest_tree_ID, verbosity)

    ! devide the sum by the numbers of elements
    call multiply_tree_with_scalar(params, hvy_block, dest_tree_ID, 1.0_rk/real(N, kind=rk), verbosity)

end subroutine
!##############################################################


!##############################################################
subroutine add_two_trees(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2, dest_tree_ID, verbosity,a,b)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    integer(kind=ik), intent(in), optional:: dest_tree_ID !< if specified result of addition will be saved here
    !< otherwise result will overwrite tree_ID1
    real(kind=rk), intent(in),optional      :: a,b ! scalar factors for addition: result = a * tree_ID1 + b *tree_ID2
    !-----------------------------------------------------------------
    logical :: verbose=.false.
    real(kind=rk) :: alpha(2)

    alpha = (/1 , 1/)
    if (present(a)) alpha(1) = a
    if (present(b)) alpha(2) = b
    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose) write(*,'("Adding trees: ",i4,",",i4)') tree_ID1, tree_ID2

    if(present(dest_tree_ID))then
        call tree_pointwise_arithmetic(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2,"+",dest_tree_ID,a=alpha(1),b=alpha(2))
    else
        call tree_pointwise_arithmetic(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2,"+",a = alpha(1),b = alpha(2))
    endif
end subroutine
!########################################################### ###

!##############################################################
subroutine substract_two_trees(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2, dest_tree_ID, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    integer(kind=ik), intent(in), optional:: dest_tree_ID !< if specified result of addition will be saved here
    !< otherwise result will overwrite tree_ID1
    !-----------------------------------------------------------------
    logical :: verbose=.false.

    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose) write(*,'("Substract trees: ",i4,",",i4)') tree_ID1, tree_ID2

    if(present(dest_tree_ID))then
        call tree_pointwise_arithmetic(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2,"-",dest_tree_ID)
    else
        call tree_pointwise_arithmetic(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2,"-")
    endif
end subroutine
!########################################################### ###



!##############################################################
! This function returns an scalar computed from the L2 scalar
! prodcut for 2 different trees
!  <f(x),g(x)> = int f(x) * g(x) dx
!!!!! Caution!!!!!
! The trees will be at finer levels after calling this routine, because
! we compute a unified grid of both trees to compute the scalarproduct.
! So you might call adapt_tree after this routine or
! all coarse_tree_2_reference_mesh.
function scalar_product_two_trees( params, hvy_block, hvy_tmp ,&
    tree_ID1, tree_ID2, verbosity)  result(sprod)
    implicit none

    !-----------------------------------------------------------------
    !> user defined parameter structure
    type (type_params), intent(inout)  :: params
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    !> heavy data array - block data
    real(kind=rk),  intent(inout)   :: hvy_block(:, :, :, :, :)
    !> heavy temp data: needed in blockxfer which is called in add_two_trees
    real(kind=rk),   intent(inout)  :: hvy_tmp(:, :, :, :, :)
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !---------------------------------------------------------------
    logical :: verbose=.false.
    character(len=cshort) :: MR
    integer(kind=ik)    :: free_tree_ID, Jmax, Bs(3), g, &
    N, k, rank, i, mpierr
    integer(kind=ik)    :: level1, level2, hvy_id1, hvy_id2, lgt_id1, lgt_id2
    integer(kind=tsize) :: treecode1, treecode2
    integer(kind=ik)    :: ix, iy, iz, iq, k1, k2, Nord, delta1, delta2, delta3

    real(kind=rk), allocatable , save ::  M(:)
    real(kind=8) :: sprod, Volume, t_elapse, t_inc(2), sprod_block, tmp
    real(kind=rk) :: x0(3), dx(3)

    integer(kind=ik) , save, allocatable :: lgt_active_ref(:,:), lgt_block_ref(:,:)
    integer(kind=ik) , save :: lgt_n_ref(2)=0_ik

    if (present(verbosity)) verbose=verbosity

    N = params%number_blocks
    rank = params%rank
    Jmax = params%Jmax
    g = params%g
    Bs= params%Bs
    Volume = product(params%domain_size(1:params%dim))

    ! depending on the predictor order, i.e. the order of the lagrange interpolation
    ! scheme. we have to use a different "FEM" Mass matrices. See:
    ! https://moodle.polymtl.ca/pluginfile.php/47880/mod_resource/content/0/week7.pdf
    ! The matrix M_ij = <phi_i,phi_j> where the phis are the lagrange basis elements
    if (params%order_predictor == "multiresolution_4th" ) then
        ! NORMALLY multiresolution_4th ONLY NEEDS 4 GHOST POINTS, BUT
        ! FOR THE WAVELET WEIGHTED SCALAR PRODUCT WE NEED 6
        ! ANYHOW IF ON CHOOSES TO HAVE ONLY 4 GHOST POINTS WE STILL COMPUTE
        ! THE SCALAR PRODUCT WITH A WARNING
        if (params%g < 5) then
            Nord = 1
            if(.not.allocated(M)) then
                allocate( M(Nord))
                M(1) = 1
                if (rank==0) write(*,*) " "
                if (rank==0) write(*,*) "--------------------------------------------"
                if (rank==0) write(*,*) "WARNING! this scalar product is imprecise,"
                if (rank==0) write(*,*) "because we do not make use of scaling functions properties,"
                if (rank==0) write(*,*) "Increase ghost nodes to n_ghosts=5 for ussage!!,"
                if (rank==0) write(*,*) "--------------------------------------------"


            endif
        else
            Nord = 6
            if (.not. allocated(M) ) then
                allocate( M(Nord))
                M  = (/ 0.8009680404299244, &
                0.13704178233326222, &
                -0.04024485728521606, &
                0.002795127767100863, &
                -7.592473960186964e-05, &
                -1.482905070349005e-07 /)
            endif
        endif
    elseif (params%order_predictor == "multiresolution_2nd" ) then
        ! linear lagrange polynomials
        Nord = 2
        if (.not. allocated(M) ) then
            allocate( M(Nord))
            M  = (/ 0.66666667, 0.16666667 /)
        endif
    else
        call abort(28754904, "The predictor method is unknown")
    endif
    !----------------------------------------------
    ! sprod = <X_i, X_j>
    t_elapse = MPI_wtime()

    t_inc(1) = t_elapse
    !=============================================
    ! Prepare the trees for pointwise arithmentic
    !=============================================
    ! The Trees can only be added when their grids are identical. At present we
    ! try to keep the finest levels of all trees. This means we refine
    ! all blocks which are not on the same level.
    if (tree_ID1 .ne. tree_ID2) then
        ! first we save the old tree structure because it will be lost after
        ! unification of both meshes
        call store_ref_meshes( lgt_block_ref, lgt_active_ref, lgt_n_ref, tree_ID1, tree_ID2)

        call refine_trees2same_lvl(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2)

        ! because all trees have the same treestructrue thier hilbertcurve is identical
        ! and therefore balance the load will try to distribute blocks with the same
        ! treecode (but on different trees) at the same rank.
        call balanceLoad_tree( params, hvy_block, tree_ID1 )

        call balanceLoad_tree( params, hvy_block, tree_ID2 )
    end if
    t_inc(2) = MPI_WTIME()

    sprod = 0.0_rk
    loop_tree1: do k1 = 1, hvy_n(tree_ID1)
        hvy_id1 = hvy_active(k1,tree_ID1)
        call hvy2lgt(lgt_id1, hvy_id1, rank, N )

        level1   = lgt_block(lgt_id1, IDX_MESH_LVL)
        treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))
        loop_tree2: do k2 = 1, hvy_n(tree_ID2) !loop over all treecodes in 2.tree to find the same block
            hvy_id2 = hvy_active(k2,tree_ID2)
            call hvy2lgt(lgt_id2, hvy_id2, rank, N )
            level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
            treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
            if (treecode1 .ne. treecode2 ) then
                cycle
            else
                call get_block_spacing_origin( params, lgt_id1, x0, dx )
                sprod_block = 0.0_rk
                if (params%dim==2) then
                    !##########################################
                    ! 2d
                    do iq = 1, params%N_eqn
                        do iy = g+1, Bs(2) + g - 1
                            do ix = g+1, Bs(1) + g - 1
                                tmp = 0.0_rk
                                do delta1 = -Nord+1, Nord-1
                                    do delta2 = -Nord+1, Nord-1
                                        tmp = tmp+ M(abs(delta1)+1) * M(abs(delta2)+1) *hvy_block(ix-delta1,iy-delta2,1,iq,hvy_id2)
                                    end do
                                end do
                                !    sprod_block = sprod_block + hvy_block(ix,iy,1,iq,hvy_id1) * &
                                !    hvy_block(ix,iy,1,iq,hvy_id2)
                                sprod_block = sprod_block + hvy_block(ix,iy,1,iq,hvy_id1) * tmp
                            end do
                        end do
                    end do
                    !##########################################
                else
                    !##########################################
                    ! 3D
                    do iq = 1, params%N_eqn
                        do iz = g+1, Bs(3) + g - 1
                            do iy = g+1, Bs(2) + g - 1
                                do ix = g+1, Bs(1) + g - 1
                                    tmp = 0.0_rk
                                    do delta1 = -Nord+1, Nord-1
                                        do delta2 = -Nord+1, Nord-1
                                            do delta3 = -Nord+1, Nord-1
                                                tmp = tmp+ &
                                                M(abs(delta1)+1) * &
                                                M(abs(delta2)+1) * &
                                                M(abs(delta3)+1) * &
                                                hvy_block(ix-delta1,iy-delta2,iz-delta3,iq,hvy_id2)
                                            end do
                                        end do
                                    end do
                                    sprod_block = sprod_block + hvy_block(ix,iy,iz,iq,hvy_id1) * tmp
                                end do
                            end do
                        end do
                    end do
                    !##########################################
                endif
                sprod =sprod + product(dx(1:params%dim))*sprod_block
                exit loop_tree2
            endif
        end do loop_tree2
    end do loop_tree1

    t_inc(1) = t_inc(2)-t_inc(1)
    t_inc(2) = MPI_wtime()-t_inc(2)
    !----------------------------------------------------
    ! sum over all Procs
    !----------------------------------------------------
    sprod_block = sprod
    call MPI_ALLREDUCE(sprod_block, sprod, 1, MPI_DOUBLE, &
    MPI_SUM,WABBIT_COMM, mpierr)
    !!! Attention we should not use prefilter of CDF44 wavelets here. Because
    !!! it will change the original data!!!
    MR = params%wavelet
    params%wavelet = "CDF40"

    if (tree_ID1 .ne. tree_ID2) then
        call coarse_tree_2_reference_mesh( params, lgt_block_ref, lgt_active_ref(:,1), lgt_n_ref(1), hvy_block, hvy_tmp, tree_ID1, verbosity=.False.)
        call coarse_tree_2_reference_mesh( params, lgt_block_ref, lgt_active_ref(:,2), lgt_n_ref(2), hvy_block, hvy_tmp, tree_ID2, verbosity=.False.)
    endif
    params%wavelet = MR

    t_elapse = MPI_WTIME() - t_elapse

    call toc( "scalra prod (prepare: refine_tree+balanceLoad_tree)", 256, t_inc(1) )
    call toc( "scalra prod (hvy_data operation)", 257, t_inc(2) )
    call toc( "scalra prod (total)", 258, t_elapse )
end function
!##############################################################





!##############################################################
! This function returns an scalar computed from the L2 scalar
! prodcut for 2 different trees
!  <f(x),g(x)> = int f(x) * g(x) dx
function scalar_product_two_trees_old( params, hvy_block, hvy_tmp, &
    tree_ID1, tree_ID2, buffer_tree_ID)  result(sprod)
    implicit none

    !-----------------------------------------------------------------
    !> user defined parameter structure
    type (type_params), intent(inout)  :: params
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    !> for the multiplication we need an additional tree as a buffer. If
    !> no tree is passed as additional argument we will create a new tree
    !> and delete it afterwards.
    !> If buffer_tree_ID is passed, we use it as a buffer and do not have
    !> to create or delete any additional tree
    integer(kind=ik) , optional, intent(in) :: buffer_tree_ID
    real(kind=rk),  intent(inout)   :: hvy_block(:, :, :, :, :)
    !> heavy temp data: needed in blockxfer which is called in add_two_trees
    real(kind=rk),   intent(inout)  :: hvy_tmp(:, :, :, :, :)
    !---------------------------------------------------------------
    integer(kind=ik) :: free_tree_ID, Jmax, Bs(3), g, &
    N, k, lgt_id, hvy_id, rank, i, mpierr
    real(kind=rk) :: sprod, Volume, t_elapse, t_inc(2)
    real(kind=rk) :: x0(3), dx(3)

    if ( present(buffer_tree_ID)) then
        free_tree_ID = buffer_tree_ID
    else
        free_tree_ID = tree_n + 1
    endif

    N = params%number_blocks
    rank = params%rank
    Jmax = params%Jmax
    g = params%g
    Bs= params%Bs
    Volume = product(params%domain_size(1:params%dim))

    !----------------------------------------------
    ! sprod = <X_i, X_j>
    t_elapse = MPI_wtime()

    !---------------------------------------------------
    ! multiply tree_ID2 * free_tree_ID -> free_tree_ID
    !---------------------------------------------------
    t_inc(1) = MPI_wtime()
    call multiply_two_trees(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2, free_tree_ID)

    t_inc(1) = MPI_wtime()-t_inc(2)
    !---------------------------------------------------
    ! adapt the mesh before summing it up
    !---------------------------------------------------
    if ( params%adapt_tree ) then
        !call adapt_tree_mesh( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        !lgt_n, lgt_sortednumlist, hvy_active, hvy_n, params%coarsening_indicator, hvy_tmp, &
        !free_tree_ID, tree_n )
    endif

    !----------------------------------------------------
    ! sum over all elements of the tree with free_tree_ID
    !-----------------------------------------------------
    t_inc(2) = MPI_wtime()
    sprod = 0.0_rk
    do k = 1, hvy_n(free_tree_ID)
        hvy_id = hvy_active(k, free_tree_ID)

        call hvy2lgt(lgt_id, hvy_id, rank, N )

        ! calculate the lattice spacing.
        ! It is needed to perform the L2 inner product
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if ( params%dim == 3 ) then
            sprod = sprod + dx(1)*dx(2)*dx(3)* sum( hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, :, hvy_id))
        else
            sprod = sprod + dx(1)*dx(2)*sum( hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, :, hvy_id))
        endif
    end do
    t_inc(2) = MPI_wtime()-t_inc(2)
    t_elapse = MPI_WTIME() - t_elapse
    !----------------------------------------------------
    ! sum over all Procs
    !----------------------------------------------------
    call MPI_ALLREDUCE(MPI_IN_PLACE, sprod, 1, MPI_DOUBLE_PRECISION, MPI_SUM,WABBIT_COMM, mpierr)

    if (.not. present(buffer_tree_ID)) then
        call delete_tree(params, free_tree_ID)
    endif
end function
!##############################################################




!##############################################################
subroutine multiply_two_trees(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2, dest_tree_ID, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout)    :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    integer(kind=ik), intent(in),optional :: dest_tree_ID !< number of the tree
    !-----------------------------------------------------------------
    logical :: verbose=.false.

    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose ) write(*,'("Multiply trees: ",i4,",",i4)') tree_ID1, tree_ID2

    if(present(dest_tree_ID))then
        call tree_pointwise_arithmetic(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2,"*",dest_tree_ID)
    else
        call tree_pointwise_arithmetic(params, hvy_block, hvy_tmp, tree_ID1, tree_ID2, "*")
    endif

end subroutine
!##############################################################



!##############################################################
!> This function compares the tree structure and block distribution of two trees.
!> If treecodes are identical but the processor ranks not, then we redistribute
!> one of the corresponding blocks, such that blocks with same treecode are on the
!> same rank.
subroutine same_block_distribution(params, hvy_block, tree_ID1, tree_ID2)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(in)      :: tree_ID1, tree_ID2 !< number of the tree
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    !-----------------------------------------------------------------
    integer(kind=ik)    :: level1, level2, lgt_id1, lgt_id2, N, fsize
    integer(kind=ik)    :: k1, k2, rank1, rank2, level_min, n_comm, Jmax
    integer(kind=tsize) :: treecode1, treecode2
    integer(kind=ik), allocatable, save :: comm_list(:,:)

    N = params%number_blocks
    fsize = params%forest_size
    Jmax = params%Jmax
    if (.not.allocated(comm_list)) allocate( comm_list( 3, params%number_procs*N) )

    !! Loop over all treecodes of both trees and check if they are identical.
    !! If the treecodes are the same but the blocks are not on the same rank,
    !! we have to send a block.
    !! The comm_list is created to save all of the block transfers.
    n_comm = 0 ! number of communications
    do k1 = 1, lgt_n(tree_ID1)
        !--------
        ! TREE 1
        !--------
        lgt_id1  = lgt_active(k1, tree_ID1)
        level1   = lgt_block(lgt_id1, IDX_MESH_LVL)

        do k2 = 1, lgt_n(tree_ID2)
            !--------
            ! TREE 2
            !--------
            lgt_id2  = lgt_active(k2, tree_ID2)
            level2   = lgt_block(lgt_id2, IDX_MESH_LVL)
            level_min= min(level1,level2)
            treecode1 = get_tc(lgt_block(lgt_id1, IDX_TC_1 : IDX_TC_2))
            treecode1 = tc_clear_until_level_b(treecode1, params%dim, level_min, Jmax)
            treecode2 = get_tc(lgt_block(lgt_id2, IDX_TC_1 : IDX_TC_2))
            treecode2 = tc_clear_until_level_b(treecode2, params%dim, level_min, Jmax)
            if (treecode1 == treecode2 )then
                ! check treestructure
                if (level1 .ne. level2) then
                    call abort(270219, "Trees need same treestructure to fetch processor block distribution")
                endif
                ! check the block distribution (leaves are on same rank)
                call lgt2proc( rank1, lgt_id1, N)
                call lgt2proc( rank2, lgt_id2, N)
                if (rank1 .ne. rank2) then
                    n_comm = n_comm + 1
                    !                write(*,*) "===============> not identical", n_comm
                    comm_list(1, n_comm) = rank1   ! sender mpirank
                    comm_list(2, n_comm) = rank2   ! receiver mpirank
                    comm_list(3, n_comm) = lgt_id1 ! block lgt_id to send
                end if
            end if
        end do ! loop over tree2
    end do ! loop over tree1
    !if (params%rank == 0) write(*,'("nr blocks redistributed: ",i9)') n_comm

    ! Transfer the blocks
    call block_xfer( params, comm_list, n_comm, hvy_block )

    ! after block transfer we have to create new lists
    call createActiveSortedLists_forest(params)

end subroutine
!##############################################################
