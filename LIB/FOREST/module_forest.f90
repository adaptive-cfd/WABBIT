
!-----------------------------------------------------------------
!> \file
!> \brief
!! Module for implementing multiple trees in a forest
!> \version 10.1.2019
!> \author P.Krah
!-----------------------------------------------------------------

module module_forest

    ! modules
    use mpi
    use module_precision
    use module_ini_files_parser_mpi, only : read_param_mpi
    use module_params
    use module_mesh

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    !**********************************************************************************************
    ! These are the important routines that are visible to WABBIT:
    !**********************************************************************************************
    PUBLIC :: add_two_trees, count_tree_hvy_n, average_trees, &
    copy_tree, multiply_two_trees, multiply_tree_with_scalar, &
    compute_tree_L2norm, delete_tree, scalar_product_two_trees, &
    same_block_distribution, prune_tree, add_pruned_to_full_tree, refine_tree
    !**********************************************************************************************

contains


    subroutine prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id)
        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, :)  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:, :)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)

        integer(kind=ik) :: k, lgt_id, hvy_id, rank, N, g ,Bs(3)

        rank = params%rank
        N = params%number_blocks
        g = params%n_ghosts
        Bs = params%Bs


        if (rank==0) write(*,'("Tree-pruning, before Nb=",i7)') lgt_n(tree_id)

        if (params%dim == 3) then
            do k = 1, hvy_n(tree_id)
                hvy_id = hvy_active(k, tree_id)
                call hvy_id_to_lgt_id( lgt_id, hvy_id, rank, N )

                ! pruning condition: all entries of the block are small
                if (.not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id) > 0.0_rk) ) then
                    ! pruning: delete the block from the tree
                    lgt_block(lgt_id, :) = -1_ik
                endif
                ! if ( (.not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id) > 0.0_rk)) .and. &
                ! (.not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 6, hvy_id) > 0.0_rk))) then
                !     ! pruning: delete the block from the tree
                !     lgt_block(lgt_id, :) = -1_ik
                ! endif
            end do
        else
            do k = 1, hvy_n(tree_id)
                hvy_id = hvy_active(k, tree_id)
                call hvy_id_to_lgt_id( lgt_id, hvy_id, rank, N )

                ! pruning condition: all entries of the block are small
                if ( .not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1, hvy_id) > 0.0_rk) ) then
                    ! pruning: delete the block from the tree
                    lgt_block(lgt_id, :) = -1_ik
                endif
                ! if ( (.not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1, hvy_id) > 0.0_rk)) .and. &
                ! (.not. any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 6, hvy_id) > 0.0_rk)) ) then
                !     ! pruning: delete the block from the tree
                !     lgt_block(lgt_id, :) = -1_ik
                ! endif
            end do
        endif

        call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
        hvy_active, hvy_n, lgt_sortednumlist, tree_n)

        if (rank==0) write(*,'("Tree-pruning, pruned to Nb=",i7)') lgt_n(tree_id)
        ! do not call update_neighbors on the pruned tree..
    end subroutine

!---------------------------------------------------------------------------------

    ! subroutine add_pruned_to_full_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    !     hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_pruned, tree_id_full)
    !     implicit none
    !     !-----------------------------------------------------------------
    !     type (type_params), intent(in) :: params   !< params structure
    !     integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
    !     integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
    !     integer(kind=ik), intent(in)      :: tree_id_pruned, tree_id_full
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
    !     Jmax = params%max_treelevel ! max treelevel
    !     rank = params%rank
    !     N = params%number_blocks
    !
    !     ! NOTES:
    !     !   We assume that the pruned tree is on JMAX (rhs level) and we assume that
    !     !   the full one is not finer, only coarser.
    !
    !     if (.not.allocated(comm_list)) allocate( comm_list( params%number_procs*N, 3 ) )
    !
    !     call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
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
    !     do k = 1, lgt_n(tree_id_pruned)
    !
    !         lgt_id1 = lgt_active(k, tree_id_pruned)
    !         level1  = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
    !
    !         ! does the block from the pruned tree exist in the full one?
    !         call does_block_exist( lgt_block(lgt_id1, 1:level1), exists, lgt_id2, &
    !         lgt_sortednumlist(:,:,tree_id_full), lgt_n(tree_id_full), tree_id_full)
    !
    !         if (exists) then
    !             ! we found the pruned trees block in the full tree : we happily copy
    !             ! its heavy data
    !
    !             ! now check on which CPU this block is currently
    !             call lgt_id_to_proc_rank( rank_pruned, lgt_id1, N)
    !             call lgt_id_to_proc_rank( rank_full, lgt_id2, N)
    !
    !             if (rank_full /= rank_pruned) then
    !                 n_comm = n_comm + 1
    !                 comm_list(n_comm, 1) = rank_pruned   ! sender mpirank
    !                 comm_list(n_comm, 2) = rank_full   ! receiver mpirank
    !                 comm_list(n_comm, 3) = lgt_id1 ! block lgt_id to send
    !                 write(*,*) "found on different rank", n_comm, rank_pruned, rank_full
    !             endif
    !         else
    !             ! it does not exist, but maybe it is coarsened by one level?
    !             ! NOTE: code will not detect if coarsened by more than one level
    !             call does_block_exist( lgt_block(lgt_id1, 1:level1-1), exists, lgt_id2, &
    !             lgt_sortednumlist(:,:,tree_id_full), lgt_n(tree_id_full), tree_id_full)
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
    !     call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    !     lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
    !
    !
    !     ! Step 2: ADDITION. now we're sure that blocks existing in both trees are on the
    !     ! same mpirank. therefore, the responsible rank can just add them together.
    !     do k = 1, hvy_n(tree_id_pruned)
    !
    !         hvy_id1  = hvy_active(k, tree_id_pruned)
    !         call hvy_id_to_lgt_id(lgt_id1, hvy_id1, params%rank, params%number_blocks)
    !         level1  = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
    !
    !         call does_block_exist(lgt_block(lgt_id1, 1:level1), exists, lgt_id2, &
    !         lgt_sortednumlist(:,:,tree_id_full), lgt_n(tree_id_full), tree_id_full)
    !
    !         ! we found the pruned trees block in the full tree : we happily copy
    !         ! its heavy data
    !         if (exists) then
    !             ! now check on which CPU this block is currently
    !             call lgt_id_to_proc_rank( rank_pruned, lgt_id1, N)
    !             call lgt_id_to_proc_rank( rank_full, lgt_id2, N)
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
    !                 call lgt_id_to_hvy_id( hvy_id1, lgt_id1, rank_pruned, N )
    !                 call lgt_id_to_hvy_id( hvy_id2, lgt_id2, rank_full  , N )
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
    !             do i = 1, lgt_n(tree_id_full)
    !                 lgt_id2 = lgt_active(i, tree_id_full)
    !                 treecode2 = treecode2int( lgt_block(lgt_id2, 1:level1) )
    !
    !                 if (treecode1 == treecode2) then
    !                     ! this is one of the sisters
    !                     call lgt_id_to_proc_rank( rank_full, lgt_id2, N)
    !                     if (params%rank == rank_full) then
    !                         call lgt_id_to_hvy_id( hvy_id2, lgt_id2, rank_full, N)
    !                         hvy_block(:,:,:,1,hvy_id2) = hvy_block(:,:,:,1,hvy_id2) + 1.0_rk
    !                     endif
    !                 endif
    !             enddo
    !         endif
    !     enddo
    !
    ! end subroutine
    subroutine add_pruned_to_full_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_pruned, tree_id_full)
        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id_pruned, tree_id_full
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, :)  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)

        integer(kind=ik) :: k, lgt_id, Jmax, hvy_id, rank, N, fsize, i
        integer(kind=ik) :: lgt_id1, lgt_id2, hvy_id1, hvy_id2
        integer(kind=ik) :: level1, level2, rank_pruned, rank_full, n_comm
        logical :: exists
        integer(kind=ik), allocatable, save :: comm_list(:,:)
        integer(kind=tsize) :: treecode1, treecode2

        fsize = params%forest_size
        Jmax = params%max_treelevel ! max treelevel
        rank = params%rank
        N = params%number_blocks

        ! NOTES:
        !   We assume that the pruned tree is on JMAX (rhs level) and we assume that
        !   the full one is finer, never coarser.

        if (.not.allocated(comm_list)) allocate( comm_list( params%number_procs*N, 3 ) )

        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
        hvy_active, hvy_n, lgt_sortednumlist, tree_n)

        ! a pruned tree has fewer entries: loop over it instead of the other one?
        ! if you find a block in the full tree -> well then that's good, copy.
        ! else: the target grid is either refined or coarsened at this position.

        ! Step 1: XFER. we look for blocks that exist in both pruned and full tree and
        ! if they are on different mpiranks, we xfer the pruned trees block to the rank
        ! holding the corresponding full trees block.

        ! Step 1a: prepare for xfer, gather all required xfers
        n_comm = 0
        do k = 1, lgt_n(tree_id_pruned)

            lgt_id1 = lgt_active(k, tree_id_pruned)
            level1  = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)

            call does_block_exist( lgt_block(lgt_id1, 1:level1), exists, lgt_id2, &
            lgt_sortednumlist(:,:,tree_id_full), lgt_n(tree_id_full), tree_id_full)

            if (exists) then
                ! we found the pruned trees block in the full tree : we happily copy
                ! its heavy data

                ! now check on which CPU this block is currently
                call lgt_id_to_proc_rank( rank_pruned, lgt_id1, N)
                call lgt_id_to_proc_rank( rank_full, lgt_id2, N)

                if (rank_full /= rank_pruned) then
                    n_comm = n_comm + 1
                    comm_list(n_comm, 1) = rank_pruned   ! sender mpirank
                    comm_list(n_comm, 2) = rank_full   ! receiver mpirank
                    comm_list(n_comm, 3) = lgt_id1 ! block lgt_id to send
                    ! write(*,*) "found on different rank", n_comm, rank_pruned, rank_full
                endif
            endif
        enddo

        ! Step 1b: actual xfer.
        call block_xfer( params, comm_list, n_comm, lgt_block, hvy_block )

        call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
        lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)


        ! Step 2: ADDITION. now we're sure that blocks existing in both trees are on the
        ! same mpirank. therefore, the responsible rank can just add them together.
        do k = 1, hvy_n(tree_id_pruned)

            hvy_id1  = hvy_active(k, tree_id_pruned)
            call hvy_id_to_lgt_id(lgt_id1, hvy_id1, params%rank, params%number_blocks)
            level1  = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)

            call does_block_exist(lgt_block(lgt_id1, 1:level1), exists, lgt_id2, &
            lgt_sortednumlist(:,:,tree_id_full), lgt_n(tree_id_full), tree_id_full)

            ! we found the pruned trees block in the full tree : we happily copy
            ! its heavy data
            if (exists) then
                ! now check on which CPU this block is currently
                call lgt_id_to_proc_rank( rank_pruned, lgt_id1, N)
                call lgt_id_to_proc_rank( rank_full, lgt_id2, N)

                if (rank_full /= rank_pruned) then
                    ! This is an erorr one should never see: the block xfer is done, and coexisting blocks should
                    ! be on the same mpirank now.
                    call abort(030719,"pruned_to_full_tree: although we xferred, we found a coexisting block on a different rank.")
                endif

                ! only responsible mpirank can perform the addition
                if (params%rank==rank_full) then
                    ! get both heavy ids
                    call lgt_id_to_hvy_id( hvy_id1, lgt_id1, rank_pruned, N )
                    call lgt_id_to_hvy_id( hvy_id2, lgt_id2, rank_full  , N )

                    ! actual addition
                    ! hvy_block(:,:,:,:,hvy_id2) = hvy_block(:,:,:,:,hvy_id2) + hvy_block(:,:,:,:,hvy_id1)
                    ! hvy_block(:,:,:,:,hvy_id2) = max(hvy_block(:,:,:,:,hvy_id2),hvy_block(:,:,:,:,hvy_id1))

                    where ( hvy_block(:,:,:,1,hvy_id2)<hvy_block(:,:,:,1,hvy_id1) )

                        hvy_block(:,:,:,1,hvy_id2) = hvy_block(:,:,:,1,hvy_id1)
                        hvy_block(:,:,:,2,hvy_id2) = hvy_block(:,:,:,2,hvy_id1)
                        hvy_block(:,:,:,3,hvy_id2) = hvy_block(:,:,:,3,hvy_id1)
                        hvy_block(:,:,:,4,hvy_id2) = hvy_block(:,:,:,4,hvy_id1)
                        hvy_block(:,:,:,5,hvy_id2) = hvy_block(:,:,:,5,hvy_id1)
                    end where

                endif

            else
                ! we did not find it. The grid has changed in the interior of the
                ! obstacle, and we can set those interior blocks to constant 1. But
                ! we need to figure out which blocks to set to 1.

                ! we can further assume that all blocks in the pruned tree are the minimum:
                ! only the interface is very important, interior blocks are dictated by
                ! gradedness. Hence: assuming that the fluid/solid interface is always
                ! on the finest level, only a refinement in the interior is possible (and
                ! not a coarsening) => look for sister blocks on higher levels

                ! find all sister blocks and add 1 to them. no xfer required.
                treecode1 = treecode2int( lgt_block(lgt_id1, 1:level1) )

                do i = 1, lgt_n(tree_id_full)
                    lgt_id2 = lgt_active(i, tree_id_full)
                    treecode2 = treecode2int( lgt_block(lgt_id2, 1:level1) )

                    if (treecode1 == treecode2) then
                        ! this is one of the sisters
                        call lgt_id_to_proc_rank( rank_full, lgt_id2, N)

                        if (params%rank == rank_full) then
                            call lgt_id_to_hvy_id( hvy_id2, lgt_id2, rank_full, N)
                            ! hvy_block(:,:,:,1,hvy_id2) = hvy_block(:,:,:,1,hvy_id2) + 1.0_rk
                            hvy_block(:,:,:,:,hvy_id2) = hvy_block(:,:,:,:,hvy_id1)
                        endif
                    endif
                enddo
            endif
        enddo

    end subroutine



    !##############################################################
    !> deletes the lgt_block data of tree with given tree_id
    !> CAUTION: active lists will be outdated!!!
    subroutine delete_tree(params, lgt_block, lgt_active, lgt_n, tree_id)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in)   :: params    !< user defined parameter structure
        integer(kind=ik), intent(inout)  :: lgt_block(:, :)!< light data array
        integer(kind=ik), intent(in)     :: lgt_active(:,:)!< list of active blocks (light data)
        integer(kind=ik), intent(in)     :: lgt_n(:)!< number of active blocks (light data)
        integer(kind=ik), intent(in)     :: tree_id!< highest tree id
        !-----------------------------------------------------------------
        integer(kind=ik)                 :: k, lgt_id

        ! loop over active list of tree
        do k = 1, lgt_n(tree_id)
            lgt_id = lgt_active(k, tree_id)
            lgt_block(lgt_id, :) = -1_ik
        end do


    end subroutine delete_tree
    !##############################################################



    !##############################################################
    !> copy hvy and lgt data of tree_id_source to tree_id_dest
    !> after copy the active and sorted lists are updated, as well
    !> as neighbors and ghost nodes
    subroutine copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_dest, tree_id_source)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in)    :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id_dest, tree_id_source !< all data from tree_id_source gets copied to destination_tree_id
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, :)  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        !-----------------------------------------------------------------
        integer(kind=ik)    :: Jmax, lgt_id_dest, lgt_id_source, hvy_id_dest, hvy_id_source, fsize
        integer(kind=ik)    :: k, N, rank
        real(kind=rk) :: t_elapse

        Jmax = params%max_treelevel ! max treelevel
        fsize= params%forest_size   ! maximal number of trees in forest
        rank = params%rank
        N = params%number_blocks

        if (tree_id_dest > fsize) call abort(0403191,"tree-copy: destination treeID ou of valid range")

        ! first we delete tree_id_dest if it is already allocated.
        ! tree_id_dest only exists if it is in the list of active trees, i.e. tree_id_dest <= tree_n
        if (tree_id_dest <= tree_n) then
            ! Caution: active lists will be outdated
            call delete_tree(params, lgt_block, lgt_active, lgt_n, tree_id_dest)
        end if


        ! Loop over the active hvy_data
        t_elapse = MPI_WTIME()
        do k = 1, hvy_n(tree_id_source)
            hvy_id_source = hvy_active(k, tree_id_source)
            call hvy_id_to_lgt_id( lgt_id_source, hvy_id_source, rank, N )

            ! first we have to find out if the hvy data belongs to the tree we want to copy
            call get_free_local_light_id( params, rank, lgt_block, lgt_id_dest)
            call lgt_id_to_hvy_id( hvy_id_dest, lgt_id_dest, rank, N )
            !--------------------
            ! Light DATA
            !--------------------
            lgt_block(lgt_id_dest, :) = lgt_block(lgt_id_source, :)  ! copy light data
            lgt_block(lgt_id_dest, Jmax + IDX_TREE_ID) = tree_id_dest  ! asign new tree_id
            !--------------------
            ! Heavy DATA
            !--------------------
            hvy_block( :, :, :, :, hvy_id_dest) = hvy_block( :, :, :, :, hvy_id_source)
        end do ! loop over source tree
        call toc( "copy_tree (copy heavy_data)", MPI_Wtime()-t_elapse )


        ! always synchronize lgt_data when you have changed lgt_block locally
        ! (i.e. when looping over lgt_ids which originate from a hvy_id)
        t_elapse = MPI_WTIME()
        call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

        call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
        lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_dest),&
        lgt_n(tree_id_dest), lgt_sortednumlist(:,:,tree_id_dest), hvy_active(:,tree_id_dest), hvy_n(tree_id_dest) )

        call toc( "copy_tree (copy synchronize hvy and lgt)", MPI_Wtime()-t_elapse )

    end subroutine
    !##############################################################

    !##############################################################
    !> multiply every block of a tree with a given value alpha
    subroutine multiply_tree_with_scalar(params, hvy_block, hvy_active, hvy_n, tree_id, &
        alpha, verbosity)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(in)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(in)   :: tree_id !< all data from tree_id2 gets copied to tree_id1
        real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        real(kind=rk), intent(in)      :: alpha !< heavy data array - block data
        integer(kind=ik), intent(in)   :: hvy_active(:, :) !< active lists
        logical, intent(in),optional   :: verbosity !< if true aditional stdout is printed
        !-----------------------------------------------------------------
        integer(kind=ik)    :: hvy_id, k
        logical :: verbose = .false.

        if (present(verbosity)) verbose=verbosity
        if (params%rank == 0 .and. verbose ) write(*,'("scalar multiplication tree: ",i3)') tree_id

        ! Loop over the active hvy_data, multiply block data by alpha
        do k = 1, hvy_n(tree_id)
            hvy_id = hvy_active(k, tree_id)
            hvy_block( :, :, :, :, hvy_id) = alpha * hvy_block( :, :, :, :, hvy_id)
        end do


    end subroutine
    !##############################################################

    !##############################################################
    !> L2 norm of all components in the statevector:
    !             norm(q_f) = sqrt( sum_ijk q_f^2 * dx_i dx_j dz_k )
    function compute_tree_L2norm(params, lgt_block, hvy_block, hvy_active, hvy_n, dF_opt, &
        tree_id_opt, verbosity) result(L2norm)

        implicit none
        !-----------------------------------------------------------------
        ! inputs:
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(in)   :: hvy_n(:)    !< number of active heavy blocks
        real(kind=rk), intent(in)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(in)   :: lgt_block(:, :)
        integer(kind=ik), intent(in)   :: hvy_active(:, :) !< active lists
        integer(kind=ik), intent(in),optional:: dF_opt ! statevector index (dF_opt = 1...,number equations) default:1
        integer(kind=ik), intent(in),optional:: tree_id_opt !< which tree is used (default:1)
        logical, intent(in),optional      :: verbosity !< if true aditional stdout is printed
        !-----------------------------------------------------------------
        ! result
        real(kind=rk) :: L2norm
        !-----------------------------------------------------------------
        integer(kind=ik)    :: lgt_id, hvy_id, ierr, tree_id = 1
        integer(kind=ik)    :: k, N, rank, g, Bs(3), dF=1
        real(kind=rk) :: norm, x0(3), dx(3)
        logical :: verbose=.false.

        if (present(dF_opt)) dF = dF_opt
        if (present(tree_id_opt)) tree_id = tree_id_opt
        if (present(verbosity)) verbose = verbosity
        rank = params%rank
        N = params%number_blocks
        Bs= params%Bs
        g = params%n_ghosts

        norm = 0.0_rk
        ! Loop over the active hvy_data
        do k = 1,hvy_n(tree_id)
            hvy_id = hvy_active(k, tree_id)
            call hvy_id_to_lgt_id(lgt_id, hvy_id, rank, N)
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
            ! ----------------------------
            ! integrate squared quantity
            ! ----------------------------
            if (params%dim == 3) then
                norm = norm + sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, &
                g+1:Bs(3)+g-1,dF , hvy_id)**2)*dx(1)*dx(2)*dx(3)
            else
                norm = norm + sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, &
                dF, hvy_id)**2)*dx(1)*dx(2)
            endif
        end do
        ! -------------------
        ! sum over all blocks
        ! -------------------
        call MPI_ALLREDUCE( norm, L2norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)

        L2norm = sqrt(L2norm)
        if (params%rank == 0 .and. verbose ) write(*,'("L2 norm: ",e12.6)') L2norm
    end function
    !##############################################################



    !##############################################################
    !> This function takes two trees and refines their treestructures
    !> until they are refined to the same max level.
    !> The resulting trees are identical in treestructure. Hence
    !> they have the same treecodes and number of blocks.
    !> Remark: You may want to balance the load after using this routine.
    subroutine refine_trees2same_lvl(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, verbosity)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
        logical, intent(in),optional      :: verbosity
        !-----------------------------------------------------------------
        integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
        integer(kind=ik)    :: k1, k2, Nblocks_2refine, level_min
        integer(kind=tsize) :: treecode1, treecode2
        logical :: verbose = .false.

        Jmax = params%max_treelevel ! max treelevel
        fsize= params%forest_size   ! maximal number of trees in forest

        if (present(verbosity)) verbose=verbosity
        if ( params%rank == 0 .and. verbose ) write(*,'("Refining trees to same level: ",i9,",",i9)') tree_id1, tree_id2

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
            do k1 = 1, lgt_n(tree_id1)
                !--------
                ! TREE 1
                !--------
                lgt_id1 = lgt_active(k1, tree_id1)
                level1  = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)

                do k2 = 1, lgt_n(tree_id2)
                    !--------
                    ! TREE 2
                    !--------
                    lgt_id2  = lgt_active(k2, tree_id2)
                    ! Skip this block if it is already tagged for refinement.
                    if ( lgt_block(lgt_id2, Jmax+IDX_REFINE_STS) == +1 ) then
                        cycle
                    endif
                    level2 = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)

                    ! The treecodes can only be compared if they have the same size
                    ! (i.e. treecode length). Therefore we have to find the minimum of both
                    ! levels.
                    level_min = min(level1, level2)
                    treecode1 = treecode2int( lgt_block(lgt_id1, 1:level_min) )
                    treecode2 = treecode2int( lgt_block(lgt_id2, 1:level_min) )

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
                            lgt_block(lgt_id2, Jmax + IDX_REFINE_STS) = +1
                        else
                            ! mark the block on tree 1 for refinement
                            lgt_block(lgt_id1, Jmax + IDX_REFINE_STS) = +1
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
                call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active(:, tree_id1), lgt_n(tree_id1), &
                lgt_sortednumlist(:,:,tree_id1), hvy_active(:,tree_id1), hvy_n(tree_id1) )

                call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active(:, tree_id2), lgt_n(tree_id2), &
                lgt_sortednumlist(:,:,tree_id2), hvy_active(:,tree_id2), hvy_n(tree_id2) )


                ! 2) refine blocks
                if ( params%dim==3 ) then
                    ! 3D:
                    call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active(:,tree_id1), hvy_n(tree_id1) )
                    call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active(:,tree_id2), hvy_n(tree_id2) )
                else
                    ! 2D:
                    call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:),&
                    hvy_active(:,tree_id1), hvy_n(tree_id1) )

                    call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:),&
                    hvy_active(:,tree_id2), hvy_n(tree_id2) )
                end if

                ! since lgt_block was synced we have to create the active lists again
                call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
                lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n )

                ! update neighbor relations
                call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:, tree_id1), &
                lgt_n(tree_id1), lgt_sortednumlist(:,:,tree_id1), hvy_active(:, tree_id1), hvy_n(tree_id1) )

                call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:, tree_id2), &
                lgt_n(tree_id2), lgt_sortednumlist(:,:,tree_id2), hvy_active(:, tree_id2), hvy_n(tree_id2) )

                ! synchronize ghosts
                call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_id1), hvy_n(tree_id1))
                call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_id2), hvy_n(tree_id2))
            endif
        end do
    end subroutine
    !##############################################################
     subroutine refine_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
        !-----------------------------------------------------------------
        integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
        integer(kind=ik)    :: k
        integer(kind=tsize) :: treecode1, treecode2
        logical :: verbose = .true.

        Jmax = params%max_treelevel ! max treelevel
        fsize= params%forest_size   ! maximal number of trees in forest

        do k = 1, lgt_n(tree_id)
            if (lgt_block( lgt_active(k,tree_id), Jmax+IDX_MESH_LVL)<12) then
                lgt_block( lgt_active(k,tree_id), Jmax+IDX_REFINE_STS) = +1
            endif
        enddo

        call respect_min_max_treelevel( params, lgt_block, lgt_active(:,tree_id), lgt_n(tree_id) )
        ! 2) refine blocks
        if ( params%dim==3 ) then
            ! 3D:
            call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active(:,tree_id), hvy_n(tree_id) )
        else
            ! 2D:
            call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:),&
            hvy_active(:,tree_id), hvy_n(tree_id) )
        end if

        ! since lgt_block was synced we have to create the active lists again
        call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
        lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n )

        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:, tree_id), &
        lgt_n(tree_id), lgt_sortednumlist(:,:,tree_id), hvy_active(:, tree_id), hvy_n(tree_id) )

        ! synchronize ghosts
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_id), hvy_n(tree_id))
    end subroutine
    !##############################################################



    !#############################################################
    !> Perform pointwise operations (+, -, /, *) for two given trees.
    !> The trees will be refined to same treestructure and then the
    !> operation will be executed on each block pointwise.
    !> Its like the Matlabs .* or ./ for trees.
    !> NOTE: if dest_tree_id is present then the result will be saved
    !>       to the given dest_tree_id. Otherwise the operation will
    !>       be saved on the first tree_id1
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
    subroutine tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, &
        operation, dest_tree_id)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
        character (len=*), intent(in)  :: operation !< which arithmetical operation (+,-,*,/) which is applied
        integer(kind=ik),optional, intent(in)::dest_tree_id !< optional for saving results to destination tree id
        !-----------------------------------------------------------------
        integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
        integer(kind=ik)    :: k1, k2, N, level_min, rank1, rank2, iq, iz, iy, ix, &
        hvy_id1, hvy_id2, Bs(3), g, lgt_id_dest, hvy_id_dest
        integer(kind=tsize) :: treecode1, treecode2
        real (kind=rk) :: t_elapse
        character (len=5) :: op !< which arithmetical operation (+,-,*,/) which is applied

        Jmax = params%max_treelevel ! max treelevel
        fsize= params%forest_size   ! maximal number of trees in forest
        N    = params%number_blocks ! number of blocks per rank
        rank = params%rank       ! proc rank
        g = params%n_ghosts         ! number of ghost nodes
        Bs= params%Bs               ! number of grid points per block

        ! decide if inplace or out of place operation
        if(present(dest_tree_id))then
            ! operation will be out of place: result saved in dest_id
            op = trim(operation)//"out"
            ! delete tree_id_dest if it is already allocated.
            ! tree_id_dest only exists if it is in the list of active trees, i.e. tree_id_dest <= tree_n
            if (dest_tree_id <= tree_n) then
                ! Caution: active lists will be outdated
                call delete_tree(params, lgt_block, lgt_active, lgt_n, dest_tree_id)
                call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
                lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n )
            end if
        else
          op = trim(operation)
        endif
        !=============================================
        ! Prepare the trees for pointwise arithmentic
        !=============================================
        ! The Trees can only be added when their grids are identical. At present we
        ! try to keep the finest levels of all trees. This means we refine
        ! all blocks which are not on the same level.

        t_elapse = MPI_WTIME()
        call refine_trees2same_lvl(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2)
        call toc( "pointwise_tree_arithmetic (refine_trees2same_lvl)", MPI_Wtime()-t_elapse )

        ! because all trees have the same treestructrue thier hilbertcurve is identical
        ! and therefore balance the load will try to distribute blocks with the same
        ! treecode (but on different trees) at the same rank.
        t_elapse = MPI_WTIME()
        call balance_load( params, lgt_block, hvy_block,  hvy_neighbor, &
        lgt_active(:, tree_id1), lgt_n(tree_id1), lgt_sortednumlist(:,:,tree_id1), &
        hvy_active(:, tree_id1), hvy_n(tree_id1), tree_id1 )

        ! since lgt_block was synced we have to create the active lists again
        call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
        lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n )

        call balance_load( params, lgt_block, hvy_block,  hvy_neighbor, &
        lgt_active(:, tree_id2), lgt_n(tree_id2), lgt_sortednumlist(:,:,tree_id2), &
        hvy_active(:, tree_id2), hvy_n(tree_id2), tree_id2 )
        call toc( "pointwise_tree_arithmetic (balancing)", MPI_Wtime()-t_elapse )

        ! since lgt_block was synced we have to create the active lists again
        call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
        lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n )

        call same_block_distribution(params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, tree_n, tree_id1, tree_id2)

        !=================================================
        ! Decide which pointwice arithmetic shell be used
        !=================================================
        !      _
        !     / \    from here we assume that trees have the same grid structure
        !    /   \   and the blocks of the same treecode are on the same rank !
        !   /     \
        !  /caution\
        ! +--------+
        select case(op)
        case("+out")
            !         #
            !         #
            !    ###########           ADDITION out of place
            !         #
            !         #
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1,tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2,tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        ! copy light data from one of the added trees
                        ! first we have to find out if the hvy data belongs to the tree we want to copy
                        call get_free_local_light_id( params, rank, lgt_block, lgt_id_dest)
                        call lgt_id_to_hvy_id( hvy_id_dest, lgt_id_dest, rank, N )
                        !--------------------
                        ! Light DATA
                        !--------------------
                        lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                        lgt_block(lgt_id_dest, Jmax + IDX_TREE_ID) = dest_tree_id  ! asign new tree_id

                        if (params%dim==2) then
                            !##########################################
                            ! 2d
                            do iq = 1, params%N_eqn
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,1,iq,hvy_id_dest) = hvy_block(ix,iy,1,iq,hvy_id1) + &
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
                                            hvy_block(ix,iy,iz,iq,hvy_id_dest) = hvy_block(ix,iy,iz,iq,hvy_id1) + &
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
        case("+")
            !         #
            !         #
            !    ###########           ADDITION
            !         #
            !         #
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1,tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2,tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        if (params%dim==2) then
                            !##########################################
                            ! 2d
                            do iq = 1, params%N_eqn
                                do iy = g+1, Bs(2) + g
                                    do ix = g+1, Bs(1) + g
                                        hvy_block(ix,iy,1,iq,hvy_id1) = hvy_block(ix,iy,1,iq,hvy_id1) + &
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
                                            hvy_block(ix,iy,iz,iq,hvy_id1) = hvy_block(ix,iy,iz,iq,hvy_id1) + &
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
        case("-")

            !
            !
            !    ###########           SUBSTRACTION
            !
            !
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1, tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                ! we want to add everything to tree1
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2,tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        !> \todo make this one faster by looping only over inner grid points
                        hvy_block(:,:,:,:,hvy_id1) = hvy_block(:,:,:,:,hvy_id1) - &
                        hvy_block(:,:,:,:,hvy_id2)
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
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1, tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                ! we want to add everything to tree1
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2, tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
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
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1,tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2,tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        ! copy light data from one of the added trees
                        ! first we have to find out if the hvy data belongs to the tree we want to copy
                        call get_free_local_light_id( params, rank, lgt_block, lgt_id_dest)
                        call lgt_id_to_hvy_id( hvy_id_dest, lgt_id_dest, rank, N )
                        !--------------------
                        ! Light DATA
                        !--------------------
                        lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                        lgt_block(lgt_id_dest, Jmax + IDX_TREE_ID) = dest_tree_id  ! asign new tree_id

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
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1,tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2,tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        ! copy light data from one of the added trees
                        ! first we have to find out if the hvy data belongs to the tree we want to copy
                        call get_free_local_light_id( params, rank, lgt_block, lgt_id_dest)
                        call lgt_id_to_hvy_id( hvy_id_dest, lgt_id_dest, rank, N )
                        !--------------------
                        ! Light DATA
                        !--------------------
                        lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                        lgt_block(lgt_id_dest, Jmax + IDX_TREE_ID) = dest_tree_id  ! asign new tree_id

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
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1,tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                ! we want to add everything to tree1
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2,tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
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
            do k1 = 1, hvy_n(tree_id1)
                hvy_id1 = hvy_active(k1,tree_id1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, rank, N )
                level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n(tree_id2)
                    hvy_id2 = hvy_active(k2,tree_id2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
                    level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        ! copy light data from one of the added trees
                        ! first we have to find out if the hvy data belongs to the tree we want to copy
                        call get_free_local_light_id( params, rank, lgt_block, lgt_id_dest)
                        call lgt_id_to_hvy_id( hvy_id_dest, lgt_id_dest, rank, N )
                        !--------------------
                        ! Light DATA
                        !--------------------
                        lgt_block(lgt_id_dest, :) = lgt_block(lgt_id1, :)  ! copy light data
                        lgt_block(lgt_id_dest, Jmax + IDX_TREE_ID) = dest_tree_id  ! asign new tree_id

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
            call abort(135,"Operation unknown")
        end select
        call toc( "pointwise_tree_arithmetic (hvy_data operation)", MPI_Wtime()-t_elapse )

        if (present(dest_tree_id)) then
          ! we have to synchronize lgt data since we were updating it locally
          call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

          call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
          lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
        endif

        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id1),&
        lgt_n(tree_id1), lgt_sortednumlist(:,:,tree_id1), hvy_active(:,tree_id1), hvy_n(tree_id1) )

        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id1), hvy_n(tree_id1) )

        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id2),&
        lgt_n(tree_id2), lgt_sortednumlist(:,:,tree_id2), hvy_active(:,tree_id2), hvy_n(tree_id2) )

        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id2), hvy_n(tree_id2) )
    end subroutine
    !##############################################################

    !##############################################################
    ! This routine sums up all trees in treeid_list and saves the sum
    ! in dest_tree_id
    subroutine sum_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, treeid_list, dest_tree_id, verbosity)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(inout) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: treeid_list(:), dest_tree_id !< number of the tree
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:, :) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
        logical, intent(in),optional      :: verbosity !< if true: additional information of processing
        !-----------------------------------------------------------------
        integer(kind=ik):: i
        logical :: verbose=.false.

        if (present(verbosity)) verbose=verbosity
        if (params%rank == 0 .and. verbose) write(*,'("Adding trees: ",i4)') treeid_list

        call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_neighbor, dest_tree_id, treeid_list(1))

        do i = 1, size(treeid_list)
          call add_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, treeid_list(i), dest_tree_id)
        enddo

    end subroutine
    !##############################################################

    !##############################################################
    !> This routine computes the average of all trees in treeid_list and
    !> saves the result in dest_tree_id
    !> result = sum(trees) / number_trees
    !> to construct a tree_id list use for example: tree_id_list= (/(id, id= 1, N)/)
    subroutine average_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, treeid_list, dest_tree_id, verbosity)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(inout) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of active trees in forest
        integer(kind=ik), intent(in)      :: treeid_list(:) !< List of tree ids you want to average
        integer(kind=ik), intent(in)      :: dest_tree_id !< tree id of the averaged tree
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:, :) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
        logical, intent(in),optional      :: verbosity !< if true: additional information of processing
        !-----------------------------------------------------------------
        integer(kind=ik):: N
        logical :: verbose=.false.

        if (present(verbosity)) verbose=verbosity
        if (params%rank == 0 .and. verbose) write(*,'("Averaging trees: ",i4)') treeid_list

        N = size(treeid_list)

        ! first sum up all trees from the given list
        call sum_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, treeid_list, dest_tree_id, verbosity)

        ! devide the sum by the numbers of elements
        call multiply_tree_with_scalar(params, hvy_block, hvy_active, hvy_n, dest_tree_id, &
        1.0_rk/N, verbosity)


    end subroutine
    !##############################################################


    !##############################################################
    subroutine add_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, dest_tree_id, verbosity)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(inout) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:, :) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
        logical, intent(in),optional      :: verbosity !< if true: additional information of processing
        integer(kind=ik), intent(in), optional:: dest_tree_id !< if specified result of addition will be saved here
                                                                !< otherwise result will overwrite tree_id1
        !-----------------------------------------------------------------
        logical :: verbose=.false.

        if (present(verbosity)) verbose=verbosity
        if (params%rank == 0 .and. verbose) write(*,'("Adding trees: ",i4,",",i4)') tree_id1, tree_id2

        if(present(dest_tree_id))then
          call tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
          hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2,"+",dest_tree_id)
        else
          call tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
          hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2,"+")
        endif
    end subroutine
    !########################################################### ###



!##############################################################
! This function returns an scalar computed from the L2 scalar
! prodcut for 2 different trees
!  <f(x),g(x)> = int f(x) * g(x) dx
function scalar_product_two_trees( params, tree_n, &
                       lgt_block,  lgt_active, lgt_n, lgt_sortednumlist, &
                       hvy_block, hvy_neighbor, hvy_active, hvy_n, hvy_tmp ,&
                       tree_id1, tree_id2, buffer_tree_id)  result(sprod)
    implicit none

    !-----------------------------------------------------------------
    !> user defined parameter structure
    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
    !> for the multiplication we need an additional tree as a buffer. If
    !> no tree is passed as additional argument we will create a new tree
    !> and delete it afterwards.
    !> If buffer_tree_id is passed, we use it as a buffer and do not have
    !> to create or delete any additional tree
    integer(kind=ik) , optional, intent(in) :: buffer_tree_id
    !> light data array
    integer(kind=ik),  intent(inout):: lgt_block(:, :)
    !> size of active lists
    integer(kind=ik),  intent(inout):: lgt_n(:), tree_n, hvy_n(:)
    !> heavy data array - block data
    real(kind=rk),  intent(inout)   :: hvy_block(:, :, :, :, :)
    !> heavy temp data: needed in blockxfer which is called in add_two_trees
    real(kind=rk),   intent(inout)  :: hvy_tmp(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), intent(inout) :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout) :: lgt_active(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout) :: hvy_active(:, :)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)       :: lgt_sortednumlist(:,:,:)
    !---------------------------------------------------------------
    integer(kind=ik) :: free_tree_id, Jmax, Bs(3), g, &
                        N_snapshots, N, k, lgt_id, hvy_id, rank, i, mpierr
    real(kind=rk) :: sprod, Volume, t_elapse, t_inc(2)
    real(kind=rk) :: x0(3), dx(3)

    if ( present(buffer_tree_id)) then
      free_tree_id = buffer_tree_id
    else
      free_tree_id = tree_n + 1
    endif

    N = params%number_blocks
    rank = params%rank
    Jmax = params%max_treelevel
    g = params%n_ghosts
    Bs= params%Bs
    Volume = product(params%domain_size(1:params%dim))
    if (params%forest_size < N_snapshots + 1 ) call abort(030319,"Forest is to small")

    !----------------------------------------------
    ! sprod = <X_i, X_j>
    t_elapse = MPI_wtime()

    !---------------------------------------------------
    ! multiply tree_id2 * free_tree_id -> free_tree_id
    !---------------------------------------------------
    t_inc(1) = MPI_wtime()
    call multiply_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, free_tree_id)
    !call write_tree_field("field1_copy_mult.h5", params, lgt_block, lgt_active, hvy_block, &
    !lgt_n, hvy_n, hvy_active, 1, free_tree_id )
    !call write_tree_field("field2_copy_mult.h5", params, lgt_block, lgt_active, hvy_block, &
    !lgt_n, hvy_n, hvy_active, 1, tree_id2 )
    !call abort(134)


    t_inc(1) = MPI_wtime()-t_inc(2)
    !---------------------------------------------------
    ! adapt the mesh before summing it up
    !---------------------------------------------------
    if ( params%adapt_mesh ) then
            !call adapt_tree_mesh( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
            !lgt_n, lgt_sortednumlist, hvy_active, hvy_n, params%coarsening_indicator, hvy_tmp, &
            !free_tree_id, tree_n )
    endif

    !----------------------------------------------------
    ! sum over all elements of the tree with free_tree_id
    !-----------------------------------------------------
    t_inc(2) = MPI_wtime()
    sprod = 0.0_rk
    do k = 1, hvy_n(free_tree_id)
      hvy_id = hvy_active(k, free_tree_id)
      call hvy_id_to_lgt_id(lgt_id, hvy_id, rank, N )
      ! calculate the lattice spacing.
      ! It is needed to perform the L2 inner product
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
      if ( params%dim == 3 ) then
        sprod = sprod + dx(1)*dx(2)*dx(3)* sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, :, hvy_id))
      else
        sprod = sprod + dx(1)*dx(2)*sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, :, hvy_id))
      endif
    end do
    t_inc(2) = MPI_wtime()-t_inc(2)
    t_elapse = MPI_WTIME() - t_elapse
    !----------------------------------------------------
    ! sum over all Procs
    !----------------------------------------------------
    call MPI_ALLREDUCE(MPI_IN_PLACE, sprod, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,WABBIT_COMM, mpierr)

    if (.not. present(buffer_tree_id)) then
       call delete_tree(params, lgt_block, lgt_active, lgt_n, free_tree_id)
    endif
end function
!##############################################################








    !##############################################################
    subroutine multiply_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, dest_tree_id, verbosity)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in)    :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:) !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:,:) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
        logical, intent(in),optional      :: verbosity !< if true: additional information of processing
        integer(kind=ik), intent(in),optional :: dest_tree_id !< number of the tree
        !-----------------------------------------------------------------
        logical :: verbose=.false.

        if (present(verbosity)) verbose=verbosity
        if (params%rank == 0 .and. verbose ) write(*,'("Multiply trees: ",i4,",",i4)') tree_id1, tree_id2

        if(present(dest_tree_id))then
          call tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
          hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2,"*",dest_tree_id)
        else
          call tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
          hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, "*")
        endif

    end subroutine
    !##############################################################



    !##############################################################
    !> This function compares the tree structure and block distribution of two trees.
    !> If treecodes are identical but the processor ranks not, then we redistribute
    !> one of the corresponding blocks, such that blocks with same treecode are on the
    !> same rank.
    subroutine same_block_distribution(params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, tree_n, tree_id1, tree_id2)

        implicit none
        !-----------------------------------------------------------------
        type (type_params), intent(in) :: params   !< params structure
        integer(kind=ik), intent(inout)   :: hvy_n(:)    !< number of active heavy blocks
        integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
        integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
        integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
        integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
        real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
        integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:, :) !< active lists
        integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:,:)
        !-----------------------------------------------------------------
        integer(kind=ik)    :: level1, level2, lgt_id1, lgt_id2, N, fsize
        integer(kind=ik)    :: k1, k2, rank1, rank2, level_min, n_comm, Jmax
        integer(kind=tsize) :: treecode1, treecode2
        integer(kind=ik), allocatable, save :: comm_list(:,:)

        N = params%number_blocks
        fsize = params%forest_size
        Jmax = params%max_treelevel
        if (.not.allocated(comm_list)) allocate( comm_list( params%number_procs*N, 3 ) )

        !! Loop over all treecodes of both trees and check if they are identical.
        !! If the treecodes are the same but the blocks are not on the same rank,
        !! we have to send a block.
        !! The comm_list is created to save all of the block transfers.
        n_comm = 0 ! number of communications
        do k1 = 1, lgt_n(tree_id1)
            !--------
            ! TREE 1
            !--------
            lgt_id1  = lgt_active(k1, tree_id1)
            level1   = lgt_block(lgt_id1, Jmax + IDX_MESH_LVL)

            do k2 = 1, lgt_n(tree_id2)
                !--------
                ! TREE 2
                !--------
                lgt_id2  = lgt_active(k2, tree_id2)
                level2   = lgt_block(lgt_id2, Jmax + IDX_MESH_LVL)
                level_min= min(level1,level2)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level_min))
                treecode2= treecode2int(lgt_block(lgt_id2, 1 : level_min))
                if (treecode1 == treecode2 )then
                    ! check treestructure
                    if (level1 .ne. level2) then
                        call abort(270219, "Trees need same treestructure to fetch processor block distribution")
                    endif
                    ! check the block distribution (leaves are on same rank)
                    call lgt_id_to_proc_rank( rank1, lgt_id1, N)
                    call lgt_id_to_proc_rank( rank2, lgt_id2, N)
                    if (rank1 .ne. rank2) then
                        n_comm = n_comm + 1
        !                write(*,*) "===============> not identical", n_comm
                        comm_list(n_comm, 1) = rank1   ! sender mpirank
                        comm_list(n_comm, 2) = rank2   ! receiver mpirank
                        comm_list(n_comm, 3) = lgt_id1 ! block lgt_id to send
                    end if
                end if
            end do ! loop over tree2
        end do ! loop over tree1
        !if (params%rank == 0) write(*,'("nr blocks redistributed: ",i9)') n_comm
        ! Transfer the blocks
        call block_xfer( params, comm_list, n_comm, lgt_block, hvy_block )
        ! after block transfer we have to create new lists
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
        lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
    end subroutine
    !##############################################################





    !##############################################################
    !> Counts the number of hvy_blocks owned by the tree with the
    !> specified tree_id
    function count_tree_hvy_n( params, lgt_block, hvy_active, hvy_n, tree_id) &
        result(tree_hvy_n)

        implicit none

        !-----------------------------------------------------------------
        integer(kind=ik),  intent(in)  :: hvy_active(:) !< list of active blocks on the current rank
        integer(kind=ik), intent(in)  :: lgt_block(:,:) !< light data array
        integer(kind=ik),  intent(in)  :: tree_id, hvy_n !< index of the tree and number of active hvy_blocks
        type (type_params), intent(in) :: params !< params structure of wabbit
        !-----------------------------------------------------------------
        integer (kind=ik) :: k, rank, lgt_id, tree_hvy_n, idxtree, N

        N = params%number_blocks
        rank = params%rank
        idxtree = params%max_treelevel + IDX_TREE_ID
        tree_hvy_n = 0

        do k = 1, hvy_n
            ! calculate lgtblock id corresponding to hvy id
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )
            if (lgt_block(lgt_id, idxtree) == tree_id) then
                tree_hvy_n = tree_hvy_n + 1
            end if
        end do
    end function count_tree_hvy_n
    !##############################################################






end module module_forest
