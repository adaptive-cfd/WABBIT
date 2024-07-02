! ********************************************************************************************
!> \brief Apply mesh coarsening: Merge tagged blocks into new, coarser blocks
! ********************************************************************************************
subroutine executeCoarsening_tree( params, hvy_block, tree_ID )
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    ! loop variables
    integer(kind=ik)                    :: k, Jmax, N, j
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: light_ids(1:8), mpirank_owners(1:8)
    integer(kind=ik), allocatable, save :: xfer_list(:,:)
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: data_rank, n_xfer, ierr, lgtID, procID, hvyID

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    Jmax = params%Jmax
    ! number of blocks to merge, 4 or 8
    N = 2**params%dim
    ! at worst every block is on a different rank
    if (.not. allocated(xfer_list)) allocate(xfer_list(3, size(lgt_block,1)))

    ! transfer counter
    n_xfer = 0


    !---------------------------------------------------------------------------
    ! first, prepare for xfer (gather information: which blocks are sent where)
    ! COLLECTIVE OPERATION
    !---------------------------------------------------------------------------
    do k = 1, lgt_n(tree_ID)
        ! Check if the block will be coarsened.
        !
        ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
        ! the active list. NOTE HERE: due to previous iterations, some light data id are already
        ! marked for xfer, (and thus given the -7 status) but they are still in active block list
        ! so: check again if this block is REALLY active (the active list is updated later)
        !
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
        ! ensureGradedness_tree removes the -1 flag if not all sister blocks share it
        lgtID = lgt_active(k, tree_ID)
        call lgt2proc(procID, lgtID, params%number_blocks)
        call lgt2hvy(hvyID, lgtID, procID, params%number_blocks)
        

        ! check if block is active: TC > 0 and block wants to be refined
        ! performance: don't construct tc and check only first int
        if ( lgt_block(lgtID, IDX_TC_1 ) >= 0 .and. lgt_block(lgtID, IDX_REFINE_STS) == -1) then
            ! find all sisters (including the block in question, so four or eight blocks)
            ! their light IDs are in "light_ids" and ordered by their last treecode-digit
            ! if this block resides on this block we can get its sisters from hvy_family
            if (procID == params%rank) then
                light_ids(1:N) = hvy_family(hvyID, 2:1+2**params%dim)
            else
                call find_sisters( params, lgtID, light_ids(1:N))
            endif

            ! figure out on which rank the sisters lie, xfer them if necessary
            do j = 1, N
                call lgt2proc( mpirank_owners(j), light_ids(j), params%number_blocks )
            enddo

            ! The merging will be done on the mpirank which holds the most of the sister blocks
            data_rank = most_common_element( mpirank_owners(1:N) )

            do j = 1, N
                if (mpirank_owners(j) /= data_rank) then
                    ! MPI xfer required. Add the xfer to the list
                    n_xfer = n_xfer + 1
                    xfer_list(1, n_xfer) = mpirank_owners(j)  ! send from this rank ..
                    xfer_list(2, n_xfer) = data_rank          ! ... to this rank
                    xfer_list(3, n_xfer) = light_ids(j)       ! transfer this block
                endif

                ! don't forget: mark all 4/8 sisters as treated here, in order not to trigger this
                ! loop again: we use the temporary status -7
                lgt_block(light_ids(j), IDX_REFINE_STS) = -7
            enddo
        endif
    enddo

    ! actual xfer
    call block_xfer( params, xfer_list, n_xfer, hvy_block )

    ! the active lists are outdates after the transfer: we need to create
    ! them or findSisters_tree will not be able to do its job
    call createActiveSortedLists_tree( params, tree_ID )

    ! actual merging
    do k = 1, lgt_n(tree_ID)
        ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
        ! the active list. NOTE HERE: due to previous loops, some light data id are already
        ! coarsened, (and thus given the -1) but they are still in active block list
        ! so: check again if this block is REALLY active (the active list is updated later)
        !
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
        ! ensureGradedness_tree removes the -1 flag and sets -7 (temporarily assigned above) flag if not all sister blocks share it
        lgtID = lgt_active(k, tree_ID)
        ! check if block is active: TC > 0 and block wants to be refined
        ! performance: don't construct tc and check only first int
        if ( lgt_block(lgtID, IDX_TC_1 ) >= 0 .and. lgt_block(lgtID, IDX_REFINE_STS) == -7) then
            ! merge the 4/8 blocks into one new block. Merging is done in two steps,
            ! first for light data (which all CPUS do redundantly, so light data is kept synched)
            ! Then only the responsible rank will perform the heavy data merging.
            call find_sisters( params, lgtID, light_ids(1:N))
            ! note the newly merged block has status 0
            call merge_blocks( params, hvy_block, light_ids(1:N) )
        endif
    enddo
end subroutine executeCoarsening_tree



! ********************************************************************************************
!> \brief Apply mesh coarsening with values currently in WD form. \n
!! Merge tagged blocks into new, coarser blocks
! ********************************************************************************************
subroutine executeCoarsening_WD_tree( params, hvy_block, tree_ID, mark_TMP_flag )
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params                       !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !< heavy data array - block data in spaghetti WD form
    integer(kind=ik), intent(in)        :: tree_ID                      !< tree_id to be coarsened
    logical, intent(in), optional       :: mark_TMP_flag                !< Set refinement of completeness to 0 or temporary flag

    ! loop variables
    integer(kind=ik)                    :: k, Jmax, N, j, rank
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: lgt_daughters(1:8), rank_daughters(1:8)
    integer(kind=tsize)                 :: treecode
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: data_rank, n_xfer, ierr, lgtID, hvyID, level_me, lgt_merge_id, digit_merge
    integer(kind=ik)                    :: nx, ny, nz, nc
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: wc
    logical                             :: markTMPflag

    integer(kind=ik)  :: iy

    markTMPflag = .false.
    if (present(mark_TMP_flag)) markTMPflag = mark_TMP_flag

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    nx = size(hvy_block,1)
    ny = size(hvy_block,2)
    nz = size(hvy_block,3)
    nc = size(hvy_block,4)

    Jmax = params%Jmax
    rank = params%rank
    ! number of blocks to merge, 4 or 8
    N = 2**params%dim
    ! at worst every block is on a different rank
    if (allocated(wc)) then
        if (size(wc, 4) < nc) deallocate(wc)
    endif
    if (.not. allocated(wc)) allocate(wc(1:nx, 1:ny, 1:nz, 1:nc) )

    ! transfer counter
    n_xfer = 0

    !---------------------------------------------------------------------------
    ! transfer all blocks that want to coarsen to Mallat format
    ! this is necessary so that move_mallat_patch_block can move the correct data
    ! elsewise we would need different logic to create the mother block out of one of the daughters
    !---------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt(lgtID, hvyID, rank, params%number_blocks)
        level_me = lgt_block( lgtID, IDX_MESH_LVL )

        if ( lgt_block(lgtID, IDX_TC_1 ) >= 0 .and. lgt_block(lgtID, IDX_REFINE_STS) == -1) then
            ! This block will be coarsened and its data needs to be transferred to Mallat for correct copying
            call spaghetti2Mallat_block(params, hvy_block(:,:,:,1:nc,hvyID), wc(:,:,:,1:nc))
            hvy_block(:,:,:,:,hvyID) = wc
        endif
    enddo

    ! setup correct patches and get indices, used for move_mallat_patch_block
    ! this is in theory only needed when we change g but if this is every done I want to avoid nasty bug finding
    call family_setup_patches(params, output_to_file=.false.)
    ! some tiny buffers depend on the number of components (nc=size(hvy_block,4))
    ! make sure they have the right size
    call xfer_ensure_correct_buffer_size(params, hvy_block)

    !---------------------------------------------------------------------------
    ! create new empty blocks on rank with most daughters if it does not exist
    !---------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        ! Check if the block will be coarsened.
        !
        ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
        ! the active list. NOTE HERE: due to previous iterations
        !
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
        ! ensureGradedness_tree removes the -1 flag if not all sister blocks share it
        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt(lgtID, hvyID, rank, params%number_blocks)
        level_me = lgt_block( lgtID, IDX_MESH_LVL )

        ! check if block is active: TC > 0 and block wants to be refined
        ! performance: don't construct tc and check only first int
        if ( lgt_block(lgtID, IDX_TC_1 ) >= 0 .and. lgt_block(lgtID, IDX_REFINE_STS) == -1) then
            ! If this block already has a mother, we do not have to create it once again
            if (hvy_family(hvyID, 1) == -1) then
                ! Get all sisters
                lgt_daughters(1:N) = hvy_family(hvyID, 2:1+2**params%dim)
    
                ! figure out on which rank the sisters lie
                do j = 1, N
                    call lgt2proc( rank_daughters(j), lgt_daughters(j), params%number_blocks )
                enddo
    
                ! The merging will be done on the mpirank which holds the most of the sister blocks
                data_rank = most_common_element( rank_daughters(1:N) )
    
                ! construct new mother block if on my rank and create light data entry for the new block
                if (data_rank == rank) then
                    ! get lgtID as first daughter block on correct rank and move values in-place
                    lgt_merge_id = -1
                    do j = 1, N
                        if (rank_daughters(j) == rank .and. lgt_merge_id == -1) then
                            lgt_merge_id = lgt_daughters(j)
                            call lgt2hvy(hvyID, lgt_merge_id, rank, params%number_blocks)
                            ! we need to compute the last digit to in-place move the patch correctly
                            treecode = get_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2 ))
                            digit_merge = tc_get_digit_at_level_b( treecode, params%dim, level_me, params%Jmax)
                            call move_mallat_patch_block(params, hvy_block, hvyID, digit_merge)
                        endif
                    enddo
                    ! call get_free_local_light_id(params, data_rank, lgt_merge_id, message="executeCoarsening_WD")
                    ! treecode = get_tc(lgt_block( lgt_daughters(1), IDX_TC_1:IDX_TC_2 ))

                    ! change meta_data of mother block
                    lgt_block( lgt_merge_id, : ) = -1
                    call set_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2), tc_clear_until_level_b(treecode, &
                        dim=params%dim, level=level_me-1, max_level=params%Jmax))
                    lgt_block( lgt_merge_id, IDX_MESH_LVL ) = level_me-1
                    if (markTMPflag) then
                        lgt_block( lgt_merge_id, IDX_REFINE_STS ) = REF_TMP_UNTREATED
                    else
                        lgt_block( lgt_merge_id, IDX_REFINE_STS ) = 0
                    endif
                    lgt_block( lgt_merge_id, IDX_TREE_ID ) = tree_ID

                    ! update sisters on my rank that they have found their mother and can be skipped
                    do j = 1, N
                        if (rank_daughters(j) == rank .and. lgt_daughters(j) /= lgt_merge_id) then
                            call lgt2hvy(hvyID, lgt_daughters(j), rank, params%number_blocks)
                            hvy_family(hvyID, 1) = lgt_merge_id
                        endif
                    enddo

                    ! write(*, '("Rank ", i0, " created a new block: ", i0, " with TC ", b64.64)') rank, lgt_merge_id, tc_clear_until_level_b(treecode, &
                    ! dim=params%dim, level=level-1, max_level=params%Jmax)
                endif
            endif
        endif
    enddo

    ! the active lists are outdated, so lets resynch
    call synchronize_lgt_data( params, refinement_status_only=.false.)
    ! update metadata but ignore neighbors - this is important as we temporarily have a non-unique grid
    ! where mothers and daughters coexist
    call updateMetadata_tree(params, tree_ID, update_neighbors=.false.)


    ! actual xfer, this works on all blocks that have a mother / daughter
    call prepare_update_family_metadata(params, tree_ID, n_xfer, size(hvy_block, 4), &
        s_M2C=.true.)
    call xfer_block_data(params, hvy_block, tree_ID, n_xfer)

    ! now the mother refinement flags have to be reset and daughter blocks to be deleted
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        level_me = lgt_block( lgtID, IDX_MESH_LVL )
        ! delete daughter blocks that wanted to refine
        if ( lgt_block(lgtID, IDX_REFINE_STS) == -1) then
            lgt_block(lgtID, :) = -1
            lgt_block(lgtID, IDX_REFINE_STS) = 0
        endif
    enddo

end subroutine executeCoarsening_WD_tree
