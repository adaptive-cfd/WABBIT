! ********************************************************************************************
!> \brief Apply mesh coarsening: Merge tagged blocks into new, coarser blocks
! ********************************************************************************************
subroutine executeCoarsening_tree( params, hvy_block, tree_ID)
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
    n_xfer = 0  ! transfer counter
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
subroutine executeCoarsening_WD_tree( params, hvy_block, tree_ID, mark_TMP_flag, no_deletion )
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params                       !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !< heavy data array - block data in spaghetti WD form
    integer(kind=ik), intent(in)        :: tree_ID                      !< tree_id to be coarsened
    logical, intent(in), optional       :: mark_TMP_flag                !< Set refinement of completeness to 0 or temporary flag
    !> do not delete but only create mother blocks and send data
    logical, intent(in), optional       :: no_deletion

    ! loop variables
    integer(kind=ik)                    :: k, Jmax, N, j, rank
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: lgt_daughters(1:8), rank_daughters(1:8)
    integer(kind=tsize)                 :: treecode
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: lgt_ID_b, hvy_ID_b, level_b, rank_b
    integer(kind=ik)                    :: data_rank, n_xfer, lgt_merge_id, digit_merge, hvy_ID
    integer(kind=ik)                    :: nx, ny, nz, nc
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: wc
    logical                             :: markTMPflag
    logical                             :: noDeletion
    real(kind=rk)                       :: t0
    integer(kind=ik), parameter :: TMP_STATUS = 17


    integer(kind=ik)  :: iy

    markTMPflag = .false.
    if (present(mark_TMP_flag)) markTMPflag = mark_TMP_flag
    noDeletion = .false.
    if (present(no_deletion)) noDeletion = no_deletion

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

    !---------------------------------------------------------------------------
    ! transfer all blocks that want to coarsen to Mallat format
    ! this is necessary so that move_mallat_patch_block can move the correct data
    ! elsewise we would need different logic to create the mother block out of one of the daughters
    !---------------------------------------------------------------------------
    t0 = MPI_Wtime()
    do k = 1, hvy_n(tree_ID)
        hvy_ID_b = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID_b, hvy_ID_b, rank, params%number_blocks)
        level_b = lgt_block( lgt_ID_b, IDX_MESH_LVL )

        if ( lgt_block(lgt_ID_b, IDX_TC_1 ) >= 0 .and. lgt_block(lgt_ID_b, IDX_REFINE_STS) == -1) then
            ! This block will be coarsened and its data needs to be transferred to Mallat for correct copying
            call spaghetti2Mallat_block(params, hvy_block(:,:,:,1:nc,hvy_ID_b), wc(:,:,:,1:nc))
            hvy_block(:,:,:,1:nc,hvy_ID_b) = wc(:,:,:,1:nc)
        endif
    enddo
    call toc( "executeCoarsening (spaghetti2Mallat)", 160, MPI_Wtime()-t0 )

    t0 = MPI_Wtime()
    ! setup correct patches and get indices, used for move_mallat_patch_block
    ! this is in theory only needed when we change g but if this is every done I want to avoid nasty bug finding
    call family_setup_patches(params, output_to_file=.false.)
    ! some tiny buffers depend on the number of components (nc=size(hvy_block,4))
    ! make sure they have the right size
    call xfer_ensure_correct_buffer_size(params, hvy_block)
    call toc( "executeCoarsening (setup buffer)", 161, MPI_Wtime()-t0 )

    !---------------------------------------------------------------------------
    ! create new empty blocks on rank with most daughters if it does not exist
    !---------------------------------------------------------------------------
    t0 = MPI_Wtime()
    do k = 1, hvy_n(tree_ID)
        ! Check if the block will be coarsened.
        ! FIRST condition: only work on light data, if block is active.
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1.

        hvy_ID_b = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID_b, hvy_ID_b, rank, params%number_blocks)
        level_b = lgt_block( lgt_ID_b, IDX_MESH_LVL )

        ! check if block is active: TC > 0 and block wants to be refined
        ! performance: don't construct tc and check only first int
        if ( lgt_block(lgt_ID_b, IDX_TC_1 ) >= 0 .and. lgt_block(lgt_ID_b, IDX_REFINE_STS) == -1) then
            ! If this block already has a mother, we do not have to create it once again
            if (hvy_family(hvy_ID_b, 1) == -1) then
                ! Get all sisters
                lgt_daughters(1:N) = hvy_family(hvy_ID_b, 2:1+2**params%dim)
    
                ! figure out on which rank the sisters lie
                do j = 1, N
                    call lgt2proc( rank_daughters(j), lgt_daughters(j), params%number_blocks )
                enddo
    
                ! The merging will be done on the mpirank which holds the most of the sister blocks
                data_rank = most_common_element( rank_daughters(1:N) )
                ! data_rank = rank_daughters(N)  ! experimental but deterministic mother choice
    
                ! construct new mother block if on my rank and create light data entry for the new block
                if (data_rank == rank) then
                    if (.not. noDeletion) then  ! create mother block from one of the sister blocks and move patch directly
                        ! get lgtID as first daughter block on correct rank and move values in-place
                        lgt_merge_id = -1
                        do j = 1, N
                            if (rank_daughters(j) == rank .and. lgt_merge_id == -1) then
                                lgt_merge_id = lgt_daughters(j)
                                call lgt2hvy(hvy_ID, lgt_merge_id, rank, params%number_blocks)
                                ! we need to compute the last digit to in-place move the patch correctly
                                treecode = get_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2 ))
                                digit_merge = tc_get_digit_at_level_b( treecode, params%dim, level_b, params%Jmax)
                                call move_mallat_patch_block(params, hvy_block, hvy_ID, digit_merge)
                            endif
                        enddo
                    else  ! create mother block as a new block
                        call get_free_local_light_id(params, data_rank, lgt_merge_id, message="executeCoarsening_WD")
                        treecode = get_tc(lgt_block( lgt_daughters(1), IDX_TC_1:IDX_TC_2 ))
                    endif

                    ! change meta_data of mother block
                    lgt_block( lgt_merge_id, : ) = -1
                    call set_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2), tc_clear_until_level_b(treecode, &
                        dim=params%dim, level=level_b-1, max_level=params%Jmax))
                    lgt_block( lgt_merge_id, IDX_MESH_LVL ) = level_b-1
                    if (markTMPflag) then
                        lgt_block( lgt_merge_id, IDX_REFINE_STS ) = REF_TMP_UNTREATED
                    else
                        lgt_block( lgt_merge_id, IDX_REFINE_STS ) = 0
                    endif
                    lgt_block( lgt_merge_id, IDX_TREE_ID ) = tree_ID

                    ! update sisters on my rank that they have found their mother and can be skipped
                    do j = 1, N
                        if (rank_daughters(j) == rank .and. lgt_daughters(j) /= lgt_merge_id) then
                            call lgt2hvy(hvy_ID_b, lgt_daughters(j), rank, params%number_blocks)
                            hvy_family(hvy_ID_b, 1) = lgt_merge_id
                        endif
                    enddo

                    ! write(*, '("Rank ", i0, " created a new block: ", i0, " with TC ", b64.64)') rank, lgt_merge_id, tc_clear_until_level_b(treecode, &
                    ! dim=params%dim, level=level-1, max_level=params%Jmax)
                endif
            endif
        endif
    enddo
    call toc( "executeCoarsening (create mothers)", 162, MPI_Wtime()-t0 )

    ! ToDo: Set mother-daughter relations inside loop so that I can skip the syncing as it is expensive
    t0 = MPI_Wtime()
    ! the active lists are outdated, so lets resynch
    call synchronize_lgt_data( params, refinement_status_only=.false.)
    ! update metadata for overlapping grid - temporarily we have a non-unique grid where mothers and daughters coexist
    call updateMetadata_tree(params, tree_ID, search_overlapping=.true., update_neighbors=.false.)
    call toc( "executeCoarsening (sync lgt + updateMetadata)", 163, MPI_Wtime()-t0 )


    ! actual xfer, this works on all blocks that have a mother / daughter and ref -1
    n_xfer = 0  ! transfer counter
    t0 = MPI_Wtime()
    call prepare_update_family_metadata(params, tree_ID, n_xfer, sync_case="D2M_ref", ncomponents=size(hvy_block, 4), s_val=-1)
    call xfer_block_data(params, hvy_block, tree_ID, n_xfer)
    call toc( "executeCoarsening (xfer_block_data)", 164, MPI_Wtime()-t0 )

    ! now the daughter blocks are to be deleted or marked that they are completed
    t0 = MPI_Wtime()
    do k = 1, lgt_n(tree_ID)
        lgt_ID_b = lgt_active(k, tree_ID)
        level_b = lgt_block( lgt_ID_b, IDX_MESH_LVL )
        call lgt2proc( rank_daughters(1), lgt_ID_b, params%number_blocks )
        if ( lgt_block(lgt_ID_b, IDX_REFINE_STS) == -1) then
            if (.not. noDeletion) then
                ! delete daughter blocks
                lgt_block(lgt_ID_b, :) = -1
            elseif (rank_daughters(1) == rank) then
                call lgt2hvy(hvy_ID_b, lgt_ID_b, rank, params%number_blocks)
                ! block has to be retransformed into spaghetti form
                call Mallat2Spaghetti_block(params, hvy_block(:,:,:,1:nc,hvy_ID_b), wc(:,:,:,1:nc))
                hvy_block(:,:,:,1:nc,hvy_ID_b) = wc(:,:,:,1:nc)
            endif
            ! mark daughter blocks as completed
            lgt_block(lgt_ID_b, IDX_REFINE_STS) = 0
        endif
    enddo
    call toc( "executeCoarsening (delete daughters or Mallat2Spaghetti)", 165, MPI_Wtime()-t0 )


end subroutine executeCoarsening_WD_tree


!> For CVS decomposition we need to update mothers from decomposed daughter values, this does exactly that
subroutine sync_D2M(params, hvy_block, tree_ID, sync_case, s_val)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params                       !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !< heavy data array - block data in spaghetti WD form
    integer(kind=ik), intent(in)        :: tree_ID                      !< tree_id to be coarsened
    character(len=*)                    :: sync_case                    !< String representing which kind of syncing we want to do
    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, dependend on sync case
    integer(kind=ik), intent(in), optional  :: s_val

    ! loop variables
    integer(kind=ik)                    :: k, sync_case_id, nc, data_rank, n_xfer, lgt_ID, hvy_ID, level_me, ref_me
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: wc

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    nc = size(hvy_block,4)
    if (allocated(wc)) then
        if (size(wc, 4) < nc) deallocate(wc)
    endif
    if (.not. allocated(wc)) allocate(wc(1:size(hvy_block,1), 1:size(hvy_block,2), 1:size(hvy_block,3), 1:nc) )

    select case(sync_case)
    case("level")
        sync_case_id = 1
    case("ref")
        sync_case_id = 2
    case default
        call abort(240809, "My language does not have so many cases, so I have no idea what you want from me.")
    end select

    !---------------------------------------------------------------------------
    ! transfer all daughters to mallat format, so that sending is smooth
    !---------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
        ref_me   = lgt_block( lgt_ID, IDX_REFINE_STS)
        ! skip some blocks
        if (sync_case_ID == 1) then
            if (level_me /= s_val) cycle
        elseif (sync_case_ID == 2) then
            if (ref_me /= s_val) cycle
        endif

        ! check if block exists and actually has a mother
        if ( lgt_block(lgt_ID, IDX_TC_1 ) >= 0 .and. .not. block_is_root(params, hvy_ID)) then
            ! This block will send or receive and its data needs to be transferred to Mallat for correct copying
            call spaghetti2Mallat_block(params, hvy_block(:,:,:,1:nc,hvy_ID), wc(:,:,:,1:nc))
            hvy_block(:,:,:,1:nc,hvy_ID) = wc(:,:,:,1:nc)
        endif
    enddo

    ! setup correct patches and get indices, used for move_mallat_patch_block
    ! this is in theory only needed when we change g but if this is every done I want to avoid nasty bug finding
    call family_setup_patches(params, output_to_file=.false.)
    ! some tiny buffers depend on the number of components (nc=size(hvy_block,4))
    ! make sure they have the right size
    call xfer_ensure_correct_buffer_size(params, hvy_block)

    ! actual xfer, this works on all blocks that are on this level and have a daughter
    n_xfer = 0  ! transfer counter
    call prepare_update_family_metadata(params, tree_ID, n_xfer, sync_case="D2M_" // sync_case, ncomponents=nc, s_val=s_val)
    call xfer_block_data(params, hvy_block, tree_ID, n_xfer, verbose_check=.true.)

    ! now the daughter blocks need to be retransformed to spaghetti
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
        ref_me   = lgt_block( lgt_ID, IDX_REFINE_STS)
        ! skip some blocks
        if (sync_case_ID == 1) then
            if (level_me /= s_val) cycle
        elseif (sync_case_ID == 2) then
            if (ref_me /= s_val) cycle
        endif

        ! check if block exists and actually has a mother
        if ( lgt_block(lgt_ID, IDX_TC_1 ) >= 0 .and. .not. block_is_root(params, hvy_ID)) then
            ! block has to be retransformed into spaghetti form
            call Mallat2Spaghetti_block(params, hvy_block(:,:,:,1:nc,hvy_ID), wc(:,:,:,1:nc))
            hvy_block(:,:,:,1:nc,hvy_ID) = wc(:,:,:,1:nc)
        endif
    enddo

end subroutine sync_D2M


!> For CVS reconstruction we need to update daughters from reconstructed mothers, this does exactly that
subroutine sync_M2D(params, hvy_block, tree_ID, sync_case, s_val)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params                       !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !< heavy data array - block data in spaghetti WD form
    integer(kind=ik), intent(in)        :: tree_ID                      !< tree_id to be coarsened
    character(len=*)                    :: sync_case                    !< String representing which kind of syncing we want to do
    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, dependend on sync case
    integer(kind=ik), intent(in), optional  :: s_val

    ! loop variables
    integer(kind=ik)                    :: k, sync_case_id, n_xfer, lgt_ID, hvy_ID, lgt_ID_m, level_m, ref_m, nc
    logical                             :: change_form
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: wc

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    nc = size(hvy_block,4)
    if (allocated(wc)) then
        if (size(wc, 4) < nc) deallocate(wc)
    endif
    if (.not. allocated(wc)) allocate(wc(1:size(hvy_block,1), 1:size(hvy_block,2), 1:size(hvy_block,3), 1:nc) )

    select case(sync_case)
    case("level")
        sync_case_id = 1
    case("ref")
        sync_case_id = 2
    case default
        call abort(240809, "My language does not have so many cases, so I have no idea what you want from me.")
    end select

    !---------------------------------------------------------------------------
    ! daughters will actually receive all values in mallat form, we need to change them to mallat to not overwrite WC
    ! block needs to check if it will receive values and if this is the case - change form
    !---------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        change_form = .false.
        
        ! check if mother exists and wants to send this block data
        if (.not. block_is_root(params, hvy_ID)) then
            lgt_ID_m = hvy_family(hvy_ID, 1)
            level_m = lgt_block( lgt_ID_m, IDX_MESH_LVL )
            ref_m = lgt_block( lgt_ID_m, IDX_REFINE_STS)
            ! for some calls we don't want to work on whole tree so skip some blocks
            if (sync_case_ID == 1) then
                if (level_m == s_val) change_form = .true.
            elseif (sync_case_ID == 2) then
                if (ref_m == s_val) change_form = .true.
            endif
        endif

        if (change_form) then
            ! block has to be transformed into Mallat form
            call Spaghetti2Mallat_block(params, hvy_block(:,:,:,1:nc,hvy_ID), wc(:,:,:,1:nc))
            hvy_block(:,:,:,1:nc,hvy_ID) = wc(:,:,:,1:nc)
        endif
    enddo

    ! setup correct patches and get indices, used for move_mallat_patch_block
    ! this is in theory only needed when we change g but if this is every done I want to avoid nasty bug finding
    call family_setup_patches(params, output_to_file=.false.)
    ! some tiny buffers depend on the number of components (nc=size(hvy_block,4))
    ! make sure they have the right size
    call xfer_ensure_correct_buffer_size(params, hvy_block)

    ! actual xfer, this works on all blocks that are on this level and have a daughter
    n_xfer = 0  ! transfer counter
    call prepare_update_family_metadata(params, tree_ID, n_xfer, sync_case="M2D_" // sync_case, ncomponents=nc, s_val=s_val)
    call xfer_block_data(params, hvy_block, tree_ID, n_xfer, verbose_check=.true.)

    !---------------------------------------------------------------------------
    ! daughters will actually receive all values in mallat form, we need to recopy them to spaghetti to continue
    ! block needs to check if it received values and if this is the case - change form
    !---------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        change_form = .false.
        
        ! check if mother exists and wants to send this block data
        if (.not. block_is_root(params, hvy_ID)) then
            lgt_ID_m = hvy_family(hvy_ID, 1)
            level_m = lgt_block( lgt_ID_m, IDX_MESH_LVL )
            ref_m = lgt_block( lgt_ID_m, IDX_REFINE_STS)
            ! for some calls we don't want to work on whole tree so skip some blocks
            if (sync_case_ID == 1) then
                if (level_m == s_val) change_form = .true.
            elseif (sync_case_ID == 2) then
                if (ref_m == s_val) change_form = .true.
            endif
        endif

        if (change_form) then
            ! block has to be retransformed into spaghetti form
            call Mallat2Spaghetti_block(params, hvy_block(:,:,:,1:nc,hvy_ID), wc(:,:,:,1:nc))
            hvy_block(:,:,:,1:nc,hvy_ID) = wc(:,:,:,1:nc)
        endif
    enddo

end subroutine sync_M2D