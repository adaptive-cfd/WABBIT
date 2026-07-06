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


!> For CVS decomposition we need to update mothers from decomposed daughter values, this does exactly that
subroutine sync_D2M(params, hvy_block, tree_ID, sync_case, s_level, s_ref, sync_debug_name)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params                       !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !< heavy data array - block data in spaghetti WD form
    integer(kind=ik), intent(in)        :: tree_ID                      !< tree_id to be coarsened
    character(len=*)                    :: sync_case                    !< String representing which kind of syncing we want to do
    !> Additional value to be considered for syncing logic, can be level or refinement status from which should be synced, dependend on sync case
    integer(kind=ik), optional, intent(in) :: s_level, s_ref
    character(len=*), optional, intent(in) :: sync_debug_name       !< name to be used in debug output files

    ! loop variables
    integer(kind=ik)                    :: k, sync_case_id, nc, ic, data_rank, n_xfer, lgt_ID, hvy_ID, level_me, ref_me
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: wc

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    nc = size(hvy_block,4)
    if (.not. allocated(wc)) allocate(wc(1:size(hvy_block,1), 1:size(hvy_block,2), 1:size(hvy_block,3), 1:1) )

    select case(sync_case)
    case("level")
        sync_case_id = 1
    case("ref")
        sync_case_id = 2
    case("level_ref")
        sync_case_id = 3
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
            if (level_me /= s_level) cycle
        elseif (sync_case_ID == 2) then
            if (ref_me /= s_ref) cycle
        elseif (sync_case_ID == 3) then
            if (level_me /= s_level .or. ref_me /= s_ref) cycle
        endif

        ! check if block exists and actually has a mother
        if ( lgt_block(lgt_ID, IDX_TC_1 ) >= 0 .and. .not. block_is_root(params, hvy_ID)) then
            ! This block will send or receive and its data needs to be transferred to Mallat for correct copying
            do ic = 1, nc
                call spaghetti2Mallat_block(params, hvy_block(:,:,:,ic:ic,hvy_ID), wc(:,:,:,1:1))
                hvy_block(:,:,:,ic,hvy_ID) = wc(:,:,:,1)
            enddo
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
    call prepare_update_family_metadata(params, tree_ID, n_xfer, sync_case="D2M_" // sync_case, ncomponents=nc, s_level=s_level, s_ref=s_ref, sync_debug_name=sync_debug_name)
    call xfer_block_data(params, hvy_block, tree_ID, n_xfer, verbose_check=.true.)

    ! now the daughter blocks need to be retransformed to spaghetti
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
        ref_me   = lgt_block( lgt_ID, IDX_REFINE_STS)
        ! skip some blocks
        if (sync_case_ID == 1) then
            if (level_me /= s_level) cycle
        elseif (sync_case_ID == 2) then
            if (ref_me /= s_ref) cycle
        elseif (sync_case_ID == 3) then
            if (level_me /= s_level .or. ref_me /= s_ref) cycle
        endif

        ! check if block exists and actually has a mother
        if ( lgt_block(lgt_ID, IDX_TC_1 ) >= 0 .and. .not. block_is_root(params, hvy_ID)) then
            ! block has to be retransformed into spaghetti form
            do ic = 1, nc
                call Mallat2Spaghetti_block(params, hvy_block(:,:,:,ic:ic,hvy_ID), wc(:,:,:,1:1))
                hvy_block(:,:,:,ic,hvy_ID) = wc(:,:,:,1)
            enddo
        endif
    enddo

end subroutine sync_D2M


!> For CVS reconstruction we need to update daughters from reconstructed mothers, this does exactly that
subroutine sync_M2D(params, hvy_block, tree_ID, sync_case, s_level, s_ref, sync_debug_name)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params                       !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !< heavy data array - block data in spaghetti WD form
    integer(kind=ik), intent(in)        :: tree_ID                      !< tree_id to be coarsened
    character(len=*)                    :: sync_case                    !< String representing which kind of syncing we want to do
    !> Additional value to be considered for syncing logic, can be level or refinement status from which should be synced, dependend on sync case
    integer(kind=ik), optional, intent(in) :: s_level, s_ref
    character(len=*), optional, intent(in) :: sync_debug_name       !< name to be used in debug output files

    ! loop variables
    integer(kind=ik)                    :: k, sync_case_id, n_xfer, lgt_ID, hvy_ID, lgt_ID_m, level_m, ref_m, nc, ic
    logical                             :: change_form
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: wc

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    nc = size(hvy_block,4)
    if (.not. allocated(wc)) allocate(wc(1:size(hvy_block,1), 1:size(hvy_block,2), 1:size(hvy_block,3), 1:1) )

    select case(sync_case)
    case("level")
        sync_case_id = 1
    case("ref")
        sync_case_id = 2
    case("level_ref")
        sync_case_id = 3
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
                if (level_m == s_level) change_form = .true.
            elseif (sync_case_ID == 2) then
                if (ref_m == s_ref) change_form = .true.
            elseif (sync_case_ID == 3) then
                if (level_m == s_level .and. ref_m == s_ref) change_form = .true.
            endif
        endif

        if (change_form) then
            ! block has to be transformed into Mallat form
            do ic = 1, nc
                call Spaghetti2Mallat_block(params, hvy_block(:,:,:,ic:ic,hvy_ID), wc(:,:,:,1:1))
                hvy_block(:,:,:,ic,hvy_ID) = wc(:,:,:,1)
            enddo
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
    call prepare_update_family_metadata(params, tree_ID, n_xfer, sync_case="M2D_" // sync_case, ncomponents=nc, s_level=s_level, s_ref=s_ref, sync_debug_name=sync_debug_name)
    call xfer_block_data(params, hvy_block, tree_ID, n_xfer)

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
                if (level_m == s_level) change_form = .true.
            elseif (sync_case_ID == 2) then
                if (ref_m == s_ref) change_form = .true.
            elseif (sync_case_ID == 3) then
                if (level_m == s_level .and. ref_m == s_ref) change_form = .true.
            endif
        endif

        if (change_form) then
            ! block has to be retransformed into spaghetti form
            do ic = 1, nc
                call Mallat2Spaghetti_block(params, hvy_block(:,:,:,ic:ic,hvy_ID), wc(:,:,:,1:1))
                hvy_block(:,:,:,ic,hvy_ID) = wc(:,:,:,1)
            enddo
        endif
    enddo

end subroutine sync_M2D