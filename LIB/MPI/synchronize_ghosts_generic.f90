!> Wrapper to synch blocks with temporary flag from finer neighbours and same-level neighbors
!> Used before wavelet decomposition
subroutine sync_TMP_from_MF(params, hvy_block, tree_ID, REF_TMP_UNTREATED, hvy_tmp, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), intent(in)   :: REF_TMP_UNTREATED             !< this block has no access to modul_mesh so we need the flag value
    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus         !< Boundary sizes in case we want to send less values

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, sync_case="MF_ref", &
        s_val=REF_TMP_UNTREATED, hvy_tmp=hvy_tmp)

end subroutine sync_TMP_from_MF



!> Wrapper to synch blocks with temporary flag from all neighbors
!> Used before wavelet decomposition
subroutine sync_TMP_from_all(params, hvy_block, tree_ID, REF_TMP_UNTREATED, hvy_tmp, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), intent(in)   :: REF_TMP_UNTREATED             !< this block has no access to modul_mesh so we need the flag value
    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus         !< Boundary sizes in case we want to send less values

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, sync_case="full_ref", &
        s_val=REF_TMP_UNTREATED, hvy_tmp=hvy_tmp)

end subroutine sync_TMP_from_all



!> Wrapper to synch from coarser neighbours and same-level neighbors
!! Used after coarse extension to update SC and WC, coarse neighbours need to be synched from hvy_tmp
subroutine sync_SCWC_from_MC(params, hvy_block, tree_ID, hvy_tmp, g_minus, g_plus, level)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout)   :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus, level

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus))    gminus = g_minus
    if (present(g_plus))      gplus = g_plus

    ! if we passed on the level then sync is level-wise
    if (present(level)) then
        call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, sync_case="MC_level", &
            s_val=level, hvy_tmp=hvy_tmp, verbose_check=.true.)
    else
        call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, sync_case="MC", &
            hvy_tmp=hvy_tmp, verbose_check=.true.)
    endif

end subroutine sync_SCWC_from_MC


!> Wrapper to synch all ghost-point patches
subroutine sync_ghosts_tree(params, hvy_block, tree_ID, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    ! set level to -1 to enable synching between all, set stati to send to all levels
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, sync_case="full")

end subroutine sync_ghosts_tree

!> Wrapper to synch all ghost-point patches, but not using the filter.
!! This is done during the RHS evaluation: the stencil size of the FD discretization is short,
!! and consequently only the first few ghost nodes are used. In the coarseExtension, those are copied
!! anyways (being near the coarse/fine interface, the appropriate filter cannot be applied).
!! Therefore, for reasons of consistency and performance, we can ignore the filter in ghost nodes
!! sync'ing. This corresponds to using non-lifted wavelets.
subroutine sync_ghosts_RHS_tree(params, hvy_block, tree_ID, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    ! set level to -1 to enable synching between all, set stati to send to all levels
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, &
    sync_case="full_RHS", ignore_Filter=.true.)

end subroutine


!> This function deals with ghost-node synching \n
!! It is a generic function with many flags, streamlining all synching process \n
!! In order to avoid confusion wrapper functions should be used everywhere in order to implement
!! specific versions. This also means that parameter changes only have to be changed in the wrappers
subroutine sync_ghosts_generic( params, hvy_block, tree_ID, sync_case, &
    g_minus, g_plus, s_val, hvy_tmp, verbose_check, ignore_Filter)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    character(len=*)          :: sync_case                     !< String representing which kind of syncing we want to do

    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag
    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, dependend on sync case
    integer(kind=ik), intent(in), optional  :: s_val
    logical, intent(in), optional  :: ignore_Filter                 !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus       !< Synch only so many ghost points

    integer(kind=ik)   :: myrank, mpisize, Bs(1:3), buffer_offset
    integer(kind=ik)   :: N, k, neighborhood, Nstages
    integer(kind=ik)   :: recver_rank, recver_hvyID, patch_size
    integer(kind=ik)   :: sender_hvyID, sender_lgtID

    integer(kind=ik) :: ijk(2,3), isend, irecv, count_send_total
    integer(kind=ik) :: bounds_type, istage, inverse, gminus, gplus
    real(kind=rk) :: t0, t1, t2
    character(len=clong) :: toc_statement

    t0 = MPI_wtime()

    if (.not. ghost_nodes_module_ready) then
        ! in order to keep the syntax clean, buffers are module-global and need to be
        ! allocated here.
        call init_ghost_nodes( params )
    endif

    ! if this mpirank has no active blocks, it has nothing to do here.
    if (hvy_n(tree_ID) == 0) return

    gminus  = params%g
    gplus   = params%g
    Bs      = params%Bs
    N       = params%number_blocks
    myrank  = params%rank
    mpisize = params%number_procs
    ! default is three stages:
    !    1. M2M, copy
    !    2. M2C, F2M, decimate, independent from coarser neighbor but needs same-level data
    !    3. M2F, C2M, interpolate, needs same-level and finer neighbor data
    Nstages = 3
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    !-----------------------------------------------------------------------
    ! set up constant arrays
    !-----------------------------------------------------------------------
    ! We frequently need to know the indices of a ghost nodes patch. Thus we save them
    ! once in a module-global array (which is faster than computing it every time with tons
    ! of IF-THEN clauses).
    ! This arrays indices are:
    ! ijkPatches([start,end], [dir], [ineighbor], [lvl_diff], [isendrecv])
    ! As g can be varied (as long as it does not exceed the maximum value params%g), it is set up
    ! each time we sync (at negligibble cost)
    call ghosts_setup_patches(params, gminus=gminus, gplus=gplus, output_to_file=.false.)
    ! some tiny buffers depend on the number of components (nc=size(hvy_block,4))
    ! make sure they have the right size
    call xfer_ensure_correct_buffer_size(params, hvy_block)

#ifdef DEV
    ! for dev check ghosts by wiping if we set all of them

    ! disabled for now
    if (sync_case == "full" .or. sync_case=="full_RHS") call reset_ghost_nodes( params, hvy_block, tree_ID)

#endif

    ! We require two stages: first, we fill all ghost nodes which are simple copy (including restriction),
    ! then in the second stage we can use interpolation and fill the remaining ones.
    do istage = 1, Nstages
        !***************************************************************************
        ! (i) stage initialization
        !***************************************************************************
        ! prepare metadata. This computes from the grid info how much data I recv and send to all other mpiranks.
        ! Also applies logic about what should be synched and saves all metadata unsorted in one array
        ! internal nodes are included in metadata but not counted
        t1 = MPI_wtime()  ! stage duration
        call prepare_ghost_synch_metadata(params, tree_ID, count_send_total, &
            istage, sync_case, ncomponents=size(hvy_block,4), s_val=s_val)
        call toc( "sync ghosts (prepare metadata)", 81, MPI_wtime()-t1 )

        !***************************************************************************
        ! (ii) sending handled by xfer_block_data
        ! If hvy_temp is present then xfer_block_data has to decide from where to grab the data
        !    - with ref in case: we decide after refinement flag present in s_val if we want to use hvy_temp
        !    - elsewise: use hvy_tmp for prediction (used in updating SC from coarser neighbours)
        !***************************************************************************
        t2 = MPI_wtime()
        if (.not. present(hvy_tmp)) then
            call xfer_block_data(params, hvy_block, tree_ID, count_send_total, ignore_Filter=ignore_Filter, &
            verbose_check=verbose_check)
        else
            if (sync_case=="full_ref" .or. sync_case=="MF_ref" .or. sync_case=="MC_ref") then
                call xfer_block_data(params, hvy_block, tree_ID, count_send_total, hvy_tmp=hvy_tmp, &
                REF_FLAG=s_val, ignore_Filter=ignore_Filter, verbose_check=verbose_check)
            else
                call xfer_block_data(params, hvy_block, tree_ID, count_send_total, hvy_tmp=hvy_tmp, &
                ignore_Filter=ignore_Filter, verbose_check=verbose_check)
            endif
        endif
        call toc( "sync ghosts (xfer_block_data)", 82, MPI_wtime()-t2 )
        write(toc_statement, '(A, i0, A)') "sync ghosts (stage ", istage, " TOTAL)"
        call toc( toc_statement, 82+istage, MPI_wtime()-t1 )

    end do ! loop over stages 1,2

    call toc( "sync ghosts (TOTAL)", 80, MPI_wtime()-t0 )

end subroutine sync_ghosts_generic



! This subroutine prepares who sends to whom. This includes:
!    - logic of different synchronization situations
!    - saving of all metadata
!    - computing of buffer sizes for metadata for both sending and receiving
! This is done strictly locally so no MPI needed here
subroutine prepare_ghost_synch_metadata(params, tree_ID, count_send, istage, sync_case, ncomponents, s_Val)

    implicit none

    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(in)    :: tree_ID             !< which tree to study

    integer(kind=ik), intent(in)    :: ncomponents         !< components can vary (for mask for example)
    integer(kind=ik), intent(out)   :: count_send          !< number of ghost patches total to be send, for looping
    !> following are variables that control the logic of where each block sends or receives
    integer(kind=ik), intent(in)    :: istage  !< current stage out of three
    character(len=*)           :: sync_case                     !< String representing which kind of syncing we want to do

    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, dependend on sync case
    integer(kind=ik), intent(in), optional  :: s_val

    ! Following are global data used but defined in module_mpi:
    !    data_recv_counter, data_send_counter
    !    meta_recv_counter, meta_send_counter
    !    meta_send_all (possibly needs renaming after this function)

    integer(kind=ik) :: k_block, myrank, N, neighborhood, ijk(2,3), inverse, ierr, lvl_diff, status, new_size, sync_id
    integer(kind=ik) :: hvyID, lgtID, level, ref
    integer(kind=ik) :: hvyID_n, lgtID_n, level_n, ref_n, rank_n
    logical :: b_send, b_recv

    myrank = params%rank
    N = params%number_blocks

    data_recv_counter(:) = 0
    data_send_counter(:) = 0
    meta_send_counter(:) = 0
    meta_recv_counter(:) = 0

    int_pos(:) = 0
    real_pos(:) = 0

    ! translate sync_case to sync_id so that we avoid string comparisons in loops
    select case(sync_case)
    case ("full")
        sync_id = 1
    case("full_RHS")
        sync_id = 11
    case("full_ref")
        sync_id = 21
    case("full_level")
        sync_id = 31
    case("MC")
        sync_id = 2
    case("MC_ref")
        sync_id = 22
    case("MC_level")
        sync_id = 32
    case("MF")
        sync_id = 3
    case("MF_ref")
        sync_id = 23
    case("MF_level")
        sync_id = 33
    case default
        call abort(240805, "No, we don't trade that here! Please ensure the sync_case is valid.")
    end select

    count_send = 0
    do k_block = 1, hvy_n(tree_ID)
        ! calculate light id
        hvyID = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgtID, hvyID, myrank, N )
        level = lgt_block( lgtID, IDX_MESH_LVL )
        ref = lgt_block( lgtID, IDX_REFINE_STS)

        ! loop over all neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists
            if ( hvy_neighbor( hvyID, neighborhood ) /= -1 ) then
                ! neighbor light data id
                lgtID_n = hvy_neighbor( hvyID, neighborhood )
                ! calculate neighbor rank
                call lgt2proc( rank_n, lgtID_n, N )
                ! neighbor heavy id
                call lgt2hvy( hvyID_n, lgtID_n, rank_n, N )

                ! define level difference: sender - receiver, so +1 means sender on higher level
                ! lvl_diff = -1 : sender to finer recver, interpolation on sender side
                ! lvl_diff =  0 : sender is same level as recver
                ! lvl_diff = +1 : sender to coarser recver, restriction is applied on sender side
                level_n = lgt_block( lgtID_n, IDX_MESH_LVL )
                lvl_diff = level - level_n
                ref_n = lgt_block( lgtID_n, IDX_REFINE_STS)

                ! Send logic, I am sender and neighbor is receiver:
                ! stage=1 && lvl_diff = 0; stage=2 && lvl_diff=+1; stage=3 && lvl_diff=-1
                ! Some special cases:
                !    REF   - only sync if the refinement value matches s_val for the neighbor
                !    Level - only send if the neighbor level matches s_val

                ! send counter. how much data will I send to my neighbors on other mpiranks?
                b_send = .false.
                ! neighbor wants to receive all patches, ids correspond to full, full_RHS, full_ref, full_level
                if (any(sync_id == (/ 1,11,21,31/)) .and. &
                    ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==+1) .or. (istage == 3 .and. lvl_diff==-1))) then
                        b_send = .true.
                ! neighbor wants to receive medium and coarse patches -> send to MF, ids correspond to MC, MC_ref, MC_level
                elseif (any(sync_id == (/ 2,22,32/)) .and. &
                    ((istage == 1 .and. lvl_diff==0) .or. (istage == 3 .and. lvl_diff==-1))) then
                        b_send = .true.
                ! neighbor wants to receive medium and fine patches -> send to MC, ids correspond to MF, MF_ref, MF_level
                elseif (any(sync_id == (/ 3,23,33/)) .and. &
                    ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==+1))) then
                        b_send = .true.
                endif
                ! special cases, first ref check and then level check
                if (any(sync_id == (/21,22,23/)) .and. b_send) b_send = s_val == ref_n ! disable sync if neighbor has wrong ref value
                if (any(sync_id == (/31,32,33/)) .and. b_send) b_send = s_val == level_n  ! disable sync if neighbor has wrong level

                if (b_send) then
                    ! choose correct size that will be send, for lvl_diff /= 0 restriction or prediction will be applied
                    if (lvl_diff==0) then
                        ijk = ijkPatches(:, :, neighborhood, SENDER)
                    else
                        ijk = ijkPatches(:, :, neighborhood, RESPRE)
                    endif

                    if (myrank /= rank_n) then
                        data_send_counter(rank_n) = data_send_counter(rank_n) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each neighborhood, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of neighborhood etc
                        ! Increase by one so number of integers can vary
                        meta_send_counter(rank_n) = meta_send_counter(rank_n) + 1
                    endif

                    ! now lets save all metadata in one array without caring for rank sorting for now
                    meta_send_all(S_META_FULL*count_send + 1) = hvyID  ! needed for same-rank sending
                    meta_send_all(S_META_FULL*count_send + 2) = ref    ! needed for hvy_tmp for adapt_tree
                    meta_send_all(S_META_FULL*count_send + 3) = hvyID_n
                    meta_send_all(S_META_FULL*count_send + 4) = rank_n
                    meta_send_all(S_META_FULL*count_send + 5) = neighborhood
                    meta_send_all(S_META_FULL*count_send + 6) = lvl_diff
                    meta_send_all(S_META_FULL*count_send + 7) = (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                    
                    count_send = count_send + 1
                endif

                ! Receive logic, I am receiver and neighbor is sender:
                ! it is defined in logic relative to receiver, so that
                ! lvl_diff = -1 : recver from finer sender, restriction on sender side
                ! lvl_diff =  0 : recver is same level as sender
                ! lvl_diff = +1 : recver from coarser sender, interpolation on sender side
                !
                ! stage=1 && lvl_diff = 0; stage=2 && lvl_diff=-1; stage=3 && lvl_diff=+1
                ! Some special cases:
                !    REF   - only sync if the refinement value matches s_val for the receiver
                !    Level - only send if the receiver level matches s_val

                ! recv counter. how much data will I recv from neighbors on other mpiranks?
                ! This is NOT the same number as before
                if (myrank /= rank_n) then  ! only receive from foreign ranks
                    b_recv = .false.
                    ! I want to receive all patches, ids correspond to full, full_RHS, full_ref, full_level
                    if (any(sync_id == (/ 1,11,21,31/)) .and. &
                        ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==-1) .or. (istage == 3 .and. lvl_diff==+1))) then
                            b_recv = .true.
                    ! I want to receive patches from medium and coarse neighbors, ids correspond to MC, MC_ref, MC_level
                    elseif (any(sync_id == (/ 2,22,32/)) .and. &
                        ((istage == 1 .and. lvl_diff==0) .or. (istage == 3 .and. lvl_diff==+1))) then
                            b_recv = .true.
                    ! I want to receive patches from medium and fine neighbors, ids correspond to MF, MF_ref, MF_level
                    elseif (any(sync_id == (/ 3,23,33/)) .and. &
                        ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==-1))) then
                            b_recv = .true.
                    endif
                    ! special cases, first ref check and then level check
                    if (any(sync_id == (/21,22,23/)) .and. b_recv) b_recv = s_val == ref ! disable sync if I have wrong ref value
                    if (any(sync_id == (/31,32,33/)) .and. b_recv) b_recv = s_val == level  ! disable sync if I have wrong level

                    if (b_recv) then
                        ijk = ijkPatches(:, :, neighborhood, RECVER)

                        data_recv_counter(rank_n) = data_recv_counter(rank_n) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each neighborhood, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of neighborhood etc
                        ! Increase by one so number of integers can vary
                        meta_recv_counter(rank_n) = meta_recv_counter(rank_n) + 1
                    endif
                endif
            endif ! neighbor exists
        end do ! loop over all possible  neighbors
    end do ! loop over all heavy active


    ! NOTE: this feature is against wabbits memory policy: we try to allocate the
    ! whole memory of the machine on startup, then work with that. however, we have to
    ! reserve portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
    ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
    ! JB: This can only trigger if we change g during the run?, and why increase by 125%? What if that is not enough?
    if (sum(data_recv_counter) + sum(meta_recv_counter)*S_META_SEND + params%number_procs > size(rData_recvBuffer, 1)) then
        ! out-of-memory case: the preallocated buffer is not large enough.
        write(*,'("rank=",i4," OOM for ghost nodes and increases its receive buffer size to 125%")') myrank
        new_size = size(rData_recvBuffer,1)*125/100
        deallocate(rData_recvBuffer)
        allocate( rData_recvBuffer(1:new_size), stat=status )
        if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
    endif

    if (sum(data_send_counter) + sum(meta_send_counter)*S_META_SEND + params%number_procs > size(rData_sendBuffer, 1)) then
        ! out-of-memory case: the preallocated buffer is not large enough.
        write(*,'("rank=",i4," OOM for ghost nodes and increases its send buffer size to 125%")') myrank
        new_size = size(rData_sendBuffer,1)*125/100
        deallocate(rData_sendBuffer)
        allocate( rData_sendBuffer(1:new_size), stat=status )
        if (status /= 0) call abort(999993, "Buffer allocation failed. Not enough memory?")
    endif
end subroutine prepare_ghost_synch_metadata