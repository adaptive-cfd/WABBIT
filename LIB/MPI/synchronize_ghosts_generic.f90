!> Wrapper to synch blocks with temporary flag from finer neighbours and same-level neighbors
!> Used before wavelet decomposition
subroutine sync_level_from_M(params, hvy_block, tree_ID, level, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), intent(in)   :: level                         !< level value to be checked against
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus         !< Boundary sizes in case we want to send less values

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, spaghetti_form=.false., sync_case="Monly_level", &
        s_val=level)

end subroutine sync_level_from_M


!> Wrapper to synch blocks with temporary flag from finer neighbours and same-level neighbors
!> Used before wavelet decomposition
subroutine sync_TMP_from_MF(params, hvy_block, tree_ID, ref_check, hvy_tmp, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), intent(in)   :: ref_check                     !< ref value to be checked against
    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus         !< Boundary sizes in case we want to send less values

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, spaghetti_form=.false., sync_case="MoF_ref", &
        s_val=ref_check, hvy_tmp=hvy_tmp)

end subroutine sync_TMP_from_MF



!> Wrapper to synch blocks with temporary flag from all neighbors
!> Used before wavelet decomposition
subroutine sync_TMP_from_all(params, hvy_block, tree_ID, ref_check, hvy_tmp, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), intent(in)   :: ref_check                     !< ref value to be checked against
    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus         !< Boundary sizes in case we want to send less values

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, spaghetti_form=.false., sync_case="full_ref", &
        s_val=ref_check, hvy_tmp=hvy_tmp)

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
        call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, spaghetti_form=.true., sync_case="MoC_level", &
            s_val=level, hvy_tmp=hvy_tmp)
    else
        call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, spaghetti_form=.true., sync_case="MoC", &
            hvy_tmp=hvy_tmp)
    endif

end subroutine sync_SCWC_from_MC


!> Wrapper to synch all ghost-point patches
subroutine sync_ghosts_tree(params, hvy_block, tree_ID, g_minus, g_plus, ignore_Filter)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus       !< set specific ghost node sizes to be synced
    logical, optional, intent(in) :: ignore_Filter                  !< ignore restriction filter, defaults to false

    integer(kind=ik) :: gminus, gplus
    logical :: do_ignore_Filter
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus
    do_ignore_Filter = .false.
    if (present(ignore_Filter)) do_ignore_Filter = ignore_Filter

    ! set level to -1 to enable synching between all, set stati to send to all levels
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, sync_case="full", spaghetti_form=.false., ignore_Filter=do_ignore_Filter)

end subroutine sync_ghosts_tree

!> Wrapper to synch all ghost-point patches, but not using the filter.
!! This is done during the RHS evaluation: the stencil size of the FD discretization is short,
!! and consequently only the first few ghost nodes are used. In the coarseExtension, those are copied
!! anyways (being near the coarse/fine interface, the appropriate filter cannot be applied).
!! Therefore, for reasons of consistency and performance, we can ignore the filter in ghost nodes
!! sync'ing. This corresponds to using non-lifted wavelets.
!
!  One might ask, why for RHS with only +-shaped stencils we still synch the diagonals (edges/corners) as well?
!  This is, because within the synching we apply the coarsening and interpolation filters at the interfaces, which need the other stencils as well!
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
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, sync_case="full_leaf", spaghetti_form=.false., ignore_Filter=.true.)

end subroutine


!> This function deals with ghost-node synching \n
!! It is a generic function with many flags, streamlining all synching process \n
!! In order to avoid confusion wrapper functions should be used everywhere in order to implement
!! specific versions. This also means that parameter changes only have to be changed in the wrappers
subroutine sync_ghosts_generic( params, hvy_block, tree_ID, sync_case, spaghetti_form, &
    g_minus, g_plus, s_val, hvy_tmp, verbose_check, ignore_Filter)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    !> String representing which kind of syncing we want to do, varies between different sync cases (betwen CMF neighbors) and restrictions (tree, level or ref)
    character(len=*)          :: sync_case
    logical, intent(in)       :: spaghetti_form                 !< If set to true, then data is in spaghetti form of SC WC

    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag
    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, used if sync case includes ref or level
    integer(kind=ik), intent(in), optional  :: s_val
    logical, intent(in), optional  :: ignore_Filter                 !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus       !< Synch only so many ghost points

    integer(kind=ik) :: count_send_total, Nstages, istage, gminus, gplus, k, hvy_id, lgt_id, i_dim
    real(kind=rk) :: t0, t1, t2, x0(1:3), dx(1:3), tolerance
    character(len=clong) :: toc_statement
    integer(kind=2) :: n_domain(1:3)

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
    ! default is three stages:
    !    1. M2M, copy
    !    2. M2C, decimate, independent from coarser neighbor but needs same-level data
    !    3. M2F, interpolate, needs same-level and finer neighbor data
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
    if (sync_case == "full") call reset_ghost_nodes( params, hvy_block, tree_ID)
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
            istage, sync_case, ncomponents=size(hvy_block,4), s_val=s_val, verbose_check=verbose_check)
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
            if (index(sync_case, "ref") > 0) then
                call xfer_block_data(params, hvy_block, tree_ID, count_send_total, hvy_tmp=hvy_tmp, &
                REF_FLAG=s_val, ignore_Filter=ignore_Filter, verbose_check=verbose_check)
            else
                call xfer_block_data(params, hvy_block, tree_ID, count_send_total, hvy_tmp=hvy_tmp, &
                ignore_Filter=ignore_Filter, verbose_check=verbose_check)
            endif
        endif
        call toc( "sync ghosts (xfer_block_data)", 82, MPI_wtime()-t2 )

        !***************************************************************************
        ! (iii) Setting of BCs
        !***************************************************************************
        t2 = MPI_wtime()
        ! For all blocks that do not have a neighbor at the borders due to non-periodic BC, we have to set the ghost patch values explicitly
        ! After stage 2 and 3 we only need to update the edges with mirrored values from ghost layers, as the faces are already set
        call bound_cond_generic(params, hvy_block, tree_ID, sync_case, spaghetti_form, s_val=s_val, edges_only=istage>=2)
        call toc( "sync ghosts (bound_cond_generic)", 83, MPI_wtime()-t2 )

        write(toc_statement, '(A, i0, A)') "sync ghosts (stage ", istage, " TOTAL)"
        call toc( toc_statement, 83+istage, MPI_wtime()-t1 )

    end do ! loop over stages 1,2,3

    call toc( "sync ghosts (TOTAL)", 80, MPI_wtime()-t0 )

end subroutine sync_ghosts_generic



! This subroutine prepares who sends to whom. This includes:
!    - logic of different synchronization situations
!    - saving of all metadata
!    - computing of buffer sizes for metadata for both sending and receiving
! This is done strictly locally so no MPI needed here
subroutine prepare_ghost_synch_metadata(params, tree_ID, count_send, istage, sync_case, ncomponents, s_Val, verbose_check)

    implicit none

    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(in)    :: tree_ID             !< which tree to study

    integer(kind=ik), intent(in)    :: ncomponents         !< components can vary (for mask for example)
    integer(kind=ik), intent(out)   :: count_send          !< number of ghost patches total to be send, for looping
    !> following are variables that control the logic of where each block sends or receives
    integer(kind=ik), intent(in)    :: istage              !< current stage out of three
    !> String representing which kind of syncing we want to do, varies between different sync cases (betwen CMF neighbors) and restrictions (tree, level or ref)
    character(len=*)                :: sync_case

    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, used if sync case includes ref or level
    integer(kind=ik), intent(in), optional  :: s_val
    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag

    ! Following are global data used but defined in module_mpi:
    !    data_recv_counter, data_send_counter
    !    meta_recv_counter, meta_send_counter
    !    meta_send_all (possibly needs renaming after this function)

    integer(kind=ik) :: k_block, myrank, N, i_n, ijk(2,3), inverse, ierr, lvl_diff, status, new_size, sync_id
    integer(kind=ik) :: hvy_ID, lgt_ID, level, ref
    integer(kind=ik) :: hvy_ID_n, lgt_ID_n, level_n, ref_n, rank_n
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
    ! at first, different neighbor restrictions are considered and set to first digit
    if (index(sync_case, "full") > 0) sync_id = 1
    if (index(sync_case, "full_leaf") > 0) sync_id = 2
    if (index(sync_case, "MoC") > 0) sync_id = 3
    if (index(sync_case, "MoF") > 0) sync_id = 4
    if (index(sync_case, "Monly") > 0) sync_id = 5

    ! now lets treat the special restrictions, set to the second digit
    if (index(sync_case, "ref") > 0) sync_id = sync_id + 10*1
    if (index(sync_case, "level") > 0) sync_id = sync_id + 10*2

    if (sync_id == 0) then
        call abort(240805, "No, we don't trade that here! Please ensure the sync_case is valid.")
    endif

    count_send = 0
    do k_block = 1, hvy_n(tree_ID)
        ! calculate light id
        hvy_ID = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_ID, hvy_ID, myrank, N )
        level = lgt_block( lgt_ID, IDX_MESH_LVL )
        ref = lgt_block( lgt_ID, IDX_REFINE_STS)

        ! for leaf-syncs we ignore non-leaf blocks
        if (mod(sync_id,10) == 2 .and. .not. block_is_leaf(params, hvy_ID)) cycle
        ! blocks with ref status REF_TMP_EMPTY do not sync anything
        if (ref == REF_TMP_EMPTY) cycle

        ! loop over all neighbors
        do i_n = 1, size(hvy_neighbor, 2)
            ! neighbor exists
            if ( hvy_neighbor( hvy_ID, i_n ) /= -1 ) then
                ! neighbor light data id
                lgt_ID_n = hvy_neighbor( hvy_ID, i_n )
                ! calculate neighbor rank
                call lgt2proc( rank_n, lgt_ID_n, N )
                ! neighbor heavy id
                call lgt2hvy( hvy_ID_n, lgt_ID_n, rank_n, N )

                ! define level difference: sender - receiver, so +1 means sender on higher level
                ! lvl_diff = -1 : sender to finer recver, interpolation on sender side
                ! lvl_diff =  0 : sender is same level as recver
                ! lvl_diff = +1 : sender to coarser recver, restriction is applied on sender side
                level_n = lgt_block( lgt_ID_n, IDX_MESH_LVL )
                lvl_diff = level - level_n
                ref_n = lgt_block( lgt_ID_n, IDX_REFINE_STS)

                ! blocks with ref status REF_TMP_EMPTY do not sync anything
                if (ref_n == REF_TMP_EMPTY) cycle

                ! Send logic, I am sender and neighbor is receiver:
                ! stage=1 && lvl_diff = 0; stage=2 && lvl_diff=+1; stage=3 && lvl_diff=-1
                ! Some special cases:
                !    REF   - only sync if the refinement value matches s_val for the neighbor
                !    Level - only send if the neighbor level matches s_val

                ! CVS grids can have medium, fine and coarse neighbors for same patches, grid is assumed to be fully updated
                ! normally, we choose the simplest one (medium, then finer, then coarser neighbor) if multiple are availbel to update the patch
                ! For sending, some blocks have to send to multiple neighbors at the same patch in special conditions:
                !    leaf blocks (without daughters) have to send to both medium and fine neighbors
                !    root blocks (without mother) have to send to both medium and coarse neighbors
                ! CVS Leaf updates assume grid is not fully updated and only update to and from other leaf blocks
                ! this is the highest lvl_diff patch (finer, then medium, then coarser neighbors)

                ! send counter. how much data will I send to my neighbors on other mpiranks?
                b_send = .false.
                ! neighbor wants to receive all patches, ids correspond to full
                if (mod(sync_id,10) == 1 .and. &
                    ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==+1) .or. (istage == 3 .and. lvl_diff==-1))) then
                        b_send = .true.
                        ! CVS: non-leaf blocks do not send to fine neighbors
                        if (lvl_diff==-1 .and. .not. block_is_leaf(params, hvy_ID, check_empty=.true.)) then
                            b_send = .false.
                        endif
                        ! CVS: non-root blocks do not send to coarse neighbors
                        if (lvl_diff==+1 .and. .not. block_is_root(params, hvy_ID, check_empty=.true.)) then
                            b_send = .false.
                        endif

                ! neighbor wants to receive all leaf-patches, ids correspond to full_leaf
                elseif (mod(sync_id,10) == 2 .and. &
                    ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==+1) .or. (istage == 3 .and. lvl_diff==-1))) then
                        b_send = .true.
                        ! leaf-wise, Fine->Medium->Coarse
                        ! check for fine neighbor
                        if (lvl_diff>=0 .and. block_has_valid_neighbor(params, hvy_ID, i_n, -1)) then
                            b_send = .false.
                        ! check for medium neighbor
                        elseif (lvl_diff==+1 .and. block_has_valid_neighbor(params, hvy_ID, i_n, 0)) then
                            b_send = .false.
                        endif

                ! neighbor wants to receive medium and coarse patches -> send to MoF, ids correspond to MoC
                elseif (mod(sync_id,10) == 3 .and. &
                    ((istage == 1 .and. lvl_diff==0) .or. (istage == 3 .and. lvl_diff==-1))) then
                        b_send = .true.
                        ! CVS: non-leaf blocks do not send to fine neighbors
                        if (lvl_diff/=0 .and. .not. block_is_leaf(params, hvy_ID, check_empty=.true.)) then
                            b_send = .false.
                        endif

                ! neighbor wants to receive medium and fine patches -> send to MoC, ids correspond to MoF
                elseif (mod(sync_id,10) == 4 .and. &
                    ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==+1))) then
                        b_send = .true.
                        ! CVS: non-root blocks do not send to coarse neighbors
                        if (lvl_diff/=0 .and. .not. block_is_root(params, hvy_ID, check_empty=.true.)) then
                            b_send = .false.
                        endif
                ! neighbor wants to receive medium -> send to M, ids correspond to M
                elseif (mod(sync_id,10) == 5 .and. &
                    ((istage == 1 .and. lvl_diff==0))) then
                        b_send = .true.
                endif

                ! special cases, first ref check and then level check, situated in second digit
                if (sync_id/10 == 1 .and. b_send) b_send = s_val == ref_n ! disable sync if neighbor has wrong ref value
                if (sync_id/10 == 2 .and. b_send) b_send = s_val == level_n  ! disable sync if neighbor has wrong level

                if (b_send) then
                    ! choose correct size that will be send, for lvl_diff /= 0 restriction or prediction will be applied
                    if (lvl_diff==0) then
                        ijk = ijkPatches(:, :, i_n, SENDER)
                    else
                        ijk = ijkPatches(:, :, i_n, RESPRE)
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
                    meta_send_all(S_META_FULL*count_send + 1) = hvy_ID  ! needed for same-rank sending
                    meta_send_all(S_META_FULL*count_send + 2) = ref    ! needed for hvy_tmp for adapt_tree
                    meta_send_all(S_META_FULL*count_send + 3) = hvy_ID_n
                    meta_send_all(S_META_FULL*count_send + 4) = rank_n
                    meta_send_all(S_META_FULL*count_send + 5) = i_n
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

                ! CVS grids can have medium, fine and coarse neighbors for same patches, grid is assumed to be fully updated
                ! normally, we choose the simplest one (medium, then finer, then coarser neighbor) if multiple are availbel to update the patch
                ! For receiving, receiving has to be restricted:
                !    non-leaf blocks (with daughters) do not receive from coarse neighbors
                !    non-root blocks (with mother) do not receive from fine neighbors
                ! CVS Leaf updates assume grid is not fully updated and only update to and from other leaf blocks
                ! this is the highest lvl_diff patch (finer, then medium, then coarser neighbors)

                ! recv counter. how much data will I recv from neighbors on other mpiranks?
                ! This is NOT the same number as before
                if (myrank /= rank_n) then  ! only receive from foreign ranks
                    b_recv = .false.

                    ! I want to receive all patches, ids correspond to full
                    if (mod(sync_id,10) == 1 .and. &
                        ((istage == 1 .and. lvl_diff==0) .or. (istage == 3 .and. lvl_diff==+1) .or. (istage == 2 .and. lvl_diff==-1))) then
                            b_recv = .true.
                            ! CVS: non-leaf blocks do not receive from coarser neighbors
                            if (lvl_diff==+1 .and. .not. block_is_leaf(params, hvy_ID, check_empty=.true.)) then
                                b_recv = .false.
                            ! CVS: leaf blocks do not receive from coarser neighbors if there is a medium one for same patch
                            elseif (lvl_diff==+1 .and. block_has_valid_neighbor(params, hvy_ID, i_n, 0)) then
                                b_recv = .false.
                            endif                                
                            ! CVS: non-root blocks do not receive fine neighbors
                            if (lvl_diff==-1 .and. .not. block_is_root(params, hvy_ID, check_empty=.true.)) then
                                b_recv = .false.
                            ! CVS: root blocks do not receive from fine neighbors if there is a medium one for same patch
                            elseif (lvl_diff==-1 .and. block_has_valid_neighbor(params, hvy_ID, i_n, 0)) then
                                b_recv = .false.
                            endif

                    ! I want to receive all leaf-patches, ids correspond to full_leaf
                    elseif (mod(sync_id,10) == 2 .and. &
                        ((istage == 1 .and. lvl_diff==0) .or. (istage == 3 .and. lvl_diff==+1) .or. (istage == 2 .and. lvl_diff==-1))) then
                            b_recv = .true.
                            ! leaf-wise, Fine->Medium->Coarse
                            ! check for fine neighbor
                            if (lvl_diff>=0 .and. block_has_valid_neighbor(params, hvy_ID, i_n, -1)) then
                                b_recv = .false.
                            ! check for medium neighbor
                            elseif (lvl_diff==+1 .and. block_has_valid_neighbor(params, hvy_ID, i_n, 0)) then
                                b_recv = .false.
                            endif

                    ! I want to receive medium and coarse patches, ids correspond to MoC
                    elseif (mod(sync_id,10) == 3 .and. &
                        ((istage == 1 .and. lvl_diff==0) .or. (istage == 3 .and. lvl_diff==+1))) then
                            b_recv = .true.
                            ! CVS: non-leaf blocks do not receive from coarser neighbors
                            if (lvl_diff/=0 .and. .not. block_is_leaf(params, hvy_ID, check_empty=.true.)) then
                                b_recv = .false.
                            ! CVS: leaf blocks do not receive from coarser neighbors if there is a medium one for same patch
                            elseif (lvl_diff/=0 .and. block_has_valid_neighbor(params, hvy_ID, i_n, 0)) then
                                b_recv = .false.
                            endif

                    ! I want to receive medium and fine patches,ids correspond to MoF
                    elseif (mod(sync_id,10) == 4 .and. &
                        ((istage == 1 .and. lvl_diff==0) .or. (istage == 2 .and. lvl_diff==-1))) then
                            b_recv = .true.
                            ! CVS: non-root blocks do not receive from fine neighbors
                            if (lvl_diff/=0 .and. .not. block_is_root(params, hvy_ID, check_empty=.true.)) then
                                b_recv = .false.
                            ! CVS: root blocks do not receive from fine neighbors if there is a medium one for same patch
                            elseif (lvl_diff/=0 .and. block_has_valid_neighbor(params, hvy_ID, i_n, 0)) then
                                b_recv = .false.
                            endif
                    ! I want to receive medium,ids correspond to M
                    elseif (mod(sync_id,10) == 5 .and. &
                        ((istage == 1 .and. lvl_diff==0))) then
                            b_recv = .true.
                    endif

                    ! special cases, first ref check and then level check
                    if (sync_id/10 == 1 .and. b_recv) b_recv = s_val == ref ! disable sync if I have wrong ref value
                    if (sync_id/10 == 2 .and. b_recv) b_recv = s_val == level  ! disable sync if I have wrong level

                    if (b_recv) then
                        ijk = ijkPatches(:, :, i_n, RECVER)

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