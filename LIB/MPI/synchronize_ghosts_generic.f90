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

    ! we set s_level to REF_TMP_UNTREATED, this value is < -1 and therefore distinctive, we use this to avoid another parameter
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, s_level=REF_TMP_UNTREATED, s_F2M=.true., s_M2M=.true., s_C2M=.false., hvy_tmp=hvy_tmp)

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

    ! we set s_level to REF_TMP_UNTREATED, this value is < -1 and therefore distinctive, we use this to avoid another parameter
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, s_level=REF_TMP_UNTREATED, s_F2M=.true., s_M2M=.true., s_C2M=.true., hvy_tmp=hvy_tmp)

end subroutine sync_TMP_from_all



!> Wrapper to synch level from coarser neighbours and same-level neighbors
!! Used after coarse extension to update SC and WC, coarse neighbours need to be synched from hvy_tmp
subroutine sync_SCWC_from_MC(params, hvy_block, tree_ID, hvy_tmp, g_minus, g_plus)
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout)   :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus

    integer(kind=ik) :: gminus, gplus
    gminus = params%g
    gplus = params%g
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

    ! set level to -1 to enable synching between all
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, s_level=-1, s_M2F=.true., s_M2M=.true., hvy_tmp=hvy_tmp, verbose_check=.true.)

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
    call sync_ghosts_generic(params, hvy_block, tree_ID, g_minus=gminus, g_plus=gplus, s_level=-1, s_M2M=.true., s_M2F=.true., s_M2C=.true.)

end subroutine sync_ghosts_tree


!> This function deals with ghost-node synching \n
!! It is a generic function with many flags, streamlining all synching process \n
!! In order to avoid confusion wrapper functions should be used everywhere in order to implement
!! specific versions. This also means that parameter changes only have to be changed in the wrappers
subroutine sync_ghosts_generic( params, hvy_block, tree_ID, g_minus, g_plus, &
    s_level, s_M2M, s_M2C, s_C2M, s_M2F, s_F2M, hvy_tmp, verbose_check)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study

    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag

    !> Level to synch, if -1 then all levels are synched, if < -1 then it is REF_TMP_UNTREATED and ref status will be checked
    integer(kind=ik), intent(in), optional  :: s_level
    logical, intent(in), optional  :: s_M2M                         !< Synch from level J   to J
    logical, intent(in), optional  :: s_M2C                         !< Synch from level J   to J-1
    logical, intent(in), optional  :: s_C2M                         !< Synch from level J-1 to J
    logical, intent(in), optional  :: s_M2F                         !< Synch from level J   to J+1
    logical, intent(in), optional  :: s_F2M                         !< Synch from level J+1 to J
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus       !< Synch only so many ghost points

    integer(kind=ik) sLevel
    logical :: SM2M, SM2C, SC2M, SM2F, SF2M

    integer(kind=ik)   :: myrank, mpisize, Bs(1:3), buffer_offset
    integer(kind=ik)   :: N, k, neighborhood, level_diff, Nstages
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

    ! initialize variables, to write 5 times .false. might be long but I tried other ways which surprisingly delivered wrong results
    sLevel = -1
    sM2M = .false.
    sM2C = .false.
    sC2M = .false.
    sM2F = .false.
    sF2M = .false.
    if (present(s_Level)) sLevel = s_Level
    if (present(s_M2M)) sM2M = s_M2M
    if (present(s_M2C)) sM2C = s_M2C
    if (present(s_C2M)) sC2M = s_C2M
    if (present(s_M2F)) sM2F = s_M2F
    if (present(s_F2M)) sF2M = s_F2M

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
    ! ijkPatches([start,end], [dir], [ineighbor], [leveldiff], [isendrecv])
    ! As g can be varied (as long as it does not exceed the maximum value params%g), it is set up
    ! each time we sync (at negligibble cost)
    call ghosts_setup_patches(params, gminus=gminus, gplus=gplus, output_to_file=.false.)
    ! some tiny buffers depend on the number of components (nc=size(hvy_block,4))
    ! make sure they have the right size
    call xfer_ensure_correct_buffer_size(params, hvy_block)

#ifdef DEV
    ! for dev check ghosts by wiping if we set all of them
    if (sLevel == -1) call reset_ghost_nodes( params, hvy_block, tree_ID, s_M2M=sM2M, s_M2C=sM2C, s_M2F=sM2F)
    ! if (present(verbose_check)) then
    !     call reset_ghost_nodes( params, hvy_block, tree_ID, s_M2M=.true., s_M2C=.true., s_M2F=.true.)
    ! endif

#endif

! Diagonal neighbors (not required for the RHS)
! 2D: 5,6,7,8
! 3D: 7-18, 19-26, 51-74
!
! 2nd IDEA: stage-free ghost nodes.
! ==> did not work yet
! It turn out, if the coarser block sends not-correctly interpolatedy points, they can be corrected on the
! receiving fine block. Only a few are affected; the ones toward the interface. The ones further away
! are directly interpolated correctly. 19 apr 2023: This can be made work, but it is tedious and does not yield an
! immense speedup neither. On my local machine, its ~5%, but on large scale parallel sims, it may be more significant.
! Idea is described in inskape notes.
!

    ! We require two stages: first, we fill all ghost nodes which are simple copy (including restriction),
    ! then in the second stage we can use interpolation and fill the remaining ones.
    do istage = 1, Nstages
        !***************************************************************************
        ! (i) stage initialization
        !***************************************************************************
        ! prepare metadata. This computes from the grid info how much data I recv and send to all other mpiranks.
        ! Also applies logic about what should be synched and saves all metadata unsorted in one array
        ! internal nodes are included in metadata but not counted
        t1 = MPI_wtime()
        t2 = MPI_wtime()  ! stage duration
        call prepare_ghost_synch_metadata(params, tree_ID, count_send_total, &
            istage, ncomponents=size(hvy_block,4), s_Level=sLevel, s_M2M = sM2M, s_M2C = sM2C, s_C2M = sC2M, s_M2F = sM2F, s_F2M = sF2M)
        call toc( "sync ghosts (prepare metadata)", 81, MPI_wtime()-t1 )

        !***************************************************************************
        ! (ii) sending handled by xfer_block_data
        ! If hvy_temp is present then xfer_block_data has to decide from where to grab the data
        !    - with sLevel < -1: we decide after refinement flag present in sLevel if we want to use hvy_temp
        !    - elsewise: use hvy_tmp for prediction (used in updating SC from coarser neighbours)
        !***************************************************************************
        t1 = MPI_wtime()
        if (.not. present(hvy_tmp)) then
            call xfer_block_data(params, hvy_block, tree_ID, count_send_total, verbose_check=verbose_check)
        else
            if (sLevel < -1) then
                call xfer_block_data(params, hvy_block, tree_ID, count_send_total, hvy_tmp=hvy_tmp, REF_FLAG=sLevel, verbose_check=verbose_check)
            else
                call xfer_block_data(params, hvy_block, tree_ID, count_send_total, hvy_tmp=hvy_tmp, verbose_check=verbose_check)
            endif
        endif
        call toc( "sync ghosts (xfer_block_data)", 82, MPI_wtime()-t1 )
        write(toc_statement, '(A, i0, A)') "sync ghosts (stage ", istage, " TOTAL)"
        call toc( toc_statement, 82+istage, MPI_wtime()-t1 )

    end do ! loop over stages 1,2

    call toc( "sync ghosts (TOTAL)", 80, MPI_wtime()-t0 )

end subroutine sync_ghosts_generic

! subroutine sync_ghosts_nostages( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!     implicit none
!
!     type (type_params), intent(in) :: params
!     !> light data array
!     integer(kind=ik), intent(in)   :: lgt_block(:, :)
!     !> heavy data array - block data
!     real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)
!     !> heavy data array - neighbor data
!     integer(kind=ik), intent(in)   :: hvy_neighbor(:,:)
!     !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)   :: hvy_active(:)
!     !> number of active blocks (heavy data)
!     integer(kind=ik), intent(in)   :: hvy_n
!
!     integer(kind=ik)   :: myrank, mpisize, ii0, ii1, Bs(1:3)
!     integer(kind=ik)   :: N, k, neighborhood, level_diff
!     integer(kind=ik)   :: recver_lgtID, recver_rank, recver_hvyID
!     integer(kind=ik)   :: sender_hvyID, sender_lgtID
!
!     integer(kind=ik) :: ijk(2,3), isend, irecv
!     integer(kind=ik) :: bounds_type, inverse
!     real(kind=rk) :: t0
!
!     t0 = MPI_wtime()
!
!     if (.not. ghost_nodes_module_ready) then
!         ! in order to keep the syntax clean, buffers are module-global and need to be
!         ! allocated here.
!         call init_ghost_nodes( params )
!     endif
!
!     ! if this mpirank has no active blocks, it has nothing to do here.
!     if (hvy_n == 0) return
!
! !
! ! New Idea (09 Apr 2023)
! !
! ! During adapt_mesh, which is the most performance-hungry part of the algorithm with biorthogonal wavelets
! ! we should sync certain levels only. This we need also for the complete wavelet transform and denoising/CVS.
! ! Synching a level should mean we syn (J,J-1). We start from Jmax. This means: sync all neighborhoods that
! ! involve level J, also on blocks of level J-1. It does not mean a block on (J-1) is completely synced, just those
! ! neighborhoods.
! !
! ! 2nd IDEA: stage-free ghost nodes.
! ! It turn out, if the coarser block sends not-correctly interpolatedy points, they can be corrected on the
! ! receiving fine block. Only a few are affected; the ones toward the interface. The ones further away
! ! are directly interpolated correctly.
!
!     Bs    = params%Bs
!     N     = params%number_blocks
!     myrank  = params%rank
!     mpisize = params%number_procs
!
!     ! call reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )
!
!         !***************************************************************************
!         ! (i) stage initialization
!         !***************************************************************************
!         int_pos(:) = 1
!         real_pos(:) = 0
!         internalNeighbor_pos = 1
!
!         ! compute, locally from the grid info, how much data I recv from and send to all
!         ! other mpiranks. this gives us the start indices of each rank in the send/recv buffers. Note
!         ! we do not count our internal nodes (.false. as last argument), as they are not put in the
!         ! buffer at any time.
!         call get_my_sendrecv_amount_with_ranks_nostages(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, &
!              data_recv_counter, data_send_counter, meta_recv_counter, meta_send_counter, &
!              count_internal=.false., ncomponents=size(hvy_block,4))
!
!
!         ! reset iMetaData_sendBuffer, but only the parts that will actually be treated.
!         do k = 1, params%number_procs
!             ii0 = sum(meta_send_counter(0:(k-1)-1)) + 1
!             ii1 = ii0 + meta_send_counter(k-1)
!             iMetaData_sendBuffer(ii0:ii1) = -99
!
!             ii0 = sum(meta_recv_counter(0:(k-1)-1)) + 1
!             ii1 = ii0 + meta_recv_counter(k-1)
!             iMetaData_recvBuffer(ii0:ii1) = -99
!         enddo
!
!
!         !***************************************************************************
!         ! (ii) prepare data for sending
!         !***************************************************************************
!         do k = 1, hvy_n
!             ! calculate light id
!             sender_hvyID = hvy_active(k)
!             call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )
!
!             ! loop over all neighbors
!             do neighborhood = 1, size(hvy_neighbor, 2)
!                 ! neighbor exists
!                 if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
!                     ! neighbor light data id
!                     recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
!                     ! calculate neighbor rank
!                     call lgt2proc( recver_rank, recver_lgtID, N )
!                     ! neighbor heavy id
!                     call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )
!                     ! define level difference: sender - receiver, so +1 means sender on higher level
!                     ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
!                     ! leveldiff =  0 : sender is same level as recver
!                     ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
!                     level_diff = lgt_block( sender_lgtID, IDX_MESH_LVL ) - lgt_block( recver_lgtID, IDX_MESH_LVL )
!
!                     if ( myrank == recver_rank ) then
!                         ! internal relation (no communication) has its own buffer (to avoid senseless copying
!                         ! from send to recv buffer)
!                         meta_send_all( internalNeighbor_pos:internalNeighbor_pos+4-1 ) = (/sender_hvyID, recver_hvyID, neighborhood, level_diff/)
!                         internalNeighbor_pos = internalNeighbor_pos + 4
!                     else
!                         ! external relation (MPI communication)
!                         call send_prepare_external( params, recver_rank, hvy_block, sender_hvyID, recver_hvyID, neighborhood, level_diff )
!                     end if ! (myrank==recver_rank)
!                 end if ! neighbor exists
!             end do ! loop over all possible  neighbors
!         end do ! loop over all heavy active
!
!         !***************************************************************************
!         ! (iii) transfer part (send/recv)
!         !***************************************************************************
!         call start_xfer_mpi( params, iMetaData_sendBuffer, rData_sendBuffer, iMetaData_recvBuffer, rData_recvBuffer, isend, irecv )
!
!
!         !***************************************************************************
!         ! (iv) Unpack received data in the ghost node layers
!         !***************************************************************************
!         ! process-internal ghost points (direct copy)
!         call unpack_ghostlayers_internal( params, hvy_block )
!
!         ! before unpacking the data we received from other ranks, we wait for the transfer
!         ! to be completed
!         call finalize_xfer_mpi(params, isend, irecv)
!
!         ! process-external ghost points (copy from buffer)
!         call unpack_ghostlayers_external( params, hvy_block )
!
!         ! call fixInterpolatedPoints_postSync( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
!
!     call toc( "WRAPPER: sync ghosts", 19, MPI_wtime()-t0 )
!
! end subroutine sync_ghosts_nostages



! This subroutine prepares who sends to whom. This includes:
!    - logic of different synchronization situations
!    - saving of all metadata
!    - computing of buffer sizes for metadata for both sending and receiving
! This is done strictly locally so no MPI needed here
subroutine prepare_ghost_synch_metadata(params, tree_ID, count_send, istage, ncomponents, &
    s_Level, s_M2M, s_M2C, s_C2M, s_M2F, s_F2M)

    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: tree_ID             !< which tree to study

    integer(kind=ik), intent(in)        :: ncomponents         !< components can vary (for mask for example)
    integer(kind=ik), intent(out)       :: count_send          !< number of ghost patches total to be send, for looping
    !> following are variables that control the logic of where each block sends or receives
    integer(kind=ik), intent(in)        :: istage  !< current stage out of three

    !> Level to synch, if -1 then all levels are synched, if < -1 then it is REF_TMP_UNTREATED and ref status will be checked
    integer(kind=ik), intent(in), optional  :: s_level
    logical, intent(in), optional  :: s_M2M                         !< Synch from level J   to J
    logical, intent(in), optional  :: s_M2C                         !< Synch from level J   to J-1
    logical, intent(in), optional  :: s_C2M                         !< Synch from level J-1 to J
    logical, intent(in), optional  :: s_M2F                         !< Synch from level J   to J+1
    logical, intent(in), optional  :: s_F2M                         !< Synch from level J+1 to J
    integer(kind=ik) sLevel
    logical :: SM2M, SM2C, SC2M, SM2F, SF2M

    ! Following are global data used but defined in module_mpi:
    !    data_recv_counter, data_send_counter
    !    meta_recv_counter, meta_send_counter
    !    meta_send_all (possibly needs renaming after this function)

    integer(kind=ik) :: k_block, sender_hvyID, sender_lgtID, sender_ref, myrank, N, neighborhood, recver_rank, recver_ref
    integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID, level, level_diff, status, new_size

    ! initialize variables, to write 5 times .false. might be long but I tried other ways which surprisingly delivered wrong results
    sLevel = -1
    sM2M = .false.
    sM2C = .false.
    sC2M = .false.
    sM2F = .false.
    sF2M = .false.
    if (present(s_Level)) sLevel = s_Level
    if (present(s_M2M)) sM2M = s_M2M
    if (present(s_M2C)) sM2C = s_M2C
    if (present(s_C2M)) sC2M = s_C2M
    if (present(s_M2F)) sM2F = s_M2F
    if (present(s_F2M)) sF2M = s_F2M

    myrank = params%rank
    N = params%number_blocks

    data_recv_counter(:) = 0
    data_send_counter(:) = 0
    meta_send_counter(:) = 0
    meta_recv_counter(:) = 0

    int_pos(:) = 0
    real_pos(:) = 0

    count_send = 0
    do k_block = 1, hvy_n(tree_ID)
        ! calculate light id
        sender_hvyID = hvy_active(k_block, tree_ID)
        call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )
        level = lgt_block( sender_lgtID, IDX_MESH_LVL )
        sender_ref = lgt_block( sender_lgtID, IDX_REFINE_STS)

        ! loop over all neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! if (skipDiagonalNeighbors) then
            !     ! Diagonal neighbors (not required for the RHS)
            !     ! 2D: 5,6,7,8
            !     ! 3D: 7-18, 19-26, 51-74
            !     if (dim==2.and.(neighborhood>=5.and.neighborhood<=8)) cycle
            !     if (dim==3.and.((neighborhood>=7.and.neighborhood<=26).or.(neighborhood>=51))) cycle
            ! endif

            ! neighbor exists
            if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
                ! neighbor light data id
                recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
                ! calculate neighbor rank
                call lgt2proc( recver_rank, recver_lgtID, N )
                ! neighbor heavy id
                call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )

                ! define level difference: sender - receiver, so +1 means sender on higher level
                ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
                ! leveldiff =  0 : sender is same level as recver
                ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
                level_diff = level - lgt_block( recver_lgtID, IDX_MESH_LVL )
                recver_ref = lgt_block( recver_lgtID, IDX_REFINE_STS)

                ! Send logic, following cases exist currently, all linked as .or.:
                ! stage=2, level_diff = +1, (sLevel=-1 and M2C) or (level=sLevel and M2C) or (level=sLevel+1 and F2M)
                !          or (sLevel<-1 and ref=sLevel and M2C) or (sLevel<-1 and ref_n=sLevel and F2M)
                ! stage=1, level_diff =  0, (sLevel=-1 and M2M) or (level=sLevel and M2M)
                !          or (sLevel<-1 and (ref=sLevel  or ref_n=sLevel) and M2M)
                ! stage=3, level_diff = -1, (sLevel=-1 and M2F) or (level=sLevel and M2F) or (level=sLevel-1 and C2M)
                !          or (sLevel<-1 and ref=sLevel and M2F) or (sLevel<-1 and ref_n=sLevel and C2M)

                ! send counter. how much data will I send to other mpiranks?
                if  ((istage==2 .and. level_diff==+1 .and. ((sLevel==-1 .and. sM2C) .or. (level==sLevel .and. sM2C) .or. (level==sLevel+1 .and. sF2M) &
                    .or. (sLevel<-1 .and. sender_ref==sLevel .and. sM2C) .or. (sLevel<-1 .and. recver_ref==sLevel .and. sF2M))) &
                .or. (istage==3 .and. level_diff==-1 .and. ((sLevel==-1 .and. sM2F) .or. (level==sLevel .and. sM2F) .or. (level==sLevel-1 .and. sC2M) &
                    .or. (sLevel<-1 .and. sender_ref==sLevel .and. sM2F) .or. (sLevel<-1 .and. recver_ref==sLevel .and. sC2M))) &
                .or. (istage==1 .and. level_diff== 0 .and. ((sLevel==-1 .and. sM2M) .or. (level==sLevel .and. sM2M) &
                    .or. (sLevel<-1 .and. (sender_ref==sLevel .or. recver_ref==sLevel) .and. sM2M)))) then
                    ! why is this RECVER and not sender? Because we adjust the data to the requirements of the
                    ! receiver before sending with interpolation or downsampling.
                    ijk = ijkPatches(:, :, neighborhood, level_diff, RECVER)

                    if (myrank /= recver_rank) then
                        data_send_counter(recver_rank) = data_send_counter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each neighborhood, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of neighborhood etc
                        ! Increase by one so number of integers can vary
                        meta_send_counter(recver_rank) = meta_send_counter(recver_rank) + 1
                    endif

                    ! now lets save all metadata in one array without caring for rank sorting for now
                    meta_send_all(S_META_FULL*count_send + 1) = sender_hvyID  ! needed for same-rank sending
                    meta_send_all(S_META_FULL*count_send + 2) = sender_ref    ! needed for hvy_tmp for adapt_tree
                    meta_send_all(S_META_FULL*count_send + 3) = recver_hvyID
                    meta_send_all(S_META_FULL*count_send + 4) = recver_rank
                    meta_send_all(S_META_FULL*count_send + 5) = neighborhood
                    meta_send_all(S_META_FULL*count_send + 6) = level_diff
                    meta_send_all(S_META_FULL*count_send + 7) = (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                    
                    count_send = count_send + 1
                endif

                ! Receive logic, following cases exist currently, all linked as .or.:
                ! stage=2, level_diff = -1, (sLevel=-1 and M2C) or (level=sLevel and F2M) or (level=sLevel-1 and M2C)
                !          or (sLevel<-1 and ref_n=sLevel and M2C) or (sLevel<-1 and ref=sLevel and F2M)
                ! stage=1, level_diff =  0, (sLevel=-1 and M2M) or (level=sLevel and M2M)
                !          or (sLevel<-1 and (ref=sLevel  or ref_n=sLevel) and M2M)
                ! stage=3, level_diff = +1, (sLevel=-1 and M2F) or (level=sLevel and C2M) or (level=sLevel+1 and M2F)
                !          or (sLevel<-1 and ref_n=sLevel and M2F) or (sLevel<-1 and ref=sLevel and C2M)

                ! recv counter. how much data will I recv from other mpiranks?
                ! This is NOT the same number as before
                if (myrank /= recver_rank) then  ! only receive from foreign ranks
                    if  ((istage==2 .and. level_diff==-1 .and. ((sLevel==-1 .and. sM2C) .or. (level==sLevel .and. sF2M) .or. (level==sLevel-1 .and. sM2C) &
                        .or. (sLevel<-1 .and. recver_ref==sLevel .and. sM2C) .or. (sLevel<-1 .and. sender_ref==sLevel .and. sF2M))) &
                    .or. (istage==3 .and. level_diff==+1 .and. ((sLevel==-1 .and. sM2F) .or. (level==sLevel .and. sC2M) .or. (level==sLevel+1 .and. sM2F) &
                        .or. (sLevel<-1 .and. recver_ref==sLevel .and. sM2F) .or. (sLevel<-1 .and. sender_ref==sLevel .and. sC2M))) &
                    .or. (istage==1 .and. level_diff== 0 .and. ((sLevel==-1 .and. sM2M) .or. (level==sLevel .and. sM2M) &
                        .or. (sLevel<-1 .and. (sender_ref==sLevel .or. recver_ref==sLevel) .and. sM2M)))) then
                        inverse = inverse_neighbor(neighborhood, dim)

                        ijk = ijkPatches(:, :, inverse, -1*level_diff, RECVER)

                        data_recv_counter(recver_rank) = data_recv_counter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each neighborhood, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of neighborhood etc
                        ! Increase by one so number of integers can vary
                        meta_recv_counter(recver_rank) = meta_recv_counter(recver_rank) + 1
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


! ! returns two lists with numbers of points I send to all other procs and how much I
! ! receive from each proc. note: strictly locally computed, NO MPI comm involved here
! subroutine get_my_sendrecv_amount_with_ranks_nostages(params, lgt_block, hvy_neighbor, hvy_active,&
!      hvy_n, data_recv_counter, data_send_counter, meta_recv_counter, meta_send_counter, &
!      count_internal, ncomponents)
!
!     implicit none
!
!     type (type_params), intent(in)      :: params
!     !> light data array
!     integer(kind=ik), intent(in)        :: lgt_block(:, :)
!     !> heavy data array - neighbor data
!     integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
!     !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_active(:)
!     !> number of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_n, ncomponents
!     integer(kind=ik), intent(inout)     :: data_recv_counter(0:), data_send_counter(0:)
!     integer(kind=ik), intent(inout)     :: meta_recv_counter(0:), meta_send_counter(0:)
!     logical, intent(in)                 :: count_internal
!
!     integer(kind=ik) :: k, sender_hvyID, sender_lgtID, myrank, N, neighborhood, recver_rank
!     integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID,level_diff, status, new_size
!
!     call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
!     N = params%number_blocks
!
!     data_recv_counter(:) = 0
!     data_send_counter(:) = 0
!     meta_recv_counter(:) = 0
!     meta_send_counter(:) = 0
!
!     do k = 1, hvy_n
!         ! calculate light id
!         sender_hvyID = hvy_active(k)
!         call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )
!
!         ! loop over all neighbors
!         do neighborhood = 1, size(hvy_neighbor, 2)
!             ! neighbor exists
!             if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
!                 ! neighbor light data id
!                 recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
!                 ! calculate neighbor rank
!                 call lgt2proc( recver_rank, recver_lgtID, N )
!                 ! neighbor heavy id
!                 call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )
!
!                 ! define level difference: sender - receiver, so +1 means sender on higher level
!                 ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
!                 ! leveldiff =  0 : sender is same level as recver
!                 ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
!                 level_diff = lgt_block( sender_lgtID, IDX_MESH_LVL ) - lgt_block( recver_lgtID, IDX_MESH_LVL )
!
!
!                 if (recver_rank /= myrank .or. count_internal) then
!                     ! it now depends on the stage if we have to sent this data
!                     ! or not.
!                     ! In stage 1, only level_diff = {+1, 0} is treated
!                     ! In stage 2, only level_diff = -1
!
!                     ! send counter. how much data will I send to other mpiranks?
!                     ! why is this RECVER and not sender? Well, complicated. The amount of data on the sender patch
!                     ! is not the same as in the receiver patch, because we interpolate or downsample. We effectively
!                     ! transfer only the data the recver wants - not the extra data.
!                     ijk = ijkPatches(:, :, neighborhood, level_diff, RECVER)
!
!                     data_send_counter(recver_rank) = data_send_counter(recver_rank) + &
!                     (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
!
!                     ! recv counter. how much data will I recv from other mpiranks?
!                     ! This is NOT the same number as before
!                     inverse = inverse_neighbor(neighborhood, dim)
!
!                     ijk = ijkPatches(:, :, inverse, -1*level_diff, RECVER)
!
!                     data_recv_counter(recver_rank) = data_recv_counter(recver_rank) + &
!                     (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
!
!                     ! counter for integer buffer: for each neighborhood, we send 6 integers as metadata
!                     ! as this is a fixed number it does not depend on the type of neighborhood etc, so
!                     ! technically one would need only one for send/recv
!                     meta_send_counter(recver_rank) = meta_send_counter(recver_rank) + 6 ! FIVE
!
!                     meta_recv_counter(recver_rank) = meta_recv_counter(recver_rank) + 6 ! FIVE
!                 endif
!
!
!
!             end if ! neighbor exists
!         end do ! loop over all possible  neighbors
!     end do ! loop over all heavy active
!
!
!
!     ! NOTE: for the int buffer, we mosly start at some index l0 and then loop unitl
!     ! we find a -99 indicating the end of the buffer. this could be avoided by using
!     ! for instead of while loops in the main routines, but I do not have time now.
!     !
!     ! In the meantime, notice we extent the amount of data by one, to copy the last -99
!     ! to the buffers
!     meta_recv_counter(:) = meta_recv_counter(:) + 1
!     meta_send_counter(:) = meta_send_counter(:) + 1
!
!
!     ! NOTE ACTUAL SEND / RECV DATA IS NEQN
!     data_recv_counter(:) = data_recv_counter(:) * ncomponents
!     data_send_counter(:) = data_send_counter(:) * ncomponents
!
!     ! NOTE: this feature is against wabbits memory policy: we try to allocate the
!     ! whole memory of the machine on startup, then work with that. however, we have to
!     ! reserver portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
!     ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
!     if (sum(data_recv_counter) > size(rData_recvBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_recvBuffer,1)*125/100
!         deallocate(rData_recvBuffer)
!         allocate( rData_recvBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
!     endif
!
!     if (sum(data_send_counter) > size(rData_sendBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_sendBuffer,1)*125/100
!         deallocate(rData_sendBuffer)
!         allocate( rData_sendBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999993, "Buffer allocation failed. Not enough memory?")
!     endif
! end subroutine
