subroutine sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, &
    syncSameLevelOnly, g_minus, g_plus  )
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in) :: params
    !> light data array
    integer(kind=ik), intent(in)   :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)   :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)   :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)   :: hvy_n
    logical, intent(in), optional  :: syncSameLevelOnly
    integer(kind=ik), optional, intent(in) :: g_minus, g_plus

    logical :: syncSameLevelOnly1=.false.

    integer(kind=ik)   :: myrank, mpisize, ii0, ii1, Bs(1:3), buffer_offset
    integer(kind=ik)   :: N, k, neighborhood, level_diff, Nstages
    integer(kind=ik)   :: recver_rank, recver_hvyID, patch_size
    integer(kind=ik)   :: sender_hvyID, sender_lgtID

    integer(kind=ik) :: ijk(2,3), isend, irecv, s(1:4), count_send_total
    integer(kind=ik) :: bounds_type, istage, inverse, gminus, gplus
    real(kind=rk) :: t0, t1

    t0 = MPI_wtime()

    if (.not. ghost_nodes_module_ready) then
        ! in order to keep the syntax clean, buffers are module-global and need to be
        ! allocated here.
        call init_ghost_nodes( params )
    endif

    ! if this mpirank has no active blocks, it has nothing to do here.
    if (hvy_n == 0) return

    syncSameLevelOnly1 = .false.
    if (present(syncSameLevelOnly)) syncSameLevelOnly1 = syncSameLevelOnly

    gminus  = params%g
    gplus   = params%g
    Bs      = params%Bs
    N       = params%number_blocks
    myrank  = params%rank
    mpisize = params%number_procs
    ! default is two stages: first, copy&decimate, then interpolate
    Nstages = 2
    ! but if we sync only sameLevel neighbors, then we don't need the second stage
    if (syncSameLevelOnly1) Nstages = 1
    ! if we sync a different number of ghost nodes
    if (present(g_minus)) gminus = g_minus
    if (present(g_plus))   gplus = g_plus

#ifdef DEV
    if (.not.syncSameLevelOnly1) call reset_ghost_nodes( params, hvy_block, hvy_active, hvy_n )
#endif

    !-----------------------------------------------------------------------
    ! set up constant arrays
    !-----------------------------------------------------------------------
    ! We frequently need to know the indices of a ghost nodes patch. Thus we save them
    ! once in a module-global array (which is faster than computing it every time with tons
    ! of IF-THEN clauses).
    ! This arrays indices are:
    ! ijkGhosts([start,end], [dir], [ineighbor], [leveldiff], [isendrecv])
    ! As g can be varied (as long as it does not exceed the maximum value params%g), it is set up
    ! each time we sync (at negligibble cost)
    call ghosts_setup_patches(params, gminus=gminus, gplus=gplus, output_to_file=.false.)
    ! some tiny buffers depend on the number of components (nc=size(hvy_block,4))
    ! make sure they have the right size
    call ghosts_ensure_correct_buffer_size(params, hvy_block)

! Diagonal neighbors (not required for the RHS)
! 2D: 5,6,7,8
! 3D: 7-18, 19-26, 51-74
!
! New Idea (09 Apr 2023)
!
! During adapt_mesh, which is the most performance-hungry part of the algorithm with biorthogonal wavelets
! we should sync certain levels only. This we need also for the complete wavelet transform and denoising/CVS.
! Synching a level should mean we syn (J,J-1). We start from Jmax. This means: sync all neighborhoods that
! involve level J, also on blocks of level J-1. It does not mean a block on (J-1) is completely synced, just those
! neighborhoods.
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
        int_pos(:) = 0
        real_pos(:) = 0

        ! prepare metadata. This computes from the grid info how much data I recv and send to all other mpiranks.
        ! Also applies logic about what should be synched and saves all metadata unsorted in one array
        ! internal nodes are included in metadata but not counted
        t1 = MPI_wtime()
        call prepare_ghost_synch_metadata(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, count_send_total, &
            istage, ncomponents=size(hvy_block,4), syncSameLevelOnly=syncSameLevelOnly1)
        call toc( "sync ghosts (prepare metadata)", MPI_wtime()-t1 )


#ifdef DEV
        ! debugging: reset buffers
        rData_sendBuffer(:) = -1
        rData_recvBuffer(:) = -1
#endif
        hvy_block(:,:,:,:,:) = 99

        !***************************************************************************
        ! (ii) prepare data for sending
        !***************************************************************************
        t1 = MPI_wtime()
        do k = 0, mpisize-1
            ! write in total size of entries at beginning of each part and pre-shift real_pos
            buffer_offset = sum(MetaData_sendCounter(0:k-1))*S_META_SEND + sum(Data_sendCounter(0:k-1)) + k + 1
            rData_sendBuffer(buffer_offset) = MetaData_sendCounter(k)
            real_pos(k) = 1 + MetaData_sendCounter(k)*S_META_SEND  ! offset real data to beginning by metadata

            ! test send/receive sizes
            ! write(*, '("Rank ", i0, " Send ", 4(i0, 1x), "Receive ", 4(i0, 1x), "Patches ", 4(i0, 1x))') myrank, Data_sendCounter, Data_recvCounter, MetaData_sendCounter
        end do

        do k = 0, count_send_total-1  ! we do MPI so lets stick to 0-based for a moment
            recver_rank = internalNeighborSyncs(S_META_FULL*k + 3)
            ! internal relation needs to prepare nothing and will be handled later
            if ( myrank /= recver_rank ) then
                ! unpack from metadata array what we need
                sender_hvyID = internalNeighborSyncs(S_META_FULL*k + 1)
                recver_hvyID = internalNeighborSyncs(S_META_FULL*k + 2)
                neighborhood = internalNeighborSyncs(S_META_FULL*k + 4)
                level_diff = internalNeighborSyncs(S_META_FULL*k + 5)
                patch_size = internalNeighborSyncs(S_META_FULL*k + 6)

                ! external relation (MPI communication)
                call send_prepare_external( params, recver_rank, hvy_block, sender_hvyID, recver_hvyID, neighborhood, level_diff, patch_size)
            end if
        end do ! loop over all patches (all heavy_n where neighborhood exists and ghost synching is applied)

        call toc( "sync ghosts (prepare data)", MPI_wtime()-t1 )
        !***************************************************************************
        ! (iii) transfer part (send/recv)
        !***************************************************************************
        t1 = MPI_wtime()
        call start_xfer_mpi( params, iMetaData_sendBuffer, rData_sendBuffer, iMetaData_recvBuffer, rData_recvBuffer, isend, irecv )
        call toc( "sync ghosts (start_xfer_mpi)", MPI_wtime()-t1 )

        !***************************************************************************
        ! (iv) Unpack received data in the ghost node layers
        !***************************************************************************
        ! process-internal ghost points (direct copy)
        t1 = MPI_wtime()
        call unpack_ghostlayers_internal( params, hvy_block, count_send_total )
        call toc( "sync ghosts (unpack internal)", MPI_wtime()-t1 )

        ! before unpacking the data we received from other ranks, we wait for the transfer
        ! to be completed
        t1 = MPI_wtime()
        call finalize_xfer_mpi(params, isend, irecv)
        call toc( "sync ghosts (finalize_xfer_mpi)", MPI_wtime()-t1 )

        ! process-external ghost points (copy from buffer)
        t1 = MPI_wtime()
        call unpack_ghostlayers_external( params, hvy_block )
        call toc( "sync ghosts (unpack external)", MPI_wtime()-t1 )

        ! edited until here
        call MPI_barrier(WABBIT_COMM, k)
        call abort(197)

    end do ! loop over stages 1,2

    call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )

end subroutine sync_ghosts

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
!              Data_recvCounter, Data_sendCounter, MetaData_recvCounter, MetaData_sendCounter, &
!              count_internal=.false., ncomponents=size(hvy_block,4))
!
!
!         ! reset iMetaData_sendBuffer, but only the parts that will actually be treated.
!         do k = 1, params%number_procs
!             ii0 = sum(MetaData_sendCounter(0:(k-1)-1)) + 1
!             ii1 = ii0 + MetaData_sendCounter(k-1)
!             iMetaData_sendBuffer(ii0:ii1) = -99
!
!             ii0 = sum(MetaData_recvCounter(0:(k-1)-1)) + 1
!             ii1 = ii0 + MetaData_recvCounter(k-1)
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
!                         internalNeighborSyncs( internalNeighbor_pos:internalNeighbor_pos+4-1 ) = (/sender_hvyID, recver_hvyID, neighborhood, level_diff/)
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
!     call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )
!
! end subroutine sync_ghosts_nostages



subroutine send_prepare_external( params, recver_rank, hvy_block, sender_hvyID, &
    recver_hvyID, neighborhood, level_diff, patch_size)
    implicit none

    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in)   :: recver_rank ! zero-based
    integer(kind=ik), intent(in)   :: sender_hvyID, recver_hvyID
    integer(kind=ik), intent(in)   :: neighborhood
    integer(kind=ik), intent(in)   :: level_diff
    integer(kind=ik), intent(in)   :: patch_size
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)

    ! Following are global data used but defined in module_mpi:
    ! res_pre_data, rDara_sendBuffer

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: buffer_offset, buffer_size, data_offset
    integer(kind=ik)   :: patch_ijk(2,3), size_buff(4), nc

    nc = size(hvy_block,4)
    if (size(res_pre_data,4) < nc) then
        size_buff(1) = size( res_pre_data, 1 )
        size_buff(2) = size( res_pre_data, 2 )
        size_buff(3) = size( res_pre_data, 3 )
        size_buff(4) = size( hvy_block, 4 )
        deallocate( res_pre_data )
        allocate( res_pre_data(size_buff(1), size_buff(2), size_buff(3), size_buff(4)) )
    endif

    ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
    ! They depend on the neighbor-relation, level difference and the bounds type.
    ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    if ( level_diff == 0 ) then
        ! simply copy the ghost node layer (no interpolation or restriction here) to a line buffer, which
        ! we will send to our neighbor mpirank
        patch_ijk = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        hvy_block( patch_ijk(1,1):patch_ijk(2,1), patch_ijk(1,2):patch_ijk(2,2), patch_ijk(1,3):patch_ijk(2,3), 1:nc, sender_hvyID) )

    else
        ! up/downsample data first, then flatten to 1D buffer
        patch_ijk = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

        call restrict_predict_data( params, res_pre_data, patch_ijk, neighborhood, level_diff, hvy_block, sender_hvyID )

        patch_ijk = ijkGhosts(:,:, neighborhood, level_diff, RESPRE)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        res_pre_data( patch_ijk(1,1):patch_ijk(2,1), patch_ijk(1,2):patch_ijk(2,2), patch_ijk(1,3):patch_ijk(2,3), 1:nc) )
    end if

    ! now append data, first lets find the positions in the array, +1 to skip count number
    buffer_offset = sum(MetaData_sendCounter(0:recver_rank-1))*4 + sum(Data_sendCounter(0:recver_rank-1)) + recver_rank + 1
    data_offset = buffer_offset + real_pos(recver_rank)

    if (buffer_size /= patch_size) then
        write(*, '("ERROR: I am confused because real buffer_size is not equivalent to theoretical one:", i0, " - ", i0)') buffer_size, patch_size
        call abort(666)
    endif

    ! set metadata, encoded in float, as ints up to 2^53 can be exactly represented with doubles this is not a problem
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 1) = recver_hvyID
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 2) = neighborhood
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 3) = level_diff
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 4) = buffer_size
#ifdef DEV
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 5) = recver_rank  ! receiver rank only for
#endif

    ! set data
    if (buffer_size>0) then
        rData_sendBuffer( data_offset:data_offset+buffer_size-1 ) = line_buffer(1:buffer_size)
    endif

    ! shift positions
    real_pos(recver_rank) = real_pos(recver_rank) + buffer_size
    int_pos(recver_rank) = int_pos(recver_rank) + 1

end subroutine send_prepare_external


subroutine unpack_ghostlayers_external( params, hvy_block )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    integer(kind=ik) :: recver_rank, myrank ! zero-based
    integer(kind=ik) :: recver_hvyID, neighborhood
    integer(kind=ik) :: level_diff, bounds_type, buffer_position, buffer_size, rank_destination
    integer(kind=ik) :: ijk1(2,3), nc, k_patches, buffer_offset, offset_data, n_patches

    myrank = params%rank
    nc = size(hvy_block,4)

    do recver_rank = 0, params%number_procs-1 ! zero-based
        ! skip my own rank, skip ranks that did not send me anything
        if (recver_rank /= myrank .and. Data_recvCounter(recver_rank) /= 0) then
            ! first get overall offset of data
            buffer_offset = sum(MetaData_sendCounter(0:recver_rank-1))*S_META_SEND + sum(Data_recvCounter(0:recver_rank-1)) + recver_rank + 1
            n_patches = int(rData_sendBuffer(buffer_offset))
            ! point to first data
            offset_data = buffer_offset + n_patches*S_META_SEND + 1

            do k_patches = 0, n_patches
                ! extract metadata
                recver_hvyID     = int(iMetaData_recvBuffer(S_META_FULL*k_patches+1))
                neighborhood     = int(iMetaData_recvBuffer(S_META_FULL*k_patches+2))
                level_diff       = int(iMetaData_recvBuffer(S_META_FULL*k_patches+3))
                buffer_size      = int(iMetaData_recvBuffer(S_META_FULL*k_patches+4))
#ifdef DEV
                rank_destination = int(iMetaData_recvBuffer(S_META_FULL*k_patches+5))
                if (rank_destination /= myrank) then
                    write(*,'("rank= ", i0, " dest= ", i0, " patch= ", i0, " recver_rank= ")') myrank, rank_destination, k_patches, recver_rank
                    call abort(7373872, "EXT this data seems to be not mine!")
                endif
#endif
                write(*, '("rank: ", i0, " hvy_ID: ", i0, " neighbourhood: ", i0, " level_diff: ", i0, " buffer_size: ", i0, " rank_dest: ", i0)') myrank, &
                    recver_hvyID, neighborhood, level_diff, buffer_size, rank_destination
 
                ! copy data to line buffer. we now need to extract this to the ghost nodes layer (2D/3D)
                line_buffer(1:buffer_size) = rData_recvBuffer( offset_data : offset_data+buffer_size-1 )


                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
                ! They depend on the neighbor-relation, level difference and the bounds type.
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                call Line2GhostLayer( params, line_buffer, ijkGhosts(:,:, neighborhood, level_diff, RECVER), hvy_block, recver_hvyID )

                ! increase buffer position marker
                offset_data = offset_data + buffer_size
            end do
        end if
    end do

end subroutine unpack_ghostlayers_external



subroutine unpack_ghostlayers_internal( params, hvy_block, count_send_total )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: count_send_total  !< total amount of patches for do loop

    integer(kind=ik) :: k_patch, recver_rank, recver_hvyID, neighborhood
    integer(kind=ik) :: sender_hvyID, level_diff
    integer(kind=ik) :: send_ijk(2,3), recv_ijk(2,3), nc, myrank

    nc = size(hvy_block,4)
    myrank  = params%rank

    do k_patch = 0, count_send_total-1
        ! check if this ghost point patch is addressed from me to me
        recver_rank = internalNeighborSyncs(S_META_FULL*k_patch + 3)
        if (recver_rank == myrank) then
            ! unpack required info:  sender_hvyID, recver_hvyID, neighborhood, level_diff
            sender_hvyID = internalNeighborSyncs(S_META_FULL*k_patch + 1)
            recver_hvyID = internalNeighborSyncs(S_META_FULL*k_patch + 2)
            neighborhood = internalNeighborSyncs(S_META_FULL*k_patch + 4)
            level_diff = internalNeighborSyncs(S_META_FULL*k_patch + 5)

            if ( level_diff == 0 ) then
                ! simply copy from sender block to receiver block (NOTE: both are on the same MPIRANK)
                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
                ! They depend on the neighbor-relation and level difference
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                send_ijk = ijkGhosts(:,:, neighborhood, level_diff, RECVER)
                recv_ijk = ijkGhosts(:,:, neighborhood, level_diff, SENDER)
    
                hvy_block( send_ijk(1,1):send_ijk(2,1), send_ijk(1,2):send_ijk(2,2), send_ijk(1,3):send_ijk(2,3), 1:nc, recver_hvyID ) = &
                hvy_block( recv_ijk(1,1):recv_ijk(2,1), recv_ijk(1,2):recv_ijk(2,2), recv_ijk(1,3):recv_ijk(2,3), 1:nc, sender_hvyID)
            else
                ! interpolation or restriction before inserting
                call restrict_predict_data( params, res_pre_data, ijkGhosts(1:2,1:3, neighborhood, level_diff, SENDER), &
                neighborhood, level_diff, hvy_block, sender_hvyID )
    
                ! copy interpolated / restricted data to ghost nodes layer
                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
                ! They depend on the neighbor-relation, level difference and the bounds type.
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    
                send_ijk = ijkGhosts(:, :, neighborhood, level_diff, RECVER)
                recv_ijk = ijkGhosts(:, :, neighborhood, level_diff, RESPRE)
    
                hvy_block( send_ijk(1,1):send_ijk(2,1), send_ijk(1,2):send_ijk(2,2), send_ijk(1,3):send_ijk(2,3), 1:nc, recver_hvyID ) = &
                res_pre_data( recv_ijk(1,1):recv_ijk(2,1), recv_ijk(1,2):recv_ijk(2,2), recv_ijk(1,3):recv_ijk(2,3), 1:nc)
            end if
        endif
    end do

end subroutine unpack_ghostlayers_internal


subroutine GhostLayer2Line( params, line_buffer, buffer_counter, hvy_data )
    implicit none

    type (type_params), intent(in)   :: params
    !> data buffer
    real(kind=rk), intent(inout)     :: line_buffer(:)
    ! buffer size
    integer(kind=ik), intent(out)    :: buffer_counter
    !> heavy block data, all data fields
    real(kind=rk), intent(inout)     :: hvy_data(:, :, :, :)

    ! loop variable
    integer(kind=ik) :: i, j, k, dF
    ! reset buffer size
    buffer_counter = 0

    ! loop over all components
    do dF = 1, size(hvy_data,4)
        do k = 1, size(hvy_data, 3) ! third dimension, note: for 2D cases k is always 1
            do j = 1, size(hvy_data, 2)
                do i = 1, size(hvy_data, 1)
                    ! increase buffer size
                    buffer_counter = buffer_counter + 1
                    ! write data buffer
                    line_buffer(buffer_counter)   = hvy_data( i, j, k, dF )
                end do
            end do
        end do
    end do

end subroutine GhostLayer2Line

!############################################################################################################

subroutine Line2GhostLayer( params, line_buffer, data_bounds, hvy_block, hvy_id )
    implicit none

    type (type_params), intent(in)  :: params
    !> data buffer
    real(kind=rk), intent(inout)    :: line_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(inout) :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)    :: hvy_id

    ! loop variable
    integer(kind=ik) :: i, j, k, dF, buffer_i

    buffer_i = 1
    ! loop over all components
    do dF = 1, size(hvy_block,4)
        do k = data_bounds(1,3), data_bounds(2,3) ! third dimension, note: for 2D cases k is always 1
            do j = data_bounds(1,2), data_bounds(2,2)
                do i = data_bounds(1,1), data_bounds(2,1)
                    ! write data buffer
                    hvy_block( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )
                    buffer_i = buffer_i + 1
                end do
            end do
        end do
    end do

end subroutine Line2GhostLayer




! subroutine AppendLineToBuffer( iMetaData_sendBuffer, rData_sendBuffer, buffer_size, recver_rank, line_buffer, &
!     hvy_id, neighborhood, level_diff )

!     implicit none

!     !> send buffers, integer and real
!     integer(kind=ik), intent(inout) :: iMetaData_sendBuffer(:)
!     real(kind=rk), intent(inout)    :: rData_sendBuffer(:)
!     ! data buffer size
!     integer(kind=ik), intent(in)    :: buffer_size
!     ! id integer
!     integer(kind=ik), intent(in)    :: recver_rank ! zero-based
!     ! restricted/predicted data buffer
!     real(kind=rk), intent(inout)    :: line_buffer(:)
!     ! data buffer intergers, receiver heavy id, neighborhood id, level difference
!     integer(kind=ik), intent(in)    :: hvy_id, neighborhood, level_diff

!     integer(kind=ik)                :: buffer_position, i0, l0

!     buffer_position = real_pos(recver_rank+1) + 1

!     i0 = sum(Data_sendCounter(0:recver_rank-1)) + buffer_position
!     l0 = sum(MetaData_sendCounter(0:recver_rank-1)) + 0 ! note as int_pos is 1-based, we here use 0-based offset

!     ! real data
!     if (buffer_size>0) then
!         rData_sendBuffer( i0:i0+buffer_size-1 ) = line_buffer(1:buffer_size)
!     endif

!     ! save position of NEXT patch in real buffer
!     real_pos(recver_rank+1) = real_pos(recver_rank+1) + buffer_size

!     ! save: neighbor id, neighborhood, level difference, buffer size
!     iMetaData_sendBuffer( l0+int_pos(recver_rank+1)   ) = hvy_id
!     iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+1 ) = neighborhood
!     iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+2 ) = level_diff
!     iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+3 ) = buffer_position
!     iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+4 ) = buffer_size
!     iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+5 ) = recver_rank ! zero-based ! FIVE

!     ! mark end of buffer with -99, will be overwritten by next element if it is not the last one
!     iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+6 ) = -99

!     int_pos(recver_rank+1) = int_pos(recver_rank+1) + 6 ! FIVE

! end subroutine AppendLineToBuffer


!############################################################################################################

subroutine start_xfer_mpi( params, iMetaData_sendBuffer, rData_sendBuffer, iMetaData_recvBuffer, rData_recvBuffer, isend, irecv )

    implicit none

    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(inout) :: iMetaData_sendBuffer(:)
    integer(kind=ik), intent(inout) :: iMetaData_recvBuffer(:)
    real(kind=rk), intent(inout)    :: rData_sendBuffer(:)
    real(kind=rk), intent(inout)    :: rData_recvBuffer(:)
    integer(kind=ik), intent(out)   :: isend, irecv

    integer(kind=ik) :: length_realBuffer, mpirank_partner
    integer(kind=ik) :: length_intBuffer, rank, ierr, tag
    integer(kind=ik) :: buffer_start

    rank = params%rank
    isend = 0
    irecv = 0

    ! these two arrays are module-global.
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! Note: internal ghost nodes do not pass by the buffer, they only require the integer
    ! buffer above. Therefore, no data (a message of zero length) will be sent to myself
    do mpirank_partner = 0, params%number_procs-1 ! zero based

        if (Data_sendCounter(mpirank_partner) > 0) then
            ! increase communication counter
            isend = isend + 1

            ! the amount of data is pre-computed in prepare_metadata
            ! hence we do know how much data we will send:
            ! 1 - amount of patches to send (metadata of metadata)
            ! metadata*num_of_integers
            ! actual data
            length_realBuffer = Data_sendCounter(mpirank_partner) + MetaData_sendCounter(mpirank_partner)*S_META_SEND + 1

            buffer_start = sum(MetaData_sendCounter(0:mpirank_partner-1))*S_META_SEND + sum(Data_sendCounter(0:mpirank_partner-1)) + mpirank_partner + 1

            ! send data
            tag = rank
            call MPI_Isend( rData_sendBuffer(buffer_start:buffer_start+length_realBuffer-1), length_realBuffer, MPI_REAL8, &
            mpirank_partner, tag, WABBIT_COMM, send_request(isend), ierr)

        end if

        if (Data_recvCounter(mpirank_partner) > 0) then
            ! increase communication counter
            irecv = irecv + 1

            ! the amount of data is pre-computed in prepare_metadata
            ! hence we do know how much data we will receive:
            ! 1 - amount of patches to receive (metadata of metadata)
            ! metadata*num_of_integers
            ! actual data
            length_realBuffer = Data_recvCounter(mpirank_partner) + MetaData_sendCounter(mpirank_partner)*S_META_SEND + 1

            buffer_start = sum(MetaData_sendCounter(0:mpirank_partner-1))*S_META_SEND + sum(Data_recvCounter(0:mpirank_partner-1)) + mpirank_partner + 1

            ! receive data
            tag = mpirank_partner
            call MPI_Irecv( rData_recvBuffer(buffer_start:buffer_start+length_realBuffer-1), length_realBuffer, MPI_REAL8, &
            mpirank_partner, MPI_ANY_TAG, WABBIT_COMM, recv_request(irecv), ierr)
        end if
    end do
end subroutine start_xfer_mpi


subroutine finalize_xfer_mpi(params, isend, irecv)
    implicit none
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in) :: isend, irecv
    integer(kind=ik) :: ierr

    ! synchronize non-blocking communications
    if (isend>0) then
        call MPI_Waitall( isend, send_request(1:isend), MPI_STATUSES_IGNORE, ierr)
    end if
    if (irecv>0) then
        call MPI_Waitall( irecv, recv_request(1:irecv), MPI_STATUSES_IGNORE, ierr)
    end if
end subroutine


! ! returns two lists with numbers of points I send to all other procs and how much I
! ! receive from each proc. note: strictly locally computed, NO MPI comm involved here
! subroutine get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active,&
!      hvy_n, Data_recvCounter, Data_sendCounter, MetaData_recvCounter, MetaData_sendCounter, istage, &
!      count_internal, ncomponents, syncSameLevelOnly)

!     implicit none

!     type (type_params), intent(in)      :: params
!     !> light data array
!     integer(kind=ik), intent(in)        :: lgt_block(:, :)
!     !> heavy data array - neighbor data
!     integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
!     !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_active(:)
!     !> number of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_n, ncomponents, istage
!     integer(kind=ik), intent(inout)     :: Data_recvCounter(0:), Data_sendCounter(0:)
!     integer(kind=ik), intent(inout)     :: MetaData_recvCounter(0:), MetaData_sendCounter(0:)
!     logical, intent(in)                 :: count_internal, syncSameLevelOnly

!     integer(kind=ik) :: k, sender_hvyID, sender_lgtID, myrank, N, neighborhood, recver_rank
!     integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID,level_diff, status, new_size

!     call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
!     N = params%number_blocks

!     Data_recvCounter(:) = 0
!     Data_sendCounter(:) = 0
!     MetaData_recvCounter(:) = 0
!     MetaData_sendCounter(:) = 0

!     do k = 1, hvy_n
!         ! calculate light id
!         sender_hvyID = hvy_active(k)
!         call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )

!         ! loop over all neighbors
!         do neighborhood = 1, size(hvy_neighbor, 2)
!             ! if (skipDiagonalNeighbors) then
!             !     ! Diagonal neighbors (not required for the RHS)
!             !     ! 2D: 5,6,7,8
!             !     ! 3D: 7-18, 19-26, 51-74
!             !     if (dim==2.and.(neighborhood>=5.and.neighborhood<=8)) cycle
!             !     if (dim==3.and.((neighborhood>=7.and.neighborhood<=26).or.(neighborhood>=51))) cycle
!             ! endif

!             ! neighbor exists
!             if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
!                 ! neighbor light data id
!                 recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
!                 ! calculate neighbor rank
!                 call lgt2proc( recver_rank, recver_lgtID, N )
!                 ! neighbor heavy id
!                 call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )

!                 ! define level difference: sender - receiver, so +1 means sender on higher level
!                 ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
!                 ! leveldiff =  0 : sender is same level as recver
!                 ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
!                 level_diff = lgt_block( sender_lgtID, IDX_MESH_LVL ) - lgt_block( recver_lgtID, IDX_MESH_LVL )


!                 if (recver_rank /= myrank .or. count_internal) then
!                     ! it now depends on the stage if we have to sent this data
!                     ! or not.
!                     ! In stage 1, only level_diff = {+1, 0} is treated
!                     ! In stage 2, only level_diff = -1

!                     ! send counter. how much data will I send to other mpiranks?
!                     if ((istage==1 .and. (level_diff==+1 .or. level_diff==0)) .or. (istage==2.and.level_diff==-1)) then
!                         if ((istage==1 .and. syncSameLevelOnly .and. level_diff==0).or.(.not.syncSameLevelOnly)) then
!                         ! why is this RECVER and not sender? Well, complicated. The amount of data on the sender patch
!                         ! is not the same as in the receiver patch, because we interpolate or downsample. We effectively
!                         ! transfer only the data the recver wants - not the extra data.
!                         ijk = ijkGhosts(:, :, neighborhood, level_diff, RECVER)

!                         Data_sendCounter(recver_rank) = Data_sendCounter(recver_rank) + &
!                         (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)

!                         ! counter for integer buffer: for each neighborhood, we send 6 integers as metadata
!                         ! as this is a fixed number it does not depend on the type of neighborhood etc, so
!                         ! technically one would need only one for send/recv
!                         MetaData_sendCounter(recver_rank) = MetaData_sendCounter(recver_rank) + 6 ! FIVE
!                         endif
!                     endif

!                     ! recv counter. how much data will I recv from other mpiranks?
!                     ! This is NOT the same number as before
!                     if ((istage==1 .and. (-1*level_diff==+1 .or. -1*level_diff==0)) .or. (istage==2.and.-1*level_diff==-1)) then
!                         if ((istage==1 .and. syncSameLevelOnly .and. level_diff==0).or.(.not.syncSameLevelOnly)) then
!                         inverse = inverse_neighbor(neighborhood, dim)

!                         ijk = ijkGhosts(:, :, inverse, -1*level_diff, RECVER)

!                         Data_recvCounter(recver_rank) = Data_recvCounter(recver_rank) + &
!                         (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)

!                         MetaData_recvCounter(recver_rank) = MetaData_recvCounter(recver_rank) + 6 ! FIVE
!                         endif
!                     endif
!                 endif
!             end if ! neighbor exists
!         end do ! loop over all possible  neighbors
!     end do ! loop over all heavy active



!     ! NOTE: for the int buffer, we mosly start at some index l0 and then loop unitl
!     ! we find a -99 indicating the end of the buffer. this could be avoided by using
!     ! for instead of while loops in the main routines, but I do not have time now.
!     !
!     ! In the meantime, notice we extent the amount of data by one, to copy the last -99
!     ! to the buffers
!     MetaData_recvCounter(:) = MetaData_recvCounter(:) + 1
!     MetaData_sendCounter(:) = MetaData_sendCounter(:) + 1


!     ! NOTE ACTUAL SEND / RECV DATA IS NEQN
!     Data_recvCounter(:) = Data_recvCounter(:) * ncomponents
!     Data_sendCounter(:) = Data_sendCounter(:) * ncomponents

!     ! NOTE: this feature is against wabbits memory policy: we try to allocate the
!     ! whole memory of the machine on startup, then work with that. however, we have to
!     ! reserve portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
!     ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
!     if (sum(Data_recvCounter) > size(rData_recvBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_recvBuffer,1)*125/100
!         deallocate(rData_recvBuffer)
!         allocate( rData_recvBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
!     endif

!     if (sum(Data_sendCounter) > size(rData_sendBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_sendBuffer,1)*125/100
!         deallocate(rData_sendBuffer)
!         allocate( rData_sendBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999993, "Buffer allocation failed. Not enough memory?")
!     endif
! end subroutine get_my_sendrecv_amount_with_ranks


! This subroutine prepares who sends to whom. This includes:
!    - logic of different synchronization situations
!    - saving of all metadata
!    - computing of buffer sizes for metadata for both sending and receiving
! This is done strictly locally so no MPI needed here
subroutine prepare_ghost_synch_metadata(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, count_send, istage, ncomponents, syncSameLevelOnly)

    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: lgt_block(:, :)     !< light data array
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)   !< heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_active(:)       !< list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n               !< number of active blocks (heavy data)

    integer(kind=ik), intent(in)        :: ncomponents         !< components can vary (for mask for example)
    integer(kind=ik), intent(out)       :: count_send          !< number of ghost patches total to be send, for looping
    !> following are variables that control the logic of where each block sends or receives
    integer(kind=ik), intent(in)        :: istage
    logical, intent(in)                 :: syncSameLevelOnly

    ! Following are global data used but defined in module_mpi:
    !    Data_recvCounter, Data_sendCounter
    !    MetaData_recvCounter, MetaData_sendCounter
    !    internalNeighborsSyncs (possibly needs renaming after this function)

    integer(kind=ik) :: k_block, sender_hvyID, sender_lgtID, myrank, N, neighborhood, recver_rank
    integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID,level_diff, status, new_size

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    N = params%number_blocks

    Data_recvCounter(:) = 0
    Data_sendCounter(:) = 0
    MetaData_sendCounter(:) = 0

    count_send = 0
    do k_block = 1, hvy_n
        ! calculate light id
        sender_hvyID = hvy_active(k_block)
        call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )

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
                level_diff = lgt_block( sender_lgtID, IDX_MESH_LVL ) - lgt_block( recver_lgtID, IDX_MESH_LVL )

                ! Send logic, following cases exist currently, all linked as .or.:
                ! stage=1, level_diff = {+1, 0}
                ! stage=2, level_diff = -1
                ! syncSameLevelOnly, level_diff = 0

                ! send counter. how much data will I send to other mpiranks?
                if  ((istage==1 .and. (level_diff==+1 .or. level_diff==0)) &
                .or. (istage==2.and.level_diff==-1) &
                .or. (istage==1 .and. syncSameLevelOnly .and. level_diff==0)) then
                    ! why is this RECVER and not sender? Because we adjust the data to the requirements of the
                    ! receiver before sending with interpolation or downsampling.
                    ijk = ijkGhosts(:, :, neighborhood, level_diff, RECVER)

                    if (myrank /= recver_rank) then
                        Data_sendCounter(recver_rank) = Data_sendCounter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each neighborhood, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of neighborhood etc
                        ! Increase by one so number of integers can vary
                        MetaData_sendCounter(recver_rank) = MetaData_sendCounter(recver_rank) + 1
                    endif

                    ! now lets save all metadata in one array without caring for rank sorting for now
                    internalNeighborSyncs(S_META_FULL*count_send + 1) = sender_hvyID  ! needed for same-rank sending
                    internalNeighborSyncs(S_META_FULL*count_send + 2) = recver_hvyID
                    internalNeighborSyncs(S_META_FULL*count_send + 3) = recver_rank
                    internalNeighborSyncs(S_META_FULL*count_send + 4) = neighborhood
                    internalNeighborSyncs(S_META_FULL*count_send + 5) = level_diff
                    internalNeighborSyncs(S_META_FULL*count_send + 6) = (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                    
                    count_send = count_send + 1
                endif

                ! Receive logic, following cases exist currently, all linked as .or.:
                ! stage=1, level_diff = {-1, 0}
                ! stage=2, level_diff = +1
                ! syncSameLevelOnly, level_diff = 0

                ! recv counter. how much data will I recv from other mpiranks?
                ! This is NOT the same number as before
                if (myrank /= recver_rank) then  ! only receive from foreign ranks
                    if  ((istage==1 .and. (-1*level_diff==+1 .or. -1*level_diff==0)) &
                    .or. (istage==2.and.-1*level_diff==-1) &
                    .or. (istage==1 .and. syncSameLevelOnly .and. level_diff==0)) then
                        inverse = inverse_neighbor(neighborhood, dim)

                        ijk = ijkGhosts(:, :, inverse, -1*level_diff, RECVER)

                        Data_recvCounter(recver_rank) = Data_recvCounter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                    endif
                endif
            endif ! neighbor exists
        end do ! loop over all possible  neighbors
    end do ! loop over all heavy active


    ! NOTE: this feature is against wabbits memory policy: we try to allocate the
    ! whole memory of the machine on startup, then work with that. however, we have to
    ! reserve portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
    ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
    ! JB: This can only trigger if we change g during the run?
    if (sum(Data_recvCounter) > size(rData_recvBuffer, 1)) then
        ! out-of-memory case: the preallocated buffer is not large enough.
        write(*,'("rank=",i4," OOM for ghost nodes and increases its receive buffer size to 125%")') myrank
        new_size = size(rData_recvBuffer,1)*125/100
        deallocate(rData_recvBuffer)
        allocate( rData_recvBuffer(1:new_size), stat=status )
        if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
    endif

    if (sum(Data_sendCounter) > size(rData_sendBuffer, 1)) then
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
!      hvy_n, Data_recvCounter, Data_sendCounter, MetaData_recvCounter, MetaData_sendCounter, &
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
!     integer(kind=ik), intent(inout)     :: Data_recvCounter(0:), Data_sendCounter(0:)
!     integer(kind=ik), intent(inout)     :: MetaData_recvCounter(0:), MetaData_sendCounter(0:)
!     logical, intent(in)                 :: count_internal
!
!     integer(kind=ik) :: k, sender_hvyID, sender_lgtID, myrank, N, neighborhood, recver_rank
!     integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID,level_diff, status, new_size
!
!     call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
!     N = params%number_blocks
!
!     Data_recvCounter(:) = 0
!     Data_sendCounter(:) = 0
!     MetaData_recvCounter(:) = 0
!     MetaData_sendCounter(:) = 0
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
!                     ijk = ijkGhosts(:, :, neighborhood, level_diff, RECVER)
!
!                     Data_sendCounter(recver_rank) = Data_sendCounter(recver_rank) + &
!                     (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
!
!                     ! recv counter. how much data will I recv from other mpiranks?
!                     ! This is NOT the same number as before
!                     inverse = inverse_neighbor(neighborhood, dim)
!
!                     ijk = ijkGhosts(:, :, inverse, -1*level_diff, RECVER)
!
!                     Data_recvCounter(recver_rank) = Data_recvCounter(recver_rank) + &
!                     (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
!
!                     ! counter for integer buffer: for each neighborhood, we send 6 integers as metadata
!                     ! as this is a fixed number it does not depend on the type of neighborhood etc, so
!                     ! technically one would need only one for send/recv
!                     MetaData_sendCounter(recver_rank) = MetaData_sendCounter(recver_rank) + 6 ! FIVE
!
!                     MetaData_recvCounter(recver_rank) = MetaData_recvCounter(recver_rank) + 6 ! FIVE
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
!     MetaData_recvCounter(:) = MetaData_recvCounter(:) + 1
!     MetaData_sendCounter(:) = MetaData_sendCounter(:) + 1
!
!
!     ! NOTE ACTUAL SEND / RECV DATA IS NEQN
!     Data_recvCounter(:) = Data_recvCounter(:) * ncomponents
!     Data_sendCounter(:) = Data_sendCounter(:) * ncomponents
!
!     ! NOTE: this feature is against wabbits memory policy: we try to allocate the
!     ! whole memory of the machine on startup, then work with that. however, we have to
!     ! reserver portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
!     ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
!     if (sum(Data_recvCounter) > size(rData_recvBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_recvBuffer,1)*125/100
!         deallocate(rData_recvBuffer)
!         allocate( rData_recvBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
!     endif
!
!     if (sum(Data_sendCounter) > size(rData_sendBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_sendBuffer,1)*125/100
!         deallocate(rData_sendBuffer)
!         allocate( rData_sendBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999993, "Buffer allocation failed. Not enough memory?")
!     endif
! end subroutine
