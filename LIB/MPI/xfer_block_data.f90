!> This function bundles actual mpi communication for currently three different communications:
!!    1. Ghost patch synchronization
!!    2. Family communication (SC between mother and daughters)
!!    3. Full block transfers
!> All transfers are bundled here so that logic can be switched between them. There are different scenarios:
!!    1. Data resides on same rank and can be copied
!!    2. Data is send / received via buffers, buffer is assumed to fit all data!
!!    3. Data is send / received directly (for whole blocks)
!> Before this step, all metadata and datasizes have to be prepared in send_counter and recv_counter
subroutine xfer_block_data(params, hvy_data, count_send_total)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_data(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: count_send_total             !< total amount of data to send from this rank

    ! Following are global data used but defined in module_mpi:
    !    data_recv_counter, data_send_counter
    !    meta_recv_counter, meta_send_counter
    !    meta_send_all (possibly needs renaming after this function)

    integer(kind=ik)   :: myrank, mpisize, Bs(1:3), buffer_offset
    integer(kind=ik)   :: k, neighborhood, level_diff
    integer(kind=ik)   :: recver_rank, recver_hvyID, patch_size
    integer(kind=ik)   :: sender_hvyID, sender_lgtID
    integer(kind=ik)   :: ijk(2,3), isend, irecv
    real(kind=rk)      :: t0

    Bs      = params%Bs
    myrank  = params%rank
    mpisize = params%number_procs

    !***************************************************************************
    ! (i) prepare data for sending
    !***************************************************************************
    t0 = MPI_wtime()
    do k = 0, mpisize-1
        ! write in total size of entries at beginning of each part and pre-shift real_pos
        buffer_offset = sum(meta_send_counter(0:k-1))*S_META_SEND + sum(data_send_counter(0:k-1)) + k + 1
        rData_sendBuffer(buffer_offset) = meta_send_counter(k)
        real_pos(k) = 1 + meta_send_counter(k)*S_META_SEND  ! offset real data to beginning by metadata

        ! ! test send/receive sizes
        ! write(*, '("Rank ", i0, " to ", i0, " Send ", 4(i0, 1x), "Receive ", 4(i0, 1x), "Patches ", 4(i0, 1x))') myrank, k, data_send_counter, data_recv_counter, meta_send_counter
    end do

    do k = 0, count_send_total-1  ! we do MPI so lets stick to 0-based for a moment
        recver_rank = meta_send_all(S_META_FULL*k + 3)
        ! internal relation needs to prepare nothing and will be handled later
        if ( myrank /= recver_rank ) then
            ! unpack from metadata array what we need
            sender_hvyID = meta_send_all(S_META_FULL*k + 1)
            recver_hvyID = meta_send_all(S_META_FULL*k + 2)
            neighborhood = meta_send_all(S_META_FULL*k + 4)
            level_diff = meta_send_all(S_META_FULL*k + 5)
            patch_size = meta_send_all(S_META_FULL*k + 6)

            ! external relation (MPI communication)
            call send_prepare_external( params, recver_rank, hvy_data, sender_hvyID, recver_hvyID, neighborhood, level_diff, patch_size)
        end if
    end do ! loop over all patches (all heavy_n where neighborhood exists and ghost synching is applied)
    call toc( "xfer_block_data (prepare data)", MPI_wtime()-t0 )

    !***************************************************************************
    ! (ii) transfer part (send/recv)
    !***************************************************************************
    t0 = MPI_wtime()
    call start_xfer_mpi( params, isend, irecv )
    call toc( "xfer_block_data (start_xfer_mpi)", MPI_wtime()-t0 )

    !***************************************************************************
    ! (iii) Unpack received data in the ghost node layers
    !***************************************************************************
    ! process-internal ghost points (direct copy)
    t0 = MPI_wtime()
    call unpack_ghostlayers_internal( params, hvy_data, count_send_total )
    call toc( "xfer_block_data (unpack internal)", MPI_wtime()-t0 )

    ! before unpacking the data we received from other ranks, we wait for the transfer
    ! to be completed
    t0 = MPI_wtime()
    call finalize_xfer_mpi(params, isend, irecv)
    call toc( "xfer_block_data (finalize_xfer_mpi)", MPI_wtime()-t0 )

    ! process-external ghost points (copy from buffer)
    t0 = MPI_wtime()
    call unpack_ghostlayers_external( params, hvy_data )
    call toc( "xfer_block_data (unpack external)", MPI_wtime()-t0 )

end subroutine xfer_block_data


subroutine send_prepare_external( params, recver_rank, hvy_data, sender_hvyID, &
    recver_hvyID, relation, level_diff, patch_size)
    implicit none

    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in)   :: recver_rank ! zero-based
    integer(kind=ik), intent(in)   :: sender_hvyID, recver_hvyID
    integer(kind=ik), intent(in)   :: relation  !< neighborhood 1:74, full block 0, family -1:-8
    integer(kind=ik), intent(in)   :: level_diff
    integer(kind=ik), intent(in)   :: patch_size
    real(kind=rk), intent(inout)   :: hvy_data(:, :, :, :, :)

    ! Following are global data used but defined in module_mpi:
    ! res_pre_data, rDara_sendBuffer

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: buffer_offset, buffer_size, data_offset
    integer(kind=ik)   :: patch_ijk(2,3), size_buff(4), nc

    nc = size(hvy_data,4)
    if (size(res_pre_data,4) < nc) then
        size_buff(1) = size( res_pre_data, 1 )
        size_buff(2) = size( res_pre_data, 2 )
        size_buff(3) = size( res_pre_data, 3 )
        size_buff(4) = size( hvy_data, 4 )
        deallocate( res_pre_data )
        allocate( res_pre_data(size_buff(1), size_buff(2), size_buff(3), size_buff(4)) )
    endif

    ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
    ! They depend on the neighbor-relation, level difference and the bounds type.
    ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    if ( level_diff == 0 .or. relation < 1 ) then
        ! simply copy the ghost node layer (no interpolation or restriction here) to a line buffer, which
        ! we will send to our neighbor mpirank
        patch_ijk = ijkPatches(:,:, relation, level_diff, SENDER)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        hvy_data( patch_ijk(1,1):patch_ijk(2,1), patch_ijk(1,2):patch_ijk(2,2), patch_ijk(1,3):patch_ijk(2,3), 1:nc, sender_hvyID) )

    else
        ! up/downsample data first for neighbor relations, then flatten to 1D buffer
        patch_ijk = ijkPatches(:,:, relation, level_diff, SENDER)

        call restrict_predict_data( params, res_pre_data, patch_ijk, relation, level_diff, hvy_data, sender_hvyID )

        patch_ijk = ijkPatches(:,:, relation, level_diff, RESPRE)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        res_pre_data( patch_ijk(1,1):patch_ijk(2,1), patch_ijk(1,2):patch_ijk(2,2), patch_ijk(1,3):patch_ijk(2,3), 1:nc) )
    end if

    ! now append data, first lets find the positions in the array, +1 to skip count number
    buffer_offset = sum(meta_send_counter(0:recver_rank-1))*S_META_SEND + sum(data_send_counter(0:recver_rank-1)) + recver_rank + 1
    data_offset = buffer_offset + real_pos(recver_rank)

    if (buffer_size /= patch_size) then
        write(*, '("ERROR: I am confused because real buffer_size is not equivalent to theoretical one:", i0, " - ", i0)') buffer_size, patch_size
        call abort(666)
    endif

    ! set metadata, encoded in float, as ints up to 2^53 can be exactly represented with doubles this is not a problem
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 1) = recver_hvyID
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 2) = relation
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 3) = level_diff
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 4) = buffer_size
#ifdef DEV
    rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 5) = recver_rank  ! receiver rank only for DEV
#endif

    ! set data
    if (buffer_size>0) then
        rData_sendBuffer( data_offset:data_offset+buffer_size-1 ) = line_buffer(1:buffer_size)
    endif

    ! shift positions
    real_pos(recver_rank) = real_pos(recver_rank) + buffer_size
    int_pos(recver_rank) = int_pos(recver_rank) + 1

end subroutine send_prepare_external



subroutine unpack_ghostlayers_external( params, hvy_data )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)

    integer(kind=ik) :: sender_rank, myrank ! zero-based
    integer(kind=ik) :: relation  !< neighborhood 1:74, full block 0, family -1:-8
    integer(kind=ik) :: recver_hvyID
    integer(kind=ik) :: level_diff, bounds_type, buffer_position, buffer_size, rank_destination
    integer(kind=ik) :: ijk1(2,3), nc, k_patches, buffer_offset, offset_data, n_patches

    myrank = params%rank
    nc = size(hvy_data,4)

    do sender_rank = 0, params%number_procs-1 ! zero-based
        ! skip my own rank, skip ranks that did not send me anything
        if (sender_rank /= myrank .and. data_recv_counter(sender_rank) /= 0) then
            ! first get overall offset of data
            buffer_offset = sum(meta_recv_counter(0:sender_rank-1))*S_META_SEND + sum(data_recv_counter(0:sender_rank-1)) + sender_rank + 1
            n_patches = int(rData_recvBuffer(buffer_offset))

            ! point to first data
            offset_data = buffer_offset + n_patches*S_META_SEND + 1

            do k_patches = 0, n_patches-1
                ! extract metadata
                recver_hvyID     = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+1))
                relation         = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+2))
                level_diff       = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+3))
                buffer_size      = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+4))
#ifdef DEV
                rank_destination = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+5))
                if (rank_destination /= myrank) then
                    write(*,'("rank= ", i0, " dest= ", i0, " patch= ", i0, " recver_rank= ", i0)') myrank, rank_destination, k_patches, sender_rank
                    call abort(7373872, "EXT this data seems to be not mine!")
                endif
#endif
 
                ! copy data to line buffer. we now need to extract this to the ghost nodes layer (2D/3D)
                line_buffer(1:buffer_size) = rData_recvBuffer( offset_data : offset_data+buffer_size-1 )


                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
                ! They depend on the neighbor-relation, level difference and the bounds type.
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                call Line2GhostLayer( params, line_buffer, ijkPatches(:,:, relation, level_diff, RECVER), hvy_data, recver_hvyID )

                ! increase buffer position marker
                offset_data = offset_data + buffer_size
            end do
        end if
    end do

end subroutine unpack_ghostlayers_external



subroutine unpack_ghostlayers_internal( params, hvy_data, count_send_total )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: count_send_total  !< total amount of patches for do loop

    integer(kind=ik) :: k_patch, recver_rank, recver_hvyID
    integer(kind=ik) :: relation  !< neighborhood 1:74, full block 0, family -1:-8
    integer(kind=ik) :: sender_hvyID, level_diff
    integer(kind=ik) :: send_ijk(2,3), recv_ijk(2,3), nc, myrank

    nc = size(hvy_data,4)
    myrank  = params%rank

    do k_patch = 0, count_send_total-1
        ! check if this ghost point patch is addressed from me to me
        recver_rank = meta_send_all(S_META_FULL*k_patch + 3)
        if (recver_rank == myrank) then
            ! unpack required info:  sender_hvyID, recver_hvyID, neighborhood, level_diff
            sender_hvyID = meta_send_all(S_META_FULL*k_patch + 1)
            recver_hvyID = meta_send_all(S_META_FULL*k_patch + 2)
            relation     = meta_send_all(S_META_FULL*k_patch + 4)
            level_diff   = meta_send_all(S_META_FULL*k_patch + 5)

            if ( level_diff == 0  .or. relation < 1 ) then
                ! simply copy from sender block to receiver block (NOTE: both are on the same MPIRANK)
                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
                ! They depend on the neighbor-relation and level difference
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                send_ijk = ijkPatches(:,:, relation, level_diff, RECVER)
                recv_ijk = ijkPatches(:,:, relation, level_diff, SENDER)
    
                hvy_data( send_ijk(1,1):send_ijk(2,1), send_ijk(1,2):send_ijk(2,2), send_ijk(1,3):send_ijk(2,3), 1:nc, recver_hvyID ) = &
                hvy_data( recv_ijk(1,1):recv_ijk(2,1), recv_ijk(1,2):recv_ijk(2,2), recv_ijk(1,3):recv_ijk(2,3), 1:nc, sender_hvyID)
            else
                ! interpolation or restriction before inserting
                call restrict_predict_data( params, res_pre_data, ijkPatches(1:2,1:3, relation, level_diff, SENDER), &
                relation, level_diff, hvy_data, sender_hvyID )
    
                ! copy interpolated / restricted data to ghost nodes layer
                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
                ! They depend on the neighbor-relation, level difference and the bounds type.
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    
                send_ijk = ijkPatches(:, :, relation, level_diff, RECVER)
                recv_ijk = ijkPatches(:, :, relation, level_diff, RESPRE)
    
                hvy_data( send_ijk(1,1):send_ijk(2,1), send_ijk(1,2):send_ijk(2,2), send_ijk(1,3):send_ijk(2,3), 1:nc, recver_hvyID ) = &
                res_pre_data( recv_ijk(1,1):recv_ijk(2,1), recv_ijk(1,2):recv_ijk(2,2), recv_ijk(1,3):recv_ijk(2,3), 1:nc)
            end if
        endif
    end do

end subroutine unpack_ghostlayers_internal

!############################################################################################################

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



subroutine Line2GhostLayer( params, line_buffer, data_bounds, hvy_data, hvy_id )
    implicit none

    type (type_params), intent(in)  :: params
    !> data buffer
    real(kind=rk), intent(inout)    :: line_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(inout) :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)    :: hvy_data(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)    :: hvy_id

    ! loop variable
    integer(kind=ik) :: i, j, k, dF, buffer_i

    buffer_i = 1
    ! loop over all components
    do dF = 1, size(hvy_data,4)
        do k = data_bounds(1,3), data_bounds(2,3) ! third dimension, note: for 2D cases k is always 1
            do j = data_bounds(1,2), data_bounds(2,2)
                do i = data_bounds(1,1), data_bounds(2,1)
                    ! write data buffer
                    hvy_data( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )
                    buffer_i = buffer_i + 1
                end do
            end do
        end do
    end do

end subroutine Line2GhostLayer

!############################################################################################################

subroutine start_xfer_mpi( params, isend, irecv)

    implicit none

    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(out)   :: isend, irecv

    ! Following are global data used but defined in module_mpi:
    ! integer(kind=ik), intent(inout) :: iMetaData_sendBuffer(:)
    ! integer(kind=ik), intent(inout) :: iMetaData_recvBuffer(:)
    ! real(kind=rk), intent(inout)    :: rData_sendBuffer(:)
    ! real(kind=rk), intent(inout)    :: rData_recvBuffer(:)

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

        if (data_send_counter(mpirank_partner) > 0) then
            ! increase communication counter
            isend = isend + 1

            ! the amount of data is pre-computed in prepare_metadata
            ! hence we do know how much data we will send:
            ! 1 - amount of patches to send (metadata of metadata)
            ! metadata*num_of_integers
            ! actual data
            length_realBuffer = data_send_counter(mpirank_partner) + meta_send_counter(mpirank_partner)*S_META_SEND + 1
            buffer_start = sum(meta_send_counter(0:mpirank_partner-1))*S_META_SEND + sum(data_send_counter(0:mpirank_partner-1)) + mpirank_partner + 1

            ! send data
            tag = rank
            call MPI_Isend( rData_sendBuffer(buffer_start:buffer_start+length_realBuffer-1), length_realBuffer, MPI_REAL8, &
            mpirank_partner, tag, WABBIT_COMM, send_request(isend), ierr)

        end if

        if (data_recv_counter(mpirank_partner) > 0) then
            ! increase communication counter
            irecv = irecv + 1

            ! the amount of data is pre-computed in prepare_metadata
            ! hence we do know how much data we will receive:
            ! 1 - amount of patches to receive (metadata of metadata)
            ! metadata*num_of_integers
            ! actual data
            length_realBuffer = data_recv_counter(mpirank_partner) + meta_recv_counter(mpirank_partner)*S_META_SEND + 1
            buffer_start = sum(meta_recv_counter(0:mpirank_partner-1))*S_META_SEND + sum(data_recv_counter(0:mpirank_partner-1)) + mpirank_partner + 1

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

!############################################################################################################

! This subroutine prepares who sends to whom. This includes:
!    - logic of different synchronization situations
!    - saving of all metadata
!    - computing of buffer sizes for metadata for both sending and receiving
! This is done strictly locally so no MPI needed here
subroutine prepare_update_family_metadata(params, lgt_block, hvy_family, hvy_active, hvy_n, count_send, ncomponents, &
        s_Level, s_M2C, s_C2M, s_M2F, s_F2M)

    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: lgt_block(:, :)     !< light data array
    integer(kind=ik), intent(in)        :: hvy_family(:,:)     !< heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_active(:)       !< list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n               !< number of active blocks (heavy data)

    integer(kind=ik), intent(in)        :: ncomponents         !< components can vary (for mask for example)
    integer(kind=ik), intent(out)       :: count_send          !< number of ghost patches total to be send, for looping

    integer(kind=ik), intent(in), optional  :: s_level              !< Level to synch, if -1 then all levels are synched
    ! specific directions, some are master ones and can work independently of s_level:
    !    1. if s_level == -1 and s_M2C then all mothers will be updated
    !    2. if s_level == -1 and s_M2F then all daughters will be updated
    logical, intent(in), optional  :: s_M2C                         !< Synch from level J   to mother J-1
    logical, intent(in), optional  :: s_C2M                         !< Synch from level J-1 to daughters J
    logical, intent(in), optional  :: s_M2F                         !< Synch from level J   to daughters J+1
    logical, intent(in), optional  :: s_F2M                         !< Synch from level J+1 to mother J
    integer(kind=ik) sLevel
    logical :: SM2C = .false., SC2M = .false., SM2F = .false., SF2M = .false.

    ! Following are global data used but defined in module_mpi:
    !    data_recv_counter, data_send_counter
    !    meta_recv_counter, meta_send_counter
    !    meta_send_all (possibly needs renaming after this function)

    integer(kind=ik) :: k_block, sender_hvyID, sender_lgtID, myrank, mylastdigit, N, family, recver_rank
    integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID, level, level_diff, status, new_size

    sLevel = -1
    if (present(s_Level)) sLevel = s_Level
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

    count_send = 0
    do k_block = 1, hvy_n
        ! calculate light id
        sender_hvyID = hvy_active(k_block)
        call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )
        level = lgt_block( sender_lgtID, IDX_MESH_LVL )
        mylastdigit = tc_get_digit_at_level_b(get_tc(lgt_block(sender_lgtID, IDX_TC_1 : IDX_TC_2)), dim=params%dim, level=level, max_level=params%Jmax)

        ! check if mother exists
        if (hvy_family(sender_hvyID, 1) /= -1) then
            ! mother light data id
            recver_lgtID = hvy_family(sender_hvyID, 1)
            ! calculate mother rank
            call lgt2proc( recver_rank, recver_lgtID, N )
            ! mother heavy id
            call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )

            ! Send logic, following cases exist currently, all linked as .or.:
            ! (sLevel=-1 and M2C) or (level=sLevel and M2C) or (level=sLevel+1 and F2M)

            ! send counter. how much data will I send to my mother?
            if  ((sLevel==-1 .and. sM2C) .or. (level==sLevel .and. sM2C) .or. (level==sLevel+1 .and. sF2M)) then

                ! why is this RECVER and not sender? Because we adjust the data to the requirements of the
                ! receiver before sending with interpolation or downsampling.
                ijk = ijkPatches(:, :, -1 - mylastdigit, +1, RECVER)

                if (myrank /= recver_rank) then
                    data_send_counter(recver_rank) = data_send_counter(recver_rank) + &
                    (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                    ! counter for integer buffer: for each family relation, we send some integers as metadata
                    ! this is a fixed number it does not depend on the type of family relation etc
                    ! Increase by one so number of integers can vary
                    meta_send_counter(recver_rank) = meta_send_counter(recver_rank) + 1
                endif

                ! now lets save all metadata in one array without caring for rank sorting for now
                meta_send_all(S_META_FULL*count_send + 1) = sender_hvyID  ! needed for same-rank sending
                meta_send_all(S_META_FULL*count_send + 2) = recver_hvyID
                meta_send_all(S_META_FULL*count_send + 3) = recver_rank
                meta_send_all(S_META_FULL*count_send + 4) = -1 - mylastdigit
                meta_send_all(S_META_FULL*count_send + 5) = level_diff
                meta_send_all(S_META_FULL*count_send + 6) = (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                
                count_send = count_send + 1
            endif

            ! Receive logic, following cases exist currently, all linked as .or.:
            ! (sLevel=-1 and M2F) or (level=sLevel and C2M) or (level=sLevel+1 and M2F)

            ! recv counter. how much data will I recv from my mother?
            ! This is NOT the same number as before
            if (myrank /= recver_rank) then  ! only receive from foreign ranks
                if  ((sLevel==-1 .and. sM2F) .or. (level==sLevel .and. sC2M) .or. (level==sLevel+1 .and. sM2F)) then
                    ijk = ijkPatches(:, :, -1 - mylastdigit, +1, RECVER)

                    data_recv_counter(recver_rank) = data_recv_counter(recver_rank) + &
                    (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                    ! counter for integer buffer: for each family relation, we send some integers as metadata
                    ! this is a fixed number it does not depend on the type of family relation etc
                    ! Increase by one so number of integers can vary
                    meta_recv_counter(recver_rank) = meta_recv_counter(recver_rank) + 1
                endif
            endif
        endif

        ! check if daughter exists
        if (hvy_family(sender_hvyID, 2+2**params%dim) /= -1) then
            ! loop over all daughters
            do family = 1, 2**params%dim
                ! daughters light data id
                recver_lgtID = hvy_family(sender_hvyID, 1+2**params%dim+family)
                ! calculate dauther rank
                call lgt2proc( recver_rank, recver_lgtID, N )
                ! daughter heavy id
                call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )

                ! Send logic, following cases exist currently, all linked as .or.:
                ! (sLevel=-1 and M2F) or (level=sLevel and M2F) or (level=sLevel-1 and C2M)

                ! send counter. how much data will I send to my mother?
                if  ((sLevel==-1 .and. sM2F) .or. (level==sLevel .and. sM2F) .or. (level==sLevel-1 .and. sC2M)) then

                    ! why is this RECVER and not sender? Because we adjust the data to the requirements of the
                    ! receiver before sending with interpolation or downsampling.
                    ijk = ijkPatches(:, :, -family, +1, RECVER)

                    if (myrank /= recver_rank) then
                        data_send_counter(recver_rank) = data_send_counter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each family relation, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of family relation etc
                        ! Increase by one so number of integers can vary
                        meta_send_counter(recver_rank) = meta_send_counter(recver_rank) + 1
                    endif

                    ! now lets save all metadata in one array without caring for rank sorting for now
                    meta_send_all(S_META_FULL*count_send + 1) = sender_hvyID  ! needed for same-rank sending
                    meta_send_all(S_META_FULL*count_send + 2) = recver_hvyID
                    meta_send_all(S_META_FULL*count_send + 3) = recver_rank
                    meta_send_all(S_META_FULL*count_send + 4) = -family
                    meta_send_all(S_META_FULL*count_send + 5) = level_diff
                    meta_send_all(S_META_FULL*count_send + 6) = (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                    
                    count_send = count_send + 1
                endif

                ! Receive logic, following cases exist currently, all linked as .or.:
                ! (sLevel=-1 and M2C) or (level=sLevel and F2M) or (level=sLevel-1 and M2C)

                ! recv counter. how much data will I recv from my mother?
                ! This is NOT the same number as before
                if (myrank /= recver_rank) then  ! only receive from foreign ranks
                    if  ((sLevel==-1 .and. sM2C) .or. (level==sLevel .and. sF2M) .or. (level==sLevel-1 .and. sM2C)) then
                        ijk = ijkPatches(:, :, -family, +1, RECVER)

                        data_recv_counter(recver_rank) = data_recv_counter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each family relation, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of family relation etc
                        ! Increase by one so number of integers can vary
                        meta_recv_counter(recver_rank) = meta_recv_counter(recver_rank) + 1
                    endif
                endif
            end do ! loop over all possible daughters
        endif ! check if daughters exist
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
end subroutine prepare_update_family_metadata