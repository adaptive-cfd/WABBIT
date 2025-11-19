!> This function bundles actual mpi communication for currently three different communications:
!!    1. Ghost patch synchronization
!!    2. Family communication (SC between mother and daughters)
!!    3. Full block transfers
!> All transfers are bundled here so that logic can be switched between them. There are different scenarios:
!!    1. Data resides on same rank and can be copied
!!    2. Data is send / received via buffers, buffer is assumed to fit all data!
!!    3. Data is send / received directly (for whole blocks)
!> Before this step, all metadata and datasizes have to be prepared in send_counter and recv_counter
subroutine xfer_block_data(params, hvy_data, tree_ID, count_send_total, verbose_check, hvy_tmp, hvy_tmp_usage, REF_FLAG, ignore_Filter)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_data(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                      !< which tree to study
    integer(kind=ik), intent(in)   :: count_send_total             !< total amount of data to send from this rank

    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag

    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    character(len=*), optional, intent(in)  :: hvy_tmp_usage       !< string indicating how to use hvy_tmp
    integer(kind=ik), intent(in), optional :: REF_FLAG             !< Flag in refinement status where we should use data from hvy_tmp
    logical, intent(in), optional  :: ignore_Filter                !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets


    ! Following are global data used but defined in module_mpi:
    !    data_recv_counter, data_send_counter
    !    meta_recv_counter, meta_send_counter
    !    meta_send_all (possibly needs renaming after this function)

    integer(kind=ik)   :: myrank, mpisize, Bs(1:3), buffer_offset, k, ijk(2,3), isend, irecv
    real(kind=rk)      :: t0
    logical :: ignoreFilter

    !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets
    ignoreFilter = .false.
    if (present(ignore_Filter)) ignoreFilter = ignore_Filter


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

        ! out-of memory case.
        ! It is more useful for users to provide a clear error message. The buffer rData_sendBuffer is simply not large enough
        ! to perform the selected transfer. It is okay to have this message checking inside the loop, because its a relatively
        ! small loop only over MPIRANKs.
        if (buffer_offset >= size(rData_sendBuffer,1)) then
            call abort(202505141, "ERROR: Out of memory! [rData_sendBuffer for MPI communication bundling overflow] Try increasing --memory. If no more memory left, possibly increase Ncpu to obtain more.")
        endif

        rData_sendBuffer(buffer_offset) = meta_send_counter(k)
        real_pos(k) = 1 + meta_send_counter(k)*S_META_SEND  ! offset real data to beginning by metadata
    end do

    call send_prepare_external(params, hvy_data, tree_ID, count_send_total, verbose_check, hvy_tmp, hvy_tmp_usage, &
    REF_FLAG, ignore_Filter=ignoreFilter)
    call toc( "xfer_block_data (prepare data)", 70, MPI_wtime()-t0 )

    !***************************************************************************
    ! (ii) transfer part (send/recv)
    !***************************************************************************
    t0 = MPI_wtime()
    call start_xfer_mpi( params, isend, irecv, verbose_check)
    call toc( "xfer_block_data (start_xfer_mpi)", 71, MPI_wtime()-t0 )

    !***************************************************************************
    ! (iii) Unpack received data in the ghost node layers
    !***************************************************************************
    ! process-internal ghost points (direct copy)
    t0 = MPI_wtime()
    call unpack_ghostlayers_internal( params, hvy_data, tree_ID, count_send_total, verbose_check, &
    hvy_tmp, hvy_tmp_usage, REF_FLAG, ignore_Filter=ignoreFilter)
    call toc( "xfer_block_data (unpack internal)", 72, MPI_wtime()-t0 )

    ! before unpacking the data we received from other ranks, we wait for the transfer
    ! to be completed
    t0 = MPI_wtime()
    call finalize_xfer_mpi(params, isend, irecv, verbose_check)
    call toc( "xfer_block_data (finalize_xfer_mpi)", 73, MPI_wtime()-t0 )

    ! process-external ghost points (copy from buffer)
    t0 = MPI_wtime()
    call unpack_ghostlayers_external( params, hvy_data, verbose_check )
    call toc( "xfer_block_data (unpack external)", 74, MPI_wtime()-t0 )

end subroutine xfer_block_data



!> \brief Prepare data to be sent with MPI \n
!! This already applies res_pre so that the recver only has to sort in the data correctly. \n
!! Fills the send buffer
subroutine send_prepare_external( params, hvy_data, tree_ID, count_send_total, verbose_check, hvy_tmp, hvy_tmp_usage, REF_FLAG, ignore_Filter )
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_data(:, :, :, :, :)
    integer(kind=ik), intent(in)   :: tree_ID                      !< which tree to study
    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag
    integer(kind=ik), intent(in)   :: count_send_total             !< total amount of data to send from this rank

    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    character(len=*), optional, intent(in)  :: hvy_tmp_usage       !< string indicating how to use hvy_tmp
    integer(kind=ik), intent(in), optional :: REF_FLAG             !< Flag in refinement status where we should use data from hvy_tmp
    logical, intent(in), optional  :: ignore_Filter                  !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets


    integer(kind=ik)   :: recver_rank, myrank, sender_hvyID, sender_ref, recver_hvyID, lvl_diff, patch_size, recver_lgtID, sender_lgtID
    integer(kind=ik)   :: relation  ! neighborhood 1:56*3, full block 0, family -1:-8

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: buffer_offset, buffer_size, data_offset
    integer(kind=ik)   :: patch_ijk(2,3), size_buff(4), nc, k_patch
    logical            :: use_hvy_TMP_rest_pred, use_hvy_tmp_same_lvl, ignoreFilter

    ignoreFilter = .false.
    if (present(ignore_Filter)) ignoreFilter = ignore_Filter



    ! Following are global data used but defined in module_mpi:
    ! res_pre_data, rDara_sendBuffer

    myrank  = params%rank


    nc = size(hvy_data,4)
    if (size(res_pre_data,4) < nc) then
        size_buff(1) = size( res_pre_data, 1 )
        size_buff(2) = size( res_pre_data, 2 )
        size_buff(3) = size( res_pre_data, 3 )
        size_buff(4) = size( hvy_data, 4 )
        deallocate( res_pre_data )
        allocate( res_pre_data(size_buff(1), size_buff(2), size_buff(3), size_buff(4)) )
    endif

    do k_patch = 0, count_send_total-1  ! we do MPI so lets stick to 0-based for a moment
        recver_rank = meta_send_all(S_META_FULL*k_patch + 4)
        ! internal relation needs to prepare nothing and will be handled later
        if ( myrank == recver_rank ) cycle

        ! unpack from metadata array what we need
        sender_hvyID = meta_send_all(S_META_FULL*k_patch + 1)
        sender_ref = meta_send_all(S_META_FULL*k_patch + 2)
        recver_hvyID = meta_send_all(S_META_FULL*k_patch + 3)
        relation = meta_send_all(S_META_FULL*k_patch + 5)
        lvl_diff = meta_send_all(S_META_FULL*k_patch + 6)
        patch_size = meta_send_all(S_META_FULL*k_patch + 7)

        ! check if to use hvy_temp, this is used in decomposition, there we need to use hvy_tmp if the sender is already decomposed
        use_hvy_TMP_rest_pred = .false.
        use_hvy_tmp_same_lvl = .false.
        if (present(hvy_tmp) .and. present(hvy_tmp_usage)) then
            if (hvy_tmp_usage == "not_flag_is_decomposed") then
                if (.not. present(REF_FLAG)) then
                    call abort(251117, "ERROR: REF_FLAG has to be provided if hvy_tmp_usage is 'not_flag_is_decomposed'")
                endif
                if (sender_ref /= REF_FLAG) then
                    use_hvy_TMP_rest_pred = .true.
                    use_hvy_tmp_same_lvl = .true.
                endif
            elseif (hvy_tmp_usage == "block_is_decomposed") then
                use_hvy_TMP_rest_pred = .true.
            else
                call abort(251117, "ERROR: hvy_tmp_usage string not recognized in send_prepare_external: " // hvy_tmp_usage)
            endif
        endif

        ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
        ! They depend on the neighbor-relation, level difference and the bounds type.
        ! The last index is 1-sender 2-receiver 3-restricted/predicted.
        if ( lvl_diff == 0 .or. relation < 1 ) then
            ! simply copy the ghost node layer (no interpolation or restriction here) to a line buffer, which
            ! we will send to our neighbor mpirank
            patch_ijk = ijkPatches(:,:, relation, SENDER)

            if (use_hvy_tmp_same_lvl) then
                call Patch2Line( params, line_buffer, buffer_size, &
                hvy_tmp( patch_ijk(1,1):patch_ijk(2,1), patch_ijk(1,2):patch_ijk(2,2), patch_ijk(1,3):patch_ijk(2,3), 1:nc, sender_hvyID) )
            else
                call Patch2Line( params, line_buffer, buffer_size, &
                hvy_data( patch_ijk(1,1):patch_ijk(2,1), patch_ijk(1,2):patch_ijk(2,2), patch_ijk(1,3):patch_ijk(2,3), 1:nc, sender_hvyID) )
            endif

        else
            ! up/downsample data first for neighbor relations, then flatten to 1D buffer
            patch_ijk = ijkPatches(:,:, relation, SENDER)

            ! for synching the SC from coarser neighbors, we always need to predict using hvy_tmp as there the original data is wavelet decomposed
            ! there we do not sync from finer neighbors, as full tree is assumed, so restriction should not occur then
            if (use_hvy_TMP_rest_pred) then
                call restrict_predict_data( params, res_pre_data, patch_ijk, relation, lvl_diff, hvy_tmp, nc, sender_hvyID, ignoreFilter )
            else
                call restrict_predict_data( params, res_pre_data, patch_ijk, relation, lvl_diff, hvy_data, nc, sender_hvyID, ignoreFilter)
            endif

            patch_ijk = ijkPatches(:,:, relation, RESPRE)

            call Patch2Line( params, line_buffer, buffer_size, &
            res_pre_data( patch_ijk(1,1):patch_ijk(2,1), patch_ijk(1,2):patch_ijk(2,2), patch_ijk(1,3):patch_ijk(2,3), 1:nc) )
        end if

        ! now append data, first lets find the positions in the array, +1 to skip count number
        buffer_offset = sum(meta_send_counter(0:recver_rank-1))*S_META_SEND + sum(data_send_counter(0:recver_rank-1)) + recver_rank + 1
        data_offset = buffer_offset + real_pos(recver_rank)

        if (buffer_size /= patch_size) then
            write(*, '("ERROR: I am confused because real buffer_size is not equivalent to theoretical one:", i0, " - ", i0)') buffer_size, patch_size
            write(*, '("relation - ", i0, " lvl_diff - ", i0)') relation, lvl_diff
            call abort(666)
        endif

        ! set metadata, encoded in float, as ints up to 2^53 can be exactly represented with doubles this is not a problem
        rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 1) = recver_hvyID
        rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 2) = relation
        rData_sendBuffer(buffer_offset + int_pos(recver_rank)*S_META_SEND + 3) = lvl_diff
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
    end do ! loop over all patches (all heavy_n where neighborhood exists and ghost synching is applied)

end subroutine send_prepare_external



!> Unpack patches which were send via MPI. \n
!! Patch by patch it reads the metadata and then assigns the data to the correct position
subroutine unpack_ghostlayers_external( params, hvy_data, verbose_check )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)

    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag

    integer(kind=ik) :: sender_rank, myrank ! zero-based
    integer(kind=ik) :: relation, i_relation  !< neighborhood 1:56*3, full block 0, family -1:-8
    integer(kind=ik) :: recver_hvyID
    integer(kind=ik) :: lvl_diff, bounds_type, buffer_position, buffer_size, rank_destination
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
                lvl_diff       = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+3))
                buffer_size      = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+4))

                ! write(*, '("Rank ", i0," Recv external ", 4(i0, 1x))') params%rank, recver_hvyID, relation, lvl_diff, buffer_size

                ! we unpack values which where thought of from sender side, now we are at receiver so we need to invert the relation
                call inverse_relation(relation, i_relation)
                lvl_diff = lvl_diff * (-1)
#ifdef DEV
                rank_destination = int(rData_recvBuffer(buffer_offset+S_META_SEND*k_patches+5))
                if (rank_destination /= myrank) then
                    write(*,'("rank= ", i0, " dest_id= ", i0, " rel= ", i0, " lvl_diff=", i0, " buff_size=", i0, " send_rank=", i0, " patch_i=", i0)') &
                        myrank, recver_hvyID, i_relation, lvl_diff, buffer_size, rank_destination, sender_rank, k_patches
                    call abort(7373872, "ERROR: this data seems to be not mine!")
                endif
#endif
 
                ! copy data to line buffer. we now need to extract this to the ghost nodes layer (2D/3D)
                line_buffer(1:buffer_size) = rData_recvBuffer( offset_data : offset_data+buffer_size-1 )


                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
                ! They depend on the neighbor-relation, level difference and the bounds type.
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                call Line2Patch( params, line_buffer, ijkPatches(:,:, i_relation, RECVER), hvy_data, recver_hvyID )

                ! increase buffer position marker
                offset_data = offset_data + buffer_size
            end do
        end if
    end do

end subroutine unpack_ghostlayers_external



!> \brief Unpack all internal patches which do not have to be sent
!! This simply copies from hvy_data, applies res_pre and then copies it back into the correct position
subroutine unpack_ghostlayers_internal( params, hvy_data, tree_ID, count_send_total, verbose_check, hvy_tmp, hvy_tmp_usage, REF_FLAG, ignore_Filter )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID                 !< which tree to study
    integer(kind=ik), intent(in)        :: count_send_total        !< total amount of patches for do loop
    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag
    !> heavy temp data array - block data of preserved values before the WD, used in adapt_tree as neighbours already might be wavelet decomposed
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    character(len=*), optional, intent(in)  :: hvy_tmp_usage       !< string indicating how to use hvy_tmp
    integer(kind=ik), intent(in), optional :: REF_FLAG             !< Flag in refinement status where we should use data from hvy_tmp
    logical, intent(in), optional  :: ignore_Filter                  !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets


    integer(kind=ik) :: k_patch, recver_rank, recver_hvyID, recver_lgtID, sender_lgtID, sender_ref
    integer(kind=ik) :: relation, i_relation  !< neighborhood 1:56*3, full block 0, family -1:-8
    integer(kind=ik) :: sender_hvyID, lvl_diff, i_lvl_diff
    integer(kind=ik) :: recv_ijk(2,3), send_ijk(2,3), buffer_ijk(2,3), nc, myrank
    logical          :: use_hvy_tmp_same_lvl, use_hvy_TMP_rest_pred, ignoreFilter

    ignoreFilter = .false.
    if (present(ignore_Filter)) ignoreFilter = ignore_Filter

    nc = size(hvy_data,4)
    myrank  = params%rank

    do k_patch = 0, count_send_total-1
        ! check if this ghost point patch is addressed from me to me
        recver_rank = meta_send_all(S_META_FULL*k_patch + 4)
        if (recver_rank == myrank) then
            ! unpack required info:  sender_hvyID, recver_hvyID, neighborhood, lvl_diff
            sender_hvyID = meta_send_all(S_META_FULL*k_patch + 1)
            sender_ref =   meta_send_all(S_META_FULL*k_patch + 2)
            recver_hvyID = meta_send_all(S_META_FULL*k_patch + 3)
            relation     = meta_send_all(S_META_FULL*k_patch + 5)
            lvl_diff   = meta_send_all(S_META_FULL*k_patch + 6)

            ! we unpack values which where thought of from sender side, receiver side needs to invert the relation
            call inverse_relation(relation, i_relation)
            i_lvl_diff = lvl_diff * (-1)

            call hvy2lgt( sender_lgtID, sender_hvyID, myrank, params%number_blocks )
            call hvy2lgt( recver_lgtID, recver_hvyID, myrank, params%number_blocks )

            ! check if to use hvy_temp, this is used in decomposition, there we need to use hvy_tmp if the sender is already decomposed
            use_hvy_TMP_rest_pred = .false.
            use_hvy_tmp_same_lvl = .false.
            if (present(hvy_tmp) .and. present(hvy_tmp_usage)) then
                if (hvy_tmp_usage == "not_flag_is_decomposed") then
                    if (.not. present(REF_FLAG)) then
                        call abort(251117, "ERROR: REF_FLAG has to be provided if hvy_tmp_usage is 'not_flag_is_decomposed'")
                    endif
                    if (sender_ref /= REF_FLAG) then
                        use_hvy_TMP_rest_pred = .true.
                        use_hvy_tmp_same_lvl = .true.
                    endif
                elseif (hvy_tmp_usage == "block_is_decomposed") then
                    use_hvy_TMP_rest_pred = .true.
                else
                    call abort(251117, "ERROR: hvy_tmp_usage string not recognized in send_prepare_external: " // hvy_tmp_usage)
                endif
            endif


            if ( lvl_diff == 0  .or. relation < 1 ) then
                ! simply copy from sender block to receiver block (NOTE: both are on the same MPIRANK)
                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
                ! They depend on the neighbor-relation and level difference
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                recv_ijk = ijkPatches(:,:, i_relation, RECVER)
                send_ijk = ijkPatches(:,:,   relation, SENDER)

                if (use_hvy_tmp_same_lvl) then
                    hvy_data( recv_ijk(1,1):recv_ijk(2,1), recv_ijk(1,2):recv_ijk(2,2), recv_ijk(1,3):recv_ijk(2,3), 1:nc, recver_hvyID ) = &
                    hvy_tmp( send_ijk(1,1):send_ijk(2,1), send_ijk(1,2):send_ijk(2,2), send_ijk(1,3):send_ijk(2,3), 1:nc, sender_hvyID)
                else
                    hvy_data( recv_ijk(1,1):recv_ijk(2,1), recv_ijk(1,2):recv_ijk(2,2), recv_ijk(1,3):recv_ijk(2,3), 1:nc, recver_hvyID ) = &
                    hvy_data( send_ijk(1,1):send_ijk(2,1), send_ijk(1,2):send_ijk(2,2), send_ijk(1,3):send_ijk(2,3), 1:nc, sender_hvyID)
                endif

                ! debug
                ! write(*, '(2(A, i0), 1(A, 6(i2, 1x)), 2(A, i0), A, 6(i2, 1x))') &
                !       "Send rel ",   relation, " lvl_d ",   lvl_diff, " ind_s ", send_ijk(:, :), &
                !     ", Recv rel ", i_relation, " lvl_d ", i_lvl_diff, " ind ",   recv_ijk(:, :)
            else
                ! interpolation or restriction before inserting
                ! when hvy_temp is given but not ref_flag the mode is to use hvy_temp for prediction, this is needed for synching the SC from coarser neighbors
                if (use_hvy_TMP_rest_pred) then
                    call restrict_predict_data( params, res_pre_data, ijkPatches(1:2,1:3, relation, SENDER), &
                    relation, lvl_diff, hvy_tmp, nc, sender_hvyID, ignoreFilter )
                else
                    call restrict_predict_data( params, res_pre_data, ijkPatches(1:2,1:3, relation, SENDER), &
                    relation, lvl_diff, hvy_data, nc, sender_hvyID, ignoreFilter )
                endif
    
                ! copy interpolated / restricted data to ghost nodes layer
                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkPatches array (see module_MPI).
                ! They depend on the neighbor-relation, level difference and the bounds type.
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    
                recv_ijk = ijkPatches(:, :, i_relation, RECVER)
                send_ijk = ijkPatches(:, :,   relation, SENDER)
                buffer_ijk = ijkPatches(:, :, relation, RESPRE)

                ! debug
                ! write(*, '(2(A, i0), 2(A, 6(i2, 1x)), 2(A, i0), A, 6(i2, 1x))') &
                !       "Send rel ",   relation, " lvl_d ",   lvl_diff, " ind_s ", send_ijk(:, :), "ind_b ", buffer_ijk(:, :), &
                !     ", Recv rel ", i_relation, " lvl_d ", i_lvl_diff, " ind ",   recv_ijk(:, :)
    
                hvy_data( recv_ijk(1,1):recv_ijk(2,1), recv_ijk(1,2):recv_ijk(2,2), recv_ijk(1,3):recv_ijk(2,3), 1:nc, recver_hvyID ) = &
                res_pre_data( buffer_ijk(1,1):buffer_ijk(2,1), buffer_ijk(1,2):buffer_ijk(2,2), buffer_ijk(1,3):buffer_ijk(2,3), 1:nc)
            end if
        endif
    end do

end subroutine unpack_ghostlayers_internal



!############################################################################################################
!> \brief Convert 2D or 3D patch from indices to 1D to insert it into one-dimensional buffer
subroutine Patch2Line( params, line_buffer, buffer_counter, hvy_data )
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

end subroutine Patch2Line



!> \brief Convert 1D line from one-dimensional buffer back to 2D or 3D patch from indices
subroutine Line2Patch( params, line_buffer, data_bounds, hvy_data, hvy_id )
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

end subroutine Line2Patch



!############################################################################################################
!> \brief Start all communications - this is non-blocking meaning we need to know how much is send and received
subroutine start_xfer_mpi( params, isend, irecv, verbose_check)

    implicit none

    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(out)   :: isend, irecv
    logical, optional, intent(in)   :: verbose_check  ! Output verbose flag

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



!> Just a small wrapper for waitall.
subroutine finalize_xfer_mpi(params, isend, irecv, verbose_check)
    implicit none
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in) :: isend, irecv
    integer(kind=ik) :: ierr
    logical, optional, intent(in)  :: verbose_check  ! Output verbose flag

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
subroutine prepare_update_family_metadata(params, tree_ID, count_send, sync_case, ncomponents, s_level, s_ref, sync_debug_name)

    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: tree_ID             !< which tree to study

    integer(kind=ik), intent(in)        :: ncomponents         !< components can vary (for mask for example)
    integer(kind=ik), intent(out)       :: count_send          !< number of ghost patches total to be send, for looping
    character(len=*)                    :: sync_case           !< String representing which kind of syncing we want to do

    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, dependend on sync case
    integer(kind=ik), intent(in), optional  :: s_level, s_ref
    character(len=*), optional, intent(in) :: sync_debug_name       !< name to be used in debug output files

    ! Following are global data used but defined in module_mpi:
    !    data_recv_counter, data_send_counter
    !    meta_recv_counter, meta_send_counter
    !    meta_send_all (possibly needs renaming after this function)

    integer(kind=ik) :: k_block, myrank, mylastdigit, N, family, ijk(2,3), inverse, ierr, lvl_diff, status, new_size, sync_id
    integer(kind=ik) :: hvy_ID, lgt_ID, level, ref
    integer(kind=ik) :: hvy_ID_n, lgt_ID_n, level_n, ref_n, rank_n
    integer(kind=tsize) :: tc
    logical :: b_send, b_recv
    character(len=clong) :: format_string, string_prepare
    integer(kind=ik)                    :: points_synched
    integer(kind=ik), allocatable, save :: points_synched_list(:)

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
    call family_setup_patches(params, output_to_file=.false.)

    myrank = params%rank
    N = params%number_blocks

    data_recv_counter(:) = 0
    data_send_counter(:) = 0
    meta_send_counter(:) = 0
    meta_recv_counter(:) = 0

    ! I searched for this for 5 hours I think, it was hidden in synchronize_ghost_generic
    int_pos(:) = 0
    real_pos(:) = 0

    ! we have two sync cases for family: update mothers from daughters or update daughters from mothers
    if (index(sync_case, "D2M") > 0) sync_id = 1
    if (index(sync_case, "M2D") > 0) sync_id = 2

    ! now lets treat the special restrictions, set to the second digit, make sure they don't override each other
    if (index(sync_case, "ref") > 0 .and. index(sync_case, "level") > 0) then
        if (.not. present(s_ref)) call abort(251114, "Sync case " // sync_case // " requires s_ref to be set!")
        if (.not. present(s_level)) call abort(251114, "Sync case " // sync_case // " requires s_level to be set!")
        sync_id = sync_id + 10*3
    elseif (index(sync_case, "ref") > 0) then
        if (.not. present(s_ref)) call abort(251114, "Sync case " // sync_case // " requires s_ref to be set!")
        sync_id = sync_id + 10*1
    elseif (index(sync_case, "level") > 0) then
        if (.not. present(s_level)) call abort(251114, "Sync case " // sync_case // " requires s_level to be set!")
        sync_id = sync_id + 10*2
    endif

    if (sync_id == 0) then
        call abort(240805, "No, we don't trade that here! Please ensure the sync_case is valid.")
    endif

    count_send = 0
    do k_block = 1, hvy_n(tree_ID)
        ! calculate light id
        hvy_ID = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_ID, hvy_ID, myrank, N )
        level = lgt_block( lgt_ID, IDX_MESH_LVL )
        tc = get_tc(lgt_block(lgt_ID, IDX_TC_1 : IDX_TC_2))
        mylastdigit = tc_get_digit_at_level_b(tc, dim=params%dim, level=level, max_level=params%Jmax)
        ref = lgt_block( lgt_ID, IDX_REFINE_STS)

        ! check if mother exists
        if (.not. block_is_root(params, hvy_ID)) then
            ! mother light data id
            lgt_ID_n = hvy_family(hvy_ID, 1)
            ! calculate mother rank
            call lgt2proc( rank_n, lgt_ID_n, N )
            ! mother heavy id
            call lgt2hvy( hvy_ID_n, lgt_ID_n, rank_n, N )
            ! mother ref
            ref_n = lgt_block( lgt_ID_n, IDX_REFINE_STS)
            level_n = lgt_block( lgt_ID_n, IDX_MESH_LVL)

            ! Send logic, following cases exist currently, all linked as .or.:
            ! sync_id=+1 -> D2M, send to mother

            ! send counter. how much data will I send to my mother?
            b_send = .false.
            if (mod(sync_id,10) == 1) b_send = .true.
            ! special cases, first ref check and then level check, situated in second digit
            if (sync_id/10 == 1 .and. b_send) b_send = s_ref == ref ! disable sync if I (sender) have wrong ref value
            if (sync_id/10 == 2 .and. b_send) b_send = s_level == level  ! disable sync if I (sender) have wrong level
            if (sync_id/10 == 3 .and. b_send) b_send = s_level == level .and. s_ref == ref  ! disable sync if I (sender) have wrong level or ref value

            if (b_send) then
                ijk = ijkPatches(:, :, -1 - mylastdigit, SENDER)

                if (myrank /= rank_n) then
                    data_send_counter(rank_n) = data_send_counter(rank_n) + &
                    (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                    ! counter for integer buffer: for each family relation, we send some integers as metadata
                    ! this is a fixed number it does not depend on the type of family relation etc
                    ! Increase by one so number of integers can vary
                    meta_send_counter(rank_n) = meta_send_counter(rank_n) + 1

                    ! write(*, '("Send to mother Size, sendID, lastDigit ", 3(i0, 1x), " tc", b32.32)') (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents, sender_hvyID, -1 - mylastdigit, lgt_block(sender_lgtID, IDX_TC_2)
                endif

                ! now lets save all metadata in one array without caring for rank sorting for now
                meta_send_all(S_META_FULL*count_send + 1) = hvy_ID  ! needed for same-rank sending
                meta_send_all(S_META_FULL*count_send + 2) = ref    ! needed for hvy_tmp for adapt_tree
                meta_send_all(S_META_FULL*count_send + 3) = hvy_ID_n
                meta_send_all(S_META_FULL*count_send + 4) = rank_n
                meta_send_all(S_META_FULL*count_send + 5) = -1 - mylastdigit
                meta_send_all(S_META_FULL*count_send + 6) = +1
                meta_send_all(S_META_FULL*count_send + 7) = (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                
                count_send = count_send + 1
            endif

            ! Receive logic, following cases exist currently, all linked as .or.:
            ! sync_id=+2 -> M2D, receive from mother

            ! recv counter. how much data will I recv from my mother?
            ! This is NOT the same number as before
            if (myrank /= rank_n) then  ! only receive from foreign ranks
                b_recv = .false.
                if (mod(sync_id,10) == 2) b_recv = .true.
                ! special cases, first ref check and then level check, situated in second digit
                if (sync_id/10 == 1 .and. b_recv) b_recv = s_ref == ref_n ! disable sync if neighbor (sender) has wrong ref value
                if (sync_id/10 == 2 .and. b_recv) b_recv = s_level == level_n  ! disable sync if neighbor (sender) has wrong level
                if (sync_id/10 == 3 .and. b_recv) b_recv = s_level == level_n .and. s_ref == ref_n  ! disable sync if neighbor (sender) has wrong level or ref value

                if (b_recv) then
                    ijk = ijkPatches(:, :, -1 - mylastdigit -8, RECVER)  ! -8 for lvl_diff=-1

                    data_recv_counter(rank_n) = data_recv_counter(rank_n) + &
                    (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                    ! counter for integer buffer: for each family relation, we send some integers as metadata
                    ! this is a fixed number it does not depend on the type of family relation etc
                    ! Increase by one so number of integers can vary
                    meta_recv_counter(rank_n) = meta_recv_counter(rank_n) + 1

                    ! write(*, '("Recv from mother Size, lastDigit ", 3(i0, 1x), " tc", b32.32)') (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents, -1 - mylastdigit, lgt_block(sender_lgtID, IDX_TC_2)
                endif
            endif
        endif

        ! check if any daughter exists
        ! in order to save blocks we directly set one of the daughters to the mother, so we should check all daughter ids
        if (.not. block_is_leaf(params, hvy_id)) then
            ! loop over all daughters
            do family = 1, 2**params%dim
                ! daughters light data id
                lgt_ID_n = hvy_family(hvy_ID, 1+2**params%dim+family)

                ! if it doesn't exist, cycle, this can be the case when a daughter has been directly assigned to the mother
                if (lgt_ID_n == -1) cycle

                ! calculate dauther rank
                call lgt2proc( rank_n, lgt_ID_n, N )
                ! daughter heavy id
                call lgt2hvy( hvy_ID_n, lgt_ID_n, rank_n, N )
                ! daughter ref
                ref_n = lgt_block( lgt_ID_n, IDX_REFINE_STS)
                level_n = lgt_block( lgt_ID_n, IDX_MESH_LVL)

                ! Send logic, following cases exist currently, all linked as .or.:
                ! sync_id=+2 -> M2D, send to daughters

                ! send counter. how much data will I send to my daughters?
                b_send = .false.
                if (mod(sync_id,10) == 2) b_send = .true.
                ! special cases, first ref check and then level check, situated in second digit
                if (sync_id/10 == 1 .and. b_send) b_send = s_ref == ref ! disable sync if I (sender) have wrong ref value
                if (sync_id/10 == 2 .and. b_send) b_send = s_level == level  ! disable sync if I (sender) have wrong level
                if (sync_id/10 == 3 .and. b_send) b_send = s_level == level .and. s_ref == ref  ! disable sync if I (sender) have wrong level or ref value

                if (b_send) then

                    ijk = ijkPatches(:, :, -family-8, SENDER)  ! -8 for lvl_diff=-1

                    if (myrank /= rank_n) then
                        data_send_counter(rank_n) = data_send_counter(rank_n) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each family relation, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of family relation etc
                        ! Increase by one so number of integers can vary
                        meta_send_counter(rank_n) = meta_send_counter(rank_n) + 1

                        ! write(*, '("Send to children Size, sendID, lastDigit ", 3(i0, 1x), " tc", b32.32)') (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents, sender_hvyID, -family, lgt_block(sender_lgtID, IDX_TC_2)
                    endif

                    ! now lets save all metadata in one array without caring for rank sorting for now
                    meta_send_all(S_META_FULL*count_send + 1) = hvy_ID  ! needed for same-rank sending
                    meta_send_all(S_META_FULL*count_send + 2) = ref    ! needed for hvy_tmp for adapt_tree
                    meta_send_all(S_META_FULL*count_send + 3) = hvy_ID_n
                    meta_send_all(S_META_FULL*count_send + 4) = rank_n
                    meta_send_all(S_META_FULL*count_send + 5) = -family-8
                    meta_send_all(S_META_FULL*count_send + 6) = -1
                    meta_send_all(S_META_FULL*count_send + 7) = (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents
                    
                    count_send = count_send + 1
                endif

                ! Receive logic, following cases exist currently, all linked as .or.:
                ! sync_id=+1 -> D2M, receive from daughters

                ! recv counter. how much data will I recv from my daughters?
                ! This is NOT the same number as before
                if (myrank /= rank_n) then  ! only receive from foreign ranks
                    b_recv = .false.
                    if (mod(sync_id,10) == 1) b_recv = .true.
                    ! special cases, first ref check and then level check, situated in second digit
                    if (sync_id/10 == 1 .and. b_recv) b_recv = s_ref == ref_n ! disable sync if neighbor (sender) has wrong ref value
                    if (sync_id/10 == 2 .and. b_recv) b_recv = s_level == level_n  ! disable sync if neighbor (sender) has wrong level
                    if (sync_id/10 == 3 .and. b_recv) b_recv = s_level == level_n .and. s_ref == ref_n  ! disable sync if neighbor (sender) has wrong level or ref value

                    if (b_recv) then
                        ijk = ijkPatches(:, :, -family, RECVER)

                        data_recv_counter(rank_n) = data_recv_counter(rank_n) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents

                        ! counter for integer buffer: for each family relation, we send some integers as metadata
                        ! this is a fixed number it does not depend on the type of family relation etc
                        ! Increase by one so number of integers can vary
                        meta_recv_counter(rank_n) = meta_recv_counter(rank_n) + 1

                        ! write(*, '("Recv from children Size, lastDigit ", 3(i0, 1x), " tc", b32.32)') (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1) * ncomponents, -family, lgt_block(sender_lgtID, IDX_TC_2)
                    endif
                endif
            end do ! loop over all possible daughters
        endif ! check if daughters exist
    end do ! loop over all heavy active

    if (params%debug_sync .and. present(sync_debug_name)) then
        if (.not. allocated(points_synched_list)) allocate(points_synched_list(1:params%number_procs))

        ! debug send volume
        call MPI_GATHER(sum(data_send_counter), 1, MPI_INTEGER, points_synched_list, 1, MPI_INTEGER, 0, WABBIT_COMM, ierr)
        if (params%rank == 0) then
            open(unit=99, file=trim("debug_sync.csv"), status="unknown", position="append")
            ! 0 is send as stage
            if (params%number_procs == 1) then
                write(99,'(A, i0, ",", i0)') trim(adjustl(sync_debug_name))//",send,", 0, points_synched_list(1)
            else
                write(format_string, '("(A, i0,",i0,"("","",i0))")') params%number_procs
                write(99,format_string) trim(adjustl(sync_debug_name))//",send,", 0, points_synched_list(1), points_synched_list(2:params%number_procs)
            endif
            close(99)
        endif
        ! debug recv volume
        call MPI_GATHER(sum(data_recv_counter), 1, MPI_INTEGER, points_synched_list, 1, MPI_INTEGER, 0, WABBIT_COMM, ierr)
        if (params%rank == 0) then
            open(unit=99, file=trim("debug_sync.csv"), status="unknown", position="append")
            ! 0 is send as stage
            if (params%number_procs == 1) then
                write(99,'(A, i0, ",", i0)') trim(adjustl(sync_debug_name))//",recv,", 0, points_synched_list(1)
            else
                write(format_string, '("(A, i0,",i0,"("","",i0))")') params%number_procs
                write(99,format_string) trim(adjustl(sync_debug_name))//",recv,", 0, points_synched_list(1), points_synched_list(2:params%number_procs)
            endif
            close(99)
        endif
    endif

    ! NOTE: this feature is against wabbits memory policy: we try to allocate the
    ! whole memory of the machine on startup, then work with that. however, we have to
    ! reserve portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
    ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
    ! JB: This can only trigger if we change g during the run?, and why increase by 125%? What if that is not enough?
    if (sum(data_recv_counter) + sum(meta_recv_counter)*S_META_SEND + params%number_procs > size(rData_recvBuffer, 1)) then
        ! out-of-memory case: the preallocated buffer is not large enough.
        write(*,'("rank=",i4," OOM for patches and increases its receive buffer size to 125%")') myrank
        new_size = size(rData_recvBuffer,1)*125/100
        deallocate(rData_recvBuffer)
        allocate( rData_recvBuffer(1:new_size), stat=status )
        if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
    endif

    if (sum(data_send_counter) + sum(meta_send_counter)*S_META_SEND + params%number_procs > size(rData_sendBuffer, 1)) then
        ! out-of-memory case: the preallocated buffer is not large enough.
        write(*,'("rank=",i4," OOM for patches and increases its send buffer size to 125%")') myrank
        new_size = size(rData_sendBuffer,1)*125/100
        deallocate(rData_sendBuffer)
        allocate( rData_sendBuffer(1:new_size), stat=status )
        if (status /= 0) call abort(999993, "Buffer allocation failed. Not enough memory?")
    endif
end subroutine prepare_update_family_metadata



!> \brief In-place move from Mallat wavelet decomposition the SC patch to a position and wipe the rest
!> When copying patches in coarsening from wavelet decomposed values into mother patches
!! we might want to do this in-place in order to create less blocks. This is exactly what
!! this function does
!  Why is this function in MPI? In order to access ijkPatches as this already has all indices available
subroutine move_mallat_patch_block(params, hvy_block, hvy_ID, digit)

    implicit none

    type (type_params), intent(in)      :: params                       !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !< heavy data array - block data
    integer(kind=ik), intent(in)        :: hvy_ID                       !< ID of the block to be concerned
    integer(kind=ik), intent(in)        :: digit                        !< to which position should the patch go

    integer(kind=ik)                    :: nc, ijk_s(2,3), ijk_r(2,3)
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: tmp_wd  ! used for data juggling

    nc = size(hvy_block, 4)
    if (allocated(tmp_wd)) then
        if (size(tmp_wd, 4) < nc) deallocate(tmp_wd)
    endif
    if (.not. allocated(tmp_wd)) allocate(tmp_wd(1:size(hvy_block, 1), 1:size(hvy_block, 2), 1:size(hvy_block, 3), 1:nc) )

    ijk_r = ijkPatches(:, :, -1 - digit -8, RECVER)  ! daughter receives from coarser level mother, -8 for lvl_diff=-1
    ijk_s = ijkPatches(:, :, -1 - digit   , SENDER)  ! mother sends to finer level daughter

    ! init tmp values
    tmp_wd(:, :, :, :) = 0.0_rk
    ! move patch
    tmp_wd(ijk_r(1,1):ijk_r(2,1), ijk_r(1,2):ijk_r(2,2), ijk_r(1,3):ijk_r(2,3), 1:nc) = &
        hvy_block(ijk_s(1,1):ijk_s(2,1), ijk_s(1,2):ijk_s(2,2), ijk_s(1,3):ijk_s(2,3), 1:nc, hvy_ID)
    hvy_block(:,:,:,1:nc, hvy_ID) = tmp_wd(:,:,:,1:nc)

end subroutine move_mallat_patch_block