
subroutine synchronize_ghosts_generic_sequence( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! MPI parameter
    integer(kind=ik)   :: myrank, mpisize
    ! grid parameter
    integer(kind=ik)   :: Bs, g, NdF
    ! loop variables
    integer(kind=ik)   :: N, k, neighborhood, level_diff
    ! id integers
    integer(kind=ik)   :: neighbor_lgt_id, neighbor_rank, hvy_id_receiver
    integer(kind=ik)   :: sender_hvy_id, sender_lgt_id

    integer(kind=ik)  :: hvyId_temp   ! just for a  consistency check
    integer(kind=ik)  :: entrySortInRound , currentSortInRound, entrySortInRound_end, iround

    integer(kind=ik) :: ijk(2,3)
    integer(kind=ik) :: bounds_type, istage, istage_buffer(1:4), rounds(1:4), inverse


    if (.not. ghost_nodes_module_ready) then
        ! in order to keep the syntax clean, buffers are module-global and need to be
        ! allocated here.
        call init_ghost_nodes( params )
    endif

    ! if this mpirank has no active blocks, it has nothing to do here.
    if (hvy_n == 0) return

    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes
    NdF   = params%number_data_fields
    N     = params%number_blocks
    myrank  = params%rank
    mpisize = params%number_procs

    ! debug check if hvy_active is sorted
    if (hvy_n>1) then
        hvyId_temp =  hvy_active(1)
        do k = 2, hvy_n
            if  (hvyId_temp> hvy_active(k))  then
                call abort(1212,' hvy_active is not sorted as assumed. Panic!')
            end if
            hvyId_temp = hvy_active(k)
        end do
    end if


    ! Stage I: send the data for entrySortInRound= 2,3,4 and effectively do the rounds 2,3,4
    !          afterwards, the ghost nodes on coarser block, including the redundant nodes, should be fine
    ! Stage II: send the data for entrySortInRound = 1 (interpolation) and do the complete sort in again 1,2,3,4
    !           the data for rounds 2,3,4 is not changed, so it is taken from the buffer for the first stage.
    do istage = 1, 2

        !***************************************************************************
        ! (i) stage initialization
        !***************************************************************************

        ! the (module-global) communication_counter is the number of neighboring relations
        ! this rank has with all other ranks (it is thus an array of mpisize)
        communication_counter(:, istage) = 0
        ! reset integer send buffer position
        int_pos(:, istage) = 2       ! TODO JR why 2? , the first filed contains the size of the XXX
        ! reset first in send buffer position
        int_send_buffer(1, :, istage) = 0
        int_send_buffer(2, :, istage) = -99

        ! compute, locally from the grid info, how much data I recv from and sent to all
        ! other mpiranks. this gives us the start indices of each rank in the send/recv buffers. Note
        ! we do not count our internal nodes (.false. as last argument), as they are not put in the
        ! buffer at any time.
        call get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, &
        recv_counter(:, istage), send_counter(:, istage), INCLUDE_REDUNDANT, .false.)


        !***************************************************************************
        ! (ii) prepare data for sending
        !***************************************************************************


        ! loop over active heavy data. NOTE: hvy_id has a linear correspondance to lgt_id,
        ! i.e.g the ordering in hvy_id and lgt_id is the same. this is very important for the
        ! secondary rule, which is that larger lgt_id wins. this works only if I treat the blocks
        ! in INCREASING lgt_id ordering.
        do k = 1, hvy_n
            ! calculate light id
            sender_hvy_id = hvy_active(k)
            call hvy_id_to_lgt_id( sender_lgt_id, sender_hvy_id, myrank, N )

            ! loop over all neighbors
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists
                if ( hvy_neighbor( sender_hvy_id, neighborhood ) /= -1 ) then

                    !  ----------------------------  determin the core ids and properties of neighbor  ------------------------------
                    ! TODO: check if info available when searching neighbor and store it in hvy_neighbor
                    ! neighbor light data id
                    neighbor_lgt_id = hvy_neighbor( sender_hvy_id, neighborhood )
                    ! calculate neighbor rank
                    call lgt_id_to_proc_rank( neighbor_rank, neighbor_lgt_id, N )
                    ! neighbor heavy id
                    call lgt_id_to_hvy_id( hvy_id_receiver, neighbor_lgt_id, neighbor_rank, N )
                    ! define level difference: sender - receiver, so +1 means sender on higher level
                    level_diff = lgt_block( sender_lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )

                    !  ----------------------------  here decide which values are taken for redundant nodes --------------------------------

                    ! here is the core of the ghost point rules
                    ! primary criterion: (very fine/historic fine) wins over (fine) wins over (same) wins over (coarse)
                    ! secondary criterion: the higher light id wins NOTE: this is an IMPLICIT rule, enforced by loop ordering ONLY.

                    ! comment: the same dominance rules within the ghos nodes are realized by the sequence of filling in the values,
                    ! first coarse then same then finer, always in the sequence of the hvy id the redundant nodes within the ghost nodes and maybe in the
                    ! redundant nodes are written several time, the one following the above rules should win
                    call set_bounds_according_to_ghost_dominance_rules( params, bounds_type, entrySortInRound, &
                    lgt_block, sender_lgt_id, neighbor_lgt_id )

                    if ( istage == 1 ) then
                        if ( level_diff == -1 ) Then
                            ! this block just receives data in this neighborhood relation, but does not send anything
                            communication_counter(neighbor_rank+1, istage) = communication_counter(neighbor_rank+1, istage) + 1
                            cycle
                        endif
                    else
                        ! in stage two leveldiff +1 and 0 are already done
                        if ( level_diff == 0 ) cycle
                        if ( level_diff == +1 ) Then
                            ! this block just receives data in this neighborhood relation, but does not send anything
                            communication_counter(neighbor_rank+1, istage) = communication_counter(neighbor_rank+1, istage) + 1
                            cycle
                        endif
                    endif

                    !----------------------------  pack describing data and node values to send ---------------------------
                    if ( myrank == neighbor_rank ) then
                        !-----------------------------------------------------------
                        ! internal relation (no communication)
                        !-----------------------------------------------------------
                        call send_prepare_internal_neighbor( neighbor_rank+1, istage, sender_hvy_id, hvy_id_receiver, neighborhood, &
                        bounds_type, level_diff, entrySortInRound )

                    else
                        !-----------------------------------------------------------
                        ! external relation (MPI communication)
                        !-----------------------------------------------------------
                        call send_prepare_external_neighbor( params, neighbor_rank+1, istage, hvy_block, communication_counter, &
                        sender_hvy_id, hvy_id_receiver, neighborhood, bounds_type, level_diff, entrySortInRound )

                    end if ! (myrank==neighbor_rank)
                end if ! neighbor exists
            end do ! loop over all possible  neighbors
        end do ! loop over all heavy active

        !***************************************************************************
        ! (iii) transfer part (send/recv)
        !***************************************************************************
        call isend_irecv_data_2( params, int_send_buffer, new_send_buffer, int_recv_buffer, &
        new_recv_buffer, communication_counter, istage )

        !***************************************************************************
        ! (iv) Unpack received data in the ghost node layers
        !***************************************************************************

        ! sort data in, ordering is important to keep dominance rules within ghost nodes.
        ! the redundant nodes owned by two blocks only should be taken care by bounds_type (include_redundant. exclude_redundant )

        if (istage == 1) Then
            entrySortInRound_end = 3
            ! We will perform these unpack rounds in the current stage, in this order...
            rounds = (/2, 3, 4, 0/)
            ! ... and take the date from those buffers
            istage_buffer = (/1, 1, 1, 0/)
        else
            entrySortInRound_end = 4
            ! We will perform these unpack rounds in the current stage, in this order...
            rounds = (/1, 2, 3, 4/)
            ! ... and take the date from those buffers
            istage_buffer = (/2, 1, 1, 1/)
        endif

        do iround = 1,  entrySortInRound_end ! rounds depend on stages, see above
            currentSortInRound = rounds(iround)

            do k = 1, mpisize
                if (k == myrank+1) then
                    !---------------------------------------------------------------
                    ! process-internal ghost points (direct copy)
                    !---------------------------------------------------------------
                    call unpack_all_ghostlayers_currentRound_internal_neighbor( params, k, istage_buffer(iround), &
                    currentSortInRound, hvy_block )

                else
                    !---------------------------------------------------------------
                    ! process-external ghost points (copy from buffer)
                    !---------------------------------------------------------------
                    call unpack_all_ghostlayers_currentRound_external_neighbor( params, k, istage_buffer(iround), &
                    currentSortInRound, hvy_block, communication_counter )

                end if  ! process-internal or external ghost points
            end do ! mpisize
        end do ! currentSortInRound
    end do ! loop over stages 1,2
end subroutine synchronize_ghosts_generic_sequence

!############################################################################################################


subroutine set_bounds_according_to_ghost_dominance_rules( params, bounds_type, entrySortInRound, &
    lgt_block, sender_lgt_id, neighbor_lgt_id )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> output of this function
    integer(kind=ik), intent(out)       :: bounds_type, entrySortInRound
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    integer(kind=ik), intent(in)        :: sender_lgt_id, neighbor_lgt_id

    integer(kind=ik)                    :: level_diff
    logical :: senderHistoricFine, recieverHistoricFine, receiverIsCoarser
    logical :: receiverIsOnSameLevel, lgtIdSenderIsHigher

    ! define level difference: sender - receiver, so +1 means sender on higher level
    level_diff = lgt_block( sender_lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )

    ! the criteria
    senderHistoricFine      = ( lgt_block( sender_lgt_id, params%max_treelevel+2)==11 )
    recieverHistoricFine    = ( lgt_block(neighbor_lgt_id, params%max_treelevel+2)==11 )
    receiverIsCoarser       = ( level_diff>0_ik )
    receiverIsOnSameLevel   = ( level_diff==0_ik )
    lgtIdSenderIsHigher     = ( neighbor_lgt_id < sender_lgt_id )

    bounds_type = EXCLUDE_REDUNDANT  ! default value, may be changed below
    ! in what round in the extraction process will this neighborhood be unpacked?
    entrySortInRound = level_diff + 2  ! now has values 1,2,3 ; is overwritten with 4 if sender is historic fine

    ! here we decide who dominates. would be simple without the historic fine
    if (senderHistoricFine) then
        ! the 4th unpack round is the last one, so setting 4 ensures that historic fine always wins
        entrySortInRound = 4
        if (recieverHistoricFine) then
            if (lgtIdSenderIsHigher)  then
                ! both are historic fine, the redundant nodes are overwritten using secondary criterion
                bounds_type = INCLUDE_REDUNDANT
            end if
        else
            ! receiver not historic fine, so sender always sends redundant nodes, no further
            ! checks on refinement level are required
            bounds_type = INCLUDE_REDUNDANT
        end if

    else  ! sender NOT historic fine,

        ! what about the neighbor/receiver, historic fine?
        if ( .not. recieverHistoricFine) then
            ! neither one is historic fine, so just do the basic rules

            ! first rule, overwrite cosarser ghost nodes
            if (receiverIsCoarser)  then ! receiver is coarser
                bounds_type = INCLUDE_REDUNDANT
            end if

            ! secondary rule: on same level decide using light id
            if (receiverIsOnSameLevel.and.lgtIdSenderIsHigher) then
                bounds_type = INCLUDE_REDUNDANT
            end if
        end if
    end if  ! else  senderHistoricFine

end subroutine


subroutine send_prepare_internal_neighbor( neighbor_rank, istage, sender_hvy_id, hvy_id_receiver, neighborhood, &
    bounds_type, level_diff, entrySortInRound )
    implicit none

    integer(kind=ik), intent(in)   :: neighbor_rank, istage
    integer(kind=ik), intent(in)   :: sender_hvy_id, hvy_id_receiver
    integer(kind=ik), intent(in)   :: neighborhood, bounds_type
    integer(kind=ik), intent(in)   :: level_diff
    integer(kind=ik), intent(in)   :: entrySortInRound

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: level_diff_indicator

    !-----------------------------------------------------------
    ! internal relation (no communication)
    !-----------------------------------------------------------
    ! pack multipe information into one number
    level_diff_indicator =  4096*sender_hvy_id + 256*bounds_type + 16*(level_diff+1) + entrySortInRound

    ! the packing has limitations: if the numbers are too large, it might fail, so check here. TODO
    if (sender_hvy_id.ne.( level_diff_indicator/4096 ) )           call abort(1212,'Packing went wrong: wrong sender_hvy_id !')
    if (modulo( level_diff_indicator/16  , 16 ) .ne. level_diff+1) call abort(1213,'Packing went wrong: wrong leveldiff !')
    if (modulo( level_diff_indicator/256 , 16 ) .ne. bounds_type)  call abort(1214,'Packing went wrong: wrong boundstype !')
    if (modulo( level_diff_indicator, 16 ) .ne. entrySortInRound)  call abort(1215,'Packing went wrong: wrong entrySortInRound !')

    ! we sort of abuse the routine AppendLineToBuffer here. In fact, we only store the integer data
    ! but do not copy the heavy data to te corresponding buffer. In that sense, we only "recall" what
    ! parameters (level_diff, entrySortInRound etc) the neighboring relation has.
    call AppendLineToBuffer( int_send_buffer, new_send_buffer, 0, neighbor_rank, line_buffer, &
    hvy_id_receiver, neighborhood, level_diff_indicator, istage )

end subroutine




subroutine send_prepare_external_neighbor( params, neighbor_rank, istage, hvy_block, communication_counter, sender_hvy_id, &
    hvy_id_receiver, neighborhood, bounds_type, level_diff, entrySortInRound )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)   :: neighbor_rank, istage
    integer(kind=ik), intent(in)   :: sender_hvy_id, hvy_id_receiver
    integer(kind=ik), intent(in)   :: neighborhood
    integer(kind=ik), intent(inout):: bounds_type
    integer(kind=ik), intent(in)   :: level_diff
    integer(kind=ik), intent(in)   :: entrySortInRound
    integer(kind=ik), intent(inout) :: communication_counter(:,:)
    !> heavy data array - block data
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :, :)

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: level_diff_indicator, buffer_size
    integer(kind=ik)   :: ijk1(2,3)

    ! count the number of communications with this mpirank. from that number, the
    ! integer buffer length can be computed while MPI exchanging data
    communication_counter(neighbor_rank, istage) = communication_counter(neighbor_rank, istage) + 1

    ! pack multipe information into one number
    level_diff_indicator = 256*bounds_type + 16*(level_diff+1) + entrySortInRound

    ! we always send INCLUDE_REDUNDANT, but possibly sort in EXCLUDE_REDUNDANT
    ! (if thats in "bounds_type" which is packed above into "level_diff_indicator")
    bounds_type = INCLUDE_REDUNDANT

    ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
    ! They depend on the neighbor-relation, level difference and the bounds type.
    ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    if ( level_diff == 0 ) then

        ! simply copy the ghost node layer (no interpolation or restriction here) to a line buffer, which
        ! we will send to our neighbor mpirank
        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, SENDER)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, sender_hvy_id) )

    else
        ! up/downsample data first, then flatten to 1D buffer
        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, SENDER)

        call restrict_predict_data( params, res_pre_data, ijk1, &
        neighborhood, level_diff, hvy_block, sender_hvy_id )

        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, RESPRE)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        res_pre_data( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) )
    end if

    ! the chunk of data is added to the MPI buffers (preparation for sending)
    call AppendLineToBuffer( int_send_buffer, new_send_buffer, buffer_size, neighbor_rank, line_buffer, &
    hvy_id_receiver, neighborhood, level_diff_indicator, istage )


end subroutine send_prepare_external_neighbor


subroutine unpack_all_ghostlayers_currentRound_external_neighbor( params, neighbor_rank, istage_buffer, &
    currentSortInRound, hvy_block, communication_counter )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: neighbor_rank, istage_buffer
    integer(kind=ik), intent(in)        :: currentSortInRound
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(inout) :: communication_counter(:,:)

    integer(kind=ik) :: l, hvy_id_receiver, neighborhood, level_diff_indicator, entrySortInRound
    integer(kind=ik) :: level_diff, bounds_type, buffer_position, buffer_size
    integer(kind=ik) :: ijk1(2,3), i0

    ! did I recv something from this rank?
    if ( (communication_counter(neighbor_rank, istage_buffer) /= 0) ) then

        l = 2  ! first field is size of data

        do while ( int_recv_buffer(l, neighbor_rank, istage_buffer) /= -99 )
            ! unpack the description of the next data chunk
            hvy_id_receiver = int_recv_buffer(l, neighbor_rank, istage_buffer)
            neighborhood = int_recv_buffer(l+1, neighbor_rank, istage_buffer)

            ! unpack & evaluate level_diff_indicator (contains multiple information, unpack it)
            level_diff_indicator = int_recv_buffer(l+2, neighbor_rank, istage_buffer)
            entrySortInRound = modulo( level_diff_indicator, 16 )

            ! check if this entry is processed in this round, otherwise cycle to next
            if (entrySortInRound /= currentSortInRound ) then
                l = l + 5  ! to read the next entry
                cycle      ! go on to next entry
            end if

            level_diff  = modulo( level_diff_indicator/16 , 16 ) - 1_ik
            bounds_type = modulo( level_diff_indicator/256, 16 )
            buffer_position = int_recv_buffer(l+3, neighbor_rank, istage_buffer)
            buffer_size     = int_recv_buffer(l+4, neighbor_rank, istage_buffer)

            ! copy data to line buffer. we now need to extract this to the ghost nodes layer (2D/3D)
            i0 = sum(recv_counter(0:neighbor_rank-1-1, istage_buffer)) + buffer_position
            line_buffer(1:buffer_size) = new_recv_buffer( i0 : i0+buffer_size-1, istage_buffer )


            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation, level difference and the bounds type.
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.

            if ( bounds_type == EXCLUDE_REDUNDANT ) then

                ! extract INCLUDE_REDUNDANT in tmp block
                call Line2GhostLayer2( params, line_buffer, ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, RECVER), tmp_block )
                ! COPY ONLY_REDUNDANT from block
                ijk1 = ijkGhosts( :, :, neighborhood, level_diff, ONLY_REDUNDANT, RECVER)

                tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) = &
                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver)

                ! copy everything to the block, INCLUDE_REDUNDANT
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :)

            else
                ! for INCLUDE_REDUNDANT, just copy
                call Line2GhostLayer( params, line_buffer, ijkGhosts(:,:, neighborhood, level_diff, bounds_type, RECVER), hvy_block, hvy_id_receiver )
            endif


            ! increase buffer postion marker
            l = l + 5
        end do
    end if

end subroutine unpack_all_ghostlayers_currentRound_external_neighbor

subroutine unpack_all_ghostlayers_currentRound_internal_neighbor( params, neighbor_rank, istage_buffer, &
    currentSortInRound, hvy_block )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: neighbor_rank, istage_buffer
    integer(kind=ik), intent(in)        :: currentSortInRound
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    integer(kind=ik) :: l, hvy_id_receiver, neighborhood, level_diff_indicator, entrySortInRound
    integer(kind=ik) :: sender_hvy_id, level_diff, bounds_type
    integer(kind=ik) :: ijk1(2,3), ijk2(2,3)



    l = 2  ! first field is size of data
    do while ( int_send_buffer(l, neighbor_rank, istage_buffer) /= -99 )
        ! unpack the description of the next data chunk
        ! required info:  sender_hvy_id, hvy_id_receiver, neighborhood, level_diff, bounds_type, entrySortInRound
        hvy_id_receiver = int_send_buffer(l, neighbor_rank, istage_buffer)
        neighborhood = int_send_buffer(l+1, neighbor_rank, istage_buffer)

        ! unpack & evaluate level_diff_indicator (contains multiple information, unpack it)
        level_diff_indicator = int_send_buffer(l+2, neighbor_rank, istage_buffer)
        entrySortInRound = modulo( level_diff_indicator, 16 )

        ! check if this entry is processed in this round, otherwise cycle to next
        if (entrySortInRound /= currentSortInRound) then
            l = l + 5  ! to read the next entry
            cycle      ! go on to next entry
        end if

        level_diff      = modulo( level_diff_indicator/16  , 16 ) - 1_ik
        bounds_type     = modulo( level_diff_indicator/256 , 16 )
        sender_hvy_id   =       ( level_diff_indicator/4096 )

        if ( level_diff == 0 ) then
            ! simply copy from sender block to receiver block (NOTE: both are on the same MPIRANK)
            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation, level difference, and the bounds type.
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.

            if (bounds_type == EXCLUDE_REDUNDANT) then
                ! step (a) into a temporary block, extract the ONLY_REDUNDANT part
                ! step (b) patch the entire INCLUDE_REDUNDANT into the block
                ! step (c) put the data form step (a) back into the block.

                ! ------- step (a) -------
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, RECVER)

                tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) = &
                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2),ijk1(1,3):ijk1(2,3), :, hvy_id_receiver)

                ! ------- step (b) -------
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, RECVER)
                ijk2 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, SENDER)

                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                hvy_block( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :, sender_hvy_id)

                ! ------- step (c) -------
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, RECVER)

                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :)

            else
                ! for INCLUDE_REDUNDANT, just copy the patch and be happy
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, RECVER)
                ijk2 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, SENDER)

                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                hvy_block( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :, sender_hvy_id)
            endif

        else  ! interpolation or restriction before inserting

            call restrict_predict_data( params, res_pre_data, ijkGhosts(1:2,1:3, neighborhood, level_diff, INCLUDE_REDUNDANT, SENDER), &
            neighborhood, level_diff, hvy_block, sender_hvy_id )

            ! copy interpolated / restricted data to ghost nodes layer
            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation, level difference and the bounds type.
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.
            if (bounds_type == EXCLUDE_REDUNDANT) then
                ! step (a) into a temporary block, extract the ONLY_REDUNDANT part
                ! step (b) patch the entire INCLUDE_REDUNDANT into the block
                ! step (c) put the data from step (a) back into the block.
                ! ------- step (a) -------
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, RECVER)

                tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) = &
                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver)

                ! ------- step (b) -------
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, RECVER)
                ijk2 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, RESPRE)

                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                res_pre_data( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :)

                ! ------- step (c) -------
                ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, RECVER)

                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :)

            else
                ijk1 = ijkGhosts(:, :, neighborhood, level_diff, INCLUDE_REDUNDANT, RECVER)
                ijk2 = ijkGhosts(:, :, neighborhood, level_diff, INCLUDE_REDUNDANT, RESPRE)

                hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                res_pre_data( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :)

            endif
        end if

        ! increase buffer postion marker
        l = l + 5
    end do

end subroutine unpack_all_ghostlayers_currentRound_internal_neighbor


!############################################################################################################

subroutine check_unique_origin(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, testOriginFlag)

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! status of the check
    logical,intent(out)                 :: testOriginFlag
    integer(kind=ik)                    :: hvy_id_k, lgt_id

    integer(kind=ik)                    :: i1, i2, iStep, j1, j2, jStep, k1, k2, kStep  , i,j,k, boundaryIndex
    integer(kind=ik)                    :: Bs, g    , levelLocal , levelOrigin , lastRedundantOrigin=-1

    integer(kind=ik)                    :: redundantOriginLgtId, local_hvy_id, localLightId, spaceDirections
    logical                             :: shouldDominate , originHistoricFine, localHistoricFine , originLghtIdHigher

    character(len=128)                  :: fileNameOrigin

    real(kind=rk), allocatable, save    :: hvy_block_test(:, :, :, :, :)

    !---------------------------------------------------------------------------
    ! Unique origin test.
    ! The idea is to fill each block on the grid with its (unique) light ID. Then
    ! we synchronize the ghost nodes, and use the "check_redundant_nodes" routine
    ! to verify that indeed on all blocks, we find the same value in the ghost
    ! nodes. Earlier tests tried comparing the function values themselves but this
    ! proved more difficult.
    !---------------------------------------------------------------------------

    if (.not. allocated(hvy_block_test)) then
        if (params%rank == 0) then
            write(*,*) "---------------ATTENTION---------------ATTENTION---------------ATTENTION---------------"
            write(*,*) " The ghost nodes test is on (check_unique_origin) and allocates a huge"
            write(*,*) " amount of memory. If this is a production run, consider disabeling it!"
            write(*,*) "---------------ATTENTION---------------ATTENTION---------------ATTENTION---------------"
        endif

        allocate( hvy_block_test (size(hvy_block,1),size(hvy_block,2),size(hvy_block,3),size(hvy_block,4),size(hvy_block,5) ) ) !its just a test, so not most time efficient..

        ! this array is global within the MODULE scope
        allocate(hvy_block_test_err(size(hvy_block,1),size(hvy_block,2),size(hvy_block,3),size(hvy_block,4),size(hvy_block,5) ) )
        hvy_block_test_err = 0.0

        ! this array is global within the MODULE scope
        allocate(hvy_block_test_val(size(hvy_block,1),size(hvy_block,2),size(hvy_block,3),size(hvy_block,4),size(hvy_block,5) ) )
        hvy_block_test_val = 0.0

        ! this array is global within the MODULE scope
        allocate(hvy_block_test_interpref(size(hvy_block,1),size(hvy_block,2),size(hvy_block,3),size(hvy_block,4),size(hvy_block,5) ) )
        hvy_block_test_interpref = 0.0
    endif

    ! Fill all blocks with their light ID, for all components and dimensions.
    ! fill all blocks, not just active ones, just to be sure. NOTE: filling all blocks with
    ! a constant is actually NOT a valid grid, as it does not respect the "finer wins over coarser"
    ! rule. One could try to respect this rule here, but that requires a logic similar to the above
    ! and this is dangerous, as we might do the same mistakes twice and conclude "hey, it works!
    ! my test is okay". So we don't do it. ...
    do hvy_id_k = 1, size( hvy_block, 5)
        call hvy_id_to_lgt_id(lgt_id, hvy_id_k, params%rank, params%number_blocks)
        hvy_block_test(:, :, :, :, hvy_id_k) = real(lgt_id, kind=rk)
    end do

    ! ... the consequence of the input field is that our routine will fix the error
    ! on the redundant nodes. However, NOTE this implies that the subsequent
    ! call to check_redundant_nodes_clean gives DIFFERENT values on interpolated
    ! points ...
    call synchronize_ghosts_generic_sequence( params, lgt_block, hvy_block_test, hvy_neighbor, hvy_active, hvy_n )

    ! .. which is why the test here CANNOT suceed on interpolated data. In the compare
    ! routine, we thus skip those points.
    testOriginFlag = .false.
    call check_redundant_nodes_clean( params, lgt_block, hvy_block_test, hvy_neighbor, hvy_active, hvy_n, testOriginFlag)


    if (testOriginFlag ) then
        ! filename is XXX.rank.dat
        call write_real5( hvy_block_test, hvy_active, hvy_n, "hvy_block_test", params%rank )
        call write_real5( hvy_block_test_err, hvy_active, hvy_n, "hvy_block_test_err", params%rank )
        call write_real5( hvy_block_test_val, hvy_active, hvy_n, "hvy_block_test_val", params%rank )
        call write_real5( hvy_block_test_interpref, hvy_active, hvy_n, "hvy_block_test_interpref", params%rank )

        call MPI_barrier(WABBIT_COMM, i1)
        ! call abort(111111,"Same origin of ghost nodes check failed - stopping.")
    endif

    ! ------------------------   check if dominace rules are fulfilled locally, globaly
    ! ------------------------   should follow by uniqueness of origin


    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    if (params%threeD_case ) then
        spaceDirections = 3
    else
        spaceDirections = 2
    end if

    do hvy_id_k = 1, hvy_n
        ! calculate light id
        local_hvy_id =  hvy_active(hvy_id_k)
        call hvy_id_to_lgt_id(localLightId, local_hvy_id  , params%rank , params%number_blocks )

        do boundaryIndex =1,spaceDirections !
            i1      = g + 1
            i2      = g + Bs
            iStep   = 1

            j1      = g + 1
            j2      = g + Bs
            jStep   = 1

            if (params%threeD_case ) then
                k1      = g+1
                k2      = g+ Bs
                kStep   = 1 !Bs -1

            else
                k1      = 1
                k2      = 1
                kStep   = 1
            end if

            select case (boundaryIndex)
            case (1)
                iStep = Bs -1 ! by this i takes the values g+1 and   g+Bs which is the redundant nodes, j, k run ov the full surface
            case (2)
                jStep = Bs -1  ! dito for j ,  in principle same
            case (3)
                kStep = Bs -1  ! dito for k ,  in principle same
            end select

            ! loop over all redundant nodes
            localHistoricFine   = (lgt_block( localLightId , params%max_treelevel+2)==11 )
            levelLocal          =  lgt_block( localLightId  , params%max_treelevel+1 )

            !                level_diff =  - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )

            ! TBD: sequence important for speed?
            do i= i1,i2,iStep
                do j = j1,j2,jStep
                    do k = k1,k2,kStep

                        redundantOriginLgtId    = int( hvy_block_test(i,j,k,1, local_hvy_id ) +0.001 , ik )  ! checking only first field, other should be the same
                        ! am i too optimistic?
                        levelOrigin             = lgt_block( redundantOriginLgtId, params%max_treelevel+1 )

                        originLghtIdHigher      = ( redundantOriginLgtId.gt.localLightId            )

                        if (.not.(redundantOriginLgtId.eq.localLightId) ) then  ! the block owns the redundant nodes, that's locally ok
                            if ( .not.(redundantOriginLgtId.eq.lastRedundantOrigin) )  then ! in many cases we get the same id over and over again..

                                originHistoricFine  =  (lgt_block( redundantOriginLgtId  , params%max_treelevel+2)==11 )
                                levelOrigin         = lgt_block( redundantOriginLgtId  , params%max_treelevel+1 )

                                ! do the check, it should only be there if it dominates the current block
                                shouldDominate = .false. ! overwritten if domination is found by one of the following conditions
                                ! is it finer? , no chekc for historic fine, since if the other is coarser, it cannot be his. fine .
                                if (levelLocal<  levelOrigin )  shouldDominate = .true.
                                ! it is historic fine but i am not
                                if (  originHistoricFine.and.(.not.localHistoricFine)) shouldDominate = .true.
                                ! both historic fine, other has higher lgt id
                                if  ( (originHistoricFine.and.localHistoricFine).and.originLghtIdHigher ) shouldDominate = .true.
                                ! none historic fine, both on same level, check if light id is higher
                                if (    (.not.originHistoricFine).and.(.not.localHistoricFine )&
                                .and.( levelLocal.eq.levelOrigin)&
                                .and.(originLghtIdHigher)            )     shouldDominate = .true.

                                ! TODO fill test in
                                if (.not.shouldDominate) then
                                    ! report error
                                    write (*,*) 'rank',  params%rank , 'hvy_id',  local_hvy_id, 'lgt_id', localLightId, 'level', levelLocal,'hF',localHistoricFine ,'i,j,k',i,j,k, &
                                    ' has origin ', redundantOriginLgtId , 'levelOrigin',  levelOrigin   , 'hF',   originHistoricFine
                                    !,originHistoricFine, localHistoricFine , originLghtIdHigher
                                    write (fileNameOrigin, "(A6,I3.3,A4)") 'origin', params%rank ,'.dat'
                                    call write_real5(hvy_block_test, hvy_active, hvy_n, fileNameOrigin, params%rank  ) ! dubug output with ghost nodes
                                    call abort(44567 ,"should not dominate, who wrote this bloody code, and this useless error message? - stopping.")
                                end if
                                ! ----------
                                lastRedundantOrigin =   redundantOriginLgtId
                            end if
                        end if
                    end do  ! k
                end do ! j
            end do ! i

        end do ! bundary index
    end do  ! active block

end subroutine check_unique_origin

subroutine GhostLayer2Line( params, line_buffer, buffer_counter, hvy_data )
    implicit none

    !> user defined parameter structure
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

    ! loop over all data fields
    do dF = 1, params%number_data_fields
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

    !> user defined parameter structure
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
    ! loop over all data fields
    do dF = 1, params%number_data_fields
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

subroutine Line2GhostLayer2( params, line_buffer, data_bounds, hvy_block )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)  :: params
    !> data buffer
    real(kind=rk), intent(inout)    :: line_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(inout) :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :)

    ! loop variable
    integer(kind=ik) :: i, j, k, dF, buffer_i

    buffer_i = 1
    ! loop over all data fields
    do dF = 1, params%number_data_fields
        do k = data_bounds(1,3), data_bounds(2,3) ! third dimension, note: for 2D cases k is always 1
            do j = data_bounds(1,2), data_bounds(2,2)
                do i = data_bounds(1,1), data_bounds(2,1)
                    ! write data buffer
                    hvy_block( i, j, k, dF ) = line_buffer( buffer_i )
                    buffer_i = buffer_i + 1
                end do
            end do
        end do
    end do

end subroutine Line2GhostLayer2


!############################################################################################################

subroutine AppendLineToBuffer( int_send_buffer, new_send_buffer, buffer_size, neighbor_rank, line_buffer, &
    hvy_id, neighborhood, level_diff, istage )

    implicit none

    !> send buffers, integer and real
    integer(kind=ik), intent(inout)        :: int_send_buffer(:,:,:)
    real(kind=rk), intent(inout)           :: new_send_buffer(:,:)
    ! data buffer size
    integer(kind=ik), intent(in)           :: buffer_size, istage
    ! id integer
    integer(kind=ik), intent(in)           :: neighbor_rank
    ! restricted/predicted data buffer
    real(kind=rk), intent(inout)           :: line_buffer(:)
    ! data buffer intergers, receiver heavy id, neighborhood id, level difference
    integer(kind=ik), intent(in)           :: hvy_id, neighborhood, level_diff

    integer(kind=ik)                       :: buffer_position, i0

    ! fill real buffer
    ! position in real buffer is stored in int buffer
    buffer_position = int_send_buffer( 1, neighbor_rank, istage ) + 1

    i0 = sum(send_counter(0:neighbor_rank-1-1, istage)) + buffer_position

    ! real data
    if (buffer_size>0) then
        new_send_buffer( i0:i0+buffer_size-1, istage  ) = line_buffer(1:buffer_size)
    endif

    ! fill int buffer
    ! sum size of single buffers on first element
    int_send_buffer(1  , neighbor_rank, istage ) = int_send_buffer(1  , neighbor_rank, istage ) + buffer_size

    ! save: neighbor id, neighborhood, level difference, buffer size
    int_send_buffer( int_pos(neighbor_rank, istage),   neighbor_rank, istage ) = hvy_id
    int_send_buffer( int_pos(neighbor_rank, istage)+1, neighbor_rank, istage ) = neighborhood
    int_send_buffer( int_pos(neighbor_rank, istage)+2, neighbor_rank, istage ) = level_diff
    int_send_buffer( int_pos(neighbor_rank, istage)+3, neighbor_rank, istage ) = buffer_position
    int_send_buffer( int_pos(neighbor_rank, istage)+4, neighbor_rank, istage ) = buffer_size
    ! mark end of buffer with -99, will be overwritten by next element if it is not the last one
    int_send_buffer( int_pos(neighbor_rank, istage)+5, neighbor_rank, istage ) = -99

    int_pos(neighbor_rank, istage) = int_pos(neighbor_rank, istage) +5
end subroutine AppendLineToBuffer


!############################################################################################################

subroutine isend_irecv_data_2( params, int_send_buffer, new_send_buffer, int_recv_buffer, new_recv_buffer, &
    communication_counter, istage )

    !---------------------------------------------------------------------------------------------
    ! modules

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> send/receive buffer, integer and real
    integer(kind=ik), intent(inout)       :: int_send_buffer(:,:,:)
    integer(kind=ik), intent(inout)       :: int_recv_buffer(:,:,:)
    real(kind=rk), intent(inout)          :: new_send_buffer(:,:)
    real(kind=rk), intent(inout)          :: new_recv_buffer(:,:)

    integer(kind=ik), intent(inout)       :: communication_counter(:,:)
    integer(kind=ik), intent(in)          :: istage

    ! process rank
    integer(kind=ik)                    :: rank
    ! MPI error variable
    integer(kind=ik)                    :: ierr

    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI request
    integer(kind=ik)                    :: send_request(params%number_procs), recv_request(params%number_procs)

    ! column number of send buffer, column number of receive buffer, real data buffer length
    integer(kind=ik)                    :: length_realBuffer, int_length, mpirank_partner

    ! loop variable
    integer(kind=ik)                    :: k, i, i0


    rank = params%rank

    ! ----------------------------------------------------------------------------------------
    ! first: integer data


    ! reset request arrays
    i = 0
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    do k = 1, params%number_procs ! one-based
        ! communication between proc rank and proc k-1
        if ( communication_counter(k, istage) > 0 ) then
            mpirank_partner = k-1 ! zero based

            ! length of integer buffer
            int_length = 5*communication_counter(k, istage) + 3

            ! increase communication counter
            i = i + 1

            ! send data
            tag = rank
            call MPI_Isend( int_send_buffer(1, k, istage), int_length, MPI_INTEGER4, &
            mpirank_partner, tag, WABBIT_COMM, send_request(i), ierr)

            ! receive data
            tag = mpirank_partner
            call MPI_Irecv( int_recv_buffer(1, k, istage), int_length, MPI_INTEGER4, &
            mpirank_partner, tag, WABBIT_COMM, recv_request(i), ierr)
        end if

    end do


    !> \todo Please check if waiting twice is really necessary
    ! synchronize non-blocking communications
    ! note: single status variable do not work with all compilers, so use MPI_STATUSES_IGNORE instead
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr)
    end if

    ! ----------------------------------------------------------------------------------------
    ! second: real data
    ! reset communication couter
    i = 0

    ! reset request arrays
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL


    do k = 1, params%number_procs
        ! communication between proc rank and proc k-1
        if ( communication_counter(k, istage) > 0 ) then
            mpirank_partner = k-1 ! zero based

            ! increase communication counter
            i = i + 1

            ! real buffer length is stored as the first entry in the integer buffer,
            ! hence we know how much data we'll receive
            length_realBuffer = recv_counter(mpirank_partner, istage)

            i0 = sum(recv_counter(0:mpirank_partner-1, istage)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! receive data
            tag = 1000*(k-1)
            call MPI_Irecv( new_recv_buffer(i0:i0+length_realBuffer-1, istage), length_realBuffer, MPI_REAL8, &
            mpirank_partner, MPI_ANY_TAG, WABBIT_COMM, recv_request(i), ierr)

            ! real buffer length is stored as the first entry in the integer buffer,
            ! hence we know how much data we'll send
            length_realBuffer = send_counter(mpirank_partner, istage)

            i0 = sum(send_counter(0:mpirank_partner-1, istage)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! send data
            tag = 1000*rank
            call MPI_Isend( new_send_buffer(i0:i0+length_realBuffer-1, istage), length_realBuffer, MPI_REAL8, &
            mpirank_partner, tag, WABBIT_COMM, send_request(i), ierr)

        end if
    end do

    ! synchronize non-blocking communications
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr)
    end if

end subroutine isend_irecv_data_2


! returns two lists with numbers of points I send ot all other procs and how much I
! receive from each proc. note: strictly locally computed, NO MPI comm involved here
subroutine get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active,&
     hvy_n, recv_list, send_list, bounds_type, count_internal)

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n, bounds_type
    integer(kind=ik), intent(inout)     :: recv_list(0:), send_list(0:)
    logical, intent(in)                 :: count_internal

    integer(kind=ik) :: k, sender_hvy_id, sender_lgt_id, myrank, N, neighborhood, neighbor_rank
    integer(kind=ik) :: ijk(2,3), inverse, ierr, hvy_id_receiver, neighbor_lgt_id,level_diff

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    N = params%number_blocks

    ! 1) we need to know how much data each mpirank gets. dry-run for counting only.
    ! I think the receiver ranks can at this occasion likewise compute the amount of data they will receive
    ! this would be good since no comm involved.
    recv_list(:) = 0
    send_list(:) = 0

    do k = 1, hvy_n
        ! calculate light id
        sender_hvy_id = hvy_active(k)
        call hvy_id_to_lgt_id( sender_lgt_id, sender_hvy_id, myrank, N )

        ! loop over all neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists
            if ( hvy_neighbor( sender_hvy_id, neighborhood ) /= -1 ) then
                ! neighbor light data id
                neighbor_lgt_id = hvy_neighbor( sender_hvy_id, neighborhood )
                ! calculate neighbor rank
                call lgt_id_to_proc_rank( neighbor_rank, neighbor_lgt_id, N )

                if (neighbor_rank /= myrank .or. count_internal) then
                    ! neighbor heavy id
                    call lgt_id_to_hvy_id( hvy_id_receiver, neighbor_lgt_id, neighbor_rank, N )
                    ! define level difference: sender - receiver, so +1 means sender on higher level
                    level_diff = lgt_block( sender_lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )

                    inverse = inverse_neighbor(neighborhood, dim)

                    ijk = ijkGhosts(:, :, inverse, -1*level_diff, bounds_type, RECVER)

                    recv_list(neighbor_rank) = recv_list(neighbor_rank) + &
                    (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)

                    ijk = ijkGhosts(:, :, neighborhood, level_diff, bounds_type, RECVER)

                    send_list(neighbor_rank) = send_list(neighbor_rank) + &
                    (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
                endif

            end if ! neighbor exists
        end do ! loop over all possible  neighbors
    end do ! loop over all heavy active

    ! NOTE ACTUAL SEND / RECV DATA IS NEQN
    recv_list(:) = recv_list(:) * params%number_data_fields
    send_list(:) = send_list(:) * params%number_data_fields
end subroutine
