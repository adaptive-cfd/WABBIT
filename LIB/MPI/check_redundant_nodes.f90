subroutine check_redundant_nodes_clean( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, stop_status )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

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
    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

    ! MPI parameter
    integer(kind=ik)                    :: myrank
    ! loop variables
    integer(kind=ik)                    :: N, k, l, dF, neighborhood, neighbor_num, level_diff
    ! id integers
    integer(kind=ik)                    :: lgt_id, neighbor_lgt_id, neighbor_rank, hvy_id
    ! type of data bounds
    ! exclude_redundant, include_redundant, only_redundant
    integer(kind=ik)                    :: data_bounds_type
    integer(kind=ik), dimension(2,3)    :: data_bounds, data_bounds2
    ! data buffer size
    integer(kind=ik)                    :: buffer_size, buffer_position
    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! number of datafields
    integer(kind=ik)                    :: NdF, id_Friend
    logical                             :: test2

 !---------------------------------------------------------------------------------------------
! variables initialization

    if (.not. ghost_nodes_module_ready) then
        call init_ghost_nodes( params )
    endif

    ! if this mpirank has no active blocks, it has nothing to do here.
    if (hvy_n == 0) return

    ! nodes test
    ! exclude_redundant, include_redundant, only_redundant
    data_bounds_type = include_redundant
    ! reset status
    stop_status = .false.

    Bs       = params%number_block_nodes
    g        = params%number_ghost_nodes
    NdF      = params%number_data_fields
    N        = params%number_blocks
    myrank   = params%rank
    neighbor_num = size(hvy_neighbor, 2)


    ! the (module-global) communication_counter is the number of neighboring relations
    ! this rank has with all other ranks (it is thus an array of number_procs)
    communication_counter(1:N_friends) = 0_ik
    ! the friends-relation is updated in every call to this routine.
    ! in the beginning all slots are free
    N_friends_used = 0
    mpirank2friend(1:params%number_procs) = -100
    friend2mpirank(1:N_friends) = -100
    ! for technical reasons, I set me as my first friend here. this solves problems
    ! if I have no blocks (and thus do not appear in the friends list)
    N_friends_used = 1
    mpirank2friend(myrank+1) = 1
    friend2mpirank(1) = myrank+1


    ! ATTENTION: if you change something here, recall to do the same in reallocate_buffers
    ! new, freshly allocated "friends" slots require consistent initialization
    ! reset integer send buffer position
    int_pos = 2
    ! reset first in send buffer position
    int_send_buffer( 1, : ) = 0
    int_send_buffer( 2, : ) = -99


    do k = 1, hvy_n
        do neighborhood = 1, neighbor_num
            ! neighbor exists
            if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                neighbor_lgt_id = hvy_neighbor( hvy_active(k), neighborhood )
                call lgt_id_to_proc_rank( neighbor_rank, neighbor_lgt_id, N )
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), myrank, N )
                call lgt_id_to_hvy_id( hvy_id, neighbor_lgt_id, neighbor_rank, N )
                ! define leveldiff: sender - receiver, so +1 means sender on higher level. sender is active block (me)
                level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )

                ! 1 = sender
                data_bounds = ijkGhosts(:,:, neighborhood, level_diff, data_bounds_type, 1)

                if ( level_diff == 0 ) then
                    !-----------------------------------------------------------
                    ! same level
                    !-----------------------------------------------------------
                    call GhostLayer2Line( params, line_buffer, buffer_size, &
                    hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), data_bounds(1,3):data_bounds(2,3), :, hvy_active(k)) )
                else
                    !-----------------------------------------------------------
                    ! different level
                    !-----------------------------------------------------------
                    ! interpoliere daten
                    call restrict_predict_data( params, res_pre_data, data_bounds, neighborhood, level_diff, data_bounds_type, hvy_block, hvy_active(k))

                    ! 3: restrict-predict
                    data_bounds2 = ijkGhosts(1:2, 1:3, neighborhood, level_diff, data_bounds_type, 3)

                    ! lese daten, verwende interpolierte daten
                    call GhostLayer2Line( params, line_buffer, buffer_size, res_pre_data( data_bounds2(1,1):data_bounds2(2,1), &
                    data_bounds2(1,2):data_bounds2(2,2), data_bounds2(1,3):data_bounds2(2,3),:) )
                end if

                call get_friend_id_for_mpirank( params, neighbor_rank, id_Friend )

                ! first: fill com matrix, count number of communication to neighboring process, needed for int buffer length
                communication_counter(id_Friend) = communication_counter(id_Friend) + 1
                ! active block send data to its neighbor block
                ! fill int/real buffer
                call AppendLineToBuffer( int_send_buffer, real_send_buffer, buffer_size, id_Friend, line_buffer, &
                hvy_id, neighborhood, level_diff )

            end if
        end do
    end do

    !***********************************************************************
    ! transfer part (send/recv)
    !***********************************************************************
    ! pretend that no communication with myself takes place, in order to skip the
    ! MPI transfer in the following routine. NOTE: you can also skip this step and just have isend_irecv_data_2
    ! transfer the data, in which case you should skip the copy part directly after isend_irecv_data_2
    communication_counter( mpirank2friend(myrank+1) ) = 0

    ! send/receive data
    call isend_irecv_data_2( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, communication_counter  )

    ! copy internal buffer (BAD! Performance penalty!)
    int_receive_buffer( 1:int_pos(mpirank2friend(myrank+1)), mpirank2friend(myrank+1) ) = &
    int_send_buffer( 1:int_pos(mpirank2friend(myrank+1)), mpirank2friend(myrank+1) )
    real_receive_buffer( 1:int_receive_buffer(1,mpirank2friend(myrank+1)), mpirank2friend(myrank+1) ) = &
    real_send_buffer( 1:int_receive_buffer(1,mpirank2friend(myrank+1)), mpirank2friend(myrank+1) )

    ! change communication_counter, equired to trigger buffer unpacking in last step
    communication_counter(mpirank2friend(myrank+1)) = 1

    !***********************************************************************
    ! Unpack received data and compare with ghost nodes data
    !***********************************************************************
    ! sortiere den real buffer ein
    ! loop over all friends
    do k = 1, N_friends_used
        if ( communication_counter(k) /= 0 ) then
            ! first element in int buffer is real buffer size
            l = 2
            ! -99 marks end of data
            do while ( int_receive_buffer(l, k) /= -99 )

                hvy_id          = int_receive_buffer(l, k)
                neighborhood    = int_receive_buffer(l+1, k)
                level_diff      = int_receive_buffer(l+2, k)
                buffer_position = int_receive_buffer(l+3, k)
                buffer_size     = int_receive_buffer(l+4, k)
                line_buffer(1:buffer_size) = real_receive_buffer( buffer_position : buffer_position-1 + buffer_size, k )

                ! data bounds (2-recv)
                data_bounds = ijkGhosts(:,:, neighborhood, level_diff, data_bounds_type, 2)

                ! compare data
                call hvy_id_to_lgt_id( lgt_id, hvy_id, myrank, N )
                call compare_hvy_data( params, line_buffer, data_bounds, hvy_block, hvy_id, stop_status, level_diff, &
                lgt_block(lgt_id, params%max_treelevel+2), treecode2int( lgt_block(lgt_id, 1:params%max_treelevel) ) )

                l = l + 5
            end do
        end if
    end do

    ! MPI sync the stop status
    test2 = stop_status
    call MPI_Allreduce(test2, stop_status, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, k )
end subroutine check_redundant_nodes_clean



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
    integer(kind=ik)   :: myrank, mpisize, irank
    ! grid parameter
    integer(kind=ik)   :: Bs, g, NdF
    ! loop variables
    integer(kind=ik)   :: N, k, dF, neighborhood, invert_neighborhood, level_diff, l, levelsToSortIn
    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: level_diff_indicator
    ! id integers
    integer(kind=ik)   :: neighbor_lgt_id, neighbor_rank, hvy_id_receiver
    integer(kind=ik)   :: sender_hvy_id, sender_lgt_id
    ! data buffer size
    integer(kind=ik)  :: buffer_size, buffer_position, data_bounds(1:2,1:3)

    integer(kind=ik)  :: hvyId_temp   ! just for a  consistency check
    integer(kind=ik)  :: entrySortInRound , currentSortInRound, entrySortInRound_start, entrySortInRound_end, iround
    integer(kind=ik)  :: ijk1(2,3), ijk2(2,3)

    ! Note each mpirank usually communicates only with a subset of all existing mpiranks.
    ! such a patner is called "friend"
    integer(kind=ik)  :: id_Friend

    integer :: bounds_type, istage, rounds(1:4)


!---------------------------------------------------------------------------------------------
! variables initialization

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

    ! stage 1: only finer blocks send to their coarser neighbors
    ! stage 2: all other relations
    do istage = 1, 2
! write(*,*) "STAGE", istage

    ! the (module-global) communication_counter is the number of neighboring relations
    ! this rank has with all other ranks (it is thus an array of mpisize)
    communication_counter(1:N_friends) = 0_ik
    ! the friends-relation is updated in every call to this routine.
    ! in the beginning all slots are free
    N_friends_used = 0
    mpirank2friend(1:mpisize) = -100
    friend2mpirank(1:N_friends) = -100

!---------------------------------------------------------------------------------------------
! main body

    ! ATTENTION: if you change something here, recall to do the same in reallocate_buffers
    ! new, freshly allocated "friends" slots require consistent initialization
    ! reset integer send buffer position
    int_pos = 2       ! TODO JR why 2? , the first filed contains the size of the XXX
    ! reset first in send buffer position
    int_send_buffer( 1, : ) = 0
    int_send_buffer( 2, : ) = -99

    ! loop over active heavy data. NOTE: hvy_id has a linear correspondance to lgt_id,
    ! i.e.g the ordering in hvy_id and lgt_id is the same. this is very important for the
    ! secondary rule, which is that larger lgt_id wins. this works only if I treat the blocks
    ! in INCREASING lgt_id ordering.
    do k = 1, hvy_n
        ! calculate light id
        sender_hvy_id = hvy_active(k)
        call hvy_id_to_lgt_id( sender_lgt_id, hvy_active(k), myrank, N )

        ! loop over all neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists
            if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                !  ----------------------------  determin the core ids and properties of neighbor  ------------------------------
                ! TODO: check if info available when searching neighbor and store it in hvy_neighbor
                ! neighbor light data id
                neighbor_lgt_id = hvy_neighbor( hvy_active(k), neighborhood )
                ! calculate neighbor rank
                call lgt_id_to_proc_rank( neighbor_rank, neighbor_lgt_id, N )
                ! neighbor heavy id
                call lgt_id_to_hvy_id( hvy_id_receiver, neighbor_lgt_id, neighbor_rank, N )
                ! define level difference: sender - receiver, so +1 means sender on higher level
                level_diff = lgt_block( sender_lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )
                !
                ! if ( istage == 1) Then
                !
                !     if (level_diff /= +1 ) cycle
                ! else
                !     if (level_diff == +1 ) cycle
                ! endif


                !  ----------------------------  here decide which values are taken for redundant nodes --------------------------------

                ! here is the core of the ghost point rules
                ! primary criterion: (very fine/historic fine) wins over (fine) wins over (same) wins over (coarse)
                ! secondary criterion: the higher light id wins NOTE: this is an IMPLICIT rule, enforced by loop ordering ONLY.

                ! comment: the same dominance rules within the ghos nodes are realized by the sequence of filling in the values,
                ! first coarse then same then finer, always in the sequence of the hvy id the redundant nodes within the ghost nodes and maybe in the
                ! redundant nodes are written several time, the one folling the above rules should win
                call set_bounds_according_to_ghost_dominance_rules( params, bounds_type, entrySortInRound, &
                     lgt_block, sender_lgt_id, neighbor_lgt_id )

                call get_friend_id_for_mpirank( params, neighbor_rank, id_Friend )

                !----------------------------  pack describing data and node values to send ---------------------------
                if ( myrank == neighbor_rank ) then
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
                    call AppendLineToBuffer( int_send_buffer, real_send_buffer, 0, id_Friend, line_buffer, &
                    hvy_id_receiver, neighborhood, level_diff_indicator )

                else
                    !-----------------------------------------------------------
                    ! external relation (MPI communication)
                    !-----------------------------------------------------------
                    ! count the number of communications with this fried. from that number, the
                    ! integer buffer length can be computed while MPI exchanging data
                    communication_counter(id_Friend) = communication_counter(id_Friend) + 1

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
                        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, 1)

                        call GhostLayer2Line( params, line_buffer, buffer_size, &
                        hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_active(k)) )

                    else
                        ! up/downsample data first, then flatten to 1D buffer
                        call restrict_predict_data( params, res_pre_data, ijkGhosts(:,:, neighborhood, level_diff, bounds_type, 1), &
                        neighborhood, level_diff, bounds_type, hvy_block, hvy_active(k) )

                        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, 3)

                        call GhostLayer2Line( params, line_buffer, buffer_size, &
                        res_pre_data( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) )
                    end if

                    ! the chunk of data is added to the MPI buffers (preparation for sending)
                    call AppendLineToBuffer( int_send_buffer, real_send_buffer, buffer_size, id_Friend, line_buffer, &
                    hvy_id_receiver, neighborhood, level_diff_indicator )

                end if ! (myrank==neighbor_rank)

            end if ! neighbor exists
        end do ! loop over all possible  neighbors
    end do ! loop over all heavy active


    !***********************************************************************
    ! transfer part (send/recv)
    !***********************************************************************
    call isend_irecv_data_2( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, communication_counter )


    !***********************************************************************
    ! Unpack received data in the ghost node layers
    !***********************************************************************
    ! sort data in, ordering is important to keep dominance rules within ghost nodes.
    ! the redundand nodes owend by two blocks only should be taken care by bounds_type (include_redundant. exclude_redundant )
    ! if (istage == 1) Then
    !     entrySortInRound_start = 1
    !     entrySortInRound_end = 1
    !
    !     rounds = (/3, 0, 0, 0/)
    ! else
    !     entrySortInRound_start = 1
    !     entrySortInRound_end = 4
    !
    !     rounds = (/1, 2, 3, 4/)
    ! endif

    if (istage == 1) Then
        entrySortInRound_start = 1
        entrySortInRound_end = 1

        rounds = (/3, 0, 0, 0/)
    else
        entrySortInRound_start = 1
        entrySortInRound_end = 4

        rounds = (/1, 2, 4, 0/)
    endif

    ! do currentSortInRound = entrySortInRound_start,  entrySortInRound_end ! coarse, same, fine, historic fine
    do iround = entrySortInRound_start,  entrySortInRound_end ! coarse, same, fine, historic fine
        currentSortInRound = rounds(iround)
        ! why now looping over mpiranks and not friends? The reason is the secondary rule,
        ! according to which the larger lgt_id wins. this works only if I treat the blocks
        ! in INCREASING lgt_id ordering. The lgt_id ordering is the same as MPIRANK ordering.
        ! Thus be careful to treat ranks in increasing order, not randomized
        do k = 1, mpisize
            id_Friend = mpirank2friend(k)

            ! skip procs which are not my friends
            if (id_Friend<0) cycle

            irank = friend2mpirank(id_Friend)
            if (irank == myrank+1) then

                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ! process-internal ghost points (direct copy)
                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                l = 2  ! first field is size of data
                do while ( int_send_buffer(l, id_Friend) /= -99 )
                    ! unpack the description of the next data chunk
                    ! required info:  sender_hvy_id, hvy_id_receiver, neighborhood, level_diff, bounds_type, entrySortInRound
                    hvy_id_receiver = int_send_buffer(l, id_Friend)
                    neighborhood = int_send_buffer(l+1, id_Friend)

                    ! unpack & evaluate level_diff_indicator (contains multiple information, unpack it)
                    level_diff_indicator = int_send_buffer(l+2, id_Friend)
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

                        if (bounds_type == exclude_redundant) then
                            ! step (a) into a temporary block, extract the ONLY_REDUNDANT part
                            ! step (b) patch the entire INCLUDE_REDUNDANT into the block
                            ! step (c) put the data form step (a) back into the block.

                            ! ------- step (a) -------
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, 2)

                            tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) = &
                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2),ijk1(1,3):ijk1(2,3), :, hvy_id_receiver)

                            ! ------- step (b) -------
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, 2)
                            ijk2 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, 1)

                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                            hvy_block( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :, sender_hvy_id)

                            ! ------- step (c) -------
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, 2)

                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                            tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :)

                        else
                            ! for INCLUDE_REDUNDANT, just copy the patch and be happy
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, 2)
                            ijk2 = ijkGhosts(:,:, neighborhood, level_diff, bounds_type, 1)

                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                            hvy_block( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :, sender_hvy_id)
                        endif

                    else  ! interpolation or restriction before inserting
                        call restrict_predict_data( params, res_pre_data, ijkGhosts(1:2,1:3, neighborhood, level_diff, INCLUDE_REDUNDANT, 1), neighborhood, level_diff, &
                        bounds_type, hvy_block, sender_hvy_id )

                        ! copy interpolated / restricted data to ghost nodes layer
                        ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
                        ! They depend on the neighbor-relation, level difference and the bounds type.
                        ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                        if (bounds_type == EXCLUDE_REDUNDANT) then
                            ! step (a) into a temporary block, extract the ONLY_REDUNDANT part
                            ! step (b) patch the entire INCLUDE_REDUNDANT into the block
                            ! step (c) put the data from step (a) back into the block.
                            ! ------- step (a) -------
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, 2)

                            tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) = &
                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver)

                            ! ------- step (b) -------
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, 2)
                            ijk2 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, 3)

                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                            res_pre_data( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :)

                            ! ------- step (c) -------
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, ONLY_REDUNDANT, 2)

                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                            tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :)

                        else
                            ijk1 = ijkGhosts(:, :, neighborhood, level_diff, INCLUDE_REDUNDANT, 2)
                            ijk2 = ijkGhosts(:, :, neighborhood, level_diff, INCLUDE_REDUNDANT, 3)

                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                            res_pre_data( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), :)

                        endif
                    end if

                    ! increase buffer postion marker
                    l = l + 5
                end do

            else

                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ! process-external ghost points (copy from buffer)
                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ! did I recv something from this rank?
                if ( (communication_counter(id_Friend) /= 0) ) then
                    l = 2  ! first field is size of data
                    do while ( int_receive_buffer(l, id_Friend) /= -99 )
                        ! unpack the description of the next data chunk
                        hvy_id_receiver = int_receive_buffer(l, id_Friend)
                        neighborhood = int_receive_buffer(l+1, id_Friend)

                        ! unpack & evaluate level_diff_indicator (contains multiple information, unpack it)
                        level_diff_indicator = int_receive_buffer(l+2, id_Friend)
                        entrySortInRound = modulo( level_diff_indicator, 16 )

                        ! check if this entry is processed in this round, otherwise cycle to next
                        if (entrySortInRound /= currentSortInRound ) then
                            l = l + 5  ! to read the next entry
                            cycle      ! go on to next entry
                        end if

                        level_diff  = modulo( level_diff_indicator/16 , 16 ) - 1_ik
                        bounds_type = modulo( level_diff_indicator/256, 16 )
                        buffer_position = int_receive_buffer(l+3, id_Friend)
                        buffer_size     = int_receive_buffer(l+4, id_Friend)

                        ! copy data to line buffer. we now need to extract this to the ghost nodes layer (2D/3D)
                        line_buffer(1:buffer_size) = real_receive_buffer( buffer_position : buffer_position-1 + buffer_size, id_Friend )

                        ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
                        ! They depend on the neighbor-relation, level difference and the bounds type.
                        ! The last index is 1-sender 2-receiver 3-restricted/predicted.

                        if ( bounds_type == EXCLUDE_REDUNDANT ) then

                            ! extract INCLUDE_REDUNDANT in tmp block
                            call Line2GhostLayer2( params, line_buffer, ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, 2), tmp_block )
                            ! COPY ONLY_REDUNDANT from block
                            ijk1 = ijkGhosts( :, :, neighborhood, level_diff, ONLY_REDUNDANT, 2)

                            tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :) = &
                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver)

                            ! copy everything to the block, INCLUDE_REDUNDANT
                            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, 2)

                            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, hvy_id_receiver ) = &
                            tmp_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :)

                        else
                            ! for INCLUDE_REDUNDANT, just copy
                             call Line2GhostLayer( params, line_buffer, ijkGhosts(:,:, neighborhood, level_diff, bounds_type, 2), hvy_block, hvy_id_receiver )
                        endif


                        ! increase buffer postion marker
                        l = l + 5
                    end do
                end if

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

! if (entrySortInRound == 2) entrySortInRound=1

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


!############################################################################################################

subroutine check_unique_origin(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)

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
    logical                             :: testOriginFlag
    integer(kind=ik)                    :: hvy_id_k, iteration , lgt_id
    real(kind=rk) :: tmp

    integer(kind=ik)                    :: i1, i2, iStep, j1, j2, jStep, k1, k2, kStep  , i,j,k, boundaryIndex
    integer(kind=ik)                    :: Bs, g    , levelLocal , levelOrigin , lastRedundantOrigin

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
        call abort(111111,"Same origin of ghost nodes check failed - stopping.")
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

subroutine add_hvy_data( params, line_buffer, data_bounds, hvy_block, hvy_synch, hvy_id )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(inout)                    :: line_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> heavy synch array
    integer(kind=1), intent(inout)                  :: hvy_synch(:, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    buffer_i = 1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all data fields
    do dF = 1, params%number_data_fields
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases kN is always 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    ! write data buffer
                    hvy_block( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id ) + line_buffer( buffer_i )

                    ! count synchronized data
                    ! note: only for first datafield
                    if (dF==1) hvy_synch( i, j, k, hvy_id ) = hvy_synch( i, j, k, hvy_id ) + 1_1

                    ! increase buffer counter
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    end do

end subroutine add_hvy_data

!############################################################################################################

subroutine compare_hvy_data( params, line_buffer, data_bounds, hvy_block, hvy_id, stop_status, level_diff, my_ref, tc )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(inout)                    :: line_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id, level_diff, my_ref
    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status
    integer(kind=tsize)::tc

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i, oddeven, bs, g

    ! error threshold
    real(kind=rk)                                   :: eps


    ! error norm
    real(kind=rk)       :: error_norm

    Bs = params%number_block_nodes
    g = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! variables initialization
    buffer_i = 1

    ! NOTE: newer versions do not compare actual data, but fill the hvy_blocks with
    ! their lgt_ids. This makes the comparison much easier, as those values differ by at least 1.0_rk
    eps = 1e-6_rk

    ! reset error norm
    error_norm = 0.0_rk

    ! the first index of the redundant points is (g+1, g+1, g+1)
    ! so if g is even, then we must compare the odd indices i,j,k on the lines
    ! of the redundant points.
    ! if g is odd, then we must compare the even ones
    ! Further note that BS is odd (always), so as odd+even=odd and odd+odd=even
    ! we can simply study the parity of g
    oddeven = mod(params%number_ghost_nodes,2)

!---------------------------------------------------------------------------------------------
! main body
    ! loop over all data fields
    do dF = 1, params%number_data_fields
        ! third dimension, note: for 2D cases k is always 1
        do k = data_bounds(1,3), data_bounds(2,3)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! first dimension
                do i = data_bounds(1,1), data_bounds(2,1)

                    if (level_diff /= -1) then
                        ! on the same or coarser level, the comparison just takes all points, no odd/even downsampling required.
                        error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i )))

                        hvy_block_test_err( i, j, k, dF, hvy_id ) = abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i ))
                        hvy_block_test_val( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id )
                        hvy_block_test_interpref( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )

                    else
                        ! if the level diff is -1, I compare with interpolated (upsampled) data. that means every EVEN
                        ! point is the result of interpolation, and not truely redundant.
                        ! Note this routine ALWAYS just compares the redundant nodes, so it will mostly be called
                        ! with a line of points (i.e. one dimension is length one)
                        !
                        ! This routine has been tested:
                        !   - old method (working version): no error found (okay)
                        !   - old method, non_uniform_mesh_correction=0; in params file -> plenty of errors (okay)
                        !   - old method, sync stage 4 deactivated: finds all occurances of "3finer blocks on corner problem" (okay)
                        !   - new method, averaging, no error found (makes sense: okay)
                        if (oddeven==0) then
                            ! even number of ghost nodes -> comparison on ODD points
                             if ( (mod(i,2)/=0) .and. (mod(j,2)/=0) .and. (mod(k,2)/=0) ) then
                                error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i )))

                                hvy_block_test_err( i, j, k, dF, hvy_id ) = abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i ))
                                hvy_block_test_val( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id )
                                hvy_block_test_interpref( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )
                            endif
                        else
                            ! odd number of ghost nodes -> comparison on EVEN points
                            if ( (mod(i,2)==0) .and. (mod(j,2)==0) .and. (mod(k,2)==0) ) then
                                error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i )))

                                hvy_block_test_err( i, j, k, dF, hvy_id ) = abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i ))
                                hvy_block_test_val( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id )
                                hvy_block_test_interpref( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )
                            endif
                        endif
                    endif
                    buffer_i = buffer_i + 1
                end do
            end do
        end do
    end do

    if (error_norm > eps)  then
        write(*,'("ERROR: difference in redundant nodes ",es12.4," level_diff=",i2, " hvy_id=",i6,1x,i6," rank=",i5 )') &
        error_norm, level_diff, nint(hvy_block( size(hvy_block,1)/2, size(hvy_block,2)/2, 1, 1, hvy_id )), hvy_id, params%rank
        write(*,*) "refinement status", my_ref, "tc=", tc
        ! stop program
        stop_status = .true.
    end if

end subroutine compare_hvy_data

!############################################################################################################

subroutine isend_irecv_data_2( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, communication_counter )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> send/receive buffer, integer and real
    integer(kind=ik), intent(inout)       :: int_send_buffer(:,:)
    integer(kind=ik), intent(inout)       :: int_receive_buffer(:,:)
    real(kind=rk), intent(inout)          :: real_send_buffer(:,:)
    real(kind=rk), intent(inout)          :: real_receive_buffer(:,:)

    integer(kind=ik), intent(inout)       :: communication_counter(:)

    ! process rank
    integer(kind=ik)                    :: rank
    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! MPI status
    !integer                             :: status(MPI_status_size)

    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI request
    integer(kind=ik)                    :: send_request(N_friends_used), recv_request(N_friends_used)

    ! column number of send buffer, column number of receive buffer, real data buffer length
    integer(kind=ik)                    :: length_realBuffer, int_length, mpirank_partner

    ! loop variable
    integer(kind=ik)                    :: k, i, ifriend


!---------------------------------------------------------------------------------------------
! variables initialization

    rank = params%rank

    ! do ifriend = 1,4
    !     write(*,*) rank, " freund:",ifriend, "is rank", &
    !     friend2mpirank(ifriend)-1, "messages:", communication_counter(ifriend)
    ! enddo
    ! call MPI_barrier(WABBIT_COMM, k)


!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! first: integer data


    ! reset request arrays
    i = 0
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! loop over com matrix
    do k = 1, N_friends_used
        ! communication between proc rank and proc k-1
        if ( communication_counter(k) > 0 ) then
            mpirank_partner = friend2mpirank(k)-1 ! zero based

            ! length of integer buffer
            int_length = 5*communication_counter(k) + 3

            ! increase communication counter
            i = i + 1

            ! send data
            tag = rank
            call MPI_Isend( int_send_buffer(1, k), int_length, MPI_INTEGER4, &
                 mpirank_partner, tag, WABBIT_COMM, send_request(i), ierr)

            ! receive data
            tag = mpirank_partner
            call MPI_Irecv( int_receive_buffer(1, k), int_length, MPI_INTEGER4, &
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

    ! loop over corresponding com matrix line
    do k = 1, N_friends_used
        ! communication between proc rank and proc k-1
        if ( communication_counter(k) > 0 ) then
            mpirank_partner = friend2mpirank(k)-1 ! zero based

            ! increase communication counter
            i = i + 1

            ! real buffer length is stored as the first entry in the integer buffer,
            ! hence we know how much data we'll receive
            length_realBuffer = int_receive_buffer(1, k)

            ! receive data
            tag = 1000*(friend2mpirank(k)-1)
            call MPI_Irecv( real_receive_buffer(1:length_realBuffer, k), length_realBuffer, MPI_REAL8, &
            mpirank_partner, MPI_ANY_TAG, WABBIT_COMM, recv_request(i), ierr)

            ! real buffer length is stored as the first entry in the integer buffer,
            ! hence we know how much data we'll receive
            length_realBuffer = int_send_buffer(1, k)

            ! send data
            tag = 1000*rank
            call MPI_Isend( real_send_buffer(1:length_realBuffer, k), length_realBuffer, MPI_REAL8, &
            mpirank_partner, tag, WABBIT_COMM, send_request(i), ierr)

        end if
    end do

    ! synchronize non-blocking communications
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
    end if

end subroutine isend_irecv_data_2

!############################################################################################################

subroutine set_synch_status( synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, &
    hvy_id, neighborhood, my_ref_status, neighbor_ref_status )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! synch stage
    integer(kind=ik), intent(in)        :: synch_stage

    ! synch status
    logical, intent(inout)    :: synch, neighbor_synch

    ! level difference
    integer(kind=ik), intent(in)        :: level_diff, my_ref_status, neighbor_ref_status

    ! heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    ! list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_id

    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! set synch stage
    ! stage 1: level +1
    ! stage 2: level 0
    ! stage 3: level -1
    ! stage 4: special
    synch = .false.
    neighbor_synch = .false.

    ! this is the zeroth stage. it corrects blocks that are on the same level, but have a different history. one is on Jmax from
    ! before, one has just gotten to Jmax via interpolation. In those cases, the former block has the status +11
    ! which indicates that its redundant nodes must overwrite the ones on the other block (which has been interpolated)
    if ((synch_stage==0) .and. (level_diff==0)) then
        if ((my_ref_status==11) .and. (neighbor_ref_status/=11)) then
            ! if a block has the +11 status, it must send data to the neighbor, if that is not +11
            synch = .true.
        elseif ((my_ref_status/=11) .and. (neighbor_ref_status==11)) then
            ! if a block is not +11 and its neighbor is, then unpack data
            neighbor_synch = .true.
        endif
    endif


    ! stage 1
    if ( (synch_stage == 1) .and. (level_diff == 1) ) then
        ! block send data
        synch = .true.
    elseif ( (synch_stage == 1) .and. (level_diff == -1) ) then
        ! neighbor send data
        neighbor_synch = .true.
    end if

    ! stage 2
    if ( (synch_stage == 2) .and. (level_diff == 0) ) then
        ! block send data
        synch = .true.
        ! neighbor send data
        neighbor_synch = .true.
    end if

    ! stage 3
    if ( (synch_stage == 3) .and. (level_diff == -1) ) then
        ! block send data
        synch = .true.
    elseif ( (synch_stage == 3) .and. (level_diff == 1) ) then
        ! neighbor send data
        neighbor_synch = .true.
    end if

    ! stage 4
    if ( (synch_stage == 4) .and. (level_diff == 0) ) then
        ! neighborhood NE
        if ( neighborhood == 5 ) then
            if ( (hvy_neighbor( hvy_id, 9) /= -1) .or. (hvy_neighbor( hvy_id, 13) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
        ! neighborhood NW
        if ( neighborhood == 6 ) then
            if ( (hvy_neighbor( hvy_id, 10) /= -1) .or. (hvy_neighbor( hvy_id, 15) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
        ! neighborhood SE
        if ( neighborhood == 7 ) then
            if ( (hvy_neighbor( hvy_id, 11) /= -1) .or. (hvy_neighbor( hvy_id, 14) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
        ! neighborhood SW
        if ( neighborhood == 8 ) then
            if ( (hvy_neighbor( hvy_id, 12) /= -1) .or. (hvy_neighbor( hvy_id, 16) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
    end if

end subroutine set_synch_status

!############################################################################################################

subroutine AppendLineToBuffer( int_send_buffer, real_send_buffer, buffer_size, id_Friend, line_buffer, &
    hvy_id, neighborhood, level_diff )

    implicit none

    !> send buffers, integer and real
    integer(kind=ik), intent(inout)        :: int_send_buffer(:,:)
    real(kind=rk), intent(inout)           :: real_send_buffer(:,:)
    ! data buffer size
    integer(kind=ik), intent(in)           :: buffer_size
    ! id integer
    integer(kind=ik), intent(in)           :: id_Friend
    ! restricted/predicted data buffer
    real(kind=rk), intent(inout)           :: line_buffer(:)
    ! data buffer intergers, receiver heavy id, neighborhood id, level difference
    integer(kind=ik), intent(in)           :: hvy_id, neighborhood, level_diff

    ! buffer position
    integer(kind=ik)                       :: buffer_position

    ! fill real buffer
    ! position in real buffer is stored in int buffer
    buffer_position = int_send_buffer( 1, id_Friend ) + 1

    ! real data
    if (buffer_size>0) then
        real_send_buffer( buffer_position : buffer_position-1 + buffer_size, id_Friend ) = line_buffer(1:buffer_size)
    endif

    ! fill int buffer
    ! sum size of single buffers on first element
    int_send_buffer(1  , id_Friend ) = int_send_buffer(1  , id_Friend ) + buffer_size

    ! save: neighbor id, neighborhood, level difference, buffer size
    int_send_buffer( int_pos(id_Friend),   id_Friend ) = hvy_id
    int_send_buffer( int_pos(id_Friend)+1, id_Friend ) = neighborhood
    int_send_buffer( int_pos(id_Friend)+2, id_Friend ) = level_diff
    int_send_buffer( int_pos(id_Friend)+3, id_Friend ) = buffer_position
    int_send_buffer( int_pos(id_Friend)+4, id_Friend ) = buffer_size
    ! mark end of buffer with -99, will be overwritten by next element if it is nt the last one
    int_send_buffer( int_pos(id_Friend)+5, id_Friend ) = -99

    int_pos(id_Friend) = int_pos(id_Friend) +5
end subroutine AppendLineToBuffer

!############################################################################################################
subroutine write_real5(data_block,hvy_active, hvy_n, fileName, rank )
    ! dump all data including ghost nodes for debugging, eg with matlab:

!    function   [data ] =  read(fileName)

!    fid=fopen(fileName, 'rb');        % Open the file.
!    [dataSize, count ] =fread(fid, 5, 'int32') ;
!    data = zeros(dataSize') ;

!    allVals = zeros(prod(dataSize),1) ;
!    allInd  = zeros(prod(dataSize),5)  ;

!    lineNum = 0 ;
!    while(1)
!     [coord, countC ] =fread(fid, 5, 'int32');
!     [val, countV ] =fread(fid, 1, 'float64') ;

!     if (countC*countV == 0 )
!         disp('all')
!         disp ( countC)
!         disp ( countV)
!         if (prod(dataSize) ~= lineNum)
!             fclose(fid) ;
!             error('wrong linenumberm better check it')
!         end
!         %!disp(lineNum)
!         break
!     end
!     lineNum = lineNum + 1 ;
!     allVals(lineNum)  = val ;
!     allInd(lineNum,:)   = coord';

!      data(coord(1), coord(2) , coord(3), coord(4),coord(5))   = val ;
!    end
!    fclose(fid) ;
!    end

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
   !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n, rank

    real(kind=rk), intent(in)       :: data_block(:, :, :, :, :)
    character(len=*), intent(in)  :: fileName
    character(len=3) :: rankname
    integer                         :: i,j,k,l,m

write(rankname,'(i3.3)') rank
    open(unit=11, file= fileName//'.'//rankname//'.dat', form='unformatted', status='replace',access='stream')

    write(*,*) "dumping hvy_n=", hvy_n, "rank=", rank
    !write(*,*) size(data_block,1), size(data_block,2), size(data_block,3), size(data_block,4), hvy_n

    write(11) size(data_block,1), size(data_block,2), size(data_block,3), size(data_block,4), hvy_n

    ! loop sequence not very quick, but i prefere this sequence
    do m = 1, hvy_n
        do i = 1, size(data_block,1)
            do j = 1, size(data_block,2)
                do k = 1, size(data_block,3)
                    do l = 1, size(data_block,4)
                        write(11) i, j, k, l, m, data_block(i, j, k, l, hvy_active(m) )
                    end do
                end do
            end do
        end do
    end do

    close(11)
end subroutine
!############################################################################################################








subroutine check_redundant_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, &
    hvy_n, stop_status, stage0, force_averaging )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

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

    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status
    ! stage0: correct blocks that are on the same level, but have a different history. one is on Jmax from
    ! before, one has just gotten to Jmax via interpolation. In those cases, the former block has the status +11
    ! which indicates that its redundant nodes must overwrite the ones on the other block (which has been interpolated)
    logical, intent(in):: stage0, force_averaging

    ! MPI parameter
    integer(kind=ik)                    :: myrank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! loop variables
    integer(kind=ik)                    :: N, k, dF, neighborhood, invert_neighborhood, neighbor_num, level_diff, l

    ! id integers
    integer(kind=ik)                    :: lgt_id, neighbor_lgt_id, neighbor_rank, hvy_id
    ! type of data bounds
    ! exclude_redundant, include_redundant, only_redundant
    integer(kind=ik)                    :: data_bounds_type
    integer(kind=ik), dimension(2,3)    :: data_bounds, data_bounds2
    ! data buffer size
    integer(kind=ik)                    :: buffer_size, buffer_position
    ! grid parameter
    integer(kind=ik)                    :: Bs, g, stage_start
    ! number of datafields
    integer(kind=ik)                    :: NdF, id_Friend
    ! type of data writing
    character(len=25)                   :: data_writing_type
    ! synch stage loop variables
    integer(kind=ik) :: synch_stage, stages
    ! synch status
    ! synch == .true. : active block sends data to neighboring block
    ! neighbor_synch == .true. : neighbor block send data to active block
    logical    :: synch, neighbor_synch, test2
write(*,*) "warning you re calling the old routine. captain."
 !---------------------------------------------------------------------------------------------
! variables initialization

    if (.not. ghost_nodes_module_ready) then
        call init_ghost_nodes( params )
    endif

    ! if this mpirank has no active blocks, it has nothing to do here.
    if (hvy_n == 0) return

    ! hack to use subroutine as redundant nodes test and for ghost nodes synchronization
    if (stop_status) then
        ! synchronization
        ! exclude_redundant, include_redundant, only_redundant
        data_bounds_type = include_redundant
        ! 'average', 'simple', 'staging', 'compare'
        data_writing_type = 'staging'

        if ( force_averaging ) then
          data_writing_type='average'
        endif

    else
        ! nodes test
        ! exclude_redundant, include_redundant, only_redundant
        data_bounds_type = include_redundant
        ! 'average', 'simple', 'staging', 'compare'
        data_writing_type = 'compare'
        ! reset status
        stop_status = .false.

    end if

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes
    NdF   = params%number_data_fields
    N     = params%number_blocks
    myrank  = params%rank
    number_procs = params%number_procs

    ! set loop number for 2D/3D case
    neighbor_num = size(hvy_neighbor, 2)


    ! the (module-global) communication_counter is the number of neighboring relations
    ! this rank has with all other ranks (it is thus an array of number_procs)
    communication_counter(1:N_friends) = 0_ik
    ! the friends-relation is updated in every call to this routine.
    ! in the beginning all slots are free
    N_friends_used = 0
    mpirank2friend(1:number_procs) = -100
    friend2mpirank(1:N_friends) = -100
    ! for technical reasons, I set me as my first friend here. this solves problems
    ! if I have no blocks (and thus do not appear in the friends list)
    N_friends_used = 1
    mpirank2friend(myrank+1) = 1
    friend2mpirank(1) = myrank+1


    ! reseting all ghost nodes to zero
    if ( (data_writing_type == 'average') .and. (data_bounds_type /= only_redundant) ) then
        do k = 1, hvy_n
            !-- x-direction
            hvy_block(1:g, :, :, :, hvy_active(k) )               = 0.0_rk
            hvy_block(Bs+g+1:Bs+2*g, :, :, :, hvy_active(k) )     = 0.0_rk
            !-- y-direction
            hvy_block(:, 1:g, :, :, hvy_active(k) )               = 0.0_rk
            hvy_block(:, Bs+g+1:Bs+2*g, :, :, hvy_active(k) )     = 0.0_rk
            !-- z-direction
            if ( params%threeD_case ) then
                hvy_block(:, :, 1:g, :, hvy_active(k) )           = 0.0_rk
                hvy_block(:, :, Bs+g+1:Bs+2*g, :, hvy_active(k) ) = 0.0_rk
            end if
        end do
    end if

    stage_start = 1
    stages = 1

    ! set number of synch stages
    if ( data_writing_type == 'staging' ) then
        ! all four stages
        stages = 4
        if (stage0) stage_start=0
    end if

!---------------------------------------------------------------------------------------------
! main body


    ! loop over all synch stages
    do synch_stage = stage_start, stages

        ! in the staging type the ghost nodes bounds depend on the stage as well
        if (data_writing_type=="staging") then
            if (synch_stage==3)  then
                data_bounds_type = exclude_redundant

            elseif (synch_stage == 0) then
                ! stage0: correct blocks that are on the same level, but have a different history. one is on Jmax from
                ! before, one has just gotten to Jmax via interpolation. In those cases, the former block has the status +11
                ! which indicates that its redundant nodes must overwrite the ones on the other block (which has been interpolated)
                data_bounds_type = only_redundant

            else
                data_bounds_type = include_redundant
            endif
        endif

        ! ATTENTION: if you change something here, recall to do the same in reallocate_buffers
        ! new, freshly allocated "friends" slots require consistent initialization
        ! reset integer send buffer position
        int_pos = 2
        ! reset first in send buffer position
        int_send_buffer( 1, : ) = 0
        int_send_buffer( 2, : ) = -99

        ! loop over active heavy data
        if (data_writing_type=="average") then
            do k = 1, hvy_n

                ! reset synch array
                ! alles auf null, knoten im block auf 1
                ! jeder später gespeicherte knoten erhöht wert um 1
                ! am ende der routine wird der wert aus dem synch array ggf. für die durchschnittsberechnung benutzt
                ! synch array hat die maximale anzahl von blöcken pro prozess alloziiert, so dass die heavy id unverändert
                ! benutzt werden kann
                ! ghost nodes layer auf 1 setzen, wenn nur die redundanten Knoten bearbeitet werden
                if (data_bounds_type == only_redundant) then
                    hvy_synch(:, :, :, hvy_active(k)) = 1
                else
                    hvy_synch(:, :, :, hvy_active(k)) = 0
                end if
                ! alles knoten im block werden auf 1 gesetzt

                ! todo: ist erstmal einfacher als nur die redundaten zu setzen, aber unnötig
                ! so gibt es aber nach der synch keine nullen mehr, kann ggf. als synch test verwendet werden?
                if ( params%threeD_case ) then
                    hvy_synch( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, hvy_active(k)) = 1
                else
                    hvy_synch( g+1:Bs+g, g+1:Bs+g, 1, hvy_active(k)) = 1
                end if

            end do
        end if

        do k = 1, hvy_n
            do neighborhood = 1, neighbor_num
                ! neighbor exists
                if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                    ! 0. ids bestimmen
                    neighbor_lgt_id = hvy_neighbor( hvy_active(k), neighborhood )
                    call lgt_id_to_proc_rank( neighbor_rank, neighbor_lgt_id, N )
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), myrank, N )
                    call lgt_id_to_hvy_id( hvy_id, neighbor_lgt_id, neighbor_rank, N )
                    ! calculate the difference between block levels
                    ! define leveldiff: sender - receiver, so +1 means sender on higher level
                    ! sender is active block (me)
                    level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )

                    ! 1. ich (aktiver block) ist der sender für seinen nachbarn
                    ! lese daten und sortiere diese in bufferform
                    ! wird auch für interne nachbarn gemacht, um gleiche routine für intern/extern zu verwenden
                    ! um diue lesbarkeit zu erhöhen werden zunächst die datengrenzen bestimmt
                    ! diese dann benutzt um die daten zu lesen
                    ! 2D/3D wird bei der datengrenzbestimmung unterschieden, so dass die tatsächliche leseroutine stark vereinfacht ist
                    ! da die interpolation bei leveldiff -1 erst bei der leseroutine stattfindet, werden als datengrenzen die für die interpolation noitwendigen bereiche angegeben
                    ! auch für restriction ist der datengrenzenbereich größer, da dann auch hier später erst die restriction stattfindet
                    !!!!!!!!!!! call calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, 'sender' )
                    data_bounds = ijkGhosts(:,:, neighborhood, level_diff, data_bounds_type, 1)

                    ! vor dem schreiben der daten muss ggf interpoliert werden
                    ! hier werden die datengrenzen ebenfalls angepasst
                    ! interpolierte daten stehen in einem extra array
                    ! dessen größe richtet sich nach dem größten möglichen interpolationsgebiet: (Bs+2*g)^3
                    ! auch die vergröberten daten werden in den interpolationbuffer geschrieben und die datengrenzen angepasst
                    if ( level_diff == 0 ) then
                        ! lese nun mit den datengrenzen die daten selbst
                        ! die gelesenen daten werden als buffervektor umsortiert
                        ! so können diese danach entweder in den buffer geschrieben werden oder an die schreiberoutine weitergegeben werden
                        ! in die lese routine werden nur die relevanten Daten (data bounds) übergeben
                        call GhostLayer2Line( params, line_buffer, buffer_size, &
                        hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), data_bounds(1,3):data_bounds(2,3), :, hvy_active(k)) )
                    else
                        ! interpoliere daten
                        call restrict_predict_data( params, res_pre_data, data_bounds, neighborhood, level_diff, data_bounds_type, hvy_block, hvy_active(k))

                        data_bounds2 = ijkGhosts(1:2, 1:3, neighborhood, level_diff, data_bounds_type, 3)
                        ! lese daten, verwende interpolierte daten
                        call GhostLayer2Line( params, line_buffer, buffer_size, res_pre_data( data_bounds2(1,1):data_bounds2(2,1), &
                                              data_bounds2(1,2):data_bounds2(2,2), data_bounds2(1,3):data_bounds2(2,3),:) )
                    end if

                    call get_friend_id_for_mpirank( params, neighbor_rank, id_Friend )

                    ! daten werden jetzt entweder in den speicher geschrieben -> schreiberoutine
                    ! oder in den send buffer geschrieben
                    ! schreiberoutine erhält die date grenzen
                    ! diese werden vorher durch erneuten calc data bounds aufruf berechnet
                    ! achtung: die nachbarschaftsbeziehung wird hier wie eine interner Kopieren ausgewertet
                    ! invertierung der nachbarschaftsbeziehung findet beim füllen des sendbuffer statt
                    if ( (myrank==neighbor_rank).and.(data_writing_type=='simple') ) then
                        ! internal neighbor and direct writing method: copy the ghost nodes as soon as possible, without passing
                        ! via the buffers first.
                        ! data bounds
                        !!!!!!!!!!!!!!call calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, 'receiver' )
                        data_bounds = ijkGhosts(:,:, neighborhood, level_diff, data_bounds_type, 2)
                        ! simply write data. No care
                        call Line2GhostLayer( params, line_buffer, data_bounds, hvy_block, hvy_id )

                    else
                        ! synch status for staging method
                        synch = .true.
                        if (data_writing_type == 'staging') then
                            call set_synch_status( synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_active(k), &
                            neighborhood, lgt_block(lgt_id,params%max_treelevel+2), lgt_block(neighbor_lgt_id,params%max_treelevel+2)  )
                        end if
                        ! first: fill com matrix, count number of communication to neighboring process, needed for int buffer length
                        communication_counter(id_Friend) = communication_counter(id_Friend) + 1

                        if (synch) then
                            ! active block send data to his neighbor block
                            ! fill int/real buffer
                            call AppendLineToBuffer( int_send_buffer, real_send_buffer, buffer_size, id_Friend, line_buffer, &
                            hvy_id, neighborhood, level_diff )
                        else
                            ! neighbor block send data to active block
                            ! write -1 to int_send buffer, placeholder
                            int_send_buffer( int_pos(id_Friend) : int_pos(id_Friend)+4  , id_Friend ) = -1
                            ! increase int buffer position
                            int_pos(id_Friend) = int_pos(id_Friend) + 5
                        end if

                    end if

                end if
            end do
        end do

        ! pretend that no communication with myself takes place, in order to skip the
        ! MPI transfer in the following routine. NOTE: you can also skip this step and just have isend_irecv_data_2
        ! transfer the data, in which case you should skip the copy part directly after isend_irecv_data_2
        communication_counter( mpirank2friend(myrank+1) ) = 0

        !***********************************************************************
        ! transfer part (send/recv)
        !***********************************************************************
        ! send/receive data
        ! note: todo, remove dummy subroutine
        ! note: new dummy subroutine sets receive buffer position accordingly to process rank
        ! note: todo: use more than non-blocking send/receive
        call isend_irecv_data_2( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, communication_counter  )

        ! fill receive buffer for internal neighbors for averaging writing type
        if ( (data_writing_type == 'average') .or. (data_writing_type == 'compare') .or. (data_writing_type == 'staging') ) then
            ! fill receive buffer
            int_receive_buffer( 1:int_pos(mpirank2friend(myrank+1))  , mpirank2friend(myrank+1) ) = &
                int_send_buffer( 1:int_pos(mpirank2friend(myrank+1))  , mpirank2friend(myrank+1) )
            real_receive_buffer( 1:int_receive_buffer(1,mpirank2friend(myrank+1)), mpirank2friend(myrank+1) ) = &
                real_send_buffer( 1:int_receive_buffer(1,mpirank2friend(myrank+1)), mpirank2friend(myrank+1) )
            ! change communication_counter, equired to trigger buffer unpacking in last step
            communication_counter(mpirank2friend(myrank+1)) = 1
        end if

        !***********************************************************************
        ! Unpack received data in the ghost node layers
        !***********************************************************************
        ! Daten einsortieren
        ! für simple, average, compare: einfach die buffer einsortieren, Reihenfolge ist egal
        ! staging: erneuter loop über alle blöcke und nachbarschaften, wenn daten notwendig, werden diese in den buffern gesucht
        if ( data_writing_type /= 'staging' ) then
            ! sortiere den real buffer ein
            ! loop over all procs
            do k = 1, N_friends_used
                if ( communication_counter(k) /= 0 ) then
                    ! neighboring proc
                    ! first element in int buffer is real buffer size
                    l = 2
                    ! -99 marks end of data
                    do while ( int_receive_buffer(l, k) /= -99 )

                        hvy_id          = int_receive_buffer(l, k)
                        neighborhood    = int_receive_buffer(l+1, k)
                        level_diff      = int_receive_buffer(l+2, k)
                        buffer_position = int_receive_buffer(l+3, k)
                        buffer_size     = int_receive_buffer(l+4, k)
                        line_buffer(1:buffer_size) = real_receive_buffer( buffer_position : buffer_position-1 + buffer_size, k )

                        ! data bounds
                        !!!!!call calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, 'receiver' )
                        data_bounds = ijkGhosts(:,:, neighborhood, level_diff, data_bounds_type, 2)
                        ! write data, hängt vom jeweiligen Fall ab
                        ! average: schreibe daten, merke Anzahl der geschriebenen Daten, Durchschnitt nach dem Einsortieren des receive buffers berechnet
                        ! simple: schreibe ghost nodes einfach in den speicher (zum Testen?!)
                        ! staging: wende staging konzept an
                        ! compare: vergleiche werte mit vorhandenen werten (nur für redundante knoten sinnvoll, als check routine)
                        select case(data_writing_type)
                            case('simple')
                                ! simply write data
                                call Line2GhostLayer( params, line_buffer, data_bounds, hvy_block, hvy_id )

                            case('average')
                                ! add data
                                call add_hvy_data( params, line_buffer, data_bounds, hvy_block, hvy_synch, hvy_id )

                            case('compare')
                                ! compare data
                                call hvy_id_to_lgt_id( lgt_id, hvy_id, myrank, N )
                                call compare_hvy_data( params, line_buffer, data_bounds, hvy_block, hvy_id, stop_status, level_diff, &
                                 lgt_block(lgt_id, params%max_treelevel+2), treecode2int( lgt_block(lgt_id, 1:params%max_treelevel) ) )

                        end select

                        ! increase buffer postion marker
                        l = l + 5

                    end do
                end if
            end do

            ! last averaging step
            if ( data_writing_type == 'average' ) then
                ! loop over active heavy data
                do k = 1, hvy_n
                    do dF = 1, NdF

                        ! calculate average for all nodes, todo: proof performance?
                        hvy_block(:, :, :, dF, hvy_active(k)) = hvy_block(:, :, :, dF, hvy_active(k)) / real( hvy_synch(:, :, :, hvy_active(k)) , kind=rk)

                    end do
                end do
            end if

        else
            ! staging type
            ! loop over active heavy data
            do k = 1, hvy_n
                ! loop over all neighbors
                do neighborhood = 1, neighbor_num
                    ! neighbor exists
                    if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                        ! invert neighborhood, needed for in buffer searching, because sender proc has invert neighborhood relation
                        invert_neighborhood = inverse_neighbor(neighborhood, dim)

                        ! 0. ids bestimmen
                        ! neighbor light data id
                        neighbor_lgt_id = hvy_neighbor( hvy_active(k), neighborhood )
                        ! calculate neighbor rank
                        call lgt_id_to_proc_rank( neighbor_rank, neighbor_lgt_id, N )
                        ! calculate light id
                        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), myrank, N )
                        ! calculate the difference between block levels
                        ! define leveldiff: sender - receiver, so +1 means sender on higher level
                        ! sender is active block (me)
                        level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_lgt_id, params%max_treelevel+1 )

                        ! set synch status
                        call set_synch_status( synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, &
                        hvy_active(k), neighborhood, lgt_block(lgt_id, params%max_treelevel+2), lgt_block(neighbor_lgt_id,params%max_treelevel+2) )
                        ! synch == .true. bedeutet, dass der aktive block seinem nachbarn daten gibt
                        ! hier sind wir aber auf der seite des empfängers, das bedeutet, neighbor_synch muss ausgewertet werden

                        if (neighbor_synch) then

                            ! search buffers for synchronized data
                            ! first element in int buffer is real buffer size
                            l = 2

                            ! -99 marks end of data
                            test2 = .false.
                            do while ( int_receive_buffer(l, mpirank2friend(neighbor_rank+1)) /= -99 )

                                ! proof heavy id and neighborhood id
                                if (  (int_receive_buffer( l,   mpirank2friend(neighbor_rank+1) ) == hvy_active(k) ) &
                                .and. (int_receive_buffer( l+1, mpirank2friend(neighbor_rank+1) ) == invert_neighborhood) ) then

                                    ! set parameter
                                    ! level diff, read from buffer because calculated level_diff is not sender-receiver
                                    level_diff      = int_receive_buffer(l+2, mpirank2friend(neighbor_rank+1))
                                    buffer_position = int_receive_buffer(l+3, mpirank2friend(neighbor_rank+1))
                                    buffer_size     = int_receive_buffer(l+4, mpirank2friend(neighbor_rank+1))
                                    line_buffer(1:buffer_size) = real_receive_buffer( buffer_position : buffer_position-1 + buffer_size, mpirank2friend(neighbor_rank+1) )

                                    ! data bounds
                                    !!!!!!!!!!!call calc_data_bounds( params, data_bounds, invert_neighborhood, level_diff, data_bounds_type, 'receiver' )
                                    data_bounds = ijkGhosts(:,:, invert_neighborhood, level_diff, data_bounds_type, 2)

                                    ! write data
                                    call Line2GhostLayer( params, line_buffer(1:buffer_size), data_bounds, hvy_block, hvy_active(k) )

                                    ! done, exit the while loop?
                                    test2=.true.
                                    exit
                                end if

                                ! increase buffer postion marker
                                l = l + 5

                            end do
                            if (test2 .eqv. .false.) call abort(777771,"not found")

                        end if

                    end if
                end do
            end do

        end if

    end do ! loop over stages

    if ( data_writing_type=='compare' ) then
        test2 = stop_status
        call MPI_Allreduce(test2, stop_status, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, k )
    endif

end subroutine check_redundant_nodes
