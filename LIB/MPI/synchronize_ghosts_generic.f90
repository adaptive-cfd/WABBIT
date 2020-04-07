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
    integer(kind=ik)   :: g, ii0, ii1
    integer(kind=ik), dimension(3) :: Bs
    ! loop variables
    integer(kind=ik)   :: N, k, neighborhood, level_diff
    ! id integers
    integer(kind=ik)   :: recver_lgtID, recver_rank, recver_hvyID
    integer(kind=ik)   :: sender_hvyID, sender_lgtID

    integer(kind=ik) :: ijk(2,3), isend, irecv
    integer(kind=ik) :: bounds_type, istage, inverse

    if (.not. ghost_nodes_module_ready) then
        ! in order to keep the syntax clean, buffers are module-global and need to be
        ! allocated here.
        call init_ghost_nodes( params )
    endif

    ! if this mpirank has no active blocks, it has nothing to do here.
    if (hvy_n == 0) return



    ! DO NOT FORGET TO REMOVE THIS
    ! call reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )



    if (size(hvy_block,4)>N_max_components) then
        call abort(160720191,"You try to ghost-sync a vector with too many components.")
    endif

    Bs    = params%Bs
    g     = params%n_ghosts
    N     = params%number_blocks
    myrank  = params%rank
    mpisize = params%number_procs



    ! We require two stages: first, we fill all ghost nodes which are simple copy (including restriction),
    ! then in the second stage we can use interpolation and fill the remaining ones.
    do istage = 1, 2

        !***************************************************************************
        ! (i) stage initialization
        !***************************************************************************
        int_pos(:) = 1
        real_pos(:) = 0

        ! compute, locally from the grid info, how much data I recv from and send to all
        ! other mpiranks. this gives us the start indices of each rank in the send/recv buffers. Note
        ! we do not count our internal nodes (.false. as last argument), as they are not put in the
        ! buffer at any time.
        call get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, &
             recv_counter, send_counter, int_recv_counter, int_send_counter, istage, &
             count_internal=.false., ncomponents=size(hvy_block,4))


        ! reset int_send_buffer, but only the parts that will actually be treated.
        do k = 1, params%number_procs
            ii0 = sum(int_send_counter(0:(k-1)-1)) + 1
            ii1 = ii0 + int_send_counter(k-1)
            int_send_buffer(ii0:ii1) = -99

            ii0 = sum(int_recv_counter(0:(k-1)-1)) + 1
            ii1 = ii0 + int_recv_counter(k-1)
            int_recv_buffer(ii0:ii1) = -99
        enddo


        !***************************************************************************
        ! (ii) prepare data for sending
        !***************************************************************************
        do k = 1, hvy_n
            ! calculate light id
            sender_hvyID = hvy_active(k)
            call hvy_id_to_lgt_id( sender_lgtID, sender_hvyID, myrank, N )

            ! loop over all neighbors
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists
                if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
                    !  ----------------------------  determin the core ids and properties of neighbor  ------------------------------
                    ! TODO: check if info available  when searching neighbor and store it in hvy_neighbor
                    ! neighbor light data id
                    recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
                    ! calculate neighbor rank
                    call lgt_id_to_proc_rank( recver_rank, recver_lgtID, N )
                    ! neighbor heavy id
                    call lgt_id_to_hvy_id( recver_hvyID, recver_lgtID, recver_rank, N )
                    ! define level difference: sender - receiver, so +1 means sender on higher level
                    level_diff = lgt_block( sender_lgtID, params%max_treelevel + IDX_MESH_LVL ) - lgt_block( recver_lgtID, params%max_treelevel + IDX_MESH_LVL )


                    ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
                    ! leveldiff =  0 : sender is same level as recver
                    ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
                    if ( istage == 1 ) then
                        ! This is preparation for sending: hence in phase 1, level_diff=-1 skips and in phase 2 [0,+1 skip]
                        if ( level_diff == -1 ) cycle
                    else
                        ! in stage two leveldiff +1 and 0 are already done
                        if ( level_diff ==  0 ) cycle
                        if ( level_diff == +1 ) cycle
                    endif

                    !----------------------------  pack describing data and node values to send ---------------------------
                    if ( myrank == recver_rank ) then
                        !-----------------------------------------------------------
                        ! internal relation (no communication)
                        !-----------------------------------------------------------
                        ! NOTE: 06/11/2018, Thomas. I checked that treating all neighbors external (ie passing by the buffers) is 10-15% slower
                        ! than treating internal separately (on irene@TGCC, using Allinea MAP profiler, 3D testcase, 144 CPUs, infiniband)
                        call send_prepare_internal_neighbor( recver_rank+1, sender_hvyID, recver_hvyID, neighborhood, &
                        level_diff )

                    else
                        !-----------------------------------------------------------
                        ! external relation (MPI communication)
                        !-----------------------------------------------------------
                        call send_prepare_external_neighbor( params, recver_rank+1, hvy_block, &
                        sender_hvyID, recver_hvyID, neighborhood, level_diff )

                    end if ! (myrank==recver_rank)
                end if ! neighbor exists
            end do ! loop over all possible  neighbors
        end do ! loop over all heavy active

        !***************************************************************************
        ! (iii) transfer part (send/recv)
        !***************************************************************************
        call isend_irecv_data_2( params, int_send_buffer, new_send_buffer, int_recv_buffer, new_recv_buffer, isend, irecv )

        ! process-internal ghost points (direct copy)
        call unpack_ghostlayers_internal_neighbor( params, myrank+1, hvy_block )

        call wait_for_xfer(params, isend, irecv)


        !***************************************************************************
        ! (iv) Unpack received data in the ghost node layers
        !***************************************************************************
        ! leveldiff = -1 : sender is finer than recver, restriction is applied on sender side
        ! leveldiff =  0 : sender is same level as recver
        ! leveldiff = +1 : sender coarser than recver, interpolation on sender side

        do k = 1, mpisize ! one-based
            if (k /= myrank+1) then
                ! process-external ghost points (copy from buffer)
                call unpack_ghostlayers_external_neighbor( params, k, hvy_block )
            end if
        end do


    end do ! loop over stages 1,2
end subroutine synchronize_ghosts_generic_sequence



subroutine send_prepare_internal_neighbor( recver_rank, sender_hvyID, recver_hvyID, neighborhood, &
    level_diff )
    implicit none

    integer(kind=ik), intent(in) :: recver_rank ! 1-based
    integer(kind=ik), intent(in) :: sender_hvyID, recver_hvyID
    integer(kind=ik), intent(in) :: neighborhood
    integer(kind=ik), intent(in) :: level_diff

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: level_diff_indicator

    !-----------------------------------------------------------
    ! internal relation (no communication)
    !-----------------------------------------------------------
    ! pack multipe information into one number
    level_diff_indicator =  4096*sender_hvyID + 16*(level_diff+1)

    ! we sort of abuse the routine AppendLineToBuffer here. In fact, we only store the integer data
    ! but do not copy the heavy data to te corresponding buffer. In that sense, we only "recall" what
    ! parameters (level_diff etc) the neighboring relation has.
    call AppendLineToBuffer( int_send_buffer, new_send_buffer, 0, recver_rank, line_buffer, &
    recver_hvyID, neighborhood, level_diff_indicator )

end subroutine




subroutine send_prepare_external_neighbor( params, recver_rank, hvy_block, sender_hvyID, &
    recver_hvyID, neighborhood, level_diff )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in)   :: recver_rank
    integer(kind=ik), intent(in)   :: sender_hvyID, recver_hvyID
    integer(kind=ik), intent(in)   :: neighborhood
    integer(kind=ik), intent(in)   :: level_diff
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: level_diff_indicator, buffer_size
    integer(kind=ik)   :: ijk1(2,3), nc

    nc = size(hvy_block,4)

    ! pack multipe information into one number
    level_diff_indicator = 16*(level_diff+1)

    ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
    ! They depend on the neighbor-relation, level difference and the bounds type.
    ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    if ( level_diff == 0 ) then
        ! simply copy the ghost node layer (no interpolation or restriction here) to a line buffer, which
        ! we will send to our neighbor mpirank
        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, sender_hvyID) )

    else
        ! up/downsample data first, then flatten to 1D buffer
        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

        call restrict_predict_data( params, res_pre_data, ijk1, neighborhood, level_diff, hvy_block, sender_hvyID )

        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, RESPRE)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        res_pre_data( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), 1:nc) )
    end if

    ! the chunk of data is added to the MPI buffers (preparation for sending)
    call AppendLineToBuffer( int_send_buffer, new_send_buffer, buffer_size, recver_rank, line_buffer, &
    recver_hvyID, neighborhood, level_diff_indicator )

end subroutine send_prepare_external_neighbor


subroutine unpack_ghostlayers_external_neighbor( params, recver_rank, hvy_block )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: recver_rank ! 1-based
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    integer(kind=ik) :: l, recver_hvyID, neighborhood, level_diff_indicator
    integer(kind=ik) :: level_diff, bounds_type, buffer_position, buffer_size, rank_destination
    integer(kind=ik) :: ijk1(2,3), i0, nc

    nc = size(hvy_block,4)

    ! did I recv something from this rank?
    if ( (recv_counter(recver_rank-1) /= 0) ) then

        ! start index of this mpirank in the int_buffer
        l = sum(int_recv_counter(0:recver_rank-1-1)) + 1

        do while ( int_recv_buffer(l) > -99 )
            ! unpack the description of the next data chunk
            recver_hvyID = int_recv_buffer(l)
            neighborhood = int_recv_buffer(l+1)

            ! unpack & evaluate level_diff_indicator (contains multiple information, unpack it)
            level_diff_indicator = int_recv_buffer(l+2)
            level_diff           = modulo( level_diff_indicator/16 , 16 ) - 1_ik

            buffer_position  = int_recv_buffer(l+3)
            buffer_size      = int_recv_buffer(l+4)
            rank_destination = int_recv_buffer(l+5) ! FIVE

            if (rank_destination /= params%rank) then
                write(*,*) "rank=", params%rank, "dest=", rank_destination, "co", int_recv_counter, &
                "l=", l, "recver_rank=", recver_rank-1
                call abort(7373872, "EXT this data seems to be not mine!")
            endif

            ! copy data to line buffer. we now need to extract this to the ghost nodes layer (2D/3D)
            i0 = sum(recv_counter(0:recver_rank-1-1)) + buffer_position
            line_buffer(1:buffer_size) = new_recv_buffer( i0 : i0+buffer_size-1 )


            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation, level difference and the bounds type.
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.
            call Line2GhostLayer( params, line_buffer, ijkGhosts(:,:, neighborhood, level_diff, RECVER), hvy_block, recver_hvyID )

            ! increase buffer postion marker
            l = l + 6 ! FIVE
        end do
    end if

end subroutine unpack_ghostlayers_external_neighbor



subroutine unpack_ghostlayers_internal_neighbor( params, recver_rank, hvy_block )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: recver_rank ! 1-based
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    integer(kind=ik) :: l, recver_hvyID, neighborhood, level_diff_indicator
    integer(kind=ik) :: sender_hvyID, level_diff, rank_destination
    integer(kind=ik) :: ijk1(2,3), ijk2(2,3), nc


    ! start index of this mpirank in the int_buffer
    l = sum(int_recv_counter(0:recver_rank-1-1)) + 1

    nc = size(hvy_block,4)

    do while ( int_recv_buffer(l) > -99 )
        ! unpack the description of the next data chunk
        ! required info:  sender_hvyID, recver_hvyID, neighborhood, level_diff
        recver_hvyID = int_recv_buffer(l)
        neighborhood = int_recv_buffer(l+1)

        ! unpack & evaluate level_diff_indicator (contains multiple information, unpack it)
        level_diff_indicator = int_recv_buffer(l+2)
        level_diff           = modulo( level_diff_indicator/16  , 16 ) - 1
        sender_hvyID         = ( level_diff_indicator/4096 )
        rank_destination     = int_recv_buffer(l+5) ! FIVE

        if (rank_destination /= params%rank) then
            write(*,*) "rank=", params%rank, "dest=", rank_destination, "co", int_recv_counter, &
            "l=", l, "recver_rank=", recver_rank-1
            call abort(7373872, "INT this data seems to be not mine!")
        endif


        if ( level_diff == 0 ) then
            ! simply copy from sender block to receiver block (NOTE: both are on the same MPIRANK)
            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation and level difference
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.

            ! ALWAYS EXCLUDE_REDUNDANT
            ! just copy the patch and be happy
            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, RECVER)
            ijk2 = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), 1:nc, recver_hvyID ) = &
            hvy_block( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), 1:nc, sender_hvyID)

        else  ! interpolation or restriction before inserting

            call restrict_predict_data( params, res_pre_data, ijkGhosts(1:2,1:3, neighborhood, level_diff, SENDER), &
            neighborhood, level_diff, hvy_block, sender_hvyID )

            ! copy interpolated / restricted data to ghost nodes layer
            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation, level difference and the bounds type.
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.

            ijk1 = ijkGhosts(:, :, neighborhood, level_diff, RECVER)
            ijk2 = ijkGhosts(:, :, neighborhood, level_diff, RESPRE)

            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), 1:nc, recver_hvyID ) = &
            res_pre_data( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), 1:nc)
        end if

        ! increase buffer position marker
        l = l + 6 ! FIVE
    end do

end subroutine unpack_ghostlayers_internal_neighbor


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




subroutine AppendLineToBuffer( int_send_buffer, new_send_buffer, buffer_size, recver_rank, line_buffer, &
    hvy_id, neighborhood, level_diff_indicator )

    implicit none

    !> send buffers, integer and real
    integer(kind=ik), intent(inout) :: int_send_buffer(:)
    real(kind=rk), intent(inout)    :: new_send_buffer(:)
    ! data buffer size
    integer(kind=ik), intent(in)    :: buffer_size
    ! id integer
    integer(kind=ik), intent(in)    :: recver_rank ! 1-based
    ! restricted/predicted data buffer
    real(kind=rk), intent(inout)    :: line_buffer(:)
    ! data buffer intergers, receiver heavy id, neighborhood id, level difference
    integer(kind=ik), intent(in)    :: hvy_id, neighborhood, level_diff_indicator

    integer(kind=ik)                :: buffer_position, i0, l0

    ! fill real buffer
    ! position in real buffer is stored in int buffer
    buffer_position = real_pos(recver_rank) + 1

    i0 = sum(send_counter(0:recver_rank-1-1)) + buffer_position
    l0 = sum(int_send_counter(0:recver_rank-1-1)) + 0 ! note as int_pos is 1-based, we here use 0-based offset

    ! real data
    if (buffer_size>0) then
        new_send_buffer( i0:i0+buffer_size-1 ) = line_buffer(1:buffer_size)
    endif

    ! save position of NEXT patch in real buffer
    real_pos(recver_rank) = real_pos(recver_rank) + buffer_size

    ! save: neighbor id, neighborhood, level difference, buffer size
    int_send_buffer( l0+int_pos(recver_rank)   ) = hvy_id
    int_send_buffer( l0+int_pos(recver_rank)+1 ) = neighborhood
    int_send_buffer( l0+int_pos(recver_rank)+2 ) = level_diff_indicator
    int_send_buffer( l0+int_pos(recver_rank)+3 ) = buffer_position
    int_send_buffer( l0+int_pos(recver_rank)+4 ) = buffer_size
    int_send_buffer( l0+int_pos(recver_rank)+5 ) = recver_rank-1 ! zero-based ! FIVE

    ! mark end of buffer with -99, will be overwritten by next element if it is not the last one
    ! int_send_buffer( l0+int_pos(recver_rank)+5 ) = -99 ! FIVE
    int_send_buffer( l0+int_pos(recver_rank)+6 ) = -99

    int_pos(recver_rank) = int_pos(recver_rank) + 6 ! FIVE
    ! int_pos(recver_rank) = int_pos(recver_rank) + 5

    ! if (l0+int_pos(recver_rank)+6 > l0+int_send_counter(recver_rank-1)-1 ) then
    !     call abort (2828, "well hell ")
    ! endif

end subroutine AppendLineToBuffer


!############################################################################################################

subroutine isend_irecv_data_2( params, int_send_buffer, new_send_buffer, int_recv_buffer, new_recv_buffer, isend, irecv )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> send/receive buffer, integer and real
    integer(kind=ik), intent(inout)       :: int_send_buffer(:)
    integer(kind=ik), intent(inout)       :: int_recv_buffer(:)
    real(kind=rk), intent(inout)          :: new_send_buffer(:)
    real(kind=rk), intent(inout)          :: new_recv_buffer(:)
    integer(kind=ik), intent(out) :: isend, irecv
    ! process rank
    integer(kind=ik) :: rank
    ! MPI error variable
    integer(kind=ik) :: ierr

    ! MPI message tag
    integer(kind=ik) :: tag
    ! column number of send buffer, column number of receive buffer, real data buffer length
    integer(kind=ik) :: length_realBuffer, mpirank_partner
    integer(kind=ik) :: length_intBuffer
    ! loop variable
    integer(kind=ik) :: k, i0, l0

    rank = params%rank
    isend = 0
    irecv = 0

    ! these two arrays are module-global.
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL


    ! for unpacking the internal ghost nodes, we require the int_recv_buffer
    ! just for my own rank to be completed: better to copy that here, now,
    ! so that we are sure to be able to unpack internal ghost nodes while we wait
    ! for the xfer
    length_intBuffer = int_recv_counter(rank)
    i0 = sum(int_recv_counter(0:rank-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner
    l0 = sum(int_send_counter(0:rank-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner
    int_recv_buffer(i0:i0+length_intBuffer-1) = int_send_buffer(l0:l0+length_intBuffer-1)



    do mpirank_partner = 0, params%number_procs-1 ! zero based

        if ( int_recv_counter(mpirank_partner) > 0 .and. mpirank_partner/=rank ) then
            ! length of integer buffer
            length_intBuffer = int_recv_counter(mpirank_partner)
            i0 = sum(int_recv_counter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! increase communication counter
            irecv = irecv + 1

            ! receive data
            tag = mpirank_partner
            call MPI_Irecv( int_recv_buffer(i0:i0+length_intBuffer-1), length_intBuffer, MPI_INTEGER4, &
            mpirank_partner, tag, WABBIT_COMM, recv_request(irecv), ierr)
        endif

        if ( int_send_counter(mpirank_partner) > 0 .and. mpirank_partner/=rank  ) then
            ! length of integer buffer
            length_intBuffer = int_send_counter(mpirank_partner)
            i0 = sum(int_send_counter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! increase communication counter
            isend = isend + 1

            ! send data
            tag = rank
            call MPI_Isend( int_send_buffer(i0:i0+length_intBuffer-1), length_intBuffer, MPI_INTEGER4, &
            mpirank_partner, tag, WABBIT_COMM, send_request(isend), ierr)
        end if

    end do

    ! Note: internal ghost nodes do not pass by the buffer, they only require the integer
    ! buffer above. Therefore, no data (a message of zero length) will be sent to myself
    do mpirank_partner = 0, params%number_procs-1 ! zero based

        if (send_counter(mpirank_partner) > 0) then

            ! increase communication counter
            isend = isend + 1

            ! the amount of data is pre-computed in get_my_sendrecv_amount_with_ranks
            ! hence we do know how much data we will receive
            length_realBuffer = send_counter(mpirank_partner)

            i0 = sum(send_counter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! send data
            tag = rank
            call MPI_Isend( new_send_buffer(i0:i0+length_realBuffer-1), length_realBuffer, MPI_REAL8, &
            mpirank_partner, tag, WABBIT_COMM, send_request(isend), ierr)

        end if

        if (recv_counter(mpirank_partner) > 0) then
            ! increase communication counter
            irecv = irecv + 1

            ! the amount of data is pre-computed in get_my_sendrecv_amount_with_ranks
            ! hence we do know how much data we will receive
            length_realBuffer = recv_counter(mpirank_partner)

            i0 = sum(recv_counter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! receive data
            tag = mpirank_partner
            call MPI_Irecv( new_recv_buffer(i0:i0+length_realBuffer-1), length_realBuffer, MPI_REAL8, &
            mpirank_partner, MPI_ANY_TAG, WABBIT_COMM, recv_request(irecv), ierr)
        end if
    end do

end subroutine isend_irecv_data_2


subroutine wait_for_xfer(params, isend, irecv)
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


! returns two lists with numbers of points I send to all other procs and how much I
! receive from each proc. note: strictly locally computed, NO MPI comm involved here
subroutine get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active,&
     hvy_n, recv_counter, send_counter, int_recv_counter, int_send_counter, istage, &
     count_internal, ncomponents)

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
    integer(kind=ik), intent(in)        :: hvy_n, ncomponents, istage
    integer(kind=ik), intent(inout)     :: recv_counter(0:), send_counter(0:)
    integer(kind=ik), intent(inout)     :: int_recv_counter(0:), int_send_counter(0:)
    logical, intent(in)                 :: count_internal

    integer(kind=ik) :: k, sender_hvyID, sender_lgtID, myrank, N, neighborhood, recver_rank
    integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID,level_diff

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    N = params%number_blocks

    recv_counter(:) = 0
    send_counter(:) = 0
    int_recv_counter(:) = 0
    int_send_counter(:) = 0

    do k = 1, hvy_n
        ! calculate light id
        sender_hvyID = hvy_active(k)
        call hvy_id_to_lgt_id( sender_lgtID, sender_hvyID, myrank, N )

        ! loop over all neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists
            if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
                ! neighbor light data id
                recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
                ! calculate neighbor rank
                call lgt_id_to_proc_rank( recver_rank, recver_lgtID, N )
                ! neighbor heavy id
                call lgt_id_to_hvy_id( recver_hvyID, recver_lgtID, recver_rank, N )

                ! define level difference: sender - receiver, so +1 means sender on higher level
                ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
                ! leveldiff =  0 : sender is same level as recver
                ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
                level_diff = lgt_block( sender_lgtID, params%max_treelevel + IDX_MESH_LVL ) - lgt_block( recver_lgtID, params%max_treelevel + IDX_MESH_LVL )


                if (recver_rank /= myrank .or. count_internal) then
                    ! it now depends on the stage if we have to sent this data
                    ! or not.
                    ! In stage 1, only level_diff = {+1, 0} is treated
                    ! In stage 2, only level_diff = -1

                    ! send counter. how much data will I send to other mpiranks?
                    if ((istage==1 .and. (level_diff==+1 .or. level_diff==0)) .or. (istage==2.and.level_diff==-1)) then
                        ! why is this RECVER and not sender? Well, complicated. The amount of data on the sender patch
                        ! is not the same as in the receiver patch, because we interpolate or downsample. We effectively
                        ! transfer only the data the recver wants - not the extra data.
                        ijk = ijkGhosts(:, :, neighborhood, level_diff, RECVER)

                        send_counter(recver_rank) = send_counter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
                    endif

                    ! recv counter. how much data will I recv from other mpiranks?
                    ! This is NOT the same number as before
                    if ((istage==1 .and. (-1*level_diff==+1 .or. -1*level_diff==0)) .or. (istage==2.and.-1*level_diff==-1)) then

                        inverse = inverse_neighbor(neighborhood, dim)

                        ijk = ijkGhosts(:, :, inverse, -1*level_diff, RECVER)

                        recv_counter(recver_rank) = recv_counter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
                    endif
                endif


                ! counter for integer buffer: for each neighborhood, we send 5 integers as metadata
                ! as this is a fixed number it does not depend on the type of neighborhood etc, so
                ! technically one would need only one for send/recv

                if ((istage==1 .and. (level_diff==+1 .or. level_diff==0)) .or. (istage==2.and.level_diff==-1)) then
                    int_send_counter(recver_rank) = int_send_counter(recver_rank) + 6 ! FIVE
                endif

                if ((istage==1 .and. (-1*level_diff==+1 .or. -1*level_diff==0)) .or. (istage==2.and.-1*level_diff==-1)) then
                    int_recv_counter(recver_rank) = int_recv_counter(recver_rank) + 6 ! FIVE
                endif

            end if ! neighbor exists
        end do ! loop over all possible  neighbors
    end do ! loop over all heavy active



    ! NOTE: for the int buffer, we mosly start at some index l0 and then loop unitl
    ! we find a -99 indicating the end of the buffer. this could be avoided by using
    ! for instead of while loops in the main routines, but I do not have time now.
    !
    ! In the meantime, notice we extent the amount of data by one, to copy the last -99
    ! to the buffers
    int_recv_counter(:) = int_recv_counter(:) + 1
    int_send_counter(:) = int_send_counter(:) + 1


    ! NOTE ACTUAL SEND / RECV DATA IS NEQN
    recv_counter(:) = recv_counter(:) * ncomponents
    send_counter(:) = send_counter(:) * ncomponents
end subroutine get_my_sendrecv_amount_with_ranks
