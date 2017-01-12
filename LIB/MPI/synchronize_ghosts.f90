! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: synchronize_ghosts.f90
! version: 0.4
! author: msr
!
! synchronize ghosts nodes
!
! input:    - params, light and heavy data
! output:   - heavy data array
!
! todo: change soubroutine, to work only on one datafield, not on all to the same time
!
! -------------------------------------------------------------------------------------------------------------------------
! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
! -------------------------------------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! 06/01/17 - use RMA to synchronize data
! ********************************************************************************************

subroutine synchronize_ghosts(  params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    ! list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: k, N, i, dF, lgt_id, hvy_id

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank, neighbor_rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! neighbor light data id
    integer(kind=ik)                    :: neighbor_light_id
    ! difference between current block and neighbor block level
    integer(kind=ik)                    :: level_diff

    ! communication lists:
    ! dim 1: list elements
    ! dim 2: columns
    !                       1   rank of sender process
    !                       2   rank of receiver process
    !                       3   sender block heavy data id
    !                       4   receiver block heavy data id
    !                       5   sender block neighborhood to receiver (dirs id)
    !                       6   difference between sender-receiver level
    ! dim 3: receiver proc rank
    integer(kind=ik), allocatable       :: com_lists(:, :, :)

    ! receiver lists: receiver_pos [position in rank and count list],
    ! receiver_rank [proc rank list], receiver_count [number of communications to receiver]
    ! receiver_N [number of neighbor procs]
    integer(kind=ik), allocatable       :: receiver_pos(:), receiver_rank(:), receiver_count(:)
    integer(kind=ik)                    :: receiver_N

    ! send/receive buffer, integer and real
    integer(kind=ik), allocatable       :: int_send_buffer(:), int_receive_buffer(:)
    real(kind=rk), allocatable          :: real_send_buffer(:), real_receive_buffer(:), proc_send_buffer(:)
    ! index of send buffer, return from proc to proc send buffer subroutine
    integer(kind=ik)                    :: buffer_i

    ! position marker in buffer arrays, length of buffer array part
    ! note: use standard integers for position variables to use them later for RMA subroutines
    integer                             :: int_pos, real_pos, int_N, real_N


    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! number of communications
    integer(kind=ik)                    :: n_com

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    integer(kind=ik), allocatable       :: com_matrix(:,:), my_com_matrix(:,:)

    ! RMA variables: windows, size variables
    integer(kind=ik)                    :: int_win, real_win
    integer(kind=MPI_ADDRESS_KIND)      :: int_length, real_length
    integer(kind=ik)                    :: int_win_start, int_win_end, real_win_start, real_win_end

    ! RMA win loop switch
    logical                             :: win_loop
    ! dummy displacement variable
    integer(kind=MPI_ADDRESS_KIND)      :: displacement

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! allocate local com_lists
    allocate( com_lists( N*16, 6, number_procs), stat=allocate_error )

    ! receiver lists
    allocate( receiver_pos( number_procs), stat=allocate_error )
    allocate( receiver_rank( number_procs), stat=allocate_error )
    allocate( receiver_count( number_procs), stat=allocate_error )

    ! allocate com matrix
    allocate( com_matrix(number_procs, number_procs), stat=allocate_error )
    allocate( my_com_matrix(number_procs, number_procs), stat=allocate_error )

    ! reset com-list, com_plan, com matrix, receiver lists
    com_lists       = -1

    com_matrix      =  0
    my_com_matrix   =  0

    receiver_pos    =  0
    receiver_rank   = -1
    receiver_count  =  0
    receiver_N      =  0

    ! reset ghost nodes for all active blocks
    ! loop over all datafields
    do dF = 2, params%number_data_fields+1
        ! loop over all active blocks
        do k = 1, hvy_n
            ! reset ghost nodes
            hvy_block(1:g, :, dF, hvy_active(k) )           = 9.0e9_rk
            hvy_block(Bs+g+1:Bs+2*g, :, dF, hvy_active(k) ) = 9.0e9_rk
            hvy_block(:, 1:g, dF, hvy_active(k) )           = 9.0e9_rk
            hvy_block(:, Bs+g+1:Bs+2*g, dF, hvy_active(k) ) = 9.0e9_rk
        end do
    end do

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! ----------------------------------------------------------------------------------------
    ! first: synchronize internal ghost nodes, create com_list for external communications

    ! loop over active heavy data
    do k = 1, hvy_n

        ! loop over all neighbors
        do i = 1, 16
            ! neighbor exists
            if ( hvy_neighbor( hvy_active(k), i ) /= -1 ) then

                ! neighbor light data id
                neighbor_light_id = hvy_neighbor( hvy_active(k), i )
                ! calculate light id
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )
                ! calculate the difference between block levels
                level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_light_id, params%max_treelevel+1 )

                ! proof if neighbor internal or external
                call lgt_id_to_proc_rank( neighbor_rank, neighbor_light_id, N )

                if ( rank == neighbor_rank ) then
                    ! calculate internal heavy id
                    call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, rank, N )

                    ! internal neighbor -> copy ghost nodes
                    call copy_ghost_nodes( params, hvy_block, hvy_active(k), hvy_id, i, level_diff )

                    ! write communications matrix
                    my_com_matrix(rank+1, rank+1) = my_com_matrix(rank+1, rank+1) + 1

                else
                    ! neighbor heavy id
                    call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, neighbor_rank, N )

                    ! check neighbor proc rank
                    if ( receiver_pos(neighbor_rank+1) == 0 ) then

                        ! first communication with neighbor proc
                        ! -------------------------------------------
                        ! set list position, increase number of neighbor procs by 1
                        receiver_N                      = receiver_N + 1
                        ! save list pos
                        receiver_pos(neighbor_rank+1)   = receiver_N
                        ! save neighbor rank
                        receiver_rank(receiver_N)       = neighbor_rank
                        ! count communications - here: first one
                        receiver_count(receiver_N)      = 1

                        ! external neighbor -> new com_lists entry (first entry)
                        com_lists( 1 , 1, neighbor_rank+1)  = rank
                        com_lists( 1 , 2, neighbor_rank+1)  = neighbor_rank
                        com_lists( 1 , 3, neighbor_rank+1)  = hvy_active(k)
                        com_lists( 1 , 4, neighbor_rank+1)  = hvy_id
                        com_lists( 1 , 5, neighbor_rank+1)  = i
                        com_lists( 1 , 6, neighbor_rank+1)  = level_diff

                    else

                        ! additional communication with neighbor proc
                        ! -------------------------------------------
                        ! count communications - +1
                        receiver_count( receiver_pos(neighbor_rank+1) )  = receiver_count( receiver_pos(neighbor_rank+1) ) + 1

                        ! external neighbor -> new com_lists entry
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 1, neighbor_rank+1)  = rank
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 2, neighbor_rank+1)  = neighbor_rank
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 3, neighbor_rank+1)  = hvy_active(k)
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 4, neighbor_rank+1)  = hvy_id
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 5, neighbor_rank+1)  = i
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 6, neighbor_rank+1)  = level_diff

                    end if

                end if

            end if
        end do

    end do

    ! write my com matrix, loop over number of receiver procs, write counted communications
    ! additional: sum communcitions
    n_com = 0
    do k = 1, receiver_N
        ! write matrix
        my_com_matrix( rank+1, receiver_rank(k)+1 ) = receiver_count( receiver_pos( receiver_rank(k)+1 ) )
        ! sum communications
        n_com = n_com + receiver_count( receiver_pos( receiver_rank(k)+1 ) )
    end do

    ! synchronize com matrix
    call MPI_Allreduce(my_com_matrix, com_matrix, number_procs*number_procs, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! save com matrix
    if ( params%debug ) then
        call write_com_matrix( com_matrix )
    end if

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts - internal and com_matrix" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - internal and com_matrix"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

    ! start time
    sub_t0 = MPI_Wtime()

    ! next steps only for more than two procs
    if ( number_procs > 1 ) then

        ! ----------------------------------------------------------------------------------------
        ! second: allocate memory for send/receive buffer
        ! buffer length:
        !                   int buffer  - number of receiver procs * sum of communications * 3
        !                   real buffer - number of receiver procs * sum of communications * (Bs+g) * g
        !                   proc buffer - number of receiver procs * sum of communications * (Bs+g) * g

        allocate( int_send_buffer( receiver_N * (n_com+1) * 3 + 1), stat=allocate_error )
        allocate( int_receive_buffer( receiver_N * (n_com+1) * 3 + 1), stat=allocate_error )

        allocate( real_send_buffer( receiver_N * params%number_data_fields *  n_com * (Bs+g) * g + 1), stat=allocate_error )
        allocate( real_receive_buffer( receiver_N * params%number_data_fields * n_com * (Bs+g) * g + 1), stat=allocate_error )

        allocate( proc_send_buffer( params%number_data_fields * n_com * (Bs+g) * g + 1), stat=allocate_error )

        real_send_buffer        = 7.0e9_rk
        real_receive_buffer     = 5.0e9_rk

        ! ----------------------------------------------------------------------------------------
        ! third: fill send buffer
        ! int buffer:  store receiver block id, neighborhood and level difference (in order of neighbor proc rank, use com matrix)
        ! real buffer: store block data (in order of neighbor proc rank, use com matrix)

        ! reset int and real buffer positions
        int_pos  = 1
        real_pos = 1

        ! loop over corresponding com matrix line
        do k = 1, number_procs

            ! other proc has to receive data
            if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

                ! first: real data
                ! ----------------
                ! number of communications: com_matrix(rank+1, k)
                i = com_matrix(rank+1, k)

                ! write real send buffer for proc k
                call create_send_buffer(params, hvy_block, com_lists( 1:i, :, k), i, proc_send_buffer, buffer_i)

                ! real buffer entry
                ! note: buffer_i - 1 => last data in buffer array
                real_send_buffer( real_pos : real_pos + buffer_i-2 ) = proc_send_buffer( 1 : buffer_i-1 )

                ! second: integer data
                ! --------------------
                ! loop over all communications to this proc
                do i = 1, com_matrix(rank+1, k)

                    ! int buffer entry: neighbor block id, neighborhood, level difference
                    int_send_buffer(int_pos)   = com_lists( i, 4, k)
                    int_send_buffer(int_pos+1) = com_lists( i, 5, k)
                    int_send_buffer(int_pos+2) = com_lists( i, 6, k)
                    ! increase int buffer position
                    int_pos = int_pos + 3

                end do

                ! third: save real buffer position
                ! ------------------------------------
                int_send_buffer(int_pos)   = real_pos
                int_send_buffer(int_pos+1) = real_pos + buffer_i-2
                int_pos = int_pos + 2

                ! fourth: increase real buffer position
                ! ------------------------------------
                real_pos = real_pos + buffer_i-1

            end if

        end do

        ! end time
        sub_t1 = MPI_Wtime()
        ! write time
        if ( params%debug ) then
            ! find free or corresponding line
            k = 1
            do while ( debug%name_comp_time(k) /= "---" )
                ! entry for current subroutine exists
                if ( debug%name_comp_time(k) == "synch. ghosts - fill send buffer" ) exit
                k = k + 1
            end do
            ! write time
            debug%name_comp_time(k) = "synch. ghosts - fill send buffer"
            debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
            debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
        end if

        ! start time
        sub_t0 = MPI_Wtime()

        ! ----------------------------------------------------------------------------------------
        ! fourth: open RMA window for int and real buffer

        ! integer window
        int_length = (int_pos-1)*ik
        call MPI_Win_create( int_send_buffer(1), int_length, ik, MPI_INFO_NULL, MPI_COMM_WORLD, int_win, ierr )

        ! real window
        real_length = (real_pos-1)*rk
        call MPI_Win_create( real_send_buffer(1), real_length, rk, MPI_INFO_NULL, MPI_COMM_WORLD, real_win, ierr )

        ! ----------------------------------------------------------------------------------------
        ! fifth: get data

        ! integer data
        ! ------------

        ! reset int_receive_buffer position
        int_pos  = 1
        real_pos = 1

        ! loop over corresponding com matrix line
        do k = 1, number_procs

            ! proc k-1 has to receive data -> so current proc also receive data from proc k-1
            if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

                ! reset int window start/end index
                int_win_start = 0
                int_win_end   = 0

                ! loop over com matrix line of sender proc to calculate start/end id in integer window
                i        = 1
                win_loop = .true.

                do while (win_loop)
                    ! proc k-1 send data
                    if ( ( com_matrix(k, i) > 0 ) .and. ( i /= k ) ) then

                        if ( (rank+1) == i ) then
                            ! entry in com matrix correspond to current proc
                            int_win_start = int_win_end + 1
                            int_win_end   = int_win_start + com_matrix(k, i) * 3 + 1
                            ! end loop
                            win_loop = .false.

                        else
                            ! other proc receive data
                            int_win_start = int_win_end + 1
                            int_win_end   = int_win_start + com_matrix(k, i) * 3 + 1

                        end if

                        ! next com matrix element
                        i = i + 1

                    else
                        ! next com matrix element
                        i = i + 1

                    end if
                end do

                ! calculate displacement
                ! displace first integers, note: displacement is special MPI integer kind, so we cant use int_win_start directly
                displacement = int_win_start-1

                ! get integer data
                !-----------------
                ! synchronize RMA
                call MPI_Win_lock(MPI_LOCK_SHARED, k-1, 0, int_win, ierr)
                ! get data
                call MPI_Get( int_receive_buffer(int_pos), (int_win_end-int_win_start)+1, MPI_INTEGER4, k-1, displacement, (int_win_end-int_win_start)+1, MPI_INTEGER4, int_win, ierr)
                ! synchronize RMA
                call MPI_Win_unlock(k-1, int_win, ierr)

                ! new integer position
                int_pos = int_pos + (int_win_end-int_win_start) + 1

                ! start and end of real data in buffer
                real_win_start = int_receive_buffer(int_pos-2)
                real_win_end   = int_receive_buffer(int_pos-1)

                ! calculate displacement
                ! displace first integers, note: displacement is special MPI integer kind, so we cant use real_win_start directly
                displacement = real_win_start-1

                ! get real data
                ! -------------
                ! synchronize RMA
                call MPI_Win_lock(MPI_LOCK_SHARED, k-1, 0, real_win, ierr)
                ! get real data
                call MPI_Get( real_receive_buffer(real_pos), (real_win_end-real_win_start)+1, MPI_REAL8, k-1, displacement, (real_win_end-real_win_start)+1, MPI_REAL8, real_win, ierr)
                ! synchronize RMA
                call MPI_Win_unlock(k-1, real_win, ierr)

                ! new real position
                real_pos = real_pos + (real_win_end-real_win_start) + 1

            end if

        end do

        ! barrier
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        ! free windows
        call MPI_Win_free( int_win, ierr)
        call MPI_Win_free( real_win, ierr)

        ! end time
        sub_t1 = MPI_Wtime()
        ! write time
        if ( params%debug ) then
            ! find free or corresponding line
            k = 1
            do while ( debug%name_comp_time(k) /= "---" )
                ! entry for current subroutine exists
                if ( debug%name_comp_time(k) == "synch. ghosts - RMA" ) exit
                k = k + 1
            end do
            ! write time
            debug%name_comp_time(k) = "synch. ghosts - RMA"
            debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
            debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
        end if

        ! start time
        sub_t0 = MPI_Wtime()

        ! ----------------------------------------------------------------------------------------
        ! sixth: write receive buffer to heavy data

        ! reset int and real position integers
        int_pos  = 1
        real_pos = 1

        ! loop over corresponding com matrix line
        do k = 1, number_procs

            ! received data from proc k-1
            if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

                ! calculate integer length
                int_N  = com_matrix(rank+1, k) * 3
                real_N = int_receive_buffer(int_pos + int_N + 1) - int_receive_buffer(int_pos + int_N) + 1

                ! read received data
                call write_receive_buffer(params, int_receive_buffer(int_pos:int_pos+int_N-1), real_receive_buffer(real_pos:real_pos+real_N-1), hvy_block)

                ! new int and real positions
                int_pos  = int_pos + int_N + 2
                real_pos = real_pos + real_N

            end if

        end do

    end if

    ! clean up
    deallocate( com_lists, stat=allocate_error )
    deallocate( receiver_pos, stat=allocate_error )
    deallocate( receiver_rank, stat=allocate_error )
    deallocate( receiver_count, stat=allocate_error )

    deallocate( com_matrix, stat=allocate_error )
    deallocate( my_com_matrix, stat=allocate_error )

    deallocate( int_send_buffer, stat=allocate_error )
    deallocate( int_receive_buffer, stat=allocate_error )
    deallocate( real_send_buffer, stat=allocate_error )
    deallocate( real_receive_buffer, stat=allocate_error )
    deallocate( proc_send_buffer, stat=allocate_error )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts - write external" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - write external"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine synchronize_ghosts
