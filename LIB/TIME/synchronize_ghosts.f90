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

    ! communication list: column
    !                       1   com_list id /= -1 => non empty line
    !                       2   rank of sender process
    !                       3   rank of receiver process
    !                       4   sender block heavy data id
    !                       5   receiver block heavy data id
    !                       6   sender block neighborhood to receiver (dirs id)
    !                       7   receiver block neighborhood to sender (dirs id)
    !                       8   difference between sender-receiver level
    integer(kind=ik), allocatable       :: com_list(:, :), my_com_list(:, :)

    ! communications plan: entrys are the number of com_list lines to fill send/receive buffer
    integer(kind=ik), allocatable       :: com_plan(:)

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

    ! allocate com_list: max linenumber -> max number of blocks * max number of neighbors * number of procs
    allocate( com_list( N*16*number_procs , 8), stat=allocate_error )
    ! allocate local com_list
    allocate( my_com_list( N*16*number_procs , 8), stat=allocate_error )
    ! allocate com_plan: max linenumber -> every proc has a neighborhood with every other proc
    !                                       * 2 (sender + receiver proc creates entry)
    allocate( com_plan( number_procs*number_procs*2), stat=allocate_error )

    ! allocate com matrix
    allocate( com_matrix(number_procs, number_procs), stat=allocate_error )
    allocate( my_com_matrix(number_procs, number_procs), stat=allocate_error )

    ! reset com-list, com_plan, com matrix
    com_list        = -1
    my_com_list     = -1
    com_plan        = -1
    com_matrix      =  0
    my_com_matrix   =  0

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
                    ! external neighbor -> new com_list entry
                    my_com_list( rank*N*16 + (k-1)*16 + i , 1)  = k+i
                    my_com_list( rank*N*16 + (k-1)*16 + i , 2)  = rank
                    my_com_list( rank*N*16 + (k-1)*16 + i , 3)  = neighbor_rank
                    my_com_list( rank*N*16 + (k-1)*16 + i , 4)  = hvy_active(k)
                    my_com_list( rank*N*16 + (k-1)*16 + i , 5)  = hvy_id
                    my_com_list( rank*N*16 + (k-1)*16 + i , 6)  = i
                    my_com_list( rank*N*16 + (k-1)*16 + i , 7)  = i
                    my_com_list( rank*N*16 + (k-1)*16 + i , 8)  = level_diff

                    ! write communications matrix
                    my_com_matrix(rank+1, neighbor_rank+1) = my_com_matrix(rank+1, neighbor_rank+1) + 1

                end if

            end if
        end do

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
            if ( debug%name_comp_time(k) == "synch. ghosts - internal" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - internal"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

    ! start time
    sub_t0 = MPI_Wtime()

    ! synchronize com_list
    call MPI_Allreduce(my_com_list, com_list, N*16*number_procs*8, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)

    k = 1
    my_com_list = -1
    ! loop over com_list
    do i = 1, size(com_list, 1)
       if ( com_list(i,1) /= -1 ) then
           ! non empty element
           my_com_list(k,:) = com_list(i,:)
           k = k + 1
       end if
    end do
    com_list = my_com_list


    ! number of communications
    n_com = k - 1

    ! ----------------------------------------------------------------------------------------
    ! second: sort com_list, create com_plan
    call sort_com_list(com_list, size(com_list,1), com_plan, size(com_plan,1), number_procs, n_com)

    ! debug writing
    if ( params%debug ) then
        call write_com_list( com_list )
    end if

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts - com. list" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - com. list"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

    ! start time
    sub_t0 = MPI_Wtime()

    ! ----------------------------------------------------------------------------------------
    ! third: start external communications
    ! synchronize ghost nodes

    i       = 1
    k       = i
    ! loop over com_plan
    do while ( com_plan(i) /= -1 )

        if ( (com_list(k, 2) == rank) .or. (com_list(k, 3) == rank) ) then

            ! proc has to send/receive data, send all datafields
            call send_receive_data( params, hvy_block, k, com_list, com_plan(i))

            ! next step in com_plan
            k = k + 2*com_plan(i)
            i = i + 2

        else

            ! nothing to do, go to next communication
            k = k + 2*com_plan(i)
            i = i + 2

        end if

    end do

    ! clean up
    deallocate( com_list, stat=allocate_error )
    deallocate( my_com_list, stat=allocate_error )
    deallocate( com_plan, stat=allocate_error )
    deallocate( com_matrix, stat=allocate_error )
    deallocate( my_com_matrix, stat=allocate_error )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts - external" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - external"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine synchronize_ghosts
