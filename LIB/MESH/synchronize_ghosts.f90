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
! -------------------------------------------------------------------------------------------------------------------------
! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
! -------------------------------------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine synchronize_ghosts( params, block_list, block_data, neighbor_list )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)        :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! neighbor list
    integer(kind=ik), intent(in)        :: neighbor_list(:)

    ! loop variables
    integer(kind=ik)                    :: k, N, i, dF

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
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

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine copy_ghost_nodes( params, block_data, sender_id, receiver_id, neighborhood, level_diff )
            use module_params
            type (type_params), intent(in)      :: params
            real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)        :: sender_id, receiver_id
            integer(kind=ik), intent(in)        :: neighborhood
            integer(kind=ik), intent(in)        :: level_diff
        end subroutine copy_ghost_nodes

        subroutine send_receive_data( params, block_data, com_id, com_list, com_number, dF )
            use module_params
            type (type_params), intent(in)      :: params
            real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)        :: com_list(:, :)
            integer(kind=ik), intent(in)        :: com_id, com_number, dF
        end subroutine send_receive_data

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! allocate com_list: max linenumber -> max number of blocks * max number of neighbors * number of procs
    allocate( com_list( N*12*number_procs , 8), stat=allocate_error )
    ! allocate local com_list
    allocate( my_com_list( N*12*number_procs , 8), stat=allocate_error )
    ! allocate com_plan: max linenumber -> every proc has a neighborhood with every other proc
    !                                       * 2 (sender + receiver proc creates entry)
    allocate( com_plan( number_procs*number_procs*2), stat=allocate_error )

    ! reset com-list, com_plan
    com_list    = -1
    my_com_list = -1
    com_plan    = -1

!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! first: synchronize internal ghost nodes, create com_list for external communications

    ! loop over heavy data
    do k = 1, N
        ! block is active
        if ( block_list(rank*N + k , 1) /= -1 ) then
            ! loop over all neighbors
            do i = 1, 16
                ! neighbor exists
                if ( neighbor_list( (k - 1)*16 + i ) /= -1 ) then

                    ! neighbor light data id
                    neighbor_light_id = neighbor_list( (k - 1)*16 + i )
                    ! calculate the difference between block levels
                    level_diff = block_list( rank*N + k, params%max_treelevel+1 ) - block_list( neighbor_light_id, params%max_treelevel+1 )

                    ! proof if neighbor internal or external
                    if ( ( neighbor_light_id > rank*N ) .and. ( neighbor_light_id < (rank+1)*N+1 ) ) then
                        ! internal neighbor -> copy ghost nodes
                        call copy_ghost_nodes( params, block_data, k, neighbor_light_id-rank*N, i, level_diff )
                    else
                        ! external neighbor -> new com_list entry
                        my_com_list( rank*N*12 + (k-1)*12 + i , 1)  = k+i
                        my_com_list( rank*N*12 + (k-1)*12 + i , 2)  = rank
                        my_com_list( rank*N*12 + (k-1)*12 + i , 3)  = (neighbor_light_id - 1) / N
                        my_com_list( rank*N*12 + (k-1)*12 + i , 4)  = k
                        my_com_list( rank*N*12 + (k-1)*12 + i , 5)  = neighbor_light_id - ( (neighbor_light_id - 1) / N )*N
                        my_com_list( rank*N*12 + (k-1)*12 + i , 6)  = i
                        my_com_list( rank*N*12 + (k-1)*12 + i , 7)  = i
                        my_com_list( rank*N*12 + (k-1)*12 + i , 8)  = level_diff

                    end if

                end if
            end do
        end if
    end do

    ! synchronize com_list
    call MPI_Allreduce(my_com_list, com_list, N*12*number_procs*8, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)

    i = 1
    ! remove empty lines, loop over com_list
    do while ( i < N*12*number_procs )
        if ( com_list(i,1) == -1 ) then
            ! empty line

            ! last line reached
            if ( i == N*12*number_procs ) exit

            ! find next non empty line
            k = i + 1
            do while ( com_list(k, 1) == -1 )
                if ( k == N*12*number_procs ) then
                    ! end of com_list reached
                    exit
                else
                    k = k + 1
                end if
            end do

            ! sort list
            com_list( i:N*12*number_procs-(k-i) , : ) = com_list( k:N*12*number_procs , : )
            com_list( N*12*number_procs-(k-i):N*12*number_procs , : ) = -1

            ! next line
            i = i + 1
        else
            ! nothing to do, go to next line
            i = i + 1
        end if
    end do

    ! count number of communications
    n_com = 1
    do while ( com_list(n_com, 1) /= -1 )
        n_com = n_com + 1
    end do
    n_com = n_com - 1

    ! ----------------------------------------------------------------------------------------
    ! second: sort com_list, create com_plan
    call sort_com_list(com_list, size(com_list,1), com_plan, size(com_plan,1), number_procs, n_com)

    ! ----------------------------------------------------------------------------------------
    ! third: start external communications
    ! synchronize ghost nodes
    i       = 1
    k       = i
    ! loop over com_plan
    do while ( com_plan(i) /= -1 )

        if ( (com_list(k, 2) == rank) .or. (com_list(k, 3) == rank) ) then

            ! proc has to send/receive data, loop over all data fields
            do dF = 2, params%number_data_fields
                call send_receive_data( params, block_data, k, com_list, com_plan(i), dF)
            end do

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

end subroutine synchronize_ghosts
