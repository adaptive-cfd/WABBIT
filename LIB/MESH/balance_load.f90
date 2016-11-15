! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: balance_load.f90
! version: 0.4
! author: msr
!
! balance the load
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine balance_load( params, block_list, block_data )

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
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)

    ! send/receive buffer, note: size is equal to block data array, because if a block want to send all his data
    real(kind=rk)                       :: buffer_data( size(block_data,1), size(block_data,2), size(block_data,3), size(block_data,4) )
    integer(kind=ik)                    :: buffer_light( params%number_blocks )

    ! light data list for working
    integer(kind=ik)                    :: my_block_list( size(block_list, 1), params%max_treelevel+2)
    ! light id start
    integer(kind=ik)                    :: my_light_start

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI status
    integer                             :: status(MPI_status_size)

    ! distribution type
    character(len=80)                   :: distribution

    ! block distribution lists
    integer(kind=ik), allocatable       :: opt_dist_list(:), dist_list(:)

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! loop variables
    integer(kind=ik)                    :: k, N, num_blocks, l, com_i, com_N

    ! com plan
    integer(kind=ik), allocatable       :: com_plan(:,:)

    ! size of data array
    integer(kind=ik)                    :: data_size

    ! free light/heavy data id
    integer(kind=ik)                    :: free_light_id, free_heavy_id

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    tag = 0

    distribution    = params%block_distribution

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! allocate block to proc lists
    allocate( opt_dist_list( number_procs ), stat=allocate_error )
    allocate( dist_list( number_procs ), stat=allocate_error )

    ! allocate com plan, maximal number of communications: every proc send every other proc something
    allocate( com_plan( number_procs*(number_procs-1), 3 ), stat=allocate_error )
    com_plan = -1

    ! number of light data
    N = size(block_list, 1)

    ! reset block count
    num_blocks = 0

    ! reset block to proc list
    dist_list = 0

    ! light data start line
    my_light_start = rank*params%number_blocks

    ! set light data list for working, only light data coresponding to proc are not zero
    my_block_list = 0
    my_block_list( my_light_start+1: my_light_start+params%number_blocks, :) = block_list( my_light_start+1: my_light_start+params%number_blocks, :)

    ! size of data array, use for readability
    data_size = size(block_data,1) * size(block_data,2) * size(block_data,3)

    ! reset send/receive buffer
    buffer_data = 9.0e9_rk

!---------------------------------------------------------------------------------------------
! main body

    select case(distribution)
        case("equal")
            ! simple uniformly distribution

            !---------------------------------------------------------------------------------
            ! first: calculate optimal block distribution
            ! count number of active blocks and current block distribution
            do k = 1, N
                ! block is active
                if ( block_list(k, 1) /= -1 ) then
                    ! count block for proc corresponing to light data
                    dist_list( (k-1) / params%number_blocks + 1 ) = dist_list( (k-1) / params%number_blocks + 1 ) + 1
                end if
            end do
            ! count blocks
            do k = 1, number_procs
                num_blocks = num_blocks + dist_list(k)
            end do

            ! optimal distribution
            opt_dist_list = (num_blocks - mod(num_blocks, number_procs))/number_procs
            ! distribute remaining blocks
            if (mod(num_blocks, number_procs) > 0) then
                opt_dist_list(1:mod(num_blocks, number_procs)) = (num_blocks - mod(num_blocks, number_procs))/number_procs + 1
            end if

            !---------------------------------------------------------------------------------
            ! second: create plan for communciations
            ! column
            !    1     sender proc
            !    2     receiver proc
            !    3     number of blocks to send

            com_i = 1
            ! loop over all procs
            do k = 1, number_procs

                ! proc want send data
                if ( (opt_dist_list(k) - dist_list(k)) < 0 ) then

                    ! find a receiver
                    do l = 1, number_procs

                        ! proc can receive data
                        if ( ( (opt_dist_list(l) - dist_list(l)) > 0 ) .and. ( l /= k ) ) then
                            ! calculate number of send/receive blocks
                            if ( ( abs(opt_dist_list(k) - dist_list(k)) - abs((opt_dist_list(l) - dist_list(l))) ) <= 0 ) then
                                ! receiver could receive all blocks from sender
                                com_N = abs(opt_dist_list(k) - dist_list(k))

                            elseif ( ( abs(opt_dist_list(k) - dist_list(k)) - abs((opt_dist_list(l) - dist_list(l))) ) > 0 ) then
                                ! receiver could not receive all blocks from sender,
                                ! so he should receive as much blocks he can take
                                com_N = abs((opt_dist_list(l) - dist_list(l)))

                            end if

                            ! create com plan
                            com_plan( com_i, 1 ) = k-1
                            com_plan( com_i, 2 ) = l-1
                            com_plan( com_i, 3 ) = com_N

                            ! change distribution list
                            dist_list(k) = dist_list(k) - com_N
                            dist_list(l) = dist_list(l) + com_N

                            ! new com was created
                            com_i = com_i + 1

                        end if

                        ! no blocks left for sending, note: a proc can have one block more than all other procs
                        if ( (opt_dist_list(k) - dist_list(k)) >= -1 ) exit

                    end do
                end if
            end do

            !---------------------------------------------------------------------------------
            ! third: send/receive data

            ! loop over all planed communications
            do k = 1, com_i-1

                if ( com_plan(k, 1) == rank ) then
                    ! proc have to send data
                    ! create send buffer, loop over all blocks
                    com_N = 1
                    l     = 1

                    do while ( com_N <= com_plan(k, 3) )

                        ! block is active
                        if ( (block_list( my_light_start + l, 1) /= -1) ) then
                            ! fill buffer, heavy data
                            buffer_data(:, :, :, com_N) = block_data(:, :, :, l)
                            ! ... light data
                            buffer_light(com_N) = my_light_start + l
                            ! count com
                            com_N = com_N + 1

                            ! delete heavy data
                            block_data(:, :, :, l) = 0.0_rk
                            ! delete light data
                            my_block_list( my_light_start + l, : ) = -1

                        end if
                        ! next block
                        l     = l + 1

                    end do

                    ! error case
                    if ( com_N-1 /= com_plan(k, 3) ) then
                        print*, "ERROR: load balancing error"
                        stop
                    end if

                    ! send data
                    call MPI_Send( buffer_light, (com_N-1), MPI_INTEGER4, com_plan(k, 2), tag, MPI_COMM_WORLD, ierr)
                    call MPI_Send( buffer_data, data_size*(com_N-1), MPI_REAL8, com_plan(k, 2), tag, MPI_COMM_WORLD, ierr)

                elseif ( com_plan(k, 2) == rank ) then
                    ! proc have to receive data
                    ! receive data
                    call MPI_Recv( buffer_light, com_plan(k, 3), MPI_INTEGER4, com_plan(k, 1), tag, MPI_COMM_WORLD, status, ierr)
                    call MPI_Recv( buffer_data, data_size*com_plan(k, 3), MPI_REAL8, com_plan(k, 1), tag, MPI_COMM_WORLD, status, ierr)

                    ! loop over all received blocks
                    do l = 1,  com_plan(k, 3)

                        ! find free light id
                        call get_free_light_id( free_heavy_id, my_block_list( my_light_start+1 : my_light_start+params%number_blocks , 1 ), params%number_blocks )
                        ! calculate light id
                        free_light_id = my_light_start + free_heavy_id

                        ! write light data
                        my_block_list( free_light_id, :) = block_list( buffer_light(l), : )

                        ! write heavy data
                        block_data(:, :, :, free_heavy_id) = buffer_data(:, :, :, l)

                    end do

                else
                    ! nothing to do
                end if

            end do

            ! synchronize light data
            block_list = 0
            call MPI_Allreduce(my_block_list, block_list, size(block_list,1)*size(block_list,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

        case("sfc1")

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: block distribution scheme is unknown"
            write(*,*) distribution
            stop

    end select

    ! clean up
    deallocate( opt_dist_list, stat=allocate_error )
    deallocate( dist_list, stat=allocate_error )
    deallocate( com_plan, stat=allocate_error )

end subroutine balance_load
