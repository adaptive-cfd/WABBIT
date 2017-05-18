!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name synchronize_external_nodes.f90
!> \version 0.5
!> \author msr
!
!> \brief synchronize external ghosts nodes
!
!> \details
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
!
!> \details
!! = log ======================================================================================
!! \n
!! 09/05/17 - create
!
! ********************************************************************************************

subroutine synchronize_external_nodes(  params, hvy_block, com_lists, com_matrix )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

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
    integer(kind=ik), intent(inout)        :: com_lists(:, :, :)

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    ! com matrix pos: position in send buffer
    integer(kind=ik), intent(inout)        :: com_matrix(:,:)

    integer(kind=ik), allocatable       :: com_matrix_pos(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, N, i, j!, dF

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! send/receive buffer, integer and real
    integer(kind=ik), allocatable       :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable          :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! length of buffer array and column number in buffer, use for readability
    integer(kind=ik)                    :: int_N, real_N, buffer_pos

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! number of communications, number of neighboring procs
    integer(kind=ik)                    :: my_n_com, n_com, n_procs

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

    ! variable for non-uniform mesh correction: remove redundant node between fine->coarse blocks
    integer(kind=ik)                    :: rmv_redundant

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    rmv_redundant = 0

    ! set non-uniform mesh correction
    if ( params%non_uniform_mesh_correction ) then
        rmv_redundant = 1
    else
        rmv_redundant = 0
    end if

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

    allocate( com_matrix_pos(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    com_matrix_pos  =  0

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! max number of communications and neighboring procs - use for buffer allocation
    ! every proc loop over com matrix line
    call max_com_num( my_n_com, n_procs, com_matrix(rank+1,:), rank )

    ! synchronize max com number, because if not:
    ! RMA buffer displacement is not fixed, so we need synchronization there
    call MPI_Allreduce(my_n_com, n_com, 1, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! for proc without neighbors: set n_procs to 1
    ! so we allocate arrays with second dimension=1
    if (n_procs==0) n_procs = 1

    ! next steps only for more than two procs
    if ( number_procs > 1 ) then

        ! ----------------------------------------------------------------------------------------
        ! second: allocate memory for send/receive buffer
        ! buffer size:
        !               number of columns: number of neighboring procs
        !               number of lines:
        !                   int buffer  - max number of communications  * 3 + 1 (length of real buffer)
        !                   real buffer - max number of communications * (Bs+g) * g * number of datafields
        !                   for 3D: real_buffer * Bs
        allocate( int_send_buffer( n_com * 3 + 1, n_procs ), stat=allocate_error )
        !call check_allocation(allocate_error)
        if ( allocate_error /= 0 ) then
            write(*,'(80("_"))')
            write(*,*) "ERROR: memory allocation fails"
            stop
        end if
        allocate( int_receive_buffer( n_com * 3 + 1, n_procs ), stat=allocate_error )
        !call check_allocation(allocate_error)
        if ( allocate_error /= 0 ) then
            write(*,'(80("_"))')
            write(*,*) "ERROR: memory allocation fails"
            stop
        end if

        if ( params%threeD_case ) then
            ! 3D:
            allocate( real_receive_buffer( n_com * (Bs+(g+rmv_redundant)) * (g+rmv_redundant) * Bs * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
            allocate( real_send_buffer( n_com * (Bs+(g+rmv_redundant)) * (g+rmv_redundant) * Bs * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
        else
            ! 2D:
            allocate( real_receive_buffer( n_com * (Bs+(g+rmv_redundant)) * (g+rmv_redundant) * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
            allocate( real_send_buffer( n_com * (Bs+(g+rmv_redundant)) * (g+rmv_redundant) * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
        end if

!        ! reset buffer for debuggung
!        if ( params%debug ) then
!            real_send_buffer        = 7.0e9_rk
!            real_receive_buffer     = 5.0e9_rk
!        end if

        ! ----------------------------------------------------------------------------------------
        ! third: fill send buffer
        ! int buffer:  store receiver block id, neighborhood and level difference (in order of neighbor proc rank, use com matrix)
        ! real buffer: store block data (in order of neighbor proc rank, use com matrix)
        ! first element of int buffer = length of real buffer (buffer_i)

        ! fill send buffer and position communication matrix
        call fill_send_buffer( params, hvy_block, com_lists, com_matrix(rank+1,:), rank, int_send_buffer, real_send_buffer )

        ! calculate position matrix: position is column in send buffer, so simply count the number of communications
        ! loop over all com_matrix elements
        do i = 1, size(com_matrix_pos,1)
            ! new line, means new proc: reset counter
            k = 1
            ! loop over communications
            do j = 1, size(com_matrix_pos,1)
                ! found external communication
                if ( (com_matrix(i,j) /= 0) .and. (i /= j) ) then
                    ! save com position
                    com_matrix_pos(i,j) = k
                    ! increase counter
                    k = k + 1

                end if
            end do
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
        ! fourth: get data for receive buffer

        ! communicate, fill receive buffer
        call fill_receive_buffer( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos  )

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
        ! fifth: write receive buffer to heavy data

        ! loop over corresponding com matrix line
        do k = 1, number_procs

            ! received data from proc k-1
            if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

                ! set buffer position and calculate length if integer/real buffer
                buffer_pos = com_matrix_pos(rank+1, k)
                int_N  = com_matrix(rank+1, k) * 3 + 1
                real_N = int_receive_buffer( 1, buffer_pos )

                ! read received data
                if ( params%threeD_case ) then
                    ! 3D:
                    call write_receive_buffer_3D(params, int_receive_buffer(2:int_N, buffer_pos), real_receive_buffer(1:real_N, buffer_pos), hvy_block )
                else
                    ! 2D:
                    call write_receive_buffer_2D(params, int_receive_buffer(2:int_N, buffer_pos), real_receive_buffer(1:real_N, buffer_pos), hvy_block(:, :, 1, :, :) )
                end if

            end if

        end do

    end if

    ! clean up
    deallocate( com_matrix_pos, stat=allocate_error )

    deallocate( int_send_buffer, stat=allocate_error )
    deallocate( int_receive_buffer, stat=allocate_error )
    deallocate( real_send_buffer, stat=allocate_error )
    deallocate( real_receive_buffer, stat=allocate_error )

end subroutine synchronize_external_nodes
