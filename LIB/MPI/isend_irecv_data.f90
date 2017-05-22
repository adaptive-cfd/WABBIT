!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name isend_irecv_data.f90
!> \version 0.4
!> \author msr
!
!> \brief transfer data with MPI non blocking subroutines
!
!>
!! input:
!!           - send buffer
!!           - com matrix and com position matrix
!!
!! output:
!!           - receive buffer
!!
!!
!! = log ======================================================================================
!! \n
!! 16/01/17 - create for v0.4
! ********************************************************************************************

subroutine isend_irecv_data( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> send/receive buffer, integer and real
    integer(kind=ik), intent(in)        :: int_send_buffer(:,:)
    integer(kind=ik), intent(out)       :: int_receive_buffer(:,:)

    real(kind=rk), intent(in)           :: real_send_buffer(:,:)
    real(kind=rk), intent(out)          :: real_receive_buffer(:,:)

    !> communications matrix: neighboring proc rank
    !> com matrix pos: position in send buffer
    integer(kind=ik), intent(in)        :: com_matrix(:,:), com_matrix_pos(:,:)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! MPI status
    !integer                             :: status(MPI_status_size)

    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI request
    integer(kind=ik)                    :: send_request(size(com_matrix,1)), recv_request(size(com_matrix,1))

    ! column number of send buffer, column number of receive buffer, real data buffer length
    integer(kind=ik)                    :: send_pos, receive_pos, real_pos

    ! loop variable
    integer(kind=ik)                    :: k, i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank            = params%rank
    number_procs    = params%number_procs

    ! set message tag
    tag = 0

!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! first: integer data

    ! reset communication counter
    i = 0

    ! reset request arrays
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! loop over corresponding com matrix line
    do k = 1, number_procs

        ! communication between proc rank and proc k-1
        if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

            ! test output
            !write(*, '( "rank ", i3, " send/receive to/from rank " , i3)') rank, k-1

            ! increase communication counter
            i = i + 1

            tag = rank+1+k

            ! receive buffer column number, read from position matrix
            receive_pos = com_matrix_pos(rank+1, k)

            ! receive data
            call MPI_Irecv( int_receive_buffer(1, receive_pos), size(int_receive_buffer,1), MPI_INTEGER4, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)

            ! send buffer column number, read position matrix
            send_pos = com_matrix_pos(rank+1, k)

            ! send data
            call MPI_Isend( int_send_buffer(1, send_pos), size(int_send_buffer,1), MPI_INTEGER4, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)

        end if

    end do

    !> \todo Please check if waiting twice is really necessary
    ! synchronize non-blocking communications
    ! note: single status variable do not work with all compilers, so use MPI_STATUSES_IGNORE instead
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
    end if

    ! ----------------------------------------------------------------------------------------
    ! second: real data

    ! reset communication couter
    i = 0

    ! reset request arrays
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! loop over corresponding com matrix line
    do k = 1, number_procs

        ! communication between proc rank and proc k-1
        if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

            ! increase communication counter
            i = i + 1

            tag = number_procs*10*(rank+1+k)

            ! receive buffer column number, read from position matrix
            receive_pos = com_matrix_pos(rank+1, k)

            ! real buffer length
            real_pos = int_receive_buffer(1, receive_pos)

            ! receive data
            call MPI_Irecv( real_receive_buffer(1, receive_pos), real_pos, MPI_REAL8, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)
            !call MPI_Irecv( real_receive_buffer(1, receive_pos), real_pos, MPI_DOUBLE_PRECISION, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)

            ! send buffer column number, read position matrix
            send_pos = com_matrix_pos(rank+1, k)

            ! real buffer length
            real_pos = int_send_buffer(1, send_pos)

            ! send data
            call MPI_Isend( real_send_buffer(1, send_pos), real_pos, MPI_REAL8, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)
            !call MPI_Isend( real_send_buffer(1, send_pos), real_pos, MPI_DOUBLE_PRECISION, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)

        end if

    end do

    ! synchronize non-blocking communications
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
    end if

end subroutine isend_irecv_data
