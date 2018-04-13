!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RMA_lock_unlock_get_data.f90
!> \version 0.4
!> \author msr
!
!> \brief get data from RMA window on other proc \n
!! use lock/unlock synchronization \n
!
!>
!! input:    
!!           - send buffer
!!           - com matrix and com position matrix
!!
!! output:   
!!           - receive buffer
!!
!! = log ======================================================================================
!! \n
!! 13/01/17 - create for v0.4
! ********************************************************************************************

subroutine RMA_lock_unlock_get_data( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )

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

    ! column number of send buffer, column number of receive buffer, real data buffer length
    integer(kind=ik)                    :: send_pos, receive_pos, real_pos

    ! RMA variables: windows, size variables
    integer(kind=ik)                    :: int_win, real_win
    integer(kind=MPI_ADDRESS_KIND)      :: int_length, real_length
    ! displacement variable
    integer(kind=MPI_ADDRESS_KIND)      :: displacement

    ! loop variable
    integer(kind=ik)                    :: k
    ! Communicator of WABBIT
    integer(kind=ik)                    :: WABBIT_COMM

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank            = params%rank
    number_procs    = params%number_procs
    WABBIT_COMM     = params%WABBIT_COMM
!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! first: create RMA windows
    ! note: size of window = size of send_buffer (2 dimensions!) * element size in byte

    ! integer window
    int_length = size(int_send_buffer,1) * size(int_send_buffer,2) * ik
    call MPI_Win_create( int_send_buffer(1,1), int_length, ik, MPI_INFO_NULL, WABBIT_COMM, int_win, ierr )

    ! real window
    real_length = size(real_send_buffer,1) * size(real_send_buffer,2) * rk
    call MPI_Win_create( real_send_buffer(1,1), real_length, rk, MPI_INFO_NULL, WABBIT_COMM, real_win, ierr )

    ! ----------------------------------------------------------------------------------------
    ! second: exchange data

    ! loop over corresponding com matrix line
    do k = 1, number_procs

        ! proc k-1 has to receive data -> so current proc also receive data from proc k-1
        if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

            ! send buffer column number, read from position matrix
            ! note: position in send buffer is written in neighboring proc line
            send_pos = com_matrix_pos(k, rank+1)

            ! receive buffer column number, read from position matrix
            ! note: position in receive buffer is written in own position matrix line
            receive_pos = com_matrix_pos(rank+1, k)

            ! set displacement
            ! note: displacement is special MPI integer kind, displace (column_pos-1) columns
            displacement = (send_pos-1)*size(int_receive_buffer,1)

            ! get integer data
            !-----------------
            ! synchronize RMA
            call MPI_Win_lock(MPI_LOCK_SHARED, k-1, 0, int_win, ierr)
            ! get data
            call MPI_Get( int_receive_buffer(1, receive_pos), size(int_receive_buffer,1), MPI_INTEGER4, k-1, displacement, size(int_receive_buffer,1), MPI_INTEGER4, int_win, ierr)
            ! synchronize RMA
            call MPI_Win_unlock(k-1, int_win, ierr)

            ! set displacement
            ! note: displacement is special MPI integer kind, displace (column_pos-1) columns
            displacement = (send_pos-1)*size(real_receive_buffer,1)

            ! real buffer length
            real_pos = int_receive_buffer(1, receive_pos)

            ! get real data
            ! -------------
            ! synchronize RMA
            call MPI_Win_lock(MPI_LOCK_SHARED, k-1, 0, real_win, ierr)
            ! get real data
            call MPI_Get( real_receive_buffer(1, receive_pos), real_pos, MPI_REAL8, k-1, displacement, real_pos, MPI_REAL8, real_win, ierr)
            ! synchronize RMA
            call MPI_Win_unlock(k-1, real_win, ierr)

        end if

    end do

    ! ----------------------------------------------------------------------------------------
    ! third: free windows
    call MPI_Win_free( int_win, ierr)
    call MPI_Win_free( real_win, ierr)

end subroutine RMA_lock_unlock_get_data
