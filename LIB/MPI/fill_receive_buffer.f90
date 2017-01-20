! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: fill_receive_buffer.f90
! version: 0.4
! author: msr
!
! fill receive buffer
!
! input:    -
! output:   - filled receive buffer
!
! = log ======================================================================================
!
! 13/01/17 - create for v0.4
! TODO: choose method with ini file
! ********************************************************************************************

subroutine fill_receive_buffer( int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! send/receive buffer, integer and real
    integer(kind=ik), intent(in)        :: int_send_buffer(:,:)
    integer(kind=ik), intent(out)       :: int_receive_buffer(:,:)

    real(kind=rk), intent(in)           :: real_send_buffer(:,:)
    real(kind=rk), intent(out)          :: real_receive_buffer(:,:)

    ! communications matrix: neighboring proc rank
    ! com matrix pos: position in send buffer
    integer(kind=ik), intent(in)        :: com_matrix(:,:), com_matrix_pos(:,:)

    ! fill method, dummy variable
    character(len=80)                   :: data_synchronization

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    data_synchronization = "RMA_lock_unlock"

!---------------------------------------------------------------------------------------------
! main body

    select case(data_synchronization)

        ! use RMA with lock/unlock synchronization, use MPI_Get
        case('RMA_lock_unlock')
            !call RMA_lock_unlock_get_data( int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )
            !call RMA_lock_unlock_put_data( int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )
            call isend_irecv_data( int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: synchronization method is unknown"
            write(*,*) data_synchronization
            stop

    end select

end subroutine fill_receive_buffer
