!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name fill_receive_buffer.f90
!> \version 0.4
!> \author msr
!
!> \brief fill receive buffer
!
!>
!! input:    -      \n
!! output:   - filled receive buffer \n
!!
!!
!! = log ======================================================================================
!! \n
!! 13/01/17 - create for v0.4
!
! ********************************************************************************************

subroutine fill_receive_buffer( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_pos )

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
    !! com matrix pos: position in send buffer
    integer(kind=ik), intent(in)        :: com_matrix(:,:)
    integer(kind=ik), intent(inout)     :: com_pos(:)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    select case(params%mpi_data_exchange)

        ! use RMA with lock/unlock synchronization, use MPI_Get
        case('RMA_lock_unlock_get')
            !call RMA_lock_unlock_get_data( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_pos )

        ! use RMA with lock/unlock synchronization, use MPI_Put
        case('RMA_lock_unlock_put')
            !call RMA_lock_unlock_put_data( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_pos )

        ! use non-blocking isend/irecv
        case('Non_blocking_Isend_Irecv')
            call isend_irecv_data( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_pos )

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: data exchange method is unknown"
            write(*,*) params%mpi_data_exchange
            call abort(99912222,"ERROR: data exchange method is unknown")

    end select

end subroutine fill_receive_buffer
