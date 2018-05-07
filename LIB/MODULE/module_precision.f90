!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_precision.f90
!> \version 0.4
!> \author msr
!
!> \brief module for data precision
!
!>
!! = log ======================================================================================
!! \n
!! 06/12/16 - create
! ********************************************************************************************

module module_precision

!---------------------------------------------------------------------------------------------
! modules
  use MPI
!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! define data precision parameters
    integer, parameter, public      :: sngl_prec=selected_real_kind(4)
    integer, parameter, public      :: dble_prec=selected_real_kind(8)

    integer, parameter, public      :: int_prec=selected_int_kind(8)

    integer, parameter, public      :: maxdigits = 16
    integer, parameter, public      :: tsize = selected_int_kind(maxdigits)

    integer, parameter, public      :: rk=dble_prec
    integer, parameter, public      :: ik=int_prec

    real(kind=rk),parameter, public :: pi  = 4 * atan(1.0_rk)

!---------------------------------------------------------------------------------------------
! variables initialization
   integer(kind=ik)                            :: WABBIT_COMM

    interface abort
      module procedure abort1, abort2, abort3
    end interface
!---------------------------------------------------------------------------------------------
! main body

contains
  
   !> \brief initialize global communicator of WABBIT MPI_World
   subroutine set_mpi_comm_global(comm)
      implicit none
        integer, intent(in) :: comm

      WABBIT_COMM=comm
  end subroutine set_mpi_comm_global


  subroutine abort1(code,msg)
    implicit none
    character(len=*), intent(in) :: msg
    integer(kind=ik), intent(in) :: code
    integer(kind=ik) :: mpierr

    write(*,*) msg
    call MPI_ABORT( WABBIT_COMM, code, mpierr)
  end subroutine

  subroutine abort2(code)
    implicit none
    integer(kind=ik), intent(in) :: code
    integer(kind=ik) :: mpierr

    call MPI_ABORT( WABBIT_COMM, code, mpierr)
  end subroutine

  subroutine abort3(msg)
    implicit none
    character(len=*), intent(in) :: msg
    integer(kind=ik) :: mpierr

    write(*,*) msg
    call MPI_ABORT( WABBIT_COMM, 666, mpierr)
  end subroutine

end module module_precision
