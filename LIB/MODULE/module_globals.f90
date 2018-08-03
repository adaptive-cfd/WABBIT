!-----------------------------------------------------------------
!> \file
!> \brief
!! Module of global parameters and functions
!> \details
!!    * index positions of lgt_data (idx_tree_number,idx_mesh_lvl,idx_refine_sts etc)
!!    * functions abort
!!    * global prints (todo)
!!    * global MPI communicator
!!
!> \version 2.8.2018
!> \author P.Krah
!-----------------------------------------------------------------

!> \brief this module contains all global parameters and functions
module module_globals

! import MPI module
 use module_precision

  implicit none
  !> global communicator for WABBIT! is capsuled with private
  !> use setter and getter functions to change it
  integer(kind=ik)       :: WABBIT_COMM


  !subroutines of this module
  interface abort
      module procedure abort1, abort2, abort3
  end interface

contains


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

end module module_globals
