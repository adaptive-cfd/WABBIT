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
  !> global communicator for WABBIT! Dont use MPI_COMM_WORLD!!!!!
  integer(kind=ik)       :: WABBIT_COMM
  ! index of light data quantities lgt_data(block_index,max_treeleve+idx_<quantity>)
  ! we define the global indices idx_<quantity> here, because it enables us to choose the order of indexing
  ! if necessary
  integer, parameter, public  :: idx_refine_sts      = 2 ! refinement status
  integer, parameter, public  :: idx_tree_nr         = 3 ! number of the tree in the forest
  integer, parameter, public  :: idx_mesh_lvl        = 1 ! current mesh level of the block
  integer, parameter, public  :: extra_lgt_fields    = 3 ! number of data fields additionaly to treecode
  !subroutines of this module
  interface abort
      module procedure abort1, abort2, abort3
  end interface

contains

!> use the abort function instead of stop!
!> this is necessary since stop, does not kill
!> the other processes
subroutine abort1(code,msg)
  implicit none
  character(len=*), intent(in) :: msg
  integer(kind=ik), intent(in) :: code
  integer(kind=ik) :: mpierr

  ! you may be tempted to use (if mprank=0) to have nicer output: do not do that.
  ! if root does not come to the error, then you will not see it. better all ranks
  ! display it.
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

  ! you may be tempted to use (if mprank=0) to have nicer output: do not do that.
  ! if root does not come to the error, then you will not see it. better all ranks
  ! display it.
  write(*,*) msg
  call MPI_ABORT( WABBIT_COMM, 666, mpierr)
end subroutine

end module module_globals
