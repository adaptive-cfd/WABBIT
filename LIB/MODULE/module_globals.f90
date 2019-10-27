!-----------------------------------------------------------------
!> \file
!> \brief
!! Module of global parameters and functions
!> \details
!!    * index positions of lgt_data (idx_tree_number,IDX_MESH_LVL,IDX_REFINE_STS etc)
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
  integer(kind=ik) :: WABBIT_COMM
  ! index of light data quantities lgt_data(block_index,max_treeleve+idx_<quantity>)
  ! we define the global indices idx_<quantity> here, because it enables us to choose the order of indexing
  ! if necessary
  integer, parameter, public  :: IDX_MESH_LVL        = 1 ! current mesh level of the block
  integer, parameter, public  :: IDX_REFINE_STS      = 2 ! refinement status
  integer, parameter, public  :: IDX_TREE_ID         = 3 ! number of the tree in the forest
  integer, parameter, public  :: EXTRA_LGT_FIELDS    = 3 ! number of data fields additionaly to treecode
  integer, parameter, public  :: tree_ID_flow = 1, tree_ID_mask = 2, tree_ID_mask_coarser = 3
  ! this parameter is a hack. in most parts of the code, a block has n_eqn component entries.
  ! universality dictates that we can also use a different number of components, for example
  ! when syn'ing the mask function (which in many cases has six entries.)
  integer, public  :: N_MAX_COMPONENTS = 6

  !subroutines of this module
  interface abort
      module procedure abort1, abort2, abort3
  end interface

contains

! helper routine to print when entering a routine. DEBUG only!
! all MPIRANKS do print
subroutine debug_header_barrier(msg)
    use mpi
    implicit none
    character(len=*), intent(in) :: msg
    integer :: mpirank, mpicode

    ! fetch my process id
    call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

    call MPI_BARRIER(WABBIT_COMM, mpicode)
    write(*,'(A,i5,1x,A)') "DEBUG_BARRIER Entry, rank=", mpirank, msg
    call MPI_BARRIER(WABBIT_COMM, mpicode)

end subroutine

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
