!> \file
!> \callgraph
subroutine error_msg(msg)
    implicit none
    character(len=*), intent(in) :: msg

    write(*,*) msg
    write(*,*) "Code stops now."
    stop
  end subroutine
