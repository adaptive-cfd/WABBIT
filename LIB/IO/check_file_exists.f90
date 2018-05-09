!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name check_file_exists.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief checks if a given file ("fname") exists. if not, code is stopped brutally
!
!>
!! input:    - name of the given file  \\
!!
!! = log ======================================================================================
!! \n
!! 22/09/17 - create
!
! ********************************************************************************************
subroutine check_file_exists(fname)
  implicit none

  character (len=*), intent(in) :: fname
  logical :: exist1

  inquire ( file=fname, exist=exist1 )
  if ( exist1 .eqv. .false.) then
    write (*,'("ERROR! file: ",A," not found")') trim(adjustl(fname))
    call abort( 191919, "File not found...."//trim(adjustl(fname)) )
  endif

end subroutine check_file_exists
