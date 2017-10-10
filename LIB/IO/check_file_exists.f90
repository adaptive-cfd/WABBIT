! checks if a given file ("fname") exists. if not, code is stopped brutally
subroutine check_file_exists(fname)
  implicit none

  character (len=*), intent(in) :: fname
  logical :: exist1

  inquire ( file=fname, exist=exist1 )
  if ( exist1 .eqv. .false.) then
    write (*,'("ERROR! file: ",A," not found")') trim(adjustl(fname))
    stop
  endif

end subroutine check_file_exists
