! overwrite and initialize file
subroutine init_empty_file( fname )
  implicit none
  character (len=*), intent(in) :: fname

  open (15, file=fname, status='replace')
  close(15)

end subroutine
