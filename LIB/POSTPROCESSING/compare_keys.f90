!> \file
! WABBIT
!> \name keyvalues.f90
!> \version 0.5
!> \author sm, engels 
!
!> \brief loads the specified *.h5 file and creates a *.key file that contains
!! min / max / mean / L2 norm of the field data. This is used for testing
!! so that we don't need to store entire fields but rather the *.key only
!! \version 10/1/18 - create commit b2719e1aa2339f4f1f83fb29bd2e4e5e81d05a2a
!*********************************************************************************************

subroutine compare_keys(help, key1, key2)
    use module_IO, only: check_file_exists
    use module_precision

    implicit none
    !> name of the file
    character(len=*), intent(in)            :: key1, key2
    !> help flag
    logical, intent(in)                     :: help

    real(kind=rk)          :: t1,t2,a1,a2,b1,b2,c1,c2,d1,d2,q1,q2, e0,e1,e2,e3,e4,e5
!-----------------------------------------------------------------------------------------------------

    if (help) then
         write(*,*) "wabbit postprocessing routine to compare keyvalues of two .key files"
         write(*,*) "./wabbit-post 2[3]D --compare-keys old.key new.key"
    else

        call check_file_exists(key1)
        call check_file_exists(key2)

        write (*,'("comparing files ",a20, a20, " for keyvalues")'), trim(adjustl(key1)), trim(adjustl(key2))

        open(59, file = key1, status = 'unknown', action='read')
        read(59,'(6(es15.8,1x))') t1,a1,b1,c1,d1,q1
        close(59)

        open(59, file = key2, status = 'unknown', action='read')
        read(59,'(6(es15.8,1x))') t2,a2,b2,c2,d2,q2
        close(59)

      write (*,'("present  : time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
      t1,a1,b1,c1,d1,q1

      write (*,'("reference: time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
      t2,a2,b2,c2,d2,q2

      ! errors:
      if (dabs(t2)>=1.0e-7_rk) then
        e0 = dabs( (t2-t1) / t2 )
      else
        e0 = dabs( (t2-t1) )
      endif

      if (dabs(a2)>=1.0e-7_rk) then
        e1 = dabs( (a2-a1) / a2 )
      else
        e1 = dabs( (a2-a1) )
      endif

      if (dabs(b2)>=1.0e-7_rk) then
        e2 = dabs( (b2-b1) / b2 )
      else
        e2 = dabs( (b2-b1) )
      endif

      if (dabs(c2)>=1.0e-7_rk) then
        e3 = dabs( (c2-c1) / c2 )
      else
        e3 = dabs( (c2-c1) )
      endif

      if (dabs(d2)>=1.0e-7_rk) then
        e4 = dabs( (d2-d1) / d2 )
      else
        e4 = dabs( (d2-d1) )
      endif

      if (dabs(q2)>=1.0e-7_rk) then
        e5 = dabs( (q2-q1) / q2 )
      else
        e5 = dabs( (q2-q1) )
      endif

      write (*,'("err(rel) : time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
      e0,e1,e2,e3,e4,e5

      if ((e1<1.0e-4_rk) .and. (e2<1.0e-4_rk) .and. (e3<1.0e-4_rk) .and. (e4<1.0e-4_rk) .and. (e0<1.0e-4_rk) .and. (e5<1.0e-4_rk)) then
        ! all cool
        write (*,*) "OKAY..."

        ! on some machines, returning an exit code (exit(1)) does not work
        ! so write your exit code in a small txt file as well. this allows unit tests
        ! on turing.
         open (15, file='return', status='replace')
         write(15,'(i1)') 0
         close(15)
        call exit(0)
      else
        ! very bad
        write (*,*) "ERROR"

        ! on some machines, returning an exit code (exit(1)) does not work
        ! so write your exit code in a small txt file as well. this allows unit tests
        ! on turing.
          open (15, file='return', status='replace')
          write(15,'(i1)') 1
          close(15)
        call exit(1)
      endif
  end if

end subroutine
