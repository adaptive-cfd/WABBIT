!> \brief loads the specified *.h5 file and creates a *.key file that contains
!! min / max / mean / L2 norm of the field data. This is used for testing
!! so that we don't need to store entire fields but rather the *.key only
!! \version 10/1/18 - create commit b2719e1aa2339f4f1f83fb29bd2e4e5e81d05a2a
!*********************************************************************************************

subroutine compare_keys(key1, key2)
    use module_helpers, only: check_file_exists
    use module_globals

    implicit none
    character(len=*), intent(in) :: key1, key2
    integer(kind=ik)             :: i
    real(kind=rk)                :: data1(1:7), data2(1:7), error(1:7)
    integer(kind=ik)             :: curves1(1:2), curves2(1:2), error_curve(1:2)
    !-----------------------------------------------------------------------------------------------------
    ! does the user need help?
    if (key1 == '--help' .or. key1 == '--h') then
        write(*,*) "wabbit postprocessing routine to compare keyvalues of two .key files"
        write(*,*) "mpi_command -n number_procs ./wabbit-post --compare-keys old.key new.key"
    else

        call check_file_exists(key1)
        call check_file_exists(key2)

        write(*,*) "Key1: ", trim(adjustl(key1))
        write(*,*) "Key2: ", trim(adjustl(key2))

        open(59, file = key1, status = 'unknown', action='read')
        read(59,'(7(es15.8,1x), 2(i12,1x))') data1, curves1
        close(59)

        open(59, file = key2, status = 'unknown', action='read')
        read(59,'(7(es15.8,1x),2(i12,1x))') data2, curves2
        close(59)

        write (*,'("present  : time=",es15.8," max=",es15.8," min=",es15.8," &
        &sum=",es15.8," sum**2=",es15.8," q=",es15.8, " grid=",es15.8," sfc_hilbert=" ,i12, " sfc_z=",i12)') &
        data1, curves1
        write (*,'("reference: time=",es15.8," max=",es15.8," min=",es15.8," &
        &sum=",es15.8," sum**2=",es15.8," q=",es15.8, " grid=",es15.8," sfc_hilbert=" ,i12, " sfc_z=",i12)') &
        data2, curves2

        ! errors:
        do i = 1, 7
            if (dabs(data2(i))>=1.0e-5_rk) then
                ! relative error
                error(i) = dabs( (data2(i)-data1(i)) / data2(i) )
            else
                ! absolute error
                error(i) = dabs( (data2(i)-data1(i)) )
            end if
        enddo
        do i = 1, 2
            error_curve(i) = abs( (curves2(i)-curves1(i)) )
        end do
        write (*,'("err(rel) : time=",es15.8," max=",es15.8," min=",es15.8," &
        &sum=",es15.8," sum**2=",es15.8," q=",es15.8, " grid=", es15.8, " sfc_hilbert=",i12, " sfc_z=", i12)') error, error_curve

        if ((maxval(error)<1.0e-4)) then !.and. maxval(error_curve)<5_ik) then
            ! all cool
            write (*,*) "okay!"

            ! on some machines, returning an exit code (exit(1)) does not work
            ! so write your exit code in a small txt file as well. this allows unit tests
            ! on turing.
            open (22, file='return', status='replace')
            write(22,'(i1)') 0
            close(22)
            call exit(0)
        else
            ! very bad
            write (*,*) "ERROR"

            ! on some machines, returning an exit code (exit(1)) does not work
            ! so write your exit code in a small txt file as well. this allows unit tests
            ! on turing.
            open (22, file='return', status='replace')
            write(22,'(i1)') 1
            close(22)
            call exit(1)
        endif
    end if

end subroutine
