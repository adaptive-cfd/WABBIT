! ********************************
! 2D AMR prototype
! --------------------------------
!
! check and make directory for given
! time step
!
! name: check_timedir.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine check_timedir(iteration, time)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(in)       :: time
    integer(kind=ik), intent(in)    :: iteration

    logical                         :: dir_e
    character(len=128)              :: name_timedir, name_iteration

    ! create directory name
    write(name_timedir, '(f6.2)') time
    write(name_iteration, '(i5.5)') iteration

    ! time dir check and writing
    inquire(file=trim(params%name_workdir) //trim(params%name_case) //'/' // trim(adjustl(trim(name_iteration))) //'_' // trim(adjustl(trim(name_timedir))), exist=dir_e)
    if (.NOT.(dir_e)) call system('mkdir ' // trim(params%name_workdir) // trim(params%name_case) //'/' // trim(adjustl(trim(name_iteration))) //'_' // trim(adjustl(trim(name_timedir))) )

end subroutine check_timedir
