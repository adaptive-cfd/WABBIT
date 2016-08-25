! ********************************
! 2D AMR prototype
! --------------------------------
! 
! check and make work directory for
! data writing
!
! name: check_workdir.f90
! date: 03.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine check_workdir()

    use module_params
    use module_blocks

    implicit none

    logical :: dir_e

    ! work dir check
    inquire(file=trim(params%name_workdir), exist=dir_e)
    if (.NOT.(dir_e)) call system('mkdir ' // trim(params%name_workdir)  )

    ! case dir check
    inquire(file=trim(params%name_workdir)//trim(params%name_case), exist=dir_e)
    if (.NOT.(dir_e)) call system('mkdir ' // trim(params%name_workdir) // trim(params%name_case)  )

end subroutine check_workdir
