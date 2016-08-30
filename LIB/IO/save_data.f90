! ********************************
! 2D AMR prototype
! --------------------------------
!
! save data main function
!
! name: save_data.f90
! date: 02.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine save_data(iteration, time)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(in) 		:: time
    integer(kind=ik), intent(in) 	:: iteration

    ! save data for all blocks
    call write_field(iteration, time)

end subroutine save_data
