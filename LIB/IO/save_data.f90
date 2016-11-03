! ********************************
! WABBIT
! --------------------------------
!
! save data main function
!
! name: save_data.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine save_data(iteration, time)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(in) 		:: time
    integer(kind=ik), intent(in) 	:: iteration

    character(len=80)               :: fname
    integer(kind=ik)                :: dF

    do dF = 1, blocks_params%number_data_fields
        call write_field(iteration, time, dF, fname)
    end do

end subroutine save_data
