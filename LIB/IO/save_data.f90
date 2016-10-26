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

    integer(kind=ik)                :: dF

    dF = blocks_params%number_data_fields

    if ( dF == 1 ) then
        call write_field(iteration, time, dF)
    else
        ! more than one field
        ! to do

    end if

end subroutine save_data
