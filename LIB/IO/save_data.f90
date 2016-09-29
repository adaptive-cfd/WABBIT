! ********************************
! WABBIT
! --------------------------------
!
! save data main function
!
! name: save_data.f90
! date: 29.09.2016
! author: msr
! version: 0.2
!
! ********************************

subroutine save_data(iteration, time, error)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(in) 		:: time, error
    integer(kind=ik), intent(in) 	:: iteration

    integer(kind=ik)                :: dF

    dF = blocks_params%number_data_fields

    if ( dF == 1 ) then
        ! exactly one field
        ! save data for all blocks
        call write_field(iteration, time, error, dF)

    else
        ! more than one field
        ! to do

    end if

end subroutine save_data
