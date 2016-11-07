! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: save_data.f90
! version: 0.4
! author: msr
!
! save data main function, call write data routine
!
! input:    - time loop parameter
!           - parameter array
!           - light data array
!           - heavy data array
! output:   -
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine save_data(iteration, time, params, block_list, block_data)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration

    ! user defined parameter structure
    type (type_params), intent(out)                 :: params
    ! light data array
    integer(kind=ik), allocatable, intent(out)      :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), allocatable, intent(out)         :: block_data(:, :, :, :)

!---------------------------------------------------------------------------------------------
! interfaces

    interface
!        subroutine write_field(block_list, number_blocks, max_treelevel)
!            use module_params
!            integer(kind=ik), allocatable, intent(out)  :: block_list(:, :)
!            integer(kind=ik), intent(in)                :: number_blocks
!            integer(kind=ik), intent(in)                :: max_treelevel
!        end subroutine write_field

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    !do dF = 1, blocks_params%number_data_fields
    !    call write_field(iteration, time, dF, fname)
    !end do

end subroutine save_data
