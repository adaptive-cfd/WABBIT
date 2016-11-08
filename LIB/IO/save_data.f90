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

subroutine save_data(iteration, time, params, block_list, block_data, neighbor_list)

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
    type (type_params), intent(in)                  :: params
    ! light data array
    integer(kind=ik), intent(in)                    :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(in)                       :: block_data(:, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)                    :: neighbor_list(:)

    ! loop variable
    integer(kind=ik)                                :: k

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine write_field(time, iteration, dF, params, block_list, block_data, neighbor_list)
            use module_params
            real(kind=rk), intent(in)                   :: time
            integer(kind=ik), intent(in)                :: iteration
            integer(kind=ik), intent(in)                :: dF
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(in)                :: block_list(:, :)
            real(kind=rk), intent(in)                   :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
        end subroutine write_field

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! real datafields start at datafield 2
    do k = 2, params%number_data_fields+1
        call write_field(time, iteration, k, params, block_list, block_data, neighbor_list)
    end do

end subroutine save_data
