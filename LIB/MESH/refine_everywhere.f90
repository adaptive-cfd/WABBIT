! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: refine_everywhere.f90
! version: 0.4
! author: msr, engels
!
! refine every block to create the wavelet safety zone
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine refine_everywhere( params, block_list, block_data )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)

    ! loop variables
    integer(kind=ik)                    :: k, N

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine respect_min_max_treelevel( block_list, max_level, min_level)
            use module_params
            integer(kind=ik), intent(inout)  :: block_list(:, :)
            integer(kind=ik), intent(in)     :: max_level, min_level
        end subroutine respect_min_max_treelevel

        subroutine refine_mesh( params, block_list, block_data )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
        end subroutine refine_mesh

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    N = size(block_list, 1)

!---------------------------------------------------------------------------------------------
! main body

    ! set status "refine" for all active blocks
    do k = 1, N
        if ( block_list(k, 1) /= -1 ) then

            block_list(k, params%max_treelevel+2 ) = 1

        end if
    end do

    ! check if block has reached maximal level
    call respect_min_max_treelevel( block_list, params%max_treelevel, params%min_treelevel )

    ! interpolate the new mesh
    call refine_mesh( params, block_list, block_data )

end subroutine refine_everywhere
