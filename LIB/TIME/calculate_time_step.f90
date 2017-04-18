!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name calculate_time_step.f90
!> \version 0.5
!> \author sm
!
!> \brief calculate time step
!
!> \details
!! input:    -  params \n
!! output:   -  time step dt
!! \n
!!
!! = log ======================================================================================
!! \n
!! 18/04/17 - create
!
! ********************************************************************************************
subroutine calculate_time_step( params, dx, dt )

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> dx
    real(kind=rk)                       :: dx
    !> time step dt
    real(kind=rk), intent(out)          :: dt

    ! loop variables
    integer(kind=ik)                    :: dF, N_dF

!---------------------------------------------------------------------------------------------
! variables initialization

    N_dF  = params%number_data_fields
!---------------------------------------------------------------------------------------------
! main body

    select case(params%time_step_method)
        case('fixed')
            dt = params%dt
        case('CFL_cond')
            ! calculate time step, loop over all data fields
            if ( params%threeD_case ) then
                do dF = 2, N_dF+1
                    dt = minval((/dt, params%CFL * dx / norm2( params%physics%u0((dF-2)*2 + 1 : (dF-2)*2 + 3 )) /))
                end do
            else
                do dF = 2, N_dF+1
                    dt = minval((/dt, params%CFL * dx / norm2( params%physics%u0((dF-2)*2 + 1 : (dF-2)*2 + 2 )) /))
                end do
            end if
    end select

end subroutine calculate_time_step
