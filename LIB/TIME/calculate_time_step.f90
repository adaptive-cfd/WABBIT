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
!>
!! input:    -  params \n
!! output:   -  time step dt \n
!!
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

    ! norm of vector u
    real(kind=rk)                       :: norm_u

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
                do dF = 1, N_dF
                    norm_u = norm2( params%physics%u0((dF-1)*2 + 1 : (dF-1)*2 + 3 ))
                    ! check for zero velocity to avoid divison by zero
                    if (norm_u < 1e-12_rk) norm_u = 9e9_rk
                    dt = minval((/dt, params%CFL * dx / norm_u /))
                end do
            else
                do dF = 1, N_dF
                    norm_u = norm2( params%physics%u0((dF-1)*2 + 1 : (dF-1)*2 + 2 ))
                    ! check for zero velocity to avoid divison by zero
                    if (norm_u < 1e-12_rk) norm_u = 9e9_rk
                    dt = minval((/dt, params%CFL * dx / norm_u /))
                end do
            end if
    end select

end subroutine calculate_time_step
