!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name array_compare.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief  soft startup funtion, is zero until time=time_release, then gently goes to one during the time period tau
!
!> \details
!! input:    
!!           - time t
!!           - time release
!!
!! output:   
!!           - value of startup conditioner
!!
!! = log ======================================================================================
!! \n
!! 24/07/17 - create
! ********************************************************************************************

real(kind=rk) function startup_conditioner(time, time_release, tau)

!---------------------------------------------------------------------------------------------
! modules
    use module_precision
!---------------------------------------------------------------------------------------------
! variables

    implicit none

    real(kind=rk), intent(in)  :: time,time_release, tau
    real(kind=rk)              :: dt
!---------------------------------------------------------------------------------------------
! main body

    dt = time-time_release

    if (time <= time_release) then
        startup_conditioner = 0.0_rk
    elseif ( ( time >time_release ).and.(time<(time_release + tau)) ) then
         startup_conditioner =  (dt**3)/(-0.5_rk*tau**3) + 3.0_rk*(dt**2)/tau**2
    else
         startup_conditioner = 1.0_rk
    endif

    return
end function startup_conditioner
