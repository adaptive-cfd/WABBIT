!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name startup_conditioner.f90
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
