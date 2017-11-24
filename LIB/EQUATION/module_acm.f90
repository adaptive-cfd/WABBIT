!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_acm.f90
!> \version 0.5
!> \author sm
!!
!! \brief module for 2D/3D acm physics
!!
!!
!! = log ======================================================================================
!! \n
!!
! ********************************************************************************************

module module_acm

!---------------------------------------------------------------------------------------------
! modules

    use module_precision

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined data structure for time independent variables
    type type_params_physics_acm
       
        ! c_0 
        real(kind=rk)                               :: c_0
        ! nu
        real(kind=rk)                               :: nu
        ! gamma_p
        real(kind=rk)                               :: gamma_p
        ! want to add forcing?
        logical                                     :: forcing
        
        ! variable names
        character(len=80), allocatable              :: names(:)

    end type type_params_physics_acm

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains


end module module_acm
