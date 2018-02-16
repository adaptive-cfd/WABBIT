!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_navier_stokes.f90
!> \version 0.4
!> \author msr
!!
!! \brief module for 2D/3D navier stokes physics
!!
!!
!! = log ======================================================================================
!! \n
!! 06/12/16 - create
! ********************************************************************************************

module module_navier_stokes

!---------------------------------------------------------------------------------------------
! modules

    use module_precision

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined data structure for time independent variables
    type type_params_physics_navier_stokes

        ! adiabatic coefficient
        real(kind=rk)                               :: gamma_

        ! specific gas constant
        real(kind=rk)                               :: Rs

        ! isochoric heat capacity
        real(kind=rk)                               :: Cv

        ! isobaric heat capacity
        real(kind=rk)                               :: Cp

        ! prandtl number
        real(kind=rk)                               :: Pr

        ! dynamic viscosity
        real(kind=rk)                               :: mu0

        ! dissipation switch
        logical                                     :: dissipation

        ! variable names
        character(len=80), allocatable              :: names(:)

    end type type_params_physics_navier_stokes

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

end module module_navier_stokes
