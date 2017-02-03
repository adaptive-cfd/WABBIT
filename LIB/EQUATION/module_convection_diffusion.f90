! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_convection_diffusion.f90
! version: 0.4
! author: msr
!
! module for 2D/3D convection diffusion physics
!
! = log ======================================================================================
!
! 06/12/16 - create
!
! ********************************************************************************************

module module_convection_diffusion

!---------------------------------------------------------------------------------------------
! modules

    use module_precision

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined data structure for time independent variables
    type type_params_convection_diffusion_physics

        ! each data field can use separat velocity and diffusion coefficient (in 2D convection-diffusion RHS)
        ! stored all velocity components in one 1D vector, same with the diffusion coefficients
        real(kind=rk), dimension(:), allocatable    :: u0
        real(kind=rk), dimension(:), allocatable    :: nu

        ! variable names
        character(len=80), allocatable              :: names(:)

    end type type_params_convection_diffusion_physics

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

end module module_convection_diffusion
