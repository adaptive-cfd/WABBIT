! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_rhs_navier_stokes.f90
! version: 0.4
! author: msr
!
! module for 2D navier stokes RHS
!
! = log ======================================================================================
!
! 06/12/16 - create
! ********************************************************************************************

module module_rhs_navier_stokes

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined data structure for time independent variables
    type type_params_rhs

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

    end type type_params_rhs

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

! data init
subroutine init_prams_rhs( params_rhs )

    implicit none

    ! user defined data structure for time independent variables
    type (type_params_rhs), intent(out)              :: params_rhs

    ! set params
    params_rhs%gamma_       = 1.4_rk
    params_rhs%Rs           = 8.3144621_rk / 0.0270_rk
    params_rhs%Cv           = params_rhs%Rs/(params_rhs%gamma_-1.0_rk)
    params_rhs%Cp           = params_rhs%cv*params_rhs%gamma_
    params_rhs%Pr           = 0.71_rk
    params_rhs%mu0          = 1e-2_rk
    params_rhs%dissipation  = .true.

end subroutine init_prams_rhs

! RHS
subroutine RHS_2D_navier_stokes( params_rhs, phi, g, Bs )

implicit none

    ! user defined data structure for time independent variables
    type (type_params_rhs), intent(out)                         :: params_rhs

    ! grid parameter
    integer(kind=ik), intent(in)                                :: g, Bs

    ! datafield
    real(kind=rk), dimension(Bs+2*g, Bs+2*g), intent(inout)     :: phi

    ! parameter
    real(kind=rk)                                               :: gamma_, gamma_1, Cv, Cp, Rs, Pr, mu0

    !---------------------------------------------------------------------------------

    ! set parameter, readability
    gamma_  = params_rhs%gamma_
    gamma_1 = gamma_ - 1.0_rk
    Rs      = params_rhs%Rs
    Cv      = params_rhs%Cv
    Cp      = params_rhs%Cp
    Pr      = params_rhs%Pr
    mu0     = params_rhs%mu0

end subroutine RHS_2D_navier_stokes

end module module_rhs_navier_stokes
