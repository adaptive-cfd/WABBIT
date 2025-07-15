!> \brief module for all operator routines
! *********************************************************************************************

module module_operators

use mpi
use module_params     ! global parameters
use module_timing
use module_treelib
use module_forestMetaData

implicit none

! I usually find it helpful to use the private keyword by itself initially, which specifies
! that everything within the module is private unless explicitly marked public.
PRIVATE

! everything is save by default
SAVE 

! Private finite difference stencil parameters
! First derivative stencils for central finite differences
real(kind=rk), parameter :: FD1_CTW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)  ! Tam & Webb (1993) optimized 7-point stencil for wave propagation
real(kind=rk), parameter :: FD1_CTWR4(-3:3) = (/-0.020843142770_rk, +0.166705904415_rk, -0.770882380518_rk, 0.0_rk, 0.770882380518_rk, -0.166705904415_rk, 0.020843142770_rk/)  ! Tam & Shen (1993) revised optimized 7-point stencil
real(kind=rk), parameter :: FD1_C2(-1:1) = (/ -0.5_rk, 0.0_rk, +0.5_rk /)  ! 2nd order central difference for first derivative
real(kind=rk), parameter :: FD1_C4(-2:2) = (/1.0_rk, -8.0_rk, 0.0_rk, 8.0_rk, -1.0_rk/) / 12.0_rk  ! 4th order central difference for first derivative
real(kind=rk), parameter :: FD1_C6(-3:3) = (/-1.0_rk, 9.0_rk, -45.0_rk, 0.0_rk, 45.0_rk, -9.0_rk, 1.0_rk/) / 60.0_rk  ! 6th order central difference for first derivative
real(kind=rk), parameter :: FD1_C8(-4:4) = (/1.0_rk/280.0_rk, -4.0_rk/105.0_rk, 1.0_rk/5.0_rk, -4.0/5.0_rk, 0.0_rk, 4.0_rk/5.0_rk, -1.0_rk/5.0_rk, 4.0_rk/105.0_rk, -1.0_rk/280.0_rk/)  ! 8th order central difference for first derivative

! Second derivative stencils for central finite differences
real(kind=rk), parameter :: FD2_C2(-1:1) = (/  1.0_rk, -2.0_rk,    1.0_rk /)  ! 2nd order central difference for second derivative
real(kind=rk), parameter :: FD2_C4(-2:2) = (/ -1.0_rk,  16.0_rk,  -30.0_rk,   16.0_rk,   -1.0_rk /) / 12.0_rk  ! 4th order central difference for second derivative
real(kind=rk), parameter :: FD2_C6(-3:3) = (/  2.0_rk, -27.0_rk,   270.0_rk, -490.0_rk,   270.0_rk,  -27.0_rk,    2.0_rk/) / 180.0_rk  ! 6th order central difference for second derivative
real(kind=rk), parameter :: FD2_C8(-4:4) = (/ -9.0_rk,  128.0_rk, -1008.0_rk, 8064.0_rk, -14350.0_rk, 8064.0_rk, -1008.0_rk, 128.0_rk, -9.0_rk /) / 5040_rk  ! 8th order central difference for second derivative

! Cross derivative coefficients (mixed partials d²/dxdy, d²/dxdz, d²/dydz)
! These are applied in a cross-shaped pattern along diagonals
! 2nd order cross derivative: 3-point stencil applied cross-wise
real(kind=rk), parameter :: FD1X_C2(-1:1) = (/ 0.25_rk, 0.0_rk, 0.25_rk /)
! 4th order cross derivative: 5-point stencil applied cross-wise
! From original implementation: 1/3 for ±1 points, -1/48 for ±2 points
real(kind=rk), parameter :: FD1X_C4(-2:2) = (/ -1.0_rk/48.0_rk, 1.0_rk/3.0_rk, 0.0_rk, 1.0_rk/3.0_rk, -1.0_rk/48.0_rk /)

! Forward finite difference stencils for first derivative (used at left boundaries)
! R0 = starting from point 0, using points to the right
real(kind=rk), parameter :: FD1_R0_1(0:1) = (/ -1.0_rk, 1.0_rk /)  ! 1st order forward difference (2-point)
real(kind=rk), parameter :: FD1_R0_2(0:2) = (/ -3.0_rk, 6.0_rk, -1.0_rk /) / 2.0_rk  ! 2nd order forward difference (3-point)
real(kind=rk), parameter :: FD1_R0_3(0:3) = (/ -11.0_rk, 18.0_rk, -9.0_rk, 2.0_rk /) / 6.0_rk  ! 3rd order forward difference (4-point)
real(kind=rk), parameter :: FD1_R0_4(0:4) = (/ -25.0_rk, 48.0_rk, -36.0_rk, 16.0_rk, -3.0_rk /) / 12.0_rk  ! 4th order forward difference (5-point)
real(kind=rk), parameter :: FD1_R0_5(0:5) = (/ -137.0_rk, 300.0_rk, -300.0_rk, 200.0_rk, -75.0_rk, 12.0_rk /) / 60.0_rk  ! 5th order forward difference (6-point)
real(kind=rk), parameter :: FD1_R0_6(0:6) = (/ -147.0_rk, 360.0_rk, -450.0_rk, 400.0_rk, -225.0_rk, 72.0_rk, -10.0_rk /) / 60.0_rk  ! 6th order forward difference (7-point)

! Biased finite difference stencils for first derivative (used near boundaries)
! R1 = starting from point -1, using mostly points to the right
real(kind=rk), parameter :: FD1_R1_3(-1:2) = (/ -2.0_rk, -3.0_rk, 6.0_rk, -1.0_rk /) / 6.0_rk  ! 3rd order biased difference (4-point, stencil from -1 to +2)
real(kind=rk), parameter :: FD1_R1_4(-1:3) = (/ -3.0_rk, -10.0_rk, 18.0_rk, -6.0_rk, 1.0_rk /) / 12.0_rk  ! 4th order biased difference (5-point, stencil from -1 to +3)

!**********************************************************************************************
! These are the important routines that are visible to WABBIT:
!**********************************************************************************************
PUBLIC :: compute_derivative, compute_vorticity, compute_vorticity_abs, compute_divergence, compute_gradient, compute_laplacian, compute_Qcriterion, componentWiseNorm_tree, componentWiseNorm_block, setup_FD1_left_stencil, setup_FD1_right_stencil, setup_FD2_stencil, setup_FD1X_stencil

contains

#include "compute_Qcriterion.f90"
#include "compute_vorticity.f90"
#include "compute_derivative.f90"
#include "compute_divergence.f90"
#include "compute_gradient.f90"
#include "compute_laplacian.f90"
#include "componentWiseNorm_tree.f90"

!**********************************************************************************************
!> \brief setup the left-side FD1 stencils
!> \details This routine sets up the left-side first derivative stencils based on the discretization order
subroutine setup_FD1_left_stencil(FD1_order_discretization_in, FD1_l_array, l_start, l_end)
    implicit none
    
    character(len=*), intent(in) :: FD1_order_discretization_in !> order of the first derivative discretization
    real(kind=rk), allocatable, intent(inout) :: FD1_l_array(:) !> left-side stencil array
    integer(kind=ik), intent(inout) :: l_start, l_end !> left-side stencil bounds
    
    if (allocated(FD1_l_array)) deallocate(FD1_l_array)
    
    select case(FD1_order_discretization_in)
        case("FD_2nd_central")
            l_start = -1
            l_end = 1
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_C2(:)
        case("FD_4th_central")
            l_start = -2
            l_end = 2
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_C4(:)
        case("FD_6th_central")
            l_start = -3
            l_end = 3
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_C6(:)
        case("FD_8th_central")
            l_start = -4
            l_end = 4
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_C8(:)
        case("FD_4th_conv_0_4")
            l_start = -4
            l_end = 0
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_R0_4(-l_start:-l_end:-1) ! left side is reversed
        case("FD_4th_conv_1_3")
            l_start = -3
            l_end = 1
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_R1_4(-l_start:-l_end:-1) ! left side is reversed
        case default
            call abort(250615, "ERROR: order of FD1 discretization not known: " // trim(FD1_order_discretization_in))
    end select
end subroutine setup_FD1_left_stencil

!**********************************************************************************************
!> \brief setup the right-side FD1 stencils
!> \details This routine sets up the right-side first derivative stencils based on the discretization order
subroutine setup_FD1_right_stencil(FD1_order_discretization_in, FD1_r_array, r_start, r_end)
    implicit none
    
    character(len=*), intent(in) :: FD1_order_discretization_in !> order of the first derivative discretization
    real(kind=rk), allocatable, intent(inout) :: FD1_r_array(:) !> right-side stencil array
    integer(kind=ik), intent(inout) :: r_start, r_end !> right-side stencil bounds
    
    if (allocated(FD1_r_array)) deallocate(FD1_r_array)
    
    select case(FD1_order_discretization_in)
        case("FD_2nd_central")
            r_start = -1
            r_end = 1
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_C2(:)
        case("FD_4th_central")
            r_start = -2
            r_end = 2
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_C4(:)
        case("FD_6th_central")
            r_start = -3
            r_end = 3
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_C6(:)
        case("FD_8th_central")
            r_start = -4
            r_end = 4
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_C8(:)
        case("FD_4th_conv_0_4")
            r_start = 0
            r_end = 4
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_R0_4(:)
        case("FD_4th_conv_1_3")
            r_start = -1
            r_end = 3
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_R1_4(:)
        case default
            call abort(250615, "ERROR: order of FD1 discretization not known: " // trim(FD1_order_discretization_in))
    end select
end subroutine setup_FD1_right_stencil

!**********************************************************************************************
!> \brief setup the FD2 stencils
!> \details This routine sets up the second derivative stencils based on the discretization order
subroutine setup_FD2_stencil(FD2_order_discretization_in, FD2_array, fd2_start, fd2_end)
    implicit none
    
    character(len=*), intent(in) :: FD2_order_discretization_in !> order of the second derivative discretization
    real(kind=rk), allocatable, intent(inout) :: FD2_array(:) !> second derivative stencil array
    integer(kind=ik), intent(inout) :: fd2_start, fd2_end !> second derivative stencil bounds
    
    if (allocated(FD2_array)) deallocate(FD2_array)
    
    select case(FD2_order_discretization_in)
        case("FD_2nd_central")
            fd2_start = -1
            fd2_end = 1
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C2(:)
        case("FD_4th_central", "FD_4th_conv_0_4", "FD_4th_conv_1_3")
            fd2_start = -2
            fd2_end = 2
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C4(:)
        case("FD_6th_central", "FD_6th_conv_0_6")
            fd2_start = -3
            fd2_end = 3
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C6(:)
        case("FD_8th_central")
            fd2_start = -4
            fd2_end = 4
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C8(:)
        case default
            call abort(250616, "ERROR: order of FD2 discretization not known: " // trim(FD2_order_discretization_in))
    end select
end subroutine setup_FD2_stencil

!**********************************************************************************************
!> \brief setup the FD1X cross-derivative stencils
!> \details This routine sets up the cross-derivative stencils based on the discretization order
subroutine setup_FD1X_stencil(FD1_order_discretization_in, FD1X_array, fd1x_start, fd1x_end)
    implicit none
    
    character(len=*), intent(in) :: FD1_order_discretization_in !> order of the first derivative discretization (used for cross derivatives)
    real(kind=rk), allocatable, intent(inout) :: FD1X_array(:) !> cross-derivative stencil array
    integer(kind=ik), intent(inout) :: fd1x_start, fd1x_end !> cross-derivative stencil bounds
    
    if (allocated(FD1X_array)) deallocate(FD1X_array)
    
    select case(FD1_order_discretization_in)
        case("FD_2nd_central")
            fd1x_start = -1
            fd1x_end = 1
            allocate(FD1X_array(fd1x_start:fd1x_end))
            FD1X_array(:) = FD1X_C2(:)
        case("FD_4th_central")
            fd1x_start = -2
            fd1x_end = 2
            allocate(FD1X_array(fd1x_start:fd1x_end))
            FD1X_array(:) = FD1X_C4(:)
        case("FD_6th_central", "FD_8th_central", "FD_4th_conv_0_4", "FD_4th_conv_1_3")
            ! For higher orders, default to 4th order cross derivatives for now
            fd1x_start = -2
            fd1x_end = 2
            allocate(FD1X_array(fd1x_start:fd1x_end))
            FD1X_array(:) = FD1X_C4(:)
        case default
            call abort(250617, "ERROR: order of FD1X (cross derivative) discretization not known: " // trim(FD1_order_discretization_in))
    end select
end subroutine setup_FD1X_stencil

end module module_operators
