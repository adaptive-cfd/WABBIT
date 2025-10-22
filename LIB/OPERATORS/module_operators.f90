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

! Composite second derivative stencils (from module_poisson)
real(kind=rk), parameter :: FD2_COMP3_1_2(-3:3) = (/ -2.0_rk, 9.0_rk, 18.0_rk, -50.0_rk, 18.0_rk, 9.0_rk, -2.0_rk /) / 36.0_rk  ! 3th order composite difference (-1,2)
real(kind=rk), parameter :: FD2_COMP4_0_4(-4:4) = (/ -75.0_rk, 544.0_rk, -1776.0_rk, 3552.0_rk, -4490.0_rk, 3552.0_rk, -1776.0_rk, 544.0_rk, -75.0_rk /) / 144.0_rk  ! 4th order composite difference (0,4)
real(kind=rk), parameter :: FD2_COMP4_2_2(-4:4) = (/ 1.0_rk, -16.0_rk, 64.0_rk, 16.0_rk, -130.0_rk, 16.0_rk, 64.0_rk, -16.0_rk, 1.0_rk /) / 144.0_rk  ! 4th order composite difference (2,2)
real(kind=rk), parameter :: FD2_COMP4_1_3(-4:4) = (/ 3.0_rk, -8.0_rk, -24.0_rk, 264.0_rk, -470.0_rk, 264.0_rk, -24.0_rk, -8.0_rk, 3.0_rk /) / 144.0_rk  ! 4th order composite difference (1,3)
real(kind=rk), parameter :: FD2_COMP5_2_3(-5:5) = (/ -6.0_rk, 105.0_rk, -590.0_rk, 1440.0_rk, 1620.0_rk, -5138.0_rk, 1620.0_rk, 1440.0_rk, -590.0_rk, 105.0_rk, -6.0_rk /) / 3600.0_rk  ! 5th order composite difference (2,3)
real(kind=rk), parameter :: FD2_COMP6_3_3(-6:6) = (/ 1.0_rk, -18.0_rk, 171.0_rk, -810.0_rk, 1935.0_rk, 828.0_rk, -4214.0_rk, 828.0_rk, 1935.0_rk, -810.0_rk, 171.0_rk, -18.0_rk, 1.0_rk /) / 3600.0_rk  ! 6th order composite difference (3,3)
real(kind=rk), parameter :: FD2_COMP6_2_4(-6:6) = (/ 2.0_rk, -40.0_rk, 217.0_rk, -520.0_rk, 270.0_rk, 4656.0_rk, -9170.0_rk, 4656.0_rk, 270.0_rk, -520.0_rk, 217.0_rk, -40.0_rk, 2.0_rk /) / 3600.0_rk  ! 6th order composite difference (2,4)
real(kind=rk), parameter :: FD2_COMP6_1_5(-6:6) = (/ 20.0_rk, 4.0_rk, -955.0_rk, 5300.0_rk, -15300.0_rk, 31560.0_rk, -41258.0_rk, 31560.0_rk, -15300.0_rk, 5300.0_rk, -955.0_rk, 4.0_rk, 20.0_rk /) / 3600.0_rk  ! 6th order composite difference (1,5)
real(kind=rk), parameter :: FD2_COMP6_0_6(-6:6) = (/ -1470.0_rk, 14184.0_rk, -63495.0_rk, 176200.0_rk, -342450.0_rk, 501840.0_rk, -569618.0_rk, 501840.0_rk, -342450.0_rk, 176200.0_rk, -63495.0_rk, 14184.0_rk, -1470.0_rk /) / 3600.0_rk  ! 6th order composite difference (0,6)
real(kind=rk), parameter :: FD2_COMP7_3_4(-7:7) = (/ -12.0_rk, 238.0_rk, -2436.0_rk, 13713.0_rk, -45612.0_rk, 83874.0_rk, 84924.0_rk, -269378.0_rk, 84924.0_rk, 83874.0_rk, -45612.0_rk, 13713.0_rk, -2436.0_rk, 238.0_rk, -12.0_rk /) / 176400.0_rk  ! 7th order composite difference (3,4)
real(kind=rk), parameter :: FD2_COMP8_3_5(-8:8) = (/ 15.0_rk, -330.0_rk, 3760.0_rk, -21966.0_rk, 74760.0_rk, -155610.0_rk, 142800.0_rk, 767730.0_rk, -1622318.0_rk, 767730.0_rk, 142800.0_rk, -155610.0_rk, 74760.0_rk, -21966.0_rk, 3760.0_rk, -330.0_rk, 15.0_rk /) / 705600.0_rk  ! 8th order composite difference (3,5)

! Cross derivative coefficients (mixed partials d²/dxdy, d²/dxdz, d²/dydz)
! These are applied in a cross-shaped pattern along diagonals
! 2nd order cross derivative: 3-point stencil applied cross-wise
real(kind=rk), parameter :: FD1X_C2(-1:1) = (/ 0.25_rk, 0.0_rk, 0.25_rk /)
! 4th order cross derivative: 5-point stencil applied cross-wise
! From original implementation: 1/3 for ±1 points, -1/48 for ±2 points
real(kind=rk), parameter :: FD1X_C4(-2:2) = (/ -1.0_rk/48.0_rk, 1.0_rk/3.0_rk, 0.0_rk, 1.0_rk/3.0_rk, -1.0_rk/48.0_rk /)

! Forward finite difference stencils for first derivative (used at left boundaries)
! R0 = starting from point 0, using points to the right
real(kind=rk), parameter :: FD1_COMP1_0_1(0:1) = (/ -1.0_rk, 1.0_rk /)  ! 1st order forward difference (2-point)
real(kind=rk), parameter :: FD1_COMP2_0_2(0:2) = (/ -3.0_rk, 6.0_rk, -1.0_rk /) / 2.0_rk  ! 2nd order forward difference (3-point)
real(kind=rk), parameter :: FD1_COMP3_0_3(0:3) = (/ -11.0_rk, 18.0_rk, -9.0_rk, 2.0_rk /) / 6.0_rk  ! 3rd order forward difference (4-point)
real(kind=rk), parameter :: FD1_COMP4_0_4(0:4) = (/ -25.0_rk, 48.0_rk, -36.0_rk, 16.0_rk, -3.0_rk /) / 12.0_rk  ! 4th order forward difference (5-point)
real(kind=rk), parameter :: FD1_COMP5_0_5(0:5) = (/ -137.0_rk, 300.0_rk, -300.0_rk, 200.0_rk, -75.0_rk, 12.0_rk /) / 60.0_rk  ! 5th order forward difference (6-point)
real(kind=rk), parameter :: FD1_COMP6_0_6(0:6) = (/ -147.0_rk, 360.0_rk, -450.0_rk, 400.0_rk, -225.0_rk, 72.0_rk, -10.0_rk /) / 60.0_rk  ! 6th order forward difference (7-point)

! Biased finite difference stencils for first derivative (used near boundaries)
! R1 = starting from point -1, using mostly points to the right
real(kind=rk), parameter :: FD1_COMP3_1_2(-1:2) = (/ -2.0_rk, -3.0_rk, 6.0_rk, -1.0_rk /) / 6.0_rk  ! 3rd order biased difference (4-point, stencil from -1 to +2)
real(kind=rk), parameter :: FD1_COMP4_1_3(-1:3) = (/ -3.0_rk, -10.0_rk, 18.0_rk, -6.0_rk, 1.0_rk /) / 12.0_rk  ! 4th order biased difference (5-point, stencil from -1 to +3)
real(kind=rk), parameter :: FD1_COMP5_2_3(-2:3) = (/  3.0_rk, -30.0_rk, -20.0_rk, 60.0_rk, -15.0_rk, 2.0_rk /) / 60.0_rk  ! 5th order biased difference (6-point, stencil from -2 to +3)
real(kind=rk), parameter :: FD1_COMP6_1_5(-1:5) = (/ -10.0_rk, -77.0_rk, 150.0_rk, -100.0_rk, 50.0_rk, -15.0_rk, 2.0_rk /) / 60.0_rk  ! 6th order biased difference (7-point, stencil from -1 to +5)
real(kind=rk), parameter :: FD1_COMP6_2_4(-2:4) = (/  2.0_rk, -24.0_rk, -35.0_rk, 80.0_rk, -30.0_rk, 8.0_rk, -1.0_rk /) / 60.0_rk  ! 6th order biased difference (7-point, stencil from -2 to +4)
real(kind=rk), parameter :: FD1_COMP7_3_4(-3:4) = (/ -4.0_rk, 42.0_rk, -252.0_rk, -105.0_rk, 420.0_rk, -126.0_rk, 28.0_rk, -3.0_rk /) / 420.0_rk  ! 7th order biased difference (8-point, stencil from -3 to +4)
real(kind=rk), parameter :: FD1_COMP8_3_5(-3:5) = (/ -5.0_rk, 60.0_rk, -420.0_rk, -378.0_rk, 1050.0_rk, -420.0_rk, 140.0_rk, -30.0_rk, 3.0_rk /) / 840.0_rk  ! 8th order biased difference (9-point, stencil from -3 to +5)

!**********************************************************************************************
! These are the important routines that are visible to WABBIT:
!**********************************************************************************************
PUBLIC :: compute_derivative, compute_vorticity, compute_vorticity_abs, compute_helicity, compute_helicity_abs, compute_divergence, compute_gradient, compute_laplacian, compute_Qcriterion, compute_dissipation, componentWiseNorm_tree, componentWiseNorm_block, setup_FD1_left_stencil, setup_FD1_right_stencil, setup_FD2_stencil, setup_FD1X_stencil

contains

#include "compute_Qcriterion.f90"
#include "compute_vorticity.f90"
#include "compute_dissipation.f90"
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
        case("FD_4th_central", "FD_4th_comp_2_2")
            l_start = -2
            l_end = 2
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_C4(:)
        case("FD_6th_central", "FD_6th_comp_3_3")
            l_start = -3
            l_end = 3
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_C6(:)
        case("FD_8th_central")
            l_start = -4
            l_end = 4
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_C8(:)
        case("FD_3th_comp_1_2")
            l_start = -2
            l_end = 1
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = FD1_COMP3_1_2(-l_start:-l_end:-1) ! left side is reversed
        case("FD_4th_comp_0_4")
            l_start = -4
            l_end = 0
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_COMP4_0_4(-l_start:-l_end:-1) ! left side is reversed
        case("FD_4th_comp_1_3")
            l_start = -3
            l_end = 1
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_COMP4_1_3(-l_start:-l_end:-1) ! left side is reversed
        case("FD_5th_comp_2_3")
            l_start = -3
            l_end = 2
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_COMP5_2_3(-l_start:-l_end:-1) ! left side is reversed
        case("FD_6th_comp_1_5")
            l_start = -5
            l_end = 1
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_COMP6_1_5(-l_start:-l_end:-1) ! left side is reversed
        case("FD_6th_comp_2_4")
            l_start = -4
            l_end = 2
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_COMP6_2_4(-l_start:-l_end:-1) ! left side is reversed
        case("FD_7th_comp_3_4")
            l_start = -4
            l_end = 3
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_COMP7_3_4(-l_start:-l_end:-1) ! left side is reversed
        case("FD_8th_comp_3_5")
            l_start = -5
            l_end = 3
            allocate(FD1_l_array(l_start:l_end))
            FD1_l_array(:) = -FD1_COMP8_3_5(-l_start:-l_end:-1) ! left side is reversed
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
        case("FD_4th_central", "FD_4th_comp_2_2")
            r_start = -2
            r_end = 2
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_C4(:)
        case("FD_6th_central", "FD_6th_comp_3_3")
            r_start = -3
            r_end = 3
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_C6(:)
        case("FD_8th_central")
            r_start = -4
            r_end = 4
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_C8(:)
        case("FD_3th_comp_1_2")
            r_start = -1
            r_end = 2
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP3_1_2(:)
        case("FD_4th_comp_0_4")
            r_start = 0
            r_end = 4
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP4_0_4(:)
        case("FD_4th_comp_1_3")
            r_start = -1
            r_end = 3
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP4_1_3(:)
        case("FD_5th_comp_2_3")
            r_start = -2
            r_end = 3
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP5_2_3(:)
        case("FD_6th_comp_2_4")
            r_start = -2
            r_end = 4
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP6_2_4(:)
        case("FD_6th_comp_1_5")
            r_start = -1
            r_end = 5
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP6_1_5(:)
        case("FD_7th_comp_3_4")
            r_start = -3
            r_end = 4
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP7_3_4(:)
        case("FD_8th_comp_3_5")
            r_start = -3
            r_end = 5
            allocate(FD1_r_array(r_start:r_end))
            FD1_r_array(:) = FD1_COMP8_3_5(:)
        case default
            call abort(250615, "ERROR: order of FD1 discretization not known: " // trim(FD1_order_discretization_in))
    end select
end subroutine setup_FD1_right_stencil

!**********************************************************************************************
!> \brief setup the FD2 stencils
!> \details This routine sets up the second derivative stencils based on the discretization order
subroutine setup_FD2_stencil(FD2_order_discretization_in, FD2_array, fd2_start, fd2_end, use_exact_discretization)
    implicit none
    
    character(len=*), intent(in) :: FD2_order_discretization_in !> order of the second derivative discretization
    real(kind=rk), allocatable, intent(inout) :: FD2_array(:) !> second derivative stencil array
    integer(kind=ik), intent(inout) :: fd2_start, fd2_end !> second derivative stencil bounds
    logical, intent(in), optional :: use_exact_discretization !> if true, use exact composite stencils; if false (default), use central stencils for composite orders
    
    logical :: use_exact
    
    ! We need the second order stencils for two reasons:
    ! 1) to compute the diffusion term
    ! 2) to solve the pressure-poisson equation
    ! For the first case, operator compatibility is not too important and we choose the best available stencil (central one)
    ! For the second case, operator compatibility is crucial and we need the correct stencils
    ! default behaviour is dictated by RHS, so for now we need to set use_exact_discretization = .true. when wanting to choose the coinciding stencil
    use_exact = .false.
    if (present(use_exact_discretization)) use_exact = use_exact_discretization
    
    if (allocated(FD2_array)) deallocate(FD2_array)
    
    select case(FD2_order_discretization_in)
        case("FD_2nd_central")
            fd2_start = -1
            fd2_end = 1
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C2(:)
        case("FD_3rd_comp_1_2")
            if (use_exact) then
                fd2_start = -3
                fd2_end = 3
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP3_1_2(:)
            ! no corresponding 3rd order central stencil, fallback to 4th order
            else
                fd2_start = -2
                fd2_end = 2
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C4(:)
            endif
        case("FD_4th_central")
            fd2_start = -2
            fd2_end = 2
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C4(:)
        case("FD_4th_comp_0_4")
            if (use_exact) then
                fd2_start = -4
                fd2_end = 4
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP4_0_4(:)
            else
                ! fall back to central 4th order stencil
                fd2_start = -2
                fd2_end = 2
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C4(:)
            endif
        case("FD_4th_comp_2_2")
            if (use_exact) then
                fd2_start = -4
                fd2_end = 4
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP4_2_2(:)
            else
                ! fall back to central 4th order stencil
                fd2_start = -2
                fd2_end = 2
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C4(:)
            endif
        case("FD_4th_comp_1_3")
            if (use_exact) then
                fd2_start = -4
                fd2_end = 4
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP4_1_3(:)
            else
                ! fall back to central 4th order stencil
                fd2_start = -2
                fd2_end = 2
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C4(:)
            endif
        case("FD_5th_comp_2_3")
            if (use_exact) then
                fd2_start = -5
                fd2_end = 5
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP5_2_3(:)
            ! no corresponding 5th order central stencil, fallback to 6th order
            else
                fd2_start = -3
                fd2_end = 3
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C6(:)
            endif
        case("FD_6th_central")
            fd2_start = -3
            fd2_end = 3
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C6(:)
        case("FD_6th_comp_3_3")
            if (use_exact) then
                fd2_start = -6
                fd2_end = 6
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP6_3_3(:)
            else
                ! fall back to central 6th order stencil
                fd2_start = -3
                fd2_end = 3
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C6(:)
            endif
        case("FD_6th_comp_2_4")
            if (use_exact) then
                fd2_start = -6
                fd2_end = 6
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP6_2_4(:)
            else
                ! fall back to central 6th order stencil
                fd2_start = -3
                fd2_end = 3
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C6(:)
            endif
        case("FD_6th_comp_1_5")
            if (use_exact) then
                fd2_start = -6
                fd2_end = 6
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP6_1_5(:)
            else
                ! fall back to central 6th order stencil
                fd2_start = -3
                fd2_end = 3
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C6(:)
            endif
        case("FD_6th_comp_0_6")
            if (use_exact) then
                fd2_start = -6
                fd2_end = 6
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP6_0_6(:)
            else
                ! fall back to central 6th order stencil
                fd2_start = -3
                fd2_end = 3
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C6(:)
            endif
        case("FD_7th_comp_3_4")
            if (use_exact) then
                fd2_start = -7
                fd2_end = 7
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP7_3_4(:)
            ! no corresponding 7th order central stencil, fallback to 8th order
            else
                fd2_start = -2
                fd2_end = 2
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C8(:)
            endif
        case("FD_8th_central")
            fd2_start = -4
            fd2_end = 4
            allocate(FD2_array(fd2_start:fd2_end))
            FD2_array(:) = FD2_C8(:)
        case("FD_8th_comp_3_5")
            if (use_exact) then
                fd2_start = -8
                fd2_end = 8
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_COMP8_3_5(:)
            else
                fd2_start = -2
                fd2_end = 2
                allocate(FD2_array(fd2_start:fd2_end))
                FD2_array(:) = FD2_C8(:)
            endif
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
        case("FD_6th_central", "FD_8th_central", "FD_4th_comp_0_4", "FD_4th_comp_1_3", "FD_4th_comp_2_2", &
             "FD_6th_comp_0_6", "FD_6th_comp_1_5", "FD_6th_comp_2_4")
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
