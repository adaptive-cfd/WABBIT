!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name compute_forcing.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief compute forcing term for artificial conmpressibility method to accelerate fluid
!
!>
!! input:    - velocity, volume integral, grid parameters \n
!! output:   - forcing term \n
!!
!!
!! = log ======================================================================================
!! \n
!! 21/07/17 - create \n
!! 19/09/17 - add 3D
!*********************************************************************************************

subroutine compute_forcing(forcing, volume_int, Lx, Ly, Lz, time)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

implicit none

    !> forcing term
    real(kind=rk), dimension(3), intent(out) :: forcing

    !> volume integral
    real(kind=rk), dimension(3), intent(in)  :: volume_int
    !> domain size
    real(kind=rk), intent(in)                :: Lx, Ly, Lz
    !> time
    real(kind=rk), intent(in)                :: time

    !> mean flow
    real(kind=rk)                            :: ux_mean, uy_mean, uz_mean
    !> startup conditioner
    real(kind=rk)                            :: startup_conditioner

!---------------------------------------------------------------------------------------------
! variables initialization
    forcing = 0.0_rk
!---------------------------------------------------------------------------------------------
! main body

    ux_mean = volume_int(1)/(Lx*Ly*Lz)
    uy_mean = volume_int(2)/(Lx*Ly*Lz)
    uz_mean = volume_int(3)/(Lx*Ly*Lz)

    forcing(1) = max(0.0_rk, 1.0_rk-ux_mean)* startup_conditioner(time, 0.0_rk, 0.5_rk)
    forcing(2) = max(0.0_rk, 0.0_rk-uy_mean)* startup_conditioner(time, 0.0_rk, 0.5_rk)
    forcing(3) = max(0.0_rk, 0.0_rk-uz_mean)* startup_conditioner(time, 0.0_rk, 0.5_rk)

end subroutine compute_forcing
