module module_insects_integration_flusi_wabbit
  ! use vars, only : xl, yl, zl, abort, &
  ! periodize_coordinate, cross, deg2rad, pi, rad2deg, root, norm2, &
  ! x0, y0, z0, startup_conditioner, rand_nbr, &
  ! itdrag, Integrals, nu


  use module_precision

  ! interp2_nonper: we need this to interpolate wing thickness and corrugation
  use module_helpers
  ! we need this to read from ini files (e.g. the wing kinematics or shape are read this way)
  use module_ini_files_parser_mpi


implicit none

!-----------------------------------------------------------------------------
  ! The derived integral quantities for fluid-structure interactions.
  type Integrals
     real(kind=rk) :: time = 0.d0
     real(kind=rk) :: EKin = 0.d0
     real(kind=rk) :: Dissip = 0.d0
     real(kind=rk) :: Divergence = 0.d0
     real(kind=rk) :: Volume = 0.d0
     real(kind=rk) :: APow = 0.d0
     real(kind=rk) :: IPow = 0.d0
     real(kind=rk) :: penalization_power = 0.d0
     real(kind=rk) :: penalization_power_x = 0.d0
     real(kind=rk) :: penalization_power_y = 0.d0
     real(kind=rk) :: penalization_power_z = 0.d0
     real(kind=rk),dimension(1:3) :: Force = 0.d0
     real(kind=rk),dimension(1:3) :: Force_unst = 0.d0
     real(kind=rk),dimension(1:3) :: Torque = 0.d0
     real(kind=rk),dimension(1:3) :: Torque_unst = 0.d0
  end type Integrals


    real(kind=rk) :: x0, y0, z0

    logical, parameter :: grid_time_dependent = .true.

    integer, parameter :: pr = rk
end module module_insects_integration_flusi_wabbit
