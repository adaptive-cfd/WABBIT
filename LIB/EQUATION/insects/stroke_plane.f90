!-------------------------------------------------------------------------------
! Stroke plane
! Input:
!       time
! Output:
!       eta_stroke: stroke plane angle
subroutine StrokePlane ( time, Insect )
  implicit none

  real(kind=rk), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  real(kind=rk) :: eta_stroke
  character(len=strlen) :: dummy

  select case (Insect%BodyMotion)
  case ("free_flight")
    eta_stroke = Insect%eta0  ! read from file
  case ("tethered")
    eta_stroke = Insect%eta0  ! read from file
  case ("command-line")
    ! for dry runs where the kinematic state is given by the command line call
    call get_command_argument(16,dummy)
    read (dummy,*) eta_stroke

    if(root) write(*,'("eta=",g12.4,"°")') eta_stroke
    eta_stroke = deg2rad(eta_stroke)

  case ("takeoff")
    eta_stroke = Insect%eta_stroke ! read from file
  case default
    eta_stroke = Insect%eta0  ! read from file
  end select

  ! save it in the insect
  Insect%eta_stroke = eta_stroke

end subroutine StrokePlane
