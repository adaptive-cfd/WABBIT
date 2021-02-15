!-----------------------
! Rotation Matrices
!-----------------------
subroutine Rx (R,angle)
  implicit none
  real(kind=rk), intent (in) :: angle
  real(kind=rk),dimension(1:3,1:3), intent(out) :: R

  R(1,:) = (/ 1.d0, 0.d0, 0.d0/)
  R(2,:) = (/ 0.d0, cos(angle), sin(angle) /)
  R(3,:) = (/ 0.d0, -sin(angle), cos(angle) /)
end subroutine


subroutine Ry (R,angle)
  implicit none
  real(kind=rk), intent (in) :: angle
  real(kind=rk),dimension(1:3,1:3), intent(out) :: R

  R(1,:) = (/ cos(angle), 0.d0, -sin(angle)/)
  R(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
  R(3,:) = (/ +sin(angle), 0.d0, cos(angle) /)
end subroutine


subroutine Rz (R,angle)
  implicit none
  real(kind=rk), intent (in) :: angle
  real(kind=rk),dimension(1:3,1:3), intent(out) :: R

  R(1,:) = (/ cos(angle), +sin(angle), 0.d0/)
  R(2,:) = (/ -sin(angle), cos(angle), 0.d0/)
  R(3,:) = (/ 0.d0, 0.d0, 1.d0  /)
end subroutine

! rotation around arbitrary axis, given by vector u (normalized!!)
! source : https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
subroutine Rarb (R, angle, u)
  implicit none
  real(kind=rk), intent (in) :: angle
  real(kind=rk), intent (in) :: u(1:3)
  real(kind=rk),dimension(1:3,1:3), intent(out) :: R

  R(1,:) = (/ cos(angle)+u(1)**2*(1.0_rk-cos(angle))       , u(1)*u(2)*(1.0_rk-cos(angle))-u(3)*sin(angle), u(1)*u(3)*(1.0_rk-cos(angle))+u(2)*sin(angle) /)
  R(2,:) = (/ u(2)*u(1)*(1.0_rk-cos(angle))+u(3)*sin(angle), cos(angle)+u(2)**2*(1.0_rk-cos(angle))       , u(2)*u(3)*(1.0_rk-cos(angle))-u(1)*sin(angle) /)
  R(3,:) = (/ u(3)*u(1)*(1.0_rk-cos(angle))-u(2)*sin(angle), u(3)*u(2)*(1.0_rk-cos(angle))+u(1)*sin(angle), cos(angle)+u(3)**2*(1.0_rk-cos(angle))        /)

  ! because wikipedia uses active rotations
  R = TRANSPOSE(R)
end subroutine



subroutine rotation_matrix_from_quaternion (ep,R)
  implicit none
  real(kind=rk), intent (in) :: ep(0:3)
  real(kind=rk),dimension(1:3,1:3), intent(out) :: R

  ! code from mod__func.f90 from Maeda
  R(1,1) = ep(1)*ep(1) +ep(0)*ep(0) -ep(2)*ep(2) -ep(3)*ep(3)
  R(1,2) = 2.0d0*(ep(1)*ep(2) -ep(3)*ep(0))
  R(1,3) = 2.0d0*(ep(1)*ep(3) +ep(2)*ep(0))
  R(2,1) = 2.0d0*(ep(1)*ep(2) +ep(3)*ep(0))
  R(2,2) = ep(2)*ep(2) +ep(0)*ep(0) -ep(1)*ep(1) -ep(3)*ep(3)
  R(2,3) = 2.0d0*(ep(2)*ep(3) -ep(1)*ep(0))
  R(3,1) = 2.0d0*(ep(1)*ep(3) -ep(2)*ep(0))
  R(3,2) = 2.0d0*(ep(2)*ep(3) +ep(1)*ep(0))
  R(3,3) = ep(3)*ep(3) +ep(0)*ep(0) -ep(1)*ep(1) -ep(2)*ep(2)

  ! note Maeda's code was from body->global
  ! but for us, M_body is global->body
  R = transpose(R)

end subroutine
