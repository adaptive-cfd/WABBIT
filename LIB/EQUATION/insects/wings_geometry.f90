
! the new routine (2/2019) creates the wings (if both wings are used, maybe just one is)
! and their solid body velocity field us. Note that us contains both contributions from
! body and wing motion.
subroutine draw_insect_wings(time, xx0, ddx, mask, mask_color, us, Insect, delete)
  implicit none

  real(kind=rk), intent(in)    :: time
  type(diptera), intent(inout) :: Insect
  real(kind=rk), intent(in)    :: xx0(1:3), ddx(1:3)
  real(kind=rk), intent(inout) :: mask(0:,0:,0:)
  real(kind=rk), intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2), intent(inout) :: mask_color(0:,0:,0:)
  logical, intent(in) :: delete

  integer :: ix, iy, iz
  real(kind=rk), dimension(1:3) :: x_glob, x_body, v_tmp
  integer(kind=2) :: c

  if (size(mask) /= size(mask_color) .or. size(us,4) /= 3) then
      write(*,*) "mask:", shape(mask), "mask_color:", shape(mask_color), "us:", shape(us)
      call abort (08021902,"Insects: arrays have wrong size..")
  endif

  if (delete) then
      where (mask_color==Insect%color_r .or. mask_color==Insect%color_l .or. &
             mask_color==Insect%color_r2 .or. mask_color==Insect%color_l2)
          mask = 0.00_rk
          us(:,:,:,1) = 0.00_rk
          us(:,:,:,2) = 0.00_rk
          us(:,:,:,3) = 0.00_rk
          mask_color = 0
      end where
  endif

  if ((dabs(Insect%time-time)>1.0d-10).and.root) then
      write(*,'("error! time=",es15.8," but Insect%time=",es15.8)') time, Insect%time
      write(*,'("Did you call Update_Insect before draw_insect_wings?")')
  endif

  ! 28/01/2019: Thomas. Discovered that this was done block based, i.e. the smoothing layer
  ! had different thickness, if some blocks happened to be at different levels (and still carry
  ! a part of the smoothing layer.) I don't know if that made sense, because the layer shrinks/expands then
  ! and because it might be discontinous. Both options are included now, default is "as before"
  ! Insect%smoothing_thickness=="local"  : smoothing_layer = c_sm * 2**-J * L/(BS-1)
  ! Insect%smoothing_thickness=="global" : smoothing_layer = c_sm * 2**-Jmax * L/(BS-1)
  ! NOTE: for FLUSI, this has no impact! Here, the grid is constant and equidistant.
  if (Insect%smoothing_thickness=="local" .or. .not. grid_time_dependent) then
      Insect%smooth = 1.0d0*maxval(ddx)
      Insect%safety = 3.5d0*Insect%smooth
  endif

  !-----------------------------------------------------------------------------
  ! Stage I: mask + us field in wing system
  !-----------------------------------------------------------------------------
  if (Insect%RightWing == "yes") then
      call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%color_r, &
      Insect%M_body, Insect%M_wing_r, Insect%x_pivot_r_b, Insect%rot_rel_wing_r_w, "R" )
  endif

  if (Insect%LeftWing == "yes") then
      call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%color_l, &
      Insect%M_body, Insect%M_wing_l, Insect%x_pivot_l_b, Insect%rot_rel_wing_l_w, "L" )
  endif

  if (Insect%RightWing2 == "yes") then
      call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%color_r2, &
      Insect%M_body, Insect%M_wing_r2, Insect%x_pivot_r2_b, Insect%rot_rel_wing_r2_w, "R" )
  endif

  if (Insect%LeftWing2 == "yes") then
      call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%color_l2, &
      Insect%M_body, Insect%M_wing_l2, Insect%x_pivot_l2_b, Insect%rot_rel_wing_l2_w, "L" )
  endif

  !-----------------------------------------------------------------------------
  ! stage II: add body motion to wing and bring us to global system
  !-----------------------------------------------------------------------------
  ! Add solid body rotation (i.e. the velocity field that originates
  ! from the body rotation and translation). Until now, the wing velocities
  ! were the only ones set plus they are in the body reference frame
  do iz = g, size(mask,3)-1-g
      x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

              c = mask_color(ix,iy,iz)
              ! skip all parts that do not belong to the wings (ie they have a different color)
              if (c==Insect%color_l .or. c==Insect%color_r .or. &
                  c==Insect%color_l2 .or. c==Insect%color_r2 ) then

                  if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                  x_body = matmul(Insect%M_body, x_glob)

                  ! add solid body rotation in the body-reference frame, if color
                  ! indicates that this part of the mask belongs to the wings
                  if (mask(ix,iy,iz) > 0.d0) then

                      ! translational part. we compute the rotational part in the body
                      ! reference frame, therefore, we must transform the body translation
                      ! velocity Insect%vc (which is in global coordinates) to the body frame
                      v_tmp = matmul(Insect%M_body, Insect%vc_body_g)

                      ! add solid body rotation to the translational velocity field. Note
                      ! that rot_body_b and x_body are in the body reference frame
                      v_tmp(1) = v_tmp(1) + Insect%rot_body_b(2)*x_body(3)-Insect%rot_body_b(3)*x_body(2)
                      v_tmp(2) = v_tmp(2) + Insect%rot_body_b(3)*x_body(1)-Insect%rot_body_b(1)*x_body(3)
                      v_tmp(3) = v_tmp(3) + Insect%rot_body_b(1)*x_body(2)-Insect%rot_body_b(2)*x_body(1)

                      ! the body motion is added to the wing motion, which is already in us
                      ! and they are also in the body refrence frame. However, us has to be
                      ! in the global reference frame, so M_body_inverse is applied
                      us(ix,iy,iz,1:3) = matmul( Insect%M_body_inv, us(ix,iy,iz,1:3)+v_tmp )
                  endif
              endif
          enddo
      enddo
  enddo

end subroutine



! Wing wrapper for different wing shapes
subroutine draw_wing(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body,&
    M_wing, x_pivot_b, rot_rel_wing_w, side)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk), intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)
  ! NOTE: for a corrugated wing, up- and downside are different, and therefore a distinction between the 
  ! left- and right wing has to be made, essentially inverting the sign of the z_wing coordinate.
  character(len=1), intent(in) :: side ! can be R or L
  character(len=strlen) :: wingshape_str
  integer(kind=2) :: wingID

  !-- wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  wingID = color_wing-1

  select case(Insect%WingShape(wingID))
  case ("pointcloud")
      call draw_wing_pointcloud(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, &
      x_pivot_b,rot_rel_wing_w)

  case ("mosquito_iams")
      call draw_wing_mosquito(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, &
      M_wing,x_pivot_b,rot_rel_wing_w)

  case ("rectangular")
      call draw_wing_rectangular(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, &
      M_wing,x_pivot_b,rot_rel_wing_w)

  case ("suzuki")
      ! this wing has a finite thickness
      call draw_wing_suzuki(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, &
      M_wing,x_pivot_b,rot_rel_wing_w)

  case ("TwoEllipses")
      call draw_wing_twoellipses(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, &
      M_wing,x_pivot_b,rot_rel_wing_w)

  case default
      ! if all other options fail, we still might load coefficients from file:
      wingshape_str = Insect%WingShape(wingID)

          ! we assume the default to be defined in fourier coefficients, the subroutine
          ! yells if it does not recongnize the wing.
          select case (Insect%wing_file_type(wingID))
          case ("fourier")
              ! ordinary fourier wing (wing planform described in polar coordinates with fourier coeffs for the radius)
              call draw_wing_fourier(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, &
          x_pivot_b,rot_rel_wing_w, side)

      case ("fourierY")
          ! fourier series for the y coordinate (used for the blade of a bristled wing)
          call draw_wing_bristled(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, &
              x_pivot_b,rot_rel_wing_w)

          case ("kleemeier")
              ! kleemeier wings is bristles with rectangular central membrane. it is separated because the rectangular
              ! membrane is bad for the fourier series.
              call draw_wing_kleemeier(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, &
              x_pivot_b,rot_rel_wing_w)

          case default
              call abort(26111901, "The wing-ini-setup has a TYPE setting that the code does not know: "//trim(adjustl(Insect%wing_file_type(wingID))))
          end select

  end select

end subroutine draw_wing

!-------------------------------------------------------------------------------

! Draws a wings that is given by a radius(theta), where the radius is given
! by a Fourier series. The Fourier coefficients are stored in the insect
! datastructure, so the function Set_Wing_Fourier_coefficients must be called
! before calling this subroutine. Fourier series is evaluated in
! Radius_Fourier
subroutine draw_wing_fourier(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, x_pivot_b, rot_rel_wing_w, side)
  implicit none

  type(diptera), intent(inout) :: Insect
  real(kind=rk), intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk), intent(inout) :: mask(0:,0:,0:)
  real(kind=rk), intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2), intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2), intent(in) :: color_wing
  real(kind=rk),intent(in) :: M_body(1:3,1:3), M_wing(1:3,1:3), x_pivot_b(1:3), rot_rel_wing_w(1:3)
  ! NOTE: for a corrugated wing, up- and downside are different, and therefore a distinction between the 
  ! left- and right wing has to be made, essentially inverting the sign of the z_wing coordinate.
  character(len=1), intent(in) :: side

  integer :: ix,iy,iz,j
  integer(kind=2) :: wingID
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3), xa(1:3), xb(1:3)
  real(kind=rk) :: R, R0, R_tmp, zz0
  real(kind=rk) :: y_tmp, x_tmp, z_tmp, s, t
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta, sign

  !-- wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  wingID = color_wing-1

  if ( ((Insect%wing_file_type(wingID)) /= "fourier") .and. ((Insect%wing_file_type(wingID)) /= "fourierY") ) call abort(26111902,"draw_wing_fourier is called with non-fourier wing...")

  if (side == "R") then
      sign = +1.0_rk
  elseif (side == "L") then
      sign = -1.0_rk
  else
      call abort(290720, "neither R nor L wing??")
  endif

  s = Insect%safety
  do iz = g, size(mask,3)-1-g
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2)+dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1)+dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
              x_body = matmul(M_body,x)
              x_wing = matmul(M_wing,x_body-x_pivot_b)

              ! bounding box check: does this point lie within the bounding box? Note Insect%wing_bounding_box
              ! is set in SET_WING_BOUNDING_BOX_FOURIER
              if ( x_wing(1) >= Insect%wing_bounding_box(1,wingID)-s &
                    .and. x_wing(1) <= Insect%wing_bounding_box(2,wingID)+s) then
                  if ( x_wing(2) >= Insect%wing_bounding_box(3,wingID)-s &
                        .and. x_wing(2) <= Insect%wing_bounding_box(4,wingID)+s) then
                      if ( x_wing(3) >= Insect%wing_bounding_box(5,wingID)-s &
                            .and. x_wing(3) <= Insect%wing_bounding_box(6,wingID)+s) then

                          !-- get normalized angle (theta)
                          theta = atan2( x_wing(2)-Insect%yc(wingID), x_wing(1)-Insect%xc(wingID) )
                          theta = ( theta + pi ) / (2.d0*pi)

                          !-- construct R by evaluating the fourier series
                          R0 = Radius_Fourier(theta,Insect,wingID)

                          !-- get smooth (radial) step function
                          R = dsqrt ( (x_wing(1)-Insect%xc(wingID))**2 + (x_wing(2)-Insect%yc(wingID))**2 )
                          R_tmp = steps(R,R0, Insect%smooth)

                          ! wing corrugation (i.e. deviation from a flat plate)
                          if ( Insect%corrugated(wingID) ) then
                              ! if the wing is corrugated, its height profile is read from ini file
                              ! and interpolated at the position on the wing
                              zz0 = interp2_nonper( x_wing(1), x_wing(2), corrugation_profile, Insect%corrugation_array_bbox(1:4) )
                          else
                              ! no corrugation - the wing is a flat surface
                              zz0 = 0.0_pr
                          endif

                          zz0 = zz0 * sign

                          ! wing thickness
                          if ( Insect%wing_thickness_distribution(wingID)=="variable") then
                              ! variable wing thickness is read from an array in the wing.ini file
                              ! and interpolated linearly at the x_wing position.
                              t = interp2_nonper( x_wing(1), x_wing(2), wing_thickness_profile, Insect%corrugation_array_bbox(1:4) )
                          else
                              ! constant thickness, read from main params.ini file
                              t = Insect%WingThickness
                          endif

                          z_tmp = steps( dabs(x_wing(3)-zz0), 0.5d0*t, Insect%smooth ) ! thickness
                          mask_tmp = z_tmp*R_tmp

                          !-----------------------------------------
                          ! set new value for mask and velocity us
                          !-----------------------------------------
                          if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0)) then

                              mask(ix,iy,iz) = mask_tmp
                              mask_color(ix,iy,iz) = color_wing

                              !------------------------------------------------
                              ! solid body rotation
                              ! Attention: the Matrix transpose(M) brings us back to the body
                              ! coordinate system, not to the inertial frame. this is done in
                              ! the main routine Draw_Insect
                              !------------------------------------------------
                              ! v_tmp is in wing system:
                              v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                              v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                              v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                              ! note we set this only if it is a part of the wing
                              ! us is now in body system (note M_wing contains the stroke plane)
                              us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                          endif
                      endif
                  endif
              endif

          enddo
      enddo
  enddo

  !-----------------------------------------------------------------------------
  ! bristles
  !-----------------------------------------------------------------------------
  ! generic fourier wings can also have bristles: they are read from an inifile
  if (Insect%bristles(wingID)) then
      ! Loop for all bristles
      do j = 1, Insect%n_bristles(wingID)
          ! start / end point (in wing coordinate system)
          xa = (/Insect%bristles_coords(wingID,j,1), Insect%bristles_coords(wingID,j,2), 0.0d0/)
          xb = (/Insect%bristles_coords(wingID,j,3), Insect%bristles_coords(wingID,j,4), 0.0d0/)
          R = Insect%bristles_coords(wingID,j,5)

          ! note input to draw_bristle in in wing coordinates
          call draw_bristle(xa, xb, R, xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, x_pivot_b, rot_rel_wing_w)
      enddo
  endif

end subroutine draw_wing_fourier


subroutine draw_wing_kleemeier(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, x_pivot_b, rot_rel_wing_w)
  implicit none

  type(diptera), intent(inout) :: Insect
  real(kind=rk), intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk), intent(inout) :: mask(0:,0:,0:)
  real(kind=rk), intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2), intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2), intent(in) :: color_wing
  real(kind=rk),intent(in) :: M_body(1:3,1:3), M_wing(1:3,1:3), x_pivot_b(1:3), rot_rel_wing_w(1:3)

  integer :: ix,iy,iz,j
  integer(kind=2) :: wingID
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3), xa(1:3), xb(1:3)
  real(kind=rk) :: R, R0, R_tmp, zz0
  real(kind=rk) :: y_tmp, x_tmp, z_tmp, s, t
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta
  real(kind=rk) :: L_membrane, c_membrane

  !-- wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  wingID = color_wing-1

  if ((Insect%wing_file_type(wingID)) /= "kleemeier") call abort(26111902,"draw_wing_kleemeier called with non-kleemeier wing...")

  s = Insect%safety
  L_membrane = Insect%L_membrane(wingID)
  c_membrane = Insect%B_membrane(wingID)

  !-----------------------------------------------------------------------------
  ! membrane (rectangle)
  !-----------------------------------------------------------------------------
  do iz = g, size(mask,3)-1-g
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

              x_body = matmul(M_body,x)
              x_wing = matmul(M_wing,x_body-x_pivot_b)

              ! spanwise length:
              if ((x_wing(2)>=-s).and.(x_wing(2)<=L_membrane+s)) then
                  ! thickness: (note left and right wing have a different orientation of the z-axis
                  ! but this does not matter since this is the same.
                  if (dabs(x_wing(3))<=0.5*Insect%WingThickness + s) then
                      ! in the x-direction, the actual wing shape plays.
                      if ((x_wing(1)>-c_membrane/2.0-s).and.(x_wing(1)<c_membrane/2.0+s)) then
                          !-- smooth length
                          if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
                              y_tmp = steps(-x_wing(2), 0.d0, Insect%smooth)
                          else
                              y_tmp = steps( x_wing(2), L_membrane, Insect%smooth)
                          endif

                          !-- smooth height
                          z_tmp = steps(dabs(x_wing(3)), 0.5d0*Insect%WingThickness, Insect%smooth) ! thickness

                          !-- smooth shape
                          if (x_wing(1)<0.d0) then
                              x_tmp = steps(-x_wing(1), c_membrane/2.0, Insect%smooth)
                          else
                              x_tmp = steps( x_wing(1), c_membrane/2.0, Insect%smooth)
                          endif

                          mask_tmp = z_tmp*y_tmp*x_tmp

                          if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0)) then
                              mask(ix,iy,iz) = mask_tmp
                              mask_color(ix,iy,iz) = color_wing
                              !------------------------------------------------
                              ! solid body rotation
                              ! Attention: the Matrix transpose(M) brings us back to the body
                              ! coordinate system, not to the inertial frame. this is done in
                              ! the main routine Draw_Insect
                              !------------------------------------------------
                              v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                              v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                              v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                              ! note we set this only if it is a part of the wing
                              us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                          endif
                      endif
                  endif
              endif
          enddo
      enddo
  enddo

  !-----------------------------------------------------------------------------
  ! bristles
  !-----------------------------------------------------------------------------
  ! generic fourier wings can also have bristles: they are read from an inifile
  if (Insect%bristles(wingID)) then
      ! Loop for all bristles
      do j = 1, Insect%n_bristles(wingID)
          ! start / end point (in wing coordinate system)
          xa = (/Insect%bristles_coords(wingID,j,1), Insect%bristles_coords(wingID,j,2), 0.0d0/)
          xb = (/Insect%bristles_coords(wingID,j,3), Insect%bristles_coords(wingID,j,4), 0.0d0/)
          R = Insect%bristles_coords(wingID,j,5)

          ! note input to draw_bristle in in wing coordinates
          call draw_bristle(xa, xb, R, xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, x_pivot_b, rot_rel_wing_w)
      enddo
  endif

end subroutine draw_wing_kleemeier


!-------------------------------------------------------------------------------
! Draws a membranous central part of a bristled wing, using the same storage spase as
! for a Fourier wing, but the algorithm is different.
subroutine draw_blade_fourier(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

  integer :: ix,iy,iz
  integer(kind=2) :: wingID
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=rk) :: R, R0, R_tmp, zz0
  real(kind=rk) :: y_tmp, x_tmp, z_tmp, s, t
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta
  real(kind=rk) :: rblade, ylte, xte, xle


  !-- wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  wingID = color_wing-1

  !-- reset the bounding box
  Insect%wing_bounding_box(1:4,wingID) = (/-1.0d0, 1.0d0, 0.0d0, 1.0d0/)

  !-- blade length
  rblade = Insect%yc(wingID)

  s = Insect%safety
  do iz = g, size(mask,3)-1-g
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2)+dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1)+dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
              x_body = matmul(M_body,x)
              x_wing = matmul(M_wing,x_body-x_pivot_b)

              ! bounding box check: does this point lie within the bounding box? Note Insect%wing_bounding_box
              ! is set in SET_WING_BOUNDING_BOX_FOURIER
              if ( x_wing(1) >= Insect%wing_bounding_box(1,wingID)-s &
                    .and. x_wing(1) <= Insect%wing_bounding_box(2,wingID)+s) then
!                  if ( x_wing(2) >= Insect%wing_bounding_box(3,wingID)-s .and. x_wing(2) <= Insect%wing_bounding_box(4,wingID)+s) then
                  if ( x_wing(2) > 0.0d0 .and. x_wing(2) < rblade ) then
                      if ( x_wing(3) >= Insect%wing_bounding_box(5,wingID)-s &
                            .and. x_wing(3) <= Insect%wing_bounding_box(6,wingID)+s) then

                          !-- calculate the polar parameter (normalized angle)
                          ylte = x_wing(2)
                          theta = dacos( 1.0d0 - 2.0d0*ylte/rblade )
                          theta = theta / (2.d0*pi)

                          !-- construct xle by evaluating the Fourier series
                          xle = Radius_Fourier(theta,Insect,wingID)

                          !-- construct xte by evaluating the Fourier series
                          xte = Radius_Fourier(1.0d0-theta,Insect,wingID)

                          !-- amplitude
                          R0 = 0.5*(xle-xte)

                          !-- get smooth rectangular function
                          R = dabs ( x_wing(1) - 0.5*(xte+xle) )
                          R_tmp = steps(R,R0, Insect%smooth)

                          ! wing corrugation (i.e. deviation from a flat plate)
                          if ( Insect%corrugated(wingID) ) then
                              ! if the wing is corrugated, its height profile is read from ini file
                              ! and interpolated at the position on the wing
                              zz0 = interp2_nonper( x_wing(1), x_wing(2), corrugation_profile, Insect%corrugation_array_bbox(1:4) )
                          else
                              ! no corrugation - the wing is a flat surface
                              zz0 = 0.0_pr
                          endif

                          ! wing thickness
                          if ( Insect%wing_thickness_distribution(wingID)=="variable") then
                              ! variable wing thickness is read from an array in the wing.ini file
                              ! and interpolated linearly at the x_wing position.
                              t = interp2_nonper( x_wing(1), x_wing(2), wing_thickness_profile, Insect%corrugation_array_bbox(1:4) )
                          else
                              ! constant thickness, read from main params.ini file
                              t = Insect%WingThickness
                          endif

                          z_tmp = steps( dabs(x_wing(3)-zz0), 0.5d0*t, Insect%smooth ) ! thickness
                          mask_tmp = z_tmp*R_tmp

                          !-----------------------------------------
                          ! set new value for mask and velocity us
                          !-----------------------------------------
                          if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0)) then
                              mask(ix,iy,iz) = mask_tmp
                              mask_color(ix,iy,iz) = color_wing
                          endif
                      endif
                  endif
              endif

          enddo
      enddo
  enddo

end subroutine draw_blade_fourier


!-------------------------------------------------------------------------------
! Draw a wing from pointcloud
subroutine draw_wing_pointcloud(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

  integer :: ix,iy,iz,npoints,matrixcols
  real(kind=rk) :: x_glob(1:3), x_wing(1:3), x(1:3), x_body(1:3)
  real(kind=rk) :: v_tmp(1:3), mask_tmp

  real(kind=rk) :: xd, yd, zd, dxinv
  real(kind=rk) :: c00, c10, c01, c11, c0, c1
  integer :: iix, iiy, iiz

  call abort(211019,"not ready: initialization has to be moved to insect_init")

  ! ----------------------------------------------------------------------------
  ! Step 0.
  ! initialization, read pointcloud from file. generate maks once, then later we
  ! move it with interpolation
  ! ----------------------------------------------------------------------------
  if ( .not. allocated(mask_wing_complete) ) then
    if (root) write(*,*) "----------------pointcloudwing initialization -------------------"
    if (root) write(*,*) "Reading pointcloud for wing from file="//trim(adjustl(Insect%pointcloudfile))

    call check_file_exists(Insect%pointcloudfile)

    call count_lines_in_ascii_file_mpi(Insect%pointcloudfile, npoints, 1)
    call count_cols_in_ascii_file_mpi(Insect%pointcloudfile, matrixcols, 1)

    if (matrixcols /= 6) call abort(1230, "pointcloudwing: input file does not have 6 cols.")

    ! allocate memory (x,y,z,nx,ny,nz)
    allocate( particle_points(1:npoints, 1:matrixcols))

    ! read the data from the ascii file. Here, we suppose it to have 6 columns (and
    ! thus do not skip any garbage cols in between, e.g. for surface colors etc)
    call read_array_from_ascii_file_mpi(Insect%pointcloudfile, particle_points(:,1:matrixcols), 1)

    if (root) write(*,*) "Internally creating the signed distance for your point cloud."
    if (root) write(*,*) "This is done redundantly on all mpiranks, so each mpirank holds a complete copy."
    if (root) write(*,*) "For the domain of this array, we choose the smallest possible."

    mask_wing_x0 = (/minval(particle_points(:,1)),minval(particle_points(:,2)),minval(particle_points(:,3))/)&
                 - dble(mask_wing_safety)*ddx
    mask_wing_xl = (/maxval(particle_points(:,1))-mask_wing_x0(1) + dble(mask_wing_safety)*ddx(1), &
                     maxval(particle_points(:,2))-mask_wing_x0(2) + dble(mask_wing_safety)*ddx(2), &
                     maxval(particle_points(:,3))-mask_wing_x0(3) + dble(mask_wing_safety)*ddx(3)/)
    mask_wing_nxyz = nint(mask_wing_xl / ddx)

    ! note each mpirank holds one complete copy of the array (this is not very efficient)
    allocate(mask_wing_complete(0:mask_wing_nxyz(1)-1,0:mask_wing_nxyz(2)-1,0:mask_wing_nxyz(3)-1))

    ! create the mask, note we do not create color or us here.
    call mask_from_pointcloud(particle_points(:,1:3), particle_points(:,4:6), &
    mask_wing_x0, ddx, mask_wing_complete, mask_wing_safety, Insect%smooth, 0.0d0*ddx(1))

    ! after generating the mask function, we do not longer need the point cloud.
    deallocate( particle_points )

    if (root) write(*,*) "Pointcloudwing: resolution is", mask_wing_nxyz
    if (root) write(*,*) "Pointcloudwing: origin is", mask_wing_x0
    if (root) write(*,*) "Pointcloudwing: domain size is", mask_wing_xl
    if (root) write(*,*) "----------------end pointcloudwing initialization ---------------"
  end if

  !-----------------------------------------------------------------------------
  ! interpolation
  !-----------------------------------------------------------------------------
  dxinv = 1.0d0 / ddx(1)

  do iz = g, size(mask,3)-1-g
    x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
    do iy = g, size(mask,2)-1-g
      x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
      do ix = g, size(mask,1)-1-g
        x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

        !-- define the various coordinate systems we are going to use
        if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot_b)


        if (x_wing(1)>mask_wing_x0(1)+2.0d0*ddx(1) .and. x_wing(1)<mask_wing_x0(1)+mask_wing_xl(1)-2.0d0*ddx(1)) then
          if (x_wing(2)>mask_wing_x0(2)+2.0d0*ddx(2) .and. x_wing(2)<mask_wing_x0(2)+mask_wing_xl(2)-2.0d0*ddx(2)) then
            if (x_wing(3)>mask_wing_x0(3)+2.0d0*ddx(3) .and. x_wing(3)<mask_wing_x0(3)+mask_wing_xl(3)-2.0d0*ddx(3)) then

              ! use 3d interpolation
              ! mask_tmp = trilinear_interp( mask_wing_x0, (/dx,dx,dx/), mask_wing_complete, x_wing, .false.)

              ! indices of cube containing the target point, lower end
              iix = floor( (x_wing(1)-mask_wing_x0(1))*dxinv )
              iiy = floor( (x_wing(2)-mask_wing_x0(2))*dxinv )
              iiz = floor( (x_wing(3)-mask_wing_x0(3))*dxinv )

              ! distance to lower point, normalized (0..1)
              xd = ( x_wing(1)-(dble(iix)*ddx(1) + mask_wing_x0(1)) ) *dxinv
              yd = ( x_wing(2)-(dble(iiy)*ddx(2) + mask_wing_x0(2)) ) *dxinv
              zd = ( x_wing(3)-(dble(iiz)*ddx(3) + mask_wing_x0(3)) ) *dxinv

              c00 = mask_wing_complete(iix,iiy  ,iiz  )*(1.d0-xd)+mask_wing_complete(iix+1,iiy  ,iiz )*xd
              c10 = mask_wing_complete(iix,iiy+1,iiz  )*(1.d0-xd)+mask_wing_complete(iix+1,iiy+1,iiz )*xd
              c01 = mask_wing_complete(iix,iiy  ,iiz+1)*(1.d0-xd)+mask_wing_complete(iix+1,iiy  ,iiz+1)*xd
              c11 = mask_wing_complete(iix,iiy+1,iiz+1)*(1.d0-xd)+mask_wing_complete(iix+1,iiy+1,iiz+1)*xd

              c0 = c00*(1.d0-yd) + c10*yd
              c1 = c01*(1.d0-yd) + c11*yd

              mask_tmp = c0*(1.d0-zd) + c1*zd

              if (mask_tmp>0.0d0) then
                ! yo
                mask(ix,iy,iz) = mask_tmp
                ! it was valid -> assign color
                mask_color(ix,iy,iz) = color_wing
              endif

            endif
          endif
        endif

      enddo
    enddo
  enddo

  ! ----------------------------------------------------------------------------
  ! add velocity field, inside the wing.
  ! ----------------------------------------------------------------------------
  do iz = g, size(mask,3)-1-g
    x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
    do iy = g, size(mask,2)-1-g
      x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
      do ix = g, size(mask,1)-1-g
        x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

        ! if this point belong to the wing we just created
        if (mask_color(ix,iy,iz)==color_wing) then
          !-- define the various coordinate systems we are going to use
          if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
          x_wing = matmul(M_wing, matmul(M_body,x)-x_pivot_b)

          !------------------------------------------------
          ! solid body rotation
          ! Attention: the Matrix transpose(M) brings us back to the body
          ! coordinate system, not to the inertial frame. this is done in
          ! the main routine Draw_Insect
          !------------------------------------------------
          v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
          v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
          v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

          ! note we set this only if it is a part of the wing
          us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
        endif
      end do
    end do
  end do
end subroutine draw_wing_pointcloud


!-------------------------------------------------------------------------------
! Suzuki's rectangular wingshape as defined in section 3.4 of my thesis (Thomas,
! "Numerical modeling of fluid-structure interaction in bio-inspired propulsion")
! This wing has finite thickness.
!-------------------------------------------------------------------------------
subroutine draw_wing_suzuki(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: color_wing
    real(kind=rk),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

    integer :: ix,iy,iz
    real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
    real(kind=rk) :: R, R0, R_tmp
    real(kind=rk) :: y_tmp, x_tmp, z_tmp, y_left, y_right
    real(kind=rk) :: v_tmp(1:3), mask_tmp, theta, x_top, x_bot

    ! wing shape (determine between which x-values (x_bot, x_top) the wing is
    ! these values depend on the spanwise direction (which is y)
    x_top = 0.0667d0
    x_bot = -0.35d0

    y_right = 1.0d0
    y_left = 0.1667d0

    do iz = g, size(mask,3)-1-g
        x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

                !-- define the various coordinate systems we are going to use
                if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                x_body = matmul(M_body,x)
                x_wing = matmul(M_wing,x_body-x_pivot_b)

                ! spanwise length:
                if ((x_wing(2)>=y_left-Insect%safety).and.(x_wing(2)<=y_right+Insect%safety)) then
                    ! thickness: (note left and right wing have a different orientation of the z-axis
                    ! but this does not matter since this is the same.
                    if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then

                        ! in the x-direction, the actual wing shape plays.
                        if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
                            !-- smooth length
                            if ( x_wing(2) < 0.5d0*(y_left+y_right) ) then
                                y_tmp = steps(-(x_wing(2)-y_left), 0.d0, Insect%smooth)
                            else
                                y_tmp = steps( (x_wing(2)-y_left), y_right-y_left, Insect%smooth)
                            endif

                            !-- smooth height
                            z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness, Insect%smooth) ! thickness

                            !-- smooth shape
                            if (x_wing(1) < 0.d0) then
                                x_tmp = steps(-x_wing(1),-x_bot, Insect%smooth)
                            else
                                x_tmp = steps( x_wing(1), x_top, Insect%smooth)
                            endif

                            mask_tmp = z_tmp*y_tmp*x_tmp

                            if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0d0)) then
                                mask(ix,iy,iz) = mask_tmp
                                mask_color(ix,iy,iz) = color_wing
                                !------------------------------------------------
                                ! solid body rotation
                                ! Attention: the Matrix transpose(M) brings us back to the body
                                ! coordinate system, not to the inertial frame. this is done in
                                ! the main routine Draw_Insect
                                !------------------------------------------------
                                v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                                v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                                v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                                ! note we set this only if it is a part of the wing
                                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                            endif

                        endif
                    endif
                endif
            enddo
        enddo
    enddo
end subroutine draw_wing_suzuki

!-------------------------------------------------------------------------------

subroutine draw_wing_rectangular(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

  integer :: ix,iy,iz
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=rk) :: R, R0, R_tmp
  real(kind=rk) :: y_tmp, x_tmp, z_tmp, y_left, y_right
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

  ! wing shape (determine between which x-values (x_bot, x_top) the wing is
  ! these values depend on the spanwise direction (which is y)
  x_top = 0.085d0 ! determinded from a calliphora wing, roughly estimated
  x_bot = -(0.294d0-x_top) ! to get the same aspect ratio as in calliphora (Engels et al., RSI2020)

  y_right = 1.0d0
  y_left = 0.0d0

  do iz = g, size(mask,3)-1-g
    x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
    do iy = g, size(mask,2)-1-g
      x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
      do ix = g, size(mask,1)-1-g
        x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

        !-- define the various coordinate systems we are going to use
        if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot_b)

        ! spanwise length:
              if ((x_wing(2)>=y_left-Insect%safety).and.(x_wing(2)<=y_right+Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then

            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
                          if ( x_wing(2) < 0.5d0*(y_left+y_right) ) then
                              y_tmp = steps(-(x_wing(2)-y_left), 0.d0, Insect%smooth)
              else
                              y_tmp = steps( (x_wing(2)-y_left), y_right-y_left, Insect%smooth)
              endif

              !-- smooth height
              z_tmp = steps(dabs(x_wing(3)), 0.5d0*Insect%WingThickness, Insect%smooth) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.d0) then
                x_tmp = steps(-x_wing(1),-x_bot, Insect%smooth)
              else
                x_tmp = steps( x_wing(1), x_top, Insect%smooth)
              endif

              mask_tmp = z_tmp*y_tmp*x_tmp

                          if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0d0)) then
                mask(ix,iy,iz) = mask_tmp
                mask_color(ix,iy,iz) = color_wing
                !------------------------------------------------
                ! solid body rotation
                ! Attention: the Matrix transpose(M) brings us back to the body
                ! coordinate system, not to the inertial frame. this is done in
                ! the main routine Draw_Insect
                !------------------------------------------------
                v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                ! note we set this only if it is a part of the wing
                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_rectangular


!-------------------------------------------------------------------------------
! Draws a wing
! here, a wing is a rigid plate of constant thickness that differs from
! a rectangular plate only in the x-direction
!
! note to save a bit of computing time, we first check the easy
! conditions (thickness and spanwise length) and then the shape
! function since this saves many evaluations of the shape.
subroutine draw_wing_twoellipses(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)


  integer :: ix,iy,iz
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=rk) :: R, R0, R_tmp,a_body
  real(kind=rk) :: y_tmp, x_tmp, z_tmp
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

  a_body = 0.5d0 * Insect%L_span

  do iz = g, size(mask,3)-1-g
    x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
    do iy = g, size(mask,2)-1-g
      x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
      do ix = g, size(mask,1)-1-g
        x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

        !-- define the various coordinate systems we are going to use
        if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot_b)

        ! spanwise length:
        if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
            ! wing shape (determine between which x-values (x_bot, x_top) the wing is
            ! these values depend on the spanwise direction (which is y)
            if ((1.d0 - ((x_wing(2)-a_body)**2)/(a_body**2)) >= 0.d0) then
              x_top =  dsqrt((Insect%b_top**2)*(1.d0-((x_wing(2)-a_body)**2)/(a_body**2)))
              x_bot = -dsqrt((Insect%b_bot**2)*(1.d0-((x_wing(2)-a_body)**2)/(a_body**2)))
            else
              x_top = 0.d0
              x_bot = 0.d0
            endif

            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
              if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
                y_tmp = steps(-x_wing(2),0.d0, Insect%smooth)
              else
                y_tmp = steps( x_wing(2),Insect%L_span, Insect%smooth)
              endif

              !-- smooth height
              z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness, Insect%smooth) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.d0) then
                x_tmp = steps(-x_wing(1),-x_bot, Insect%smooth)
              else
                x_tmp = steps( x_wing(1), x_top, Insect%smooth)
              endif

              mask_tmp = z_tmp*y_tmp*x_tmp

              if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0)) then
                mask(ix,iy,iz) = mask_tmp
                mask_color(ix,iy,iz) = color_wing
                !------------------------------------------------
                ! solid body rotation
                ! Attention: the Matrix transpose(M) brings us back to the body
                ! coordinate system, not to the inertial frame. this is done in
                ! the main routine Draw_Insect
                !------------------------------------------------
                v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                ! note we set this only if it is a part of the wing
                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_twoellipses

!-------------------------------------------------------------------------------
! Draws a wing of a mosquito. it is a simple ellipse shape, as presented in
! [1] Iams "Flight stability of mosquitos: A reduced model" SIAM J. Appl. Math. 74(5) 1535--1550 (2014)
subroutine draw_wing_mosquito(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

  integer :: ix,iy,iz
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=rk) :: R, R0, R_tmp,a_wing
  real(kind=rk) :: y_tmp, x_tmp, z_tmp
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

  a_wing = 0.5d0 * Insect%L_span
  Insect%b_top = 0.1474d0
  Insect%b_bot = 0.1474d0

  do iz = g, size(mask,3)-1-g
    x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
    do iy = g, size(mask,2)-1-g
      x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
      do ix = g, size(mask,1)-1-g
        x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

        !-- define the various coordinate systems we are going to use
        if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot_b)

        ! spanwise length:
        if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
            ! wing shape (determine between which x-values (x_bot, x_top) the wing is
            ! these values depend on the spanwise direction (which is y)
            if ((1.d0 - ((x_wing(2)-a_wing)**2)/(a_wing**2)) >= 0.d0) then
              x_top =  dsqrt((Insect%b_top**2)*(1.d0-((x_wing(2)-a_wing)**2)/(a_wing**2)))
              x_bot = -dsqrt((Insect%b_bot**2)*(1.d0-((x_wing(2)-a_wing)**2)/(a_wing**2)))
            else
              x_top = 0.d0
              x_bot = 0.d0
            endif

            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
              if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
                y_tmp = steps(-x_wing(2),0.d0, Insect%smooth)
              else
                y_tmp = steps( x_wing(2),Insect%L_span, Insect%smooth)
              endif

              !-- smooth height
              z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness, Insect%smooth) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.d0) then
                x_tmp = steps(-x_wing(1),-x_bot, Insect%smooth)
              else
                x_tmp = steps( x_wing(1), x_top, Insect%smooth)
              endif

              mask_tmp = z_tmp*y_tmp*x_tmp

              if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0)) then
                mask(ix,iy,iz) = mask_tmp
                mask_color(ix,iy,iz) = color_wing
                !------------------------------------------------
                ! solid body rotation
                ! Attention: the Matrix transpose(M) brings us back to the body
                ! coordinate system, not to the inertial frame. this is done in
                ! the main routine Draw_Insect
                !------------------------------------------------
                v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                ! note we set this only if it is a part of the wing
                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_mosquito


!-------------------------------------------------------------------------------
! Bristled wing 
!-------------------------------------------------------------------------------
subroutine draw_wing_bristled(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

  integer :: ix,iy,iz,j
  integer(kind=2) :: wingID
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3),xa(1:3),xb(1:3)
  real(kind=rk) :: R,s
  real(kind=rk) :: v_tmp(1:3)

  !-- wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  wingID = color_wing-1

  ! Draw the membranous blade using Fourier series
  call draw_blade_fourier(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,&
       x_pivot_b,rot_rel_wing_w)

  ! Set the solid velocity
  s = Insect%safety
  do iz = g, size(mask,3)-1-g
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2)+dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1)+dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
              x_body = matmul(M_body,x)
              x_wing = matmul(M_wing,x_body-x_pivot_b)

              ! bounding box check: does this point lie within the bounding box? Note Insect%wing_bounding_box
              ! is set in SET_WING_BOUNDING_BOX_FOURIER
              if ( x_wing(1) >= Insect%wing_bounding_box(1,wingID)-s &
                    .and. x_wing(1) <= Insect%wing_bounding_box(2,wingID)+s) then
                  if ( x_wing(2) >= Insect%wing_bounding_box(3,wingID)-s &
                        .and. x_wing(2) <= Insect%wing_bounding_box(4,wingID)+s) then
                      if ( x_wing(3) >= Insect%wing_bounding_box(5,wingID)-s &
                            .and. x_wing(3) <= Insect%wing_bounding_box(6,wingID)+s) then
                          !-----------------------------------------
                          ! set new value for solid velocity us
                          !-----------------------------------------
                          if ( (mask(ix,iy,iz)>0.0).and.(mask_color(ix,iy,iz)==color_wing) ) then
                              !------------------------------------------------
                              ! solid body rotation
                              ! Attention: the Matrix transpose(M) brings us back to the body
                              ! coordinate system, not to the inertial frame. this is done in
                              ! the main routine Draw_Insect
                              !------------------------------------------------
                              v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                              v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                              v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)
                              ! note we set this only if it is a part of the wing
                              us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                          endif
                      endif
                  endif
              endif

          enddo
      enddo
  enddo

  !-----------------------------------------------------------------------------
  ! bristles
  !-----------------------------------------------------------------------------
  ! generic fourier wings can also have bristles: they are read from an inifile
  if (Insect%bristles(wingID)) then
      ! Loop for all bristles
      do j = 1, Insect%n_bristles(wingID)
          ! start / end point (in wing coordinate system)
          xa = (/Insect%bristles_coords(wingID,j,1), Insect%bristles_coords(wingID,j,2), 0.0d0/)
          xb = (/Insect%bristles_coords(wingID,j,3), Insect%bristles_coords(wingID,j,4), 0.0d0/)
          R = Insect%bristles_coords(wingID,j,5)

          ! note input to draw_bristle in in wing coordinates
          call draw_bristle(xa, xb, R, xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, x_pivot_b, rot_rel_wing_w)
      enddo
  endif

  !-----------------------------------------------------------------------------
  ! effective membrane
  !-----------------------------------------------------------------------------
  if (Insect%bristles_simplex(wingID)) then
      ! Loop for all bristles
      do j = 1, Insect%n_bristles(wingID)-1
          ! draw a membrane triangular element
          call draw_triangle(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, x_pivot_b, rot_rel_wing_w, &
                            Insect%bristles_coords(wingID,j,1), Insect%bristles_coords(wingID,j,2), &
                            Insect%bristles_coords(wingID,j,3), Insect%bristles_coords(wingID,j,4), &
                            Insect%bristles_coords(wingID,j+1,1), Insect%bristles_coords(wingID,j+1,2))
          ! draw a membrane triangular element
          call draw_triangle(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_body, M_wing, x_pivot_b, rot_rel_wing_w, &
                            Insect%bristles_coords(wingID,j,3), Insect%bristles_coords(wingID,j,4), &
                            Insect%bristles_coords(wingID,j+1,3), Insect%bristles_coords(wingID,j+1,4), &
                            Insect%bristles_coords(wingID,j+1,1), Insect%bristles_coords(wingID,j+1,2))

      enddo
  endif

end subroutine draw_wing_bristled




!-------------------------------------------------------------------------------
! evaluates the fourier series given in the ai, bi
! NOTE: angle theta is NORMALIZED!  theta = ( theta + pi ) / (2.d0*pi)
!-------------------------------------------------------------------------------
real(kind=rk) function Radius_Fourier( theta, Insect, wingID )
  implicit none
  integer :: i,j, n_radius
  real(kind=rk) :: R0, theta2, dphi, area
  type(diptera),intent(inout)::Insect
  real(kind=rk), intent(in) :: theta
  integer(kind=2), intent(in) :: wingID ! wing id number

  n_radius = 25000
  dphi = (2.d0*pi) / (dble(n_radius-1))


  ! evaluate the entire R(theta) once with very fine resolution, so when
  ! calling it for the second time we only need linear interpolation.
  if (.not.Insect%wings_radius_table_ready(wingID)) then
    !---------------------------------------------------------------------------
    ! fill radius table
    !---------------------------------------------------------------------------
    Insect%R0_table(:,wingID) = 0.d0
    ! loop over all thetas and compute the radius for all of them, store it
    ! in the table Insect%R0_table
    do j = 1, n_radius
      R0 = Insect%a0_wings(wingID) / 2.d0
      theta2 = dble(j-1) * dphi
      ! evaluate Fourier series
      do i = 1, Insect%nfft_wings(wingID)
        R0 = R0 + Insect%ai_wings(i,wingID)*dcos(2.d0*pi*dble(i)*theta2) &
                + Insect%bi_wings(i,wingID)*dsin(2.d0*pi*dble(i)*theta2)
      enddo
      Insect%R0_table(j,wingID)=R0
    enddo

    ! skip setup on next call
    Insect%wings_radius_table_ready(wingID) = .true.
  endif

  ! linear interpolation, if already stored the radius
  j = floor( theta / dphi ) + 1
  Radius_Fourier = Insect%R0_table(j,wingID) + ((theta-dble(j-1)*dphi) / dphi) &
                 * (Insect%R0_table(j+1,wingID)-Insect%R0_table(j,wingID))
end function



!---------------------------------------------------------------------------
! compute wing surface area for Fourier wings
!---------------------------------------------------------------------------
! from the Fourier series, we can directly compute the wing Area (surface)
! the aera is the double integral A = \int(0,2pi) \int(0,R(theta)) r dr dtheta
subroutine compute_wing_surface( Insect, wingID, area )
  implicit none
  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(out) :: area
  integer(kind=2), intent(in) :: wingID ! wing id number

  real(kind=rk) :: dr, dtheta, r
  real(kind=rk) :: R0, theta

  ! the numerical precision is a tad exagerated, but since this is cheap, we do
  ! not really care.
  dr = 1.0d-4
  dtheta = 2.d0*pi/1000.d0
  theta = 0.d0
  r = 0.d0
  area = 0.d0

  ! this method currently works only for Fourier wings (i.e. the shape is described
  ! as Fourier coefficients). We return zero for non-Fourier wings.
  if ( Insect%nfft_wings(wingID) /= 0 ) then
    ! solve the double integral
    do while ( theta < 2.d0*pi )
        R0 = Radius_Fourier(theta,Insect,wingID)
        r = 0.d0
        do while (r < R0)
            area = area + r*dtheta*dr
            r = r+dr
        end do
        theta = theta + dtheta;
    end do
  endif
end subroutine


!-------------------------------------------------------------------------------
! Here all hard-coded fourier series coefficients for different wings shapes are
! collected. This routine is only called once per time step, and it doesn't do
! anything when called for the second time.
! In the first call, the arrays ai and bi, that hold the Fourier
! coefficients, are allocated. Then they are filled with the values corresponding
! to Insect%WingShape. If the routine is called with an unkown wing shape, it
! stops the code. This prevents errors for wings that are NOT given by Fourier
! series.
!-------------------------------------------------------------------------------
subroutine Setup_Wing_Fourier_coefficients(Insect, wingID)
  implicit none
  real(kind=rk) :: xroot, yroot
  type(diptera),intent(inout)::Insect
  integer(kind=2), intent(in) :: wingID ! wing id number
  character(len=strlen) :: wingshape_str

  if (Insect%wingsetup_done(wingID)) then
    ! the second call is just a return statement
    return
  endif

  Insect%a0_wings(wingID) = 0.d0
  Insect%ai_wings(:,wingID) = 0.d0
  Insect%bi_wings(:,wingID) = 0.d0

  !-----------------------------------------
  ! hard-coded Fourier coefficients for R(theta)
  !-----------------------------------------
  select case (Insect%WingShape(wingID))
  case ('pieris-brassicae1')
    ! butterfly (p. brassicae) wingshape, extracted from frame #633 of flight
    ! recording "160909_flight38"
    ! we have 29 points on the wing
    Insect%nfft_wings(wingID) = 30
    Insect%a0_wings(wingID) = 0.8262378
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.0333388,0.0088668,0.0910448,0.0051568,-0.0252108,-0.0125628,&
    -0.0080008,-0.0054058,0.0110598,0.0090758,-0.0076998,-0.0071098,&
    0.0018768,0.0003048,0.0033128,0.0049088,-0.0021308,-0.0051658,&
    -0.0003278,0.0022188,0.0015628,0.0008548,-0.0009428,-0.0027778,&
    -0.0005728,0.0021988,0.0010718,-0.0007128,-0.0010658,-0.0008818 /)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/-0.0974478,0.0697448,0.0277868,-0.0461138,-0.0059608,0.0080508,&
    -0.0059118,0.0113048,0.0038478,-0.0157198,-0.0070308,0.0034218,&
    0.0030868,0.0040238,0.0017568,-0.0027368,-0.0050618,-0.0008578,&
    0.0033778,0.0009348,-0.0012898,-0.0011198,-0.0016878,0.0000398,&
    0.0029238,0.0009108,-0.0020888,-0.0012848,-0.0003728,0.0004738  /)
    Insect%yc(wingID) = 0.3184928
    Insect%xc(wingID) = -0.2459908

  case ('drosophila')
    !********************************************
    ! Drosophila wing from Jan Gruber's png file
    !********************************************
    Insect%nfft_wings(wingID) = 40
    Insect%a0_wings(wingID) = 0.5140278
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
        (/0.1276258,-0.1189758,-0.0389458,0.0525938,0.0151538,-0.0247938,&
          -0.0039188,0.0104848,-0.0030638,-0.0064578,0.0042208,0.0043248,&
          -0.0026878,-0.0021458,0.0017688,0.0006398,-0.0013538,-0.0002038,&
          0.0009738,0.0002508,-0.0003548,-0.0003668,-0.0002798,0.0000568,&
          0.0003358,0.0001408,-0.0002208,0.0000028,0.0004348,0.0001218,&
          -0.0006458,-0.0003498,0.0007168,0.0003288,-0.0007078,-0.0001368,&
          0.0007828,0.0001458,-0.0007078,-0.0001358/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
        (/-0.1072518,-0.0449318,0.0296558,0.0265668,-0.0043988,-0.0113218,&
          -0.0003278,0.0075028,0.0013598,-0.0057338,-0.0021228,0.0036178,&
          0.0013328,-0.0024128,-0.0007688,0.0011478,0.0003158,-0.0005528,&
          0.0000458,0.0003768,0.0002558,0.0000168,-0.0006018,-0.0006338,&
          0.0001718,0.0007758,0.0001328,-0.0005888,-0.0001088,0.0006298,&
          0.0000318,-0.0008668,-0.0000478,0.0009048,0.0001198,-0.0008248,&
          -0.0000788,0.0007028,-0.0000118,-0.0006608/)

    ! wing root point
    xroot =+0.1122
    yroot =-0.0157
    ! center of circle
    Insect%xc(wingID) =-0.1206 + xroot
    Insect%yc(wingID) = 0.3619 + yroot
  case ('drosophila_mutated')
    !********************************************
    ! mutated Drosophila wing from Jan Gruber's png file
    !********************************************
    Insect%nfft_wings(wingID) = 70
    Insect%a0_wings(wingID) = 0.4812548
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
           (/0.1593968, -0.1056828, -0.0551518, 0.0508748, 0.0244538, -0.0264738,&
            -0.0080828, 0.0181228, 0.0023648, -0.0134578, -0.0037068, 0.0064508,&
            0.0028748, -0.0014258, -0.0006028, -0.0008898, -0.0020408, 0.0009218,&
            0.0029938, 0.0002768, -0.0026968, -0.0011518, 0.0017798, 0.0016538,&
            -0.0006098, -0.0012998, -0.0001918, 0.0003478, 0.0001408, 0.0003098,&
            0.0001078, -0.0005568, -0.0005998, 0.0006128, 0.0009078, -0.0003798,&
            -0.0009268, 0.0002128, 0.0009098, -0.0000598, -0.0010668, -0.0003428,&
            0.0009228, 0.0007688, -0.0003568, -0.0010458, -0.0004378, 0.0008738,&
            0.0009478, -0.0004108, -0.0012248, -0.0000638, 0.0013148, 0.0004978,&
            -0.0010638, -0.0007148, 0.0006338, 0.0007438, -0.0003278, -0.0006078,&
            0.0001838, 0.0003768, -0.0001698, -0.0002148, 0.0001318, 0.0001628,&
            -0.0000878, 0.0000068, 0.0001478, -0.0001128/)

    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
          (/-0.1132588, -0.0556428, 0.0272098, 0.0221478, -0.0063798, -0.0059078,&
              0.0043788, 0.0043208, -0.0003308, -0.0026598, -0.0013158, 0.0025178,&
              0.0022438, -0.0023798, -0.0037048, 0.0001528, 0.0031218, 0.0022248,&
              -0.0007428, -0.0027298, -0.0018298, 0.0014538, 0.0028888, 0.0000648,&
              -0.0023508, -0.0009418, 0.0017848, 0.0016578, -0.0008058, -0.0017348,&
              -0.0001368, 0.0011138, 0.0004218, -0.0005918, -0.0002798, 0.0002388,&
              0.0002148, 0.0001408, 0.0000218, -0.0005138, -0.0003458, 0.0008208,&
              0.0009888, -0.0007468, -0.0015298, 0.0002728, 0.0015588, 0.0002758,&
              -0.0012498, -0.0006908,0.0008718, 0.0008848, -0.0003038, -0.0008048,&
              -0.0001538, 0.0005418, 0.0003658, -0.0001988, -0.0003938, 0.0000048,&
              0.0003008, 0.0000538, -0.0002748, -0.0000598, 0.0002898, 0.0001398,&
              -0.0002108, -0.0001888, 0.0001838, 0.0001888 /)

    ! wing root point
    xroot =+0.1122
    yroot =-0.0157
    ! center of circle
    Insect%xc(wingID) =-0.1206 + xroot
    Insect%yc(wingID) = 0.3619 + yroot
  case ('drosophila_sandberg')
    !********************************************
    !  Drosophila wing from Ramamurti & Sandberg ( JEB 210, 881-896, 2007)
    !********************************************
    Insect%nfft_wings(wingID) = 24
    Insect%a0_wings(wingID) = 0.4995578
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.0164168,-0.1621518,0.0030938,0.0601108,-0.0083988,-0.0199988,&
    0.0049048,0.0047878,-0.0005648,-0.0001108,-0.0008638,-0.0006928,&
    0.0006608,0.0001978,0.0001558,0.0006878,-0.0007498,-0.0008018,&
    0.0003878,0.0007028,0.0000408,-0.0001108,-0.0001068,-0.0003958 &
    /)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/-0.2083518,-0.0106488,0.0878308,-0.0018168,-0.0338278,0.0045768,&
    0.0113778,-0.0020678,-0.0026928,0.0002758,-0.0000838,-0.0001298,&
    0.0004118,0.0005638,-0.0001018,-0.0006918,-0.0002268,0.0005238,&
    0.0004008,-0.0001818,-0.0003038,-0.0000068,-0.0001218,0.0002008 &
    /)
    Insect%xc(wingID) =-0.0235498
    Insect%yc(wingID) = 0.1531398
  case ('drosophila_maeda')
    !********************************************
    !  Drosophila wing from Maeda and Liu, similar to Liu and Aono, BB2009
    !********************************************
    Insect%nfft_wings(wingID) = 25
    !Insect%a0_wings(wingID) = 0.591294836514357
    !ai = (/0.11389995408864588, -0.08814321795213981, -0.03495210456149335,&
    !0.024972085605453047, 0.009422293191002384, -0.01680813499169695,&
    !-0.006006435254421029, 0.012157932943676907, 0.00492283934032996,&
    !-0.009882103857127606, -0.005421102356676356, 0.007230876076797827,&
    !0.005272314598249222, -0.004519437431722127, -0.004658072133773225,&
    !0.0030795046767766853, 0.003970792618725898, -0.0016315879319092456,&
    !-0.002415442110272326, 0.0011118187761994598, 0.001811261693911865,&
    !-2.6496695842951815E-4, -0.0012472769174353662, -1.7427507835680091E-4,&
    !0.0010049640224536927/)
    !bi = (/0.0961275426181888, 0.049085916171592914, -0.022051083533094627,&
    !-0.014004783021121204, 0.012955446778711292, 0.006539648525493488,&
    !-0.011873438993933363, -0.00691719567010525, 0.008479044683798266,&
    !0.0045388280405204194, -0.008252172088956379, -0.005091347100627815,&
    !0.004626409662755484, 0.004445034936616318, -0.0030708884306814804,&
    !-0.004428808427471962, 0.0014113707529017868, 0.003061279043478891,&
    !-8.658653756413232E-4, -0.002153349816945423, 3.317570161883452E-4,&
    !0.001573518502682025, 2.14583094242007E-4, -0.0011299834277813852,&
    !-5.172854674801216E-4/)
    Insect%a0_wings(wingID) = 0.585432698694358
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.113400475583443, -0.0862823485047213, -0.0346234482214816,&
    0.0237625254732323,0.00902498439287132,-0.0158926757445186,&
    -0.00549384372979449,0.0114928668063701,0.00431222381497978,&
    -0.00951270119733201,-0.00484045133879639,0.00706223174320460,&
    0.00473736389439926,-0.00449539769983697,-0.00418487169011745,&
    0.00320520052884641,0.00355631891057573,-0.00183155403463614,&
    -0.00191680264797099,0.00144768631289857,0.00135580122365068,&
    -0.000579638217642394,-0.000818378434108882,0.000132570375969864,&
    0.000683325977327827/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/0.0939265226506824,0.0486063180327962,-0.0206591129298861,&
    -0.0136085709392758,0.0118575265347540,0.00604510770670991,&
    -0.0110263907282936,-0.00636979352727611,0.00786779216718321,&
    0.00390493804324433,-0.00797763174406198,-0.00450123591642554,&
    0.00445099872769504,0.00387237248613979,-0.00305464314668877,&
    -0.00398381251524846,0.00144450353105449,0.00257445316700965,&
    -0.00104247508055041,-0.00167946127380679,0.000577428923826108,&
    0.00114016779684690,-2.63209684213992e-05,-0.000753899380930065,&
    -0.000294894042986087/)

    !Insect%xc(wingID) = 0.0 ! original mesh
    Insect%xc(wingID) = 0.0473 ! shifted towards t.e. to 1/4 of the root chord ("+" sign here)
    !Insect%xc(wingID) = -0.0728 ! shifted towards l.e., to 0.2cmean from the l.e. (Liu and Aono BB 2009)
    !Insect%yc(wingID) = 0.7
    !Insect%yc(wingID) = 0.712 ! measured using kinematics snapshots
    Insect%yc(wingID) = 0.702 ! According to Maeda's email, Jun 21, 2014
  case ('drosophila_sun')
    !********************************************
    !  Drosophila virilis wing from Chen and Sun, Acta Mech Sin 2014
    !********************************************
    Insect%nfft_wings(wingID) = 25
    Insect%a0_wings(wingID) = 0.5427789795180327
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.10879717599747996, -0.11445383382313232, -0.02023898255134523,&
    0.04903268079573884, 0.0012813019346402, -0.02397317767942499,&
    0.0013575396713610029, 0.0108149787395804, -0.001514114743464855,&
    -0.005364275911068656, 3.6505751634048205E-4, 0.002640180169907162,&
    -3.2673259786225535E-4, -0.0014323857426683313, 2.431115176929324E-4,&
    5.392319229992534E-4, -4.5833334881866856E-4, -1.3216432233072333E-4,&
    6.563502263270568E-4, 2.0750829321817808E-4, -4.807960800434886E-4,&
    -2.9006261005712504E-4, 2.7578746591965946E-4, 2.7519915193569374E-4,&
    -3.0570604954113513E-4/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/-0.09385487982296374, -0.010821846776797858, 0.030052970821579587,&
    0.005312859230387492, -0.006054695188204192, -0.0015421303479118118,&
    -0.002533264559802815, -0.0014806147599133366, 0.003640199794653037,&
    0.0020416212413267134, -0.0024948946435721206, -7.83017244422372E-4,&
    0.0021574389122894035, 2.6950667683726845E-4, -0.00131044444112179,&
    6.404390762251693E-5, 2.513250728448789E-4, -4.7634735716375334E-4,&
    -1.5949800516545527E-5, 5.001276053841919E-4, 8.445613796483002E-5,&
    -5.510759077970704E-4, -3.3722093938416713E-4, 3.524656540450335E-4,&
    2.9924999100355387E-4/)
    Insect%xc(wingID) = 0.0
    Insect%yc(wingID) = 0.399446382250523
  case ('bumblebee')
    !********************************************
    !  Bumblebee
    !  http://www.entomology.umn.edu/museum/links/coursefiles/JPEG%20images/Hymenoptera%20web%20jpeg/Bombus-wings.jpg
    !********************************************
    Insect%nfft_wings(wingID) = 25
    Insect%a0_wings(wingID) = 0.594557593733011d0
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/-0.0128037920989526d0,-0.106777418654552d0,0.0380851289321982d0,&
    0.0330548081081197d0,-0.0178496286355627d0,-0.00328588543359649d0,&
    0.0108246924336137d0,-0.00489302388329943d0,-0.00708808441172961d0,&
    0.00518244772891516d0,0.00445979960844562d0,-0.000108072056165527d0,&
    0.00204437845603716d0,0.00147176382618797d0,-0.00229559463098105d0,&
    -0.000514633972391526d0,0.00134150515430486d0,-0.000149860228261824d0,&
    9.01456938813568d-05,0.00150639712261487d0,0.000914624010720407d0,&
    -0.000737650894315551d0,-0.000843849923321745d0,-0.000354971670482499d0,&
    -0.000382956472432449d0/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/-0.0158061138788171d0,0.0308243584184200d0,-0.00903330410923372d0,&
    -0.0185758334697500d0,-0.000924452934252486d0,-0.00242101213359519d0,&
    -0.00204549530064489d0,0.00291468131401423d0,-0.000140755032337495d0,&
    -0.00135036427128534d0,0.00141285439042451d0,-0.000334215276598231d0,&
    -0.00161521722061879d0,-0.000164055684312904d0,-0.000256278551727569d0,&
    -0.000740258481681094d0,0.000847498161852221d0,0.00157442110960973d0,&
    -0.000559835622451578d0,-0.000617498559228280d0,0.00115413452523474d0,&
    0.000322564770099778d0,-0.000917375185844477d0,4.44819399488798d-05,&
    0.000710028654602170d0/)
    Insect%xc(wingID) = -0.1d0
    Insect%yc(wingID) = 0.501549263807117d0

  case ('b_ignitus')
    !********************************************
    !  Bumblebee B. ignitus
    !  Digitized from images taken at Liu Lab
    !********************************************
    Insect%nfft_wings(wingID) = 25
    Insect%a0_wings(wingID) = 0.536472532931637d0
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/-0.0447167394708177d0,-0.106357727795917d0,0.0504418160417239d0,&
    0.0217275689429364d0,-0.0259085955164794d0,0.00272535910748833d0,&
    0.00925289824790763d0,-0.00453010629382665d0,-0.000726647749565597d0,&
    0.00258280099999843d0,-0.00193033529765617d0,-0.00121090519402499d0,&
    0.00149872968653121d0,0.000716207684720514d0,-0.000205317764190544d0,&
    0.000120507537444963d0,-0.000381477942805165d0,-0.000364957961985063d0,&
    -6.70598716926467d-05,0.000166365788794039d0,0.000332993591840758d0,&
    -0.000225912231239784d0,-0.000554023819155716d0,0.000352735383706648d0,&
    0.000650085631908143d0/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/-0.0580660125663764d0,0.0271775529659247d0,0.0178916506228727d0,&
    -0.0196983386855655d0,-0.00865040473524334d0,0.0112078637630294d0,&
    0.00505882127179290d0,-0.00516874871678530d0,-0.000418585234573997d0,&
    0.00248996756589669d0,-0.00248081765717699d0,-0.00165307115885468d0,&
    0.00236884835642553d0,0.000920860396041608d0,-0.00160449459432319d0,&
    7.96078949775159d-05,0.000716588388745441d0,0.000306756717543478d0,&
    0.000310638954298390d0,-0.000523512353114016d0,-0.000773372382092419d0,&
    1.97258594500968d-05,0.000261943571939630d0,0.000262003935722642d0,&
    0.000278542046262820d0/)
    Insect%xc(wingID) = -0.13d0
    Insect%yc(wingID) = 0.434820393790595d0

  case ('paratuposa_flatwing')
    !********************************************
    !  Paratuposa wing simplified as a flat plate
    !  Digitized from 3D model, Moscow University entomology lab
    !********************************************
    Insect%nfft_wings(wingID) = 15
    Insect%a0_wings(wingID) = 0.694542662069373d0
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/-0.134736163793655d0,0.00530251847896251d0,-0.0345113221312334d0,&
    0.00564308389276391d0,-0.0151286715792430d0,0.00702004741152472d0,&
    0.00144649560886655d0,-0.00185566405410384d0,-0.00275905041561011d0,&
    0.00217239130911607d0,-0.00106500370428430d0,-0.000750733476326611d0,&
    0.00280149738648434d0,-0.00182306390466332d0,-0.000432849087666278d0/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.0931417992571149d0,0.0598684008118805d0,-0.0297533040215975d0,&
    0.00720759394160553d0,-0.00591704536021243d0,-0.0103773128578859d0,&
    0.00480389428474622d0,-0.000497969127742664d0,-0.000267837077569406d0,&
    -0.00102010912896721d0,0.00218717420988224d0,-0.00290492852428728d0,&
    0.000860757518054172d0,0.00147755983849289d0,-0.000638118479966807d0/)
    Insect%xc(wingID) = -0.3d0
    Insect%yc(wingID) = 0.7d0

  case ('paratuposa_flatelytra')
    !********************************************
    !  Paratuposa elytra simplified as a flat plate
    !  Digitized from 3D model, Moscow University entomology lab
    !********************************************
    Insect%nfft_wings(wingID) = 21
    Insect%a0_wings(wingID) = 0.565146320110115d0 * 0.62d0
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.0977588377084158d0,-0.134895186363991d0,-0.0415626398480136d0,&
    0.0462811349514129d0,0.0116178197112501d0,-0.0173654394721538d0,&
    -0.00143307308642128d0,0.00463364466230196d0,-0.00245823642621685d0,&
    -0.000295909979786882d0,0.00293439703102658d0,-0.000757888513431659d0,&
    -0.00232308360478113d0,0.000590866168893430d0,0.00134477744241136d0,&
    -0.000179797491049942d0,-0.000579366975013497d0,0.000298632558506859d0,&
    0.000121897474918553d0,-0.000129056828603867d0,6.65231631883479d-5/) &
    * 0.62d0
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.0873846736908491d0,0.0261657542571480d0,-0.0361196299337545d0,&
    -0.0169805639649468d0,0.0120703518387247d0,0.00786957912240014d0,&
    -0.00255959633597845d0,-0.00262639690803298d0,-0.000425744504894171d0,&
    -0.000897607832134308d0,0.00110940805947636d0,0.00119153774474213d0,&
    -0.00114698647381269d0,-0.000899734725370745d0,0.000630784528183092d0,&
    0.000278761021018954d0,-0.000290963412817785d0,6.98199270677237d-5,&
    0.000195892295619984d0,-0.000190195259814270d0,-4.82888469073568d-5/) &
    * 0.62d0
    Insect%xc(wingID) = -0.02d0 * 0.62d0
    Insect%yc(wingID) = 0.65d0 * 0.62d0

  case ('flapper_sane')
    !********************************************
    !  Mechanical model from Sane and Dickinson, JEB 205, 2002
    !  'The aerodynamic effects...'
    !********************************************
    Insect%nfft_wings(wingID) = 25
    Insect%a0_wings(wingID) = 0.5379588906565078
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/0.135338653455782,-0.06793162622123261,-0.0398235167675977,&
    0.006442194893963269,0.0012783260416583853,-0.007014398516674715,&
    0.0017710765408983137,0.006401601802033519,-2.970619204124993E-4,&
    -0.0038483478773981405,-6.180958756568494E-4,8.015784831786756E-4,&
    -6.957513357109226E-4,-1.4028929172227943E-4,0.0013484885717868547,&
    4.827827498543977E-4,-9.747844462919694E-4,-5.838504331939134E-4,&
    2.72834004831554E-4,2.8152492682871664E-5,-1.2802199282558645E-4,&
    4.117887216124469E-4,3.364169982438278E-4,-3.33258003686823E-4,&
    -3.5615733035757616E-4/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/2.686408368800394E-4,0.01649582345310688,0.01288513083639708,&
    0.004711436946785864,-0.0035725088809005073,-0.00898640397179334,&
    -0.003856509905612652,0.004536524572892801,0.004849677692836578,&
    2.9194421255236984E-4,-7.512780802871473E-4,7.12685261783966E-4,&
    -1.5519932673320404E-4,-0.0012695469974603026,2.2861692091158138E-4,&
    0.0016461316319681953,5.257476721137781E-4,-7.686482830046961E-4,&
    -3.108879176661735E-4,2.2437540206568518E-4,-2.578427217327782E-4,&
    -2.5120263516966855E-4,4.1693453021778877E-4,3.9290173948150096E-4,&
    -1.9762601237675826E-4/)
    Insect%xc(wingID) = 0.0
    Insect%yc(wingID) = 0.6
  case ('flapper_dickinsonII')
    !********************************************
    ! Digitized from Dickinson et al 1999 Science, figure 1A, drawing
    ! of the mechanical robot
    !********************************************
    Insect%nfft_wings(wingID) = 20
    Insect%a0_wings(wingID) = 0.6442788
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/0.0482978,-0.1208378,0.0061008,0.0356718,-0.0148328,-0.0109958,&
    0.0110268,0.0018538,-0.0061998,0.0015458,0.0025508,-0.0017538,&
    -0.0002578,0.0015018,-0.0003158,-0.0006048,0.0007168,-0.0001568,&
    -0.0005018,0.0004118/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/-0.0521708,0.0051828,0.0369428,-0.0002868,-0.0177448,0.0023218,&
    0.0081378,-0.0036288,-0.0038168,0.0031348,0.0011858,-0.0023828,&
    -0.0001638,0.0016098,-0.0004768,-0.0007188,0.0007228,0.0002278,&
    -0.0005798,0.0001228/)
    Insect%yc(wingID) = 0.5282438
    Insect%xc(wingID) = -0.1184548

  case ('robofly_dickinson')
    !********************************************
    ! Digitized from the hand drawn figure M. Dickinson sent via email, which
    ! contained the exact location of the pivot point. He also sent a CAD drawing
    ! which looks slightly different, and had no pivot point marked.
    !********************************************
    Insect%nfft_wings(wingID) = 28
    Insect%a0_wings(wingID) = 0.5313628
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/-0.0245658,-0.0842918,0.0218028,0.0105418,-0.0095288,0.0012928,&
    0.0021928,0.0000328,-0.0007648,-0.0015808,0.0013808,0.0013068,&
    -0.0010748,0.0002408,-0.0000378,-0.0010888,0.0008248,0.0004708,&
    -0.0003988,0.0002658,-0.0003178,-0.0004218,0.0002768,0.0000818,&
    0.0000318,0.0001228,-0.0001918,-0.0000558/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
     (/-0.0905448,0.0278058,0.0392558,-0.0125248,-0.0159598,0.0048268,&
    0.0038898,-0.0028828,0.0012618,0.0012998,-0.0019058,0.0003118,&
    0.0003198,-0.0004298,0.0006388,-0.0000648,-0.0002308,0.0002518,&
    -0.0003948,0.0000928,0.0004478,-0.0003078,-0.0000888,0.0001638,&
    -0.0002348,0.0001398,0.0001398,-0.0002358/)
    Insect%yc(wingID) = 0.4645238
    Insect%xc(wingID) = -0.0716018

  case ('hawkmoth1')
    ! this wingshape is digitized from figure 1 from Kim et al. "Hovering and forward flight of the hawkmoth
    ! M. sexta: trim search and 6DOF dynamic stability characterization (Bioinspir Biomim. 10 (2015) 056012)"
    ! its area is about 0.30, which is lower than most references found. in their paper, they state A=0.3788,
    ! which raises question if fig 1 is to-scale drawing or sketch.
    Insect%nfft_wings(wingID) = 28
    Insect%a0_wings(wingID) = 0.5860758
    Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.0219308,-0.1252418,0.0154668,0.0356038,-0.0203008,-0.0061968,&
    0.0178288,0.0002728,-0.0089908,0.0022758,0.0022948,-0.0046148,&
    -0.0008808,0.0032598,-0.0003708,-0.0019528,0.0006858,0.0008268,&
    -0.0008358,-0.0000718,0.0008938,0.0000348,-0.0004598,0.0004428,&
    0.0003158,-0.0003108,-0.0000658,0.0002798/)
    Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
    (/0.0062418,0.0452798,0.0303808,-0.0184998,-0.0179088,0.0068018,&
    0.0030268,-0.0058108,0.0017748,0.0040588,-0.0033678,-0.0030638,&
    0.0021898,0.0007208,-0.0015538,0.0007128,0.0016948,-0.0003828,&
    -0.0005898,0.0006388,0.0002888,-0.0005258,0.0000808,0.0002248,&
    -0.0004308,-0.0002758,0.0002298,-0.0000548/)
    Insect%yc(wingID) = 0.4171918
    Insect%xc(wingID) = -0.0395258
    case ('hawkmoth2')
        ! this wingshape is digitized from https://en.wikipedia.org/wiki/Manduca_sexta#/media/File:Manduca_sexta_female_sjh.JPG
        ! it has a greater aerea (A=0.40), but the original image is tricky since it is rotated. we therefore used a bit of modeling
        ! for this wing shape.
        Insect%nfft_wings(wingID) = 18
        Insect%a0_wings(wingID) = 0.6617728
        Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = &
        (/-0.0837648,-0.0802108,0.0703808,0.0069808,-0.0183478,0.0156518,&
        -0.0000308,-0.0153718,-0.0011538,0.0032378,-0.0005008,0.0020798,&
        0.0019888,-0.0009568,-0.0016378,-0.0010208,0.0005658,0.0009028/)
        Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = &
        (/0.0086968,0.0763208,0.0216658,-0.0322558,-0.0125988,0.0042128,&
        -0.0066278,-0.0040288,0.0093858,0.0045358,-0.0043238,0.0006298,&
        0.0010848,-0.0028958,0.0007268,0.0022578,-0.0013068,-0.0003538/)
        Insect%yc(wingID) = 0.3946798
        Insect%xc(wingID) = -0.2157968

    case default

        ! if all other options fail, we still might load coefficients from file:
        wingshape_str = Insect%WingShape(wingID)
        if (index(wingshape_str,"from_file::") /= 0) then
            !-------------------------------------------------------------------------
            ! wing shape is read from ini-file
            !-------------------------------------------------------------------------
            call Setup_Wing_from_inifile(Insect, wingID, trim(adjustl(wingshape_str( 12:len_trim(wingshape_str) ))))

        else
            ! now we theres an error...
            write (*,*) "Insect module: trying to set up fourier descriptors for wing&
            & shape but the type Insect%WingShape is unknown! :: "// Insect%WingShape

            call abort(554329, "Insect module: trying to set up fourier descriptors for wing&
            & shape but the type Insect%WingShape is unknown!")
        end if
    end select

  ! skip this routine in the future and let other routines know that the wing
  ! shape is ready to be used.
  Insect%wingsetup_done(wingID) = .true.

  ! for many cases, it is important that Lspan and Lchord are known, but that is
  ! tedious for Fourier shapes, as the use cannot see it from the cooefficients.
  ! Therefore, we compute the max / min of x / y here and store the result
  call set_wing_bounding_box_fourier( Insect, wingID )

  ! this is the old defaut value:
  if (maxval(Insect%corrugation_array_bbox) == 0.0_rk) then
      Insect%corrugation_array_bbox(1:4) = Insect%wing_bounding_box(1:4,wingID)
  endif

  if (root) then
    write(*,'(30("-"))')
    write(*,'("Insect module: Setup_Wing_Fourier_coefficients")')
    write(*,'("Wing shape is ",A)') trim(adjustl(Insect%WingShape(wingID)))
    write(*,'("nfft_wings=",i3)') Insect%nfft_wings(wingID)
    write(*,'(30("-"))')
  endif

end subroutine Setup_Wing_Fourier_coefficients



!-------------------------------------------------------------------------------
! Initialize the wing from an ini file
!-------------------------------------------------------------------------------
! We read:
!     - Fourier coefficients for the radius wing shape
!     - Wing thickness profile (constant or variable)
!     - Wing corrugation profile (flat or corrugated)
!-------------------------------------------------------------------------------
subroutine Setup_Wing_from_inifile( Insect, wingID, fname )
  implicit none
  type(diptera),intent(inout) :: Insect
  character(len=*), intent(in) :: fname

  type(inifile) :: ifile
  real(kind=rk), allocatable :: tmparray(:,:)
  character(len=strlen) :: type_str
  integer :: a,b
  integer(kind=2), intent(in) :: wingID ! wing id number
  real(kind=rk) :: init_thickness

  if (root) then
    write(*,'(80("-"))')
    write(*,'("Reading wing shape from file ",A)') fname
    write(*,'(80("-"))')
  endif
  ! instead of the hard-coded values above, read fourier coefficients for wings from
  ! an ini-file
  call read_ini_file_mpi(ifile, fname, .true.  )


  ! check if this file seem to be valid:
  call read_param_mpi(ifile, "Wing", "type", Insect%wing_file_type(wingID), "none")

  if (Insect%wing_file_type(wingID) /= "fourier" .and. Insect%wing_file_type(wingID) /= "fourierY" .and. Insect%wing_file_type(wingID) /= "kleemeier" ) then
    call abort(6652, "ini file for wing does not seem to be correct type...")
  endif

  if (Insect%wing_file_type(wingID) == "kleemeier" ) then
      call read_param_mpi(ifile, "Wing", "B_membrane", Insect%B_membrane(wingID), 8.6_rk/130.0_rk)
      call read_param_mpi(ifile, "Wing", "L_membrane", Insect%L_membrane(wingID), 100.0_rk/130.0_rk)
  endif
  !-----------------------------------------------------------------------------
  ! Read fourier coeffs for wing radius
  !-----------------------------------------------------------------------------
  if (Insect%wing_file_type(wingID) == "fourier" .or. Insect%wing_file_type(wingID) == "fourierY") then
      call read_param_mpi( ifile, "Wing", "a0_wings", Insect%a0_wings(wingID), 0.d0)

      ! NOTE: Annoyingly, the fujitsu SXF90 compiler cannot handle allocatable arrays
      ! as arguments. so we have to split the routine in one part that returns the size
      ! of the array, then let the caller allocate, then read the matrix. very tedious.
      ! fetch size of matrix
      call param_matrix_size_mpi( ifile, "Wing", "ai_wings", a, b)
      ! allocate matrix
      allocate( tmparray(1:a,1:b) )
      ! read matrix
      call param_matrix_read_mpi( ifile, "Wing", "ai_wings", tmparray)
      Insect%nfft_wings = size(tmparray,2)
      Insect%ai_wings(1:Insect%nfft_wings(wingID),wingID) = tmparray(1,:)
      deallocate(tmparray)


      call param_matrix_size_mpi( ifile, "Wing", "bi_wings", a, b)
      ! allocate matrix
      allocate( tmparray(1:a,1:b) )
      ! read matrix
      call param_matrix_read_mpi( ifile, "Wing", "bi_wings", tmparray)
      Insect%nfft_wings = size(tmparray,2)
      Insect%bi_wings(1:Insect%nfft_wings(wingID),wingID) = tmparray(1,:)
      deallocate(tmparray)

      ! wing mid-point (of course in wing system..)
      call read_param_mpi(ifile,"Wing","x0w",Insect%xc(wingID), 0.d0)
      call read_param_mpi(ifile,"Wing","y0w",Insect%yc(wingID), 0.d0)
  endif

  !-----------------------------------------------------------------------------
  ! wing thickness
  !-----------------------------------------------------------------------------
  call read_param_mpi(ifile,"Wing","wing_thickness_distribution",Insect%wing_thickness_distribution(wingID), "constant")
  if ( Insect%wing_thickness_distribution(wingID) == "constant") then

      if (root) write(*,*) "Wing thickness is constant along the wing"

      ! wing thickness (NOTE: overwrites settings in other params file)
      if ( (Insect%WingThickness>0.0d0) .and. (Insect%WingThickness<1.0d0) ) then
         init_thickness = Insect%WingThickness ! Use existing value if it is reasonable
      else
         init_thickness = 0.05d0 ! This is the defauls value otherwise, because we may not know dx here
      endif
      
      call read_param_mpi(ifile,"Wing","wing_thickness_value",Insect%WingThickness, init_thickness)

  elseif ( Insect%wing_thickness_distribution(wingID) == "variable") then

      if (root) write(*,*) "Wing thickness is variable, i.e. t = t(x,y)"

      ! read matrix from ini file, see comments on SXF90 compiler
      call param_matrix_size_mpi(ifile,"Wing","wing_thickness_profile",a,b)
      call Allocate_Arrays(Insect,"wing_thickness_profile",a,b)
      call param_matrix_read_mpi(ifile,"Wing","wing_thickness_profile",wing_thickness_profile)

  else
      call abort(77623, " Insect wing thickness distribution is unknown (must be constant or variable)")
  endif


  !-----------------------------------------------------------------------------
  ! bristles
  !-----------------------------------------------------------------------------
  call read_param_mpi(ifile, "Wing","bristles", Insect%bristles(wingID), .false.)
  call read_param_mpi(ifile, "Wing","bristles_simplex", Insect%bristles_simplex(wingID), .false.)
  if (Insect%bristles(wingID)) then
      call param_matrix_size_mpi( ifile, "Wing", "bristles_coords", a, b)

      ! number of bristles on this wing
      Insect%n_bristles(wingID) = a

      ! check memory allocation
      if (.not. allocated(Insect%bristles_coords)) then
          allocate(Insect%bristles_coords(1:4, 1:a, 1:b)) ! four wings...
      elseif ( size(Insect%bristles_coords, 2) .ne. a ) then
          call abort(76237, " Unequal number of bristles not yet supported. Modify wing shape files.")
      endif

      call param_matrix_read_mpi( ifile, "Wing", "bristles_coords", Insect%bristles_coords(wingID,:,:))
  endif


  !-----------------------------------------------------------------------------
  ! wing corrugation
  !-----------------------------------------------------------------------------
  call read_param_mpi(ifile,"Wing","corrugated",Insect%corrugated(wingID), .false.)
  if (Insect%corrugated(wingID)) then
      if (root) write(*,*) "wing is corrugated, z=z(x,y)"
      ! read matrix from ini file, see comments on SXF90 compiler
      call param_matrix_size_mpi(ifile,"Wing","corrugation_profile",a,b)
      call Allocate_Arrays(Insect,"corrugation_profile",a,b)
      call param_matrix_read_mpi(ifile,"Wing","corrugation_profile",corrugation_profile)

  else
      if (root) write(*,*) "wing is flat (non-corrugated), z==0"
  endif

  call read_param_mpi(ifile,"Wing","corrugation_array_bbox",Insect%corrugation_array_bbox, (/0.0_rk,0.0_rk,0.0_rk,0.0_rk/))

end subroutine Setup_Wing_from_inifile


!-------------------------------------------------------------------------------
! Setup extends of wing for reduction of computational time
!-------------------------------------------------------------------------------
! for many cases, it is important that Lspan and Lchord are known, but that is
! tedious for Fourier shapes, as the use cannot see it from the cooefficients.
! Therefore, we compute the max / min of x / y / z  here and store the result
! NOTE: This code is executed only once.
!-------------------------------------------------------------------------------
subroutine set_wing_bounding_box_fourier( Insect, wingID )
  implicit none
  type(diptera),intent(inout) :: Insect
  real(kind=rk) :: theta, xmin,xmax, ymin, ymax, R, x, y, theta_prime, tmp
  integer(kind=2), intent(in) :: wingID ! wing id number

  theta = 0.d0
  xmin = 999.d9
  ymin = 999.d9
  xmax = -999.d9
  ymax = -999.d9

  ! construct the wing border by looping over the angle theta, look for smallest and largest x,y values
  ! note flusi uses a normalized angle between [0,1)
  do while ( theta < 1.d0 )
    ! note this is normalized angle
    R = Radius_Fourier( theta, Insect, wingID )

    theta_prime = 2.d0*pi*theta - pi
    x = Insect%xc(wingID) + R * cos( theta_prime )
    y = Insect%yc(wingID) + R * sin( theta_prime )

    ! NOTE: A word on theta_prime: it is the angle described with the positve x-axis
    ! and indeed rotates in positve z-direction. That means in the plane
    !
    !  ^ x_wing
    !  |
    !  o-------> y_wing
    !
    ! It rotates CLOCKWISE, starting from the x_wing axis (zero is the x-axis, vertical)

    xmin = min( xmin, x )
    ymin = min( ymin, y )

    xmax = max( xmax, x )
    ymax = max( ymax, y )

    theta = theta + 1.0d-3
  end do

  Insect%wing_bounding_box(1:4,wingID) = (/xmin, xmax, ymin, ymax/)

  ! the bounding box in z-direction depends on the wing thicnkess (constant or not)
  ! and the corrugation
  if ( Insect%wing_thickness_distribution(wingID) == "constant" ) then
    if ( Insect%corrugated(wingID) ) then
      Insect%wing_bounding_box(5,wingID) = minval(corrugation_profile) - Insect%WingThickness / 2.0_pr
      Insect%wing_bounding_box(6,wingID) = maxval(corrugation_profile) + Insect%WingThickness / 2.0_pr
    else
      ! constant thickness, no corrugation is the classical flat case:
      Insect%wing_bounding_box(5,wingID) = -Insect%WingThickness / 2.0_pr
      Insect%wing_bounding_box(6,wingID) = +Insect%WingThickness / 2.0_pr
    endif
  else
    if ( Insect%corrugated(wingID) ) then
      ! minimum of lower surface
      Insect%wing_bounding_box(5,wingID) = minval(corrugation_profile-wing_thickness_profile/2.0_pr)
      ! maximum of upper surface
      Insect%wing_bounding_box(6,wingID) = maxval(corrugation_profile+wing_thickness_profile/2.0_pr)
    else
      ! bounding box is +- largest thickness  simply
      Insect%wing_bounding_box(5,wingID) = -maxval(wing_thickness_profile / 2.0_pr)
      Insect%wing_bounding_box(6,wingID) =  maxval(wing_thickness_profile / 2.0_pr)
    endif
  end if

  ! Mirror the bounding box for left wings
  if ( (wingID == 1) .or. (wingID == 3) ) then
    tmp = Insect%wing_bounding_box(6,wingID)
    Insect%wing_bounding_box(6,wingID) = -Insect%wing_bounding_box(5,wingID)
    Insect%wing_bounding_box(5,wingID) = -tmp
  end if

  if (root) then
    write(*,'("Effective (=the real surface) wing lengths are:")')
    write(*,'("Lspan=",es15.8,"Lchord=",es15.8)') ymax-ymin, xmax-xmin
    write(*,'("Bounding box is:")')
    write(*,'("xwmin=",es15.8," xwmax=",es15.8)') Insect%wing_bounding_box(1:2,wingID)
    write(*,'("ywmin=",es15.8," ywmax=",es15.8)') Insect%wing_bounding_box(3:4,wingID)
    write(*,'("zwmin=",es15.8," zwmax=",es15.8)') Insect%wing_bounding_box(5:6,wingID)
  endif
end subroutine set_wing_bounding_box_fourier




subroutine draw_bristle(x1w, x2w, R0, xx0, ddx, mask, mask_color, us, Insect, color_val, M_body, M_wing, x_pivot_b, rot_rel_wing_w)
    implicit none

    real(kind=rk), dimension(1:3), intent(in):: x1w, x2w
    real(kind=rk),intent(in)::R0
    type(diptera),intent(inout)::Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: color_val
    real(kind=rk),intent(in) :: M_body(1:3,1:3), M_wing(1:3,1:3), x_pivot_b(1:3), rot_rel_wing_w(1:3)

    real(kind=rk),dimension(1:3) ::  cb, rb, ab, u, vp, x1, x2, x_wing, x_body, v_tmp, x
    real(kind=rk),dimension(1:3) :: x_glob, e_x, tmp, e_r, e_3
    real(kind=rk) :: ceta1, ceta2, ceta3, R, RR0, clength, safety, t
    integer :: ix,iy,iz

    integer, dimension(1:3) :: lbounds, ubounds
    integer :: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: Nsafety

    safety = Insect%safety
    Nsafety = nint(safety / minval(ddx))

    ! bounds of the current patch of data
    lbounds = g
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) - 1 - g

    RR0 = R0 + safety


    !---------------------------------------------------------------------------
    ! transform coordinates to global system. they are defined in the wing system
    !---------------------------------------------------------------------------
    x1 = matmul( transpose(M_wing), x1w) + x_pivot_b
    x1 = matmul( transpose(M_body), x1) + Insect%xc_body_g

    x2 = matmul( transpose(M_wing), x2w) + x_pivot_b
    x2 = matmul( transpose(M_body), x2) + Insect%xc_body_g


    !---------------------------------------------------------------------------
    ! define bristle (cylinder) coordinate system
    !---------------------------------------------------------------------------
    ! unit vector in bristle (cylinder) axis direction and bristle (cylinder) length
    e_x = x2 - x1
    clength = norm2(e_x)
    e_x = e_x / clength

    ! radial unit vector
    ! use a vector perpendicular to e_x, since it is a azimuthal symmetry
    ! it does not really matter which one. however, we must be sure that the vector
    ! we use and the e_x vector are not colinear -- their cross product is the zero vector, if that is the case
    e_r = (/0.d0, 0.d0, 0.d0/)
    do while ( norm2(e_r) <= 1.0d-12 )
        e_r = cross( (/rand_nbr(),rand_nbr(),rand_nbr()/), e_x)
    enddo
    e_r = e_r / norm2(e_r)

    ! third (also radial) unit vector, simply the cross product of the others
    e_3 = cross(e_x,e_r)
    e_3 = e_3 / norm2(e_3)


    !---------------------------------------------------------------------------
    ! bounding box of the vicinity of the bristle (cylinder).
    !---------------------------------------------------------------------------
    t = minval( (/x1(1)+RR0*e_r(1), x1(1)-RR0*e_r(1), x1(1)+RR0*e_3(1), x1(1)-RR0*e_3(1), &
                  x2(1)+RR0*e_r(1), x2(1)-RR0*e_r(1), x2(1)+RR0*e_3(1), x2(1)-RR0*e_3(1) /) )
    xmin = nint( (t-xx0(1)) / ddx(1) ) - Nsafety

    t = maxval( (/x1(1)+RR0*e_r(1), x1(1)-RR0*e_r(1), x1(1)+RR0*e_3(1), x1(1)-RR0*e_3(1), &
                  x2(1)+RR0*e_r(1), x2(1)-RR0*e_r(1), x2(1)+RR0*e_3(1), x2(1)-RR0*e_3(1) /) )
    xmax = nint( (t-xx0(1)) / ddx(1) ) + Nsafety

    t = minval( (/x1(2)+RR0*e_r(2), x1(2)-RR0*e_r(2), x1(2)+RR0*e_3(2), x1(2)-RR0*e_3(2), &
                  x2(2)+RR0*e_r(2), x2(2)-RR0*e_r(2), x2(2)+RR0*e_3(2), x2(2)-RR0*e_3(2) /) )
    ymin = nint( (t-xx0(2)) / ddx(2) ) - Nsafety

    t = maxval( (/x1(2)+RR0*e_r(2), x1(2)-RR0*e_r(2), x1(2)+RR0*e_3(2), x1(2)-RR0*e_3(2), &
                  x2(2)+RR0*e_r(2), x2(2)-RR0*e_r(2), x2(2)+RR0*e_3(2), x2(2)-RR0*e_3(2) /) )
    ymax = nint( (t-xx0(2)) / ddx(2) ) + Nsafety

    t = minval( (/x1(3)+RR0*e_r(3), x1(3)-RR0*e_r(3), x1(3)+RR0*e_3(3), x1(3)-RR0*e_3(3), &
                  x2(3)+RR0*e_r(3), x2(3)-RR0*e_r(3), x2(3)+RR0*e_3(3), x2(3)-RR0*e_3(3) /) )
    zmin = nint( (t-xx0(3)) / ddx(3) ) - Nsafety

    t = maxval( (/x1(3)+RR0*e_r(3), x1(3)-RR0*e_r(3), x1(3)+RR0*e_3(3), x1(3)-RR0*e_3(3), &
                  x2(3)+RR0*e_r(3), x2(3)-RR0*e_r(3), x2(3)+RR0*e_3(3), x2(3)-RR0*e_3(3) /) )
    zmax = nint( (t-xx0(3)) / ddx(3) ) + Nsafety


    ! first we draw the cylinder, then the endpoint spheres
    do iz = max(zmin, lbounds(3)), min(zmax, ubounds(3))
        x_glob(3) = xx0(3) + dble(iz)*ddx(3)

        do iy = max(ymin, lbounds(2)), min(ymax, ubounds(2))
            x_glob(2) = xx0(2) + dble(iy)*ddx(2)

            do ix = max(xmin, lbounds(1)), min(xmax, ubounds(1))
                x_glob(1) = xx0(1) + dble(ix)*ddx(1)
                ! if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! cb is the distance to the cylinder mid-point
                cb = 0.5d0*(x1+x2) - x_glob
                ! rb is the length of the clinder
                rb = x1 - x2

                ! this is a spherical bounding box, centered around the mid-point
                if ( sum(cb**2) < 0.25*sum(rb**2) ) then ! the 0.25 is from the 0.5 squared
                    ab = x_glob - x1
                    u = x2 - x1

                    vp = cross(ab, u)
                    R = sqrt( sum(vp**2) / sum(u**2) )

                    if (R <= R0+safety) then
                        t = steps(R, R0, Insect%smooth)
                        if (t >= mask(ix,iy,iz)) then

                            mask(ix,iy,iz) = t
                            mask_color(ix,iy,iz) = color_val

                            x_body = matmul(M_body, x_glob - Insect%xc_body_g)
                            x_wing = matmul(M_wing, x_body - x_pivot_b)

                            !---------------------------------------------------
                            ! solid body rotation
                            ! Attention: the Matrix transpose(M) brings us back to the body
                            ! coordinate system, not to the inertial frame. this is done in
                            ! the main routine Draw_Insect
                            !---------------------------------------------------
                            v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                            v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                            v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                            ! note we set this only if it is a part of the wing
                            ! note velocity is to be set in BODY coordinate system.
                            us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                        endif
                    endif
                endif
            enddo
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! endpoint sphere x2
    !---------------------------------------------------------------------------
    Nsafety = nint( (R0+Insect%safety) / minval(ddx))

    ! bounds of the current patch of data
    lbounds = g
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) - 1 - g

    ! bounding box of the vicinity of the sphere.
    xmin = nint( (x2(1)-xx0(1)) / ddx(1) ) - (Nsafety)
    xmax = nint( (x2(1)-xx0(1)) / ddx(1) ) + (Nsafety)
    ymin = nint( (x2(2)-xx0(2)) / ddx(2) ) - (Nsafety)
    ymax = nint( (x2(2)-xx0(2)) / ddx(2) ) + (Nsafety)
    zmin = nint( (x2(3)-xx0(3)) / ddx(3) ) - (Nsafety)
    zmax = nint( (x2(3)-xx0(3)) / ddx(3) ) + (Nsafety)


    do iz = max(zmin,lbounds(3)), min(zmax,ubounds(3))
        x(3) = xx0(3) + dble(iz)*ddx(3) - x2(3)

        do iy = max(ymin,lbounds(2)), min(ymax,ubounds(2))
            x(2) = xx0(2) + dble(iy)*ddx(2) - x2(2)

            do ix = max(xmin,lbounds(1)), min(xmax,ubounds(1))
                x(1) = xx0(1) + dble(ix)*ddx(1) - x2(1)

                if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                ! the bounding box check is incorporated in the loop bounds - no if clause!
                ! compute radius
                R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
                if ( R <= R0+Insect%safety ) then
                    t = steps(R, R0, Insect%smooth)
                    if ( t >= mask(ix,iy,iz) ) then
                        ! set new value
                        mask(ix,iy,iz) = t
                        mask_color(ix,iy,iz) = color_val

                        x_body = matmul(M_body, x + x2 - Insect%xc_body_g)
                        x_wing = matmul(M_wing, x_body - x_pivot_b)

                        !---------------------------------------------------
                        ! solid body rotation
                        ! Attention: the Matrix transpose(M) brings us back to the body
                        ! coordinate system, not to the inertial frame. this is done in
                        ! the main routine Draw_Insect
                        !---------------------------------------------------
                        v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                        v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                        v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                        ! note we set this only if it is a part of the wing
                        ! note velocity is to be set in BODY coordinate system.
                        us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                    endif
                endif
            enddo
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! endpoint sphere x1
    !---------------------------------------------------------------------------
    Nsafety = nint( (R0+Insect%safety) / minval(ddx))

    ! bounds of the current patch of data
    lbounds = g
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) - 1 - g

    ! bounding box of the vicinity of the sphere.
    xmin = nint( (x1(1)-xx0(1)) / ddx(1) ) - (Nsafety)
    xmax = nint( (x1(1)-xx0(1)) / ddx(1) ) + (Nsafety)
    ymin = nint( (x1(2)-xx0(2)) / ddx(2) ) - (Nsafety)
    ymax = nint( (x1(2)-xx0(2)) / ddx(2) ) + (Nsafety)
    zmin = nint( (x1(3)-xx0(3)) / ddx(3) ) - (Nsafety)
    zmax = nint( (x1(3)-xx0(3)) / ddx(3) ) + (Nsafety)


    do iz = max(zmin,lbounds(3)), min(zmax,ubounds(3))
        x(3) = xx0(3) + dble(iz)*ddx(3) - x1(3)

        do iy = max(ymin,lbounds(2)), min(ymax,ubounds(2))
            x(2) = xx0(2) + dble(iy)*ddx(2) - x1(2)

            do ix = max(xmin,lbounds(1)), min(xmax,ubounds(1))
                x(1) = xx0(1) + dble(ix)*ddx(1) - x1(1)

                if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                ! the bounding box check is incorporated in the loop bounds - no if clause!
                ! compute radius
                R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
                if ( R <= R0+Insect%safety ) then
                    t = steps(R, R0, Insect%smooth)
                    if ( t >= mask(ix,iy,iz) ) then
                        ! set new value
                        mask(ix,iy,iz) = t
                        mask_color(ix,iy,iz) = color_val

                        x_body = matmul(M_body, x + x1 - Insect%xc_body_g)
                        x_wing = matmul(M_wing, x_body - x_pivot_b)

                        !---------------------------------------------------
                        ! solid body rotation
                        ! Attention: the Matrix transpose(M) brings us back to the body
                        ! coordinate system, not to the inertial frame. this is done in
                        ! the main routine Draw_Insect
                        !---------------------------------------------------
                        v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                        v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                        v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                        ! note we set this only if it is a part of the wing
                        ! note velocity is to be set in BODY coordinate system.
                        us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                    endif
                endif
            enddo
        enddo
    enddo

end subroutine


!-------------------------------------------------------------------------------
! Draw a triangle determined by points x1,y1 x2,y2 x3,y3 
!-------------------------------------------------------------------------------
subroutine draw_triangle(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_body,M_wing,x_pivot_b,rot_rel_wing_w, &
                        x1,y1,x2,y2,x3,y3)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing
  real(kind=rk),intent(in) :: M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)
  real(kind=rk),intent(in) :: x1,y1,x2,y2,x3,y3

  integer :: ix,iy,iz
  integer :: v1,v2,v3
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=rk) :: v_tmp(1:3),mask_tmp,theta,xc,yc,xt,yt,xmax,xmin,ymax,ymin

  !-- bounding box
  xmin = min(x1,min(x2,x3))
  xmax = max(x1,max(x2,x3))
  ymin = min(y1,min(y2,y3))
  ymax = max(y1,max(y2,y3))

  !-- centroid
  xc = (x1+x2+x3)/3
  yc = (y1+y2+y3)/3

  do iz = g, size(mask,3)-1-g
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

              x_body = matmul(M_body,x)
              x_wing = matmul(M_wing,x_body-x_pivot_b)

              !-- test point
              xt = x_wing(1)
              yt = x_wing(2)

              !-- if inside bounding box
              if ( (xt>=xmin-Insect%safety) .and. (xt<=xmax+Insect%safety) .and. &
                   (yt>=ymin-Insect%safety) .and. (yt<=ymax+Insect%safety) ) then

                  if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
                      !-- determine which side
                      v1 = f_same_side_point(x1,y1,x2,y2,xc,yc,xt,yt);
                      v2 = f_same_side_point(x2,y2,x3,y3,xc,yc,xt,yt);
                      v3 = f_same_side_point(x3,y3,x1,y1,xc,yc,xt,yt);

                      !-- if inside, draw
                      if ( (v1==1) .and. (v2==1) .and. (v3==1) ) then

                           !-- smooth height
                           mask_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness, Insect%smooth) ! thickness

                           if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0d0)) then
                               mask(ix,iy,iz) = mask_tmp
                               mask_color(ix,iy,iz) = color_wing
                               !------------------------------------------------
                               ! solid body rotation
                               ! Attention: the Matrix transpose(M) brings us back to the body
                               ! coordinate system, not to the inertial frame. this is done in
                               ! the main routine Draw_Insect
                               !------------------------------------------------
                               v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                               v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                               v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                               ! note we set this only if it is a part of the wing
                               us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
                           endif
                      endif
                  endif 
              endif
          enddo
      enddo
  enddo
end subroutine draw_triangle


!-------------------------------------------------------------------------------
! Given a line through two points x1,y1 and x2,y2 
! determine if points x3,y3 and x4,y4 are on the same side
! Return 1 if true
! Return 0 if false
!-------------------------------------------------------------------------------
function f_same_side_point(x1,y1,x2,y2,x3,y3,x4,y4)
    implicit none

    real(kind=rk), intent(in) :: x1,y1,x2,y2,x3,y3,x4,y4
    real(kind=rk) :: rr,dx21,m1,b1,b3,b4
    integer :: f_same_side_point

    ! Initialize the result
    rr = 0
    
    ! Denominator in the slope of the line
    dx21 = (x2-x1) 

    ! If the slope is small enough
    if (dabs(dx21) > 1.0d-5) then
        m1 = (y2-y1)/dx21
        b1 = y1-m1*x1
        b3 = y3-m1*x3
        b4 = y4-m1*x4
        if ( ( (b1>=b3) .and. (b1>=b4) ) .or. ( (b1<=b3) .and. (b1<=b4) ) ) then
            rr = 1
        endif
    ! If the slope is close to 90deg, swap x and y
    else
        m1 = dx21/(y2-y1)
        b1 = x1-m1*y1 
        b3 = x3-m1*y3
        b4 = x4-m1*y4
        if ( ( (b1>=b3) .and. (b1>=b4) ) .or. ( (b1<=b3) .and. (b1<=b4) ) ) then
            rr = 1
        endif
    endif
    f_same_side_point = rr
end function f_same_side_point
