
! the new routine (2/2019) creates the wings (if both wings are used, maybe just one is)
! and their solid body velocity field us. Note that us contains both contributions from
! body and wing motion.
subroutine draw_insect_wings(time, xx0, ddx, mask, mask_color, us, Insect)
  implicit none

  real(kind=rk), intent(in)    :: time
  type(diptera), intent(inout) :: Insect
  real(kind=rk), intent(in)    :: xx0(1:3), ddx(1:3)
  real(kind=rk), intent(inout) :: mask(0:,0:,0:)
  real(kind=rk), intent(inout) :: us(0:,0:,0:,1:)
  real(kind=rk), intent(inout) :: mask_color(0:,0:,0:)

  integer :: ix, iy, iz
  real(kind=rk), dimension(1:3) :: x_glob, x_body, v_tmp
  integer(kind=2) :: c

  if (size(mask) /= size(mask_color) .or. size(us,4) /= 3) then
    write(*,*) "mask:", shape(mask), "mask_color:", shape(mask_color), "us:", shape(us)
    call abort (08021902,"Insects: arrays have wrong size..")
  endif


  ! sometimes we have the geometry type insect but it has no wings (for example for fractal_tree), we then want to skip the rest
  if ((.not.(Insect%Wings(1)%used)).and. &
      (.not.(Insect%Wings(2)%used)).and. &
      (.not.(Insect%Wings(3)%used)).and. &
      (.not.(Insect%Wings(4)%used)) ) return

  if ((dabs(Insect%time-time)>1.0d-10) .and. root) then
    write(*,'("error! time=",es15.8," but Insect%time=",es15.8)') time, Insect%time
    write(*,'("Did you call Update_Insect before draw_insect_wings?")')
  endif

  !-----------------------------------------------------------------------------
  ! Stage I: mask + us field in BODY system
  !-- wingID: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  !-----------------------------------------------------------------------------

  ! NOTE: For an unknown reason, a different orderung (L,R,L2,R2) instead of (R,L,R2,L2) 
  ! changes the results slightly and the unit tests fail.

  if (Insect%Wings(2)%used) then
    call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%Wings(2)%color, 2_2, &
    Insect%M_g2b, Insect%Wings(2)%M_b2w, Insect%Wings(2)%x_pivot_b, &
    Insect%Wings(2)%rot_rel_wing_w, Insect%Wings(2)%side )
  endif

  if (Insect%Wings(1)%used) then
    call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%Wings(1)%color, 1_2, &
    Insect%M_g2b, Insect%Wings(1)%M_b2w, Insect%Wings(1)%x_pivot_b, &
    Insect%Wings(1)%rot_rel_wing_w, Insect%Wings(1)%side )
  endif

  if (Insect%Wings(4)%used) then
    call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%Wings(4)%color, 4_2, &
    Insect%M_g2b, Insect%Wings(4)%M_b2w, Insect%Wings(4)%x_pivot_b, &
    Insect%Wings(4)%rot_rel_wing_w, Insect%Wings(4)%side )
  endif

  if (Insect%Wings(3)%used) then
    call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%Wings(3)%color, 3_2, &
    Insect%M_g2b, Insect%Wings(3)%M_b2w, Insect%Wings(3)%x_pivot_b, &
    Insect%Wings(3)%rot_rel_wing_w, Insect%Wings(3)%side )
  endif

  !-----------------------------------------------------------------------------
  ! stage II: add body motion to wing and bring us to global system
  !-----------------------------------------------------------------------------
  ! Add solid body rotation (i.e. the velocity field that originates
  ! from the body rotation and translation). Until now, the wing velocities
  ! were the only ones set plus they are in the body reference frame
  do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
    x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
    do iy = g, size(mask,2)-1-g
        x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
        do ix = g, size(mask,1)-1-g
            x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

            c = mask_color(ix,iy,iz)
            ! skip all parts that do not belong to the wings (ie they have a different color)
            ! real comparison should usually be done with a tolerance, but since we only ever set color values and do no arithmetics, this is fine
            if (c==Insect%color_l .or. c==Insect%color_r .or. &
                c==Insect%color_l2 .or. c==Insect%color_r2 ) then

                ! disabled
                ! if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                x_body = matmul(Insect%M_g2b, x_glob)

                ! add solid body rotation in the body-reference frame, if color
                ! indicates that this part of the mask belongs to the wings
                if (mask(ix,iy,iz) > 0.0_rk) then

                    ! translational part. we compute the rotational part in the body
                    ! reference frame, therefore, we must transform the body translation
                    ! velocity Insect%vc (which is in global coordinates) to the body frame
                    v_tmp = matmul(Insect%M_g2b, Insect%vc_body_g)

                    ! add solid body rotation to the translational velocity field. Note
                    ! that rot_body_b and x_body are in the body reference frame
                    v_tmp(1) = v_tmp(1) + Insect%rot_body_b(2)*x_body(3)-Insect%rot_body_b(3)*x_body(2)
                    v_tmp(2) = v_tmp(2) + Insect%rot_body_b(3)*x_body(1)-Insect%rot_body_b(1)*x_body(3)
                    v_tmp(3) = v_tmp(3) + Insect%rot_body_b(1)*x_body(2)-Insect%rot_body_b(2)*x_body(1)

                    ! the body motion is added to the wing motion, which is already in us
                    ! and they are also in the body refrence frame. However, us has to be
                    ! in the global reference frame, so M_b2g is applied
                    us(ix,iy,iz,1:3) = matmul( Insect%M_b2g, us(ix,iy,iz,1:3)+v_tmp )
                endif
            endif
        enddo
    enddo
  enddo

end subroutine



! Wing wrapper for different wing shapes
subroutine draw_wing(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b,&
    M_b2w, x_pivot_b, rot_rel_wing_w, side)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing, wingID
  real(kind=rk), intent(in) :: M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)
  ! NOTE: for a corrugated wing, up- and downside are different, and therefore a distinction between the
  ! left- and right wing has to be made, essentially inverting the sign of the z_wing coordinate.
  character(len=1), intent(in) :: side ! can be R or L

  select case(Insect%Wings(wingID)%WingShape)
  case ("rectangular")
      call draw_wing_rectangular(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, &
      M_b2w,x_pivot_b, rot_rel_wing_w)

  case ("suzuki_butterfly")
      call draw_wing_rectangular_suzuki_butterfly(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, &
      M_b2w,x_pivot_b, rot_rel_wing_w)

  case ("suzuki")
      ! this wing has a finite thickness
      call draw_wing_suzuki(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, &
      M_b2w,x_pivot_b, rot_rel_wing_w)

  case ("TwoEllipses")
      call draw_wing_twoellipses(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, &
      M_b2w,x_pivot_b, rot_rel_wing_w)

  case default
      ! if all other options fail, we still might load coefficients from file:
      ! we assume the default to be defined in fourier coefficients, the subroutine
      ! yells if it does not recongnize the wing.
      select case (Insect%Wings(wingID)%WingShapeType)
      case ("fourier", "linear")
          ! ordinary fourier wing (wing planform described in polar coordinates with fourier coeffs for the radius)
          call draw_wing_fourier(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, M_b2w, &
          x_pivot_b, rot_rel_wing_w, side)

      case ("fourierY")
          ! fourier series for the y coordinate (used for the blade of a bristled wing)
          call draw_wing_bristled(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, M_b2w, &
          x_pivot_b, rot_rel_wing_w, side)

      case ("polygon")
          ! polygon can be an arbitrily shaped wing with coordinates given in the ini-file
          call draw_wing_polygon(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, M_b2w, &
          x_pivot_b, rot_rel_wing_w)

      case default
          call abort(26111901, "The wing-ini-setup has a TYPE setting that the code does not know: "//trim(adjustl(Insect%Wings(wingID)%WingShapeType)))
      end select

  end select

end subroutine draw_wing

!-------------------------------------------------------------------------------

! Draws a wings that is given by a radius(theta), where the radius is given
! by a Fourier series. The Fourier coefficients are stored in the insect
! datastructure, so the function Set_Wing_Fourier_coefficients must be called
! before calling this subroutine. Fourier series is evaluated in
! Radius_Fourier
subroutine draw_wing_fourier(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, M_b2w, x_pivot_b, rot_rel_wing_w, side)
  implicit none

  type(diptera), intent(inout) :: Insect
  real(kind=rk), intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk), intent(inout) :: mask(0:,0:,0:)
  real(kind=rk), intent(inout) :: us(0:,0:,0:,1:)
  real(kind=rk), intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2), intent(in) :: color_wing, wingID  !< wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  real(kind=rk),intent(in) :: M_g2b(1:3,1:3), M_b2w(1:3,1:3), x_pivot_b(1:3), rot_rel_wing_w(1:3)
  ! NOTE: for a corrugated wing, up- and downside are different, and therefore a distinction between the
  ! left- and right wing has to be made, essentially inverting the sign of the z_wing coordinate.
  character(len=1), intent(in) :: side

  integer :: ix,iy,iz,j
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3), xa(1:3), xb(1:3)
  real(kind=rk) :: R, R0, R_tmp, zz0
  real(kind=rk) :: y_tmp, x_tmp, z_tmp, s, t, D
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta, wsign
  logical :: variable_wing_thickness

  if ( ((Insect%Wings(wingID)%WingShapeType)/="linear") .and. ((Insect%Wings(wingID)%WingShapeType)/="fourier") .and. ((Insect%Wings(wingID)%WingShapeType)/="fourierY") ) then
      call abort(26111902,"draw_wing_fourier is called with a wing that is neither linear/fourier/fourierY ...")
  endif

  if (side == "R") then
      wsign = +1.0_rk
  elseif (side == "L") then
      wsign = -1.0_rk
  else
      call abort(290720, "neither R nor L wing??")
  endif

  variable_wing_thickness = .false.
  if (Insect%wings(wingID)%wing_thickness_distribution == "variable") then
      variable_wing_thickness = .true.
  endif

  s = Insect%safety
  do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2)+dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1)+dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              ! disabled
              ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
              x_body = matmul(M_g2b,x)
              x_wing = matmul(M_b2w,x_body-x_pivot_b)

              ! bounding box check: does this point lie within the bounding box? Note Insect%wing_bounding_box
              ! is set in SET_WING_BOUNDING_BOX_FOURIER
              if ( x_wing(1) >= Insect%Wings(wingID)%wing_bounding_box(1)-s &
                    .and. x_wing(1) <= Insect%Wings(wingID)%wing_bounding_box(2)+s) then
                  if ( x_wing(2) >= Insect%Wings(wingID)%wing_bounding_box(3)-s &
                        .and. x_wing(2) <= Insect%Wings(wingID)%wing_bounding_box(4)+s) then
                      if ( x_wing(3) >= Insect%Wings(wingID)%wing_bounding_box(5)-s &
                            .and. x_wing(3) <= Insect%Wings(wingID)%wing_bounding_box(6)+s) then

                          !-- get normalized angle (theta)
                          theta = atan2( x_wing(2)-Insect%Wings(wingID)%yc, x_wing(1)-Insect%Wings(wingID)%xc )
                          ! note flusi uses an angle between [0, 2*pi)
                          theta = theta + pi

                          !-- construct R by evaluating the fourier series
                          R0 = Radius_Fourier(theta, Insect, wingID)

                          !-- get smooth (radial) step function
                          R = dsqrt ( (x_wing(1)-Insect%Wings(wingID)%xc)**2 + (x_wing(2)-Insect%Wings(wingID)%yc)**2 )
                          R_tmp = step(R,R0, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)

                          ! wing corrugation (i.e. deviation from a flat plate)
                          if ( Insect%Wings(wingID)%corrugated ) then
                              ! if the wing is corrugated, its height profile is read from ini file
                              ! and interpolated at the position on the wing
                              zz0 = interp2_nonper( x_wing(1), x_wing(2), Insect%Wings(wingID)%corrugation_profile, &
                              Insect%Wings(wingID)%corrugation_array_bbox )
                          else
                              ! no corrugation - the wing is a flat surface
                              zz0 = 0.0_rk
                          endif

                          zz0 = zz0 * wsign

                          ! wing thickness
                          if ( variable_wing_thickness ) then
                              ! variable wing thickness is read from an array in the wing.ini file
                              ! and interpolated linearly at the x_wing position.
                              t = interp2_nonper( x_wing(1), x_wing(2), Insect%Wings(wingID)%wing_thickness_profile, &
                              Insect%Wings(wingID)%corrugation_array_bbox )
                          else
                              ! constant thickness, read from main params.ini file
                              t = Insect%Wings(wingID)%WingThickness
                          endif

                          ! wing damage: apply a mask [0,1] (including possibly smoothing)
                          ! to the wing
                          if ( Insect%Wings(wingID)%damaged ) then
                              D = interp2_nonper( x_wing(1), x_wing(2), Insect%Wings(wingID)%damage_mask, Insect%Wings(wingID)%corrugation_array_bbox )
                          else
                              D = 1.0_rk
                          endif

                          z_tmp = step( dabs(x_wing(3)-zz0), 0.5_rk*t, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int )
                          ! mask function approximated as product of 1D mask functions:
                          mask_tmp = z_tmp*R_tmp*D

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
                              ! us is now in body system (note M_b2w contains the stroke plane)
                              us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
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
  if (Insect%Wings(wingID)%bristles) then
      ! Loop for all bristles
      do j = 1, Insect%Wings(wingID)%n_bristles

        !   ! wing corrugation (i.e. deviation from a flat plate)
        !   if ( Insect%Wings(wingID)%corrugated ) then
        !       ! if the wing is corrugated, its height profile is read from ini file
        !       ! and interpolated at the position on the wing
        !       zz0 = interp2_nonper( Insect%Wings(wingID)%bristles_coords(j,1), Insect%Wings(wingID)%bristles_coords(j,2), &
        !       Insect%Wings(wingID)%corrugation_profile, Insect%Wings(wingID)%corrugation_array_bbox, size(Insect%Wings(wingID)%corrugation_profile,1), &
        !       size(Insect%Wings(wingID)%corrugation_profile, 2) )
        !   else
        !       ! no corrugation - the wing is a flat surface
        !       zz0 = 0.0_rk
        !   endif

        !   zz0 = zz0 * wsign

          ! start / end point (in wing coordinate system)
          xa = (/Insect%Wings(wingID)%bristles_coords(j,1), Insect%Wings(wingID)%bristles_coords(j,2), 0.0_rk/)
          xb = (/Insect%Wings(wingID)%bristles_coords(j,3), Insect%Wings(wingID)%bristles_coords(j,4), 0.0_rk/)
          R  = Insect%Wings(wingID)%bristles_coords(j,5)

          ! note input to draw_bristle in in wing coordinates
          call draw_bristle(xa, xb, R, xx0, ddx, mask, mask_color, us, Insect, color_wing, M_g2b, M_b2w, x_pivot_b, rot_rel_wing_w)
      enddo
  endif

end subroutine draw_wing_fourier


!-------------------------------------------------------------------------------
! Draws a membranous central part of a bristled wing, using the same storage space as
! for a Fourier wing, but the algorithm is different.
subroutine draw_blade_fourier(xx0, ddx, mask, mask_color, us,Insect,color_wing,wingID,M_g2b,M_b2w,x_pivot_b,rot_rel_wing_w,side)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing, wingID  !< wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  real(kind=rk),intent(in)::M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)
  character(len=1), intent(in) :: side

  integer :: ix,iy,iz
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=rk) :: R, R0, R_tmp, zz0
  real(kind=rk) :: y_tmp, x_tmp, z_tmp, s, t
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta
  real(kind=rk) :: rblade, ylte, xte, xle, wsign

  ! NOTE: prescribed wing deformation is untested work in progress! -TE 02/2026
  integer :: i, j, k, a, b, c
  real(kind=rk) :: tt, t1, t2, def1, def2
  !KVN-2025<<<<<

  !-- reset the bounding box
  Insect%Wings(wingID)%wing_bounding_box(1:4) = (/-1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk/)

  !-- blade length
  rblade = Insect%Wings(wingID)%yc

  if (side == "R") then
      wsign = +1.0_rk
  elseif (side == "L") then
      wsign = -1.0_rk
  else
      call abort(290720, "neither R nor L wing??")
  endif

  s = Insect%safety
  do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2)+dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1)+dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              ! disabled
              ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
              x_body = matmul(M_g2b, x)
              x_wing = matmul(M_b2w, x_body-x_pivot_b)

              ! bounding box check: does this point lie within the bounding box? Note Insect%wing_bounding_box
              ! is set in SET_WING_BOUNDING_BOX_FOURIER
              if ( x_wing(1) >= Insect%Wings(wingID)%wing_bounding_box(1)-s .and. x_wing(1) <= Insect%Wings(wingID)%wing_bounding_box(2)+s) then
                  if ( x_wing(2) > 0.0_rk .and. x_wing(2) < rblade ) then
                      if ( x_wing(3) >= Insect%Wings(wingID)%wing_bounding_box(5)-s .and. x_wing(3) <= Insect%Wings(wingID)%wing_bounding_box(6)+s) then
                            !-- calculate the polar parameter (normalized angle)
                            ylte = x_wing(2)
                            theta = dacos( 1.0_rk - 2.0_rk*ylte/rblade )
                            theta = theta / (2.0_rk*pi)
                            
                            !-- construct xle by evaluating the Fourier series
                            xle = Radius_Fourier( 2.0_rk*pi*theta,Insect,wingID)
                                
                            !-- construct xte by evaluating the Fourier series
                            xte = Radius_Fourier( 2.0_rk*pi*(1.0_rk-theta),Insect,wingID)
                            
                            !-- amplitude
                            R0 = 0.5_rk * (xle-xte)
                            
                            !-- get smooth rectangular function
                            R = dabs ( x_wing(1) - 0.5*(xte+xle) )
                            R_tmp = step(R, R0, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            
                            ! wing corrugation (i.e. deviation from a flat plate)
                            if ( Insect%Wings(wingID)%corrugated ) then
                                ! if the wing is corrugated, its height profile is read from ini file
                                ! and interpolated at the position on the wing
                                zz0 = interp2_nonper( x_wing(1), x_wing(2), Insect%Wings(wingID)%corrugation_profile, Insect%Wings(wingID)%corrugation_array_bbox )
                            else
                              ! no corrugation - the wing is a flat surface
                              zz0 = 0.0_rk
                            endif
                            
                            ! !KVN-2025>>>>>
                            ! ! wing deformation
                            ! ! NOTE: prescribed wing deformation is untested work in progress! -TE 02/2026
                            ! if ( Insect%deformable(wingID) ) then
                            !   a = Insect%deformation_a(wingID)
                            !   b = Insect%deformation_b(wingID)
                            !   c = Insect%deformation_c(wingID)
                            !   tt = mod(Insect%time, Insect%deformations(a*c,1,wingID))
                            !   do k = 1, c-1
                            !     t1 = Insect%deformations(a*k-1,1,wingID)
                            !     t2 = Insect%deformations(a*k+1,1,wingID)
                            !     !if (Insect%time >= t1 .AND. Insect%time <= t2) then
                            !     if (tt >= t1 .AND. tt <= t2) then
                            !       exit
                            !     endif
                            !   enddo
                            !   !if (k == c-1) then
                            !   !       t1 = Insect%deformations(a*c,1,wingID)
                            !   !       k = mod(k,c)
                            !   !       t1 = t1 + Insect%deformations(a*k-1,1,wingID)
                            !   !       t2 = t1 + Insect%deformations(a*k+1,1,wingID)
                            !   !end if
                            !   do j = 1, b
                            !     do i = 1, a
                            !       if (k>=c-1) then
                            !         Insect%deformation_profile(i,j,wingID) = Insect%deformations(a*k-1,j+1,wingID)
                            !       else
                            !         def1 = Insect%deformations(a*k-1,j+1,wingID)
                            !         def2 = Insect%deformations(a*k+1,j+1,wingID)
                            !         Insect%deformation_profile(i,j,wingID) = def1 + (tt-t1)/(t2-t1)*(def2-def1)
                            !         !Insect%deformation_profile(i,j,wingID) = def1 + (Insect%time-t1)/(t2-t1)*(def2-def1)
                            !         !Insect%deformation_profile(i,j,wingID) = (Insect%deformations(a*k-1,j+1,wingID)+Insect%deformations(a*k+1,j+1,wingID))/2_rk
                            !       endif
                            !     enddo
                            !   enddo
                            !   zz0 = interp2_nonper( x_wing(1), x_wing(2), Insect%deformation_profile(:,:,wingID), Insect%deformation_array_bbox(1:4,wingID), a, b )
                            !   Insect%wing_bounding_box(5,wingID) = min(Insect%wing_bounding_box(5,wingID), zz0 - Insect%WingThickness / 2.0_pr)
                            !   Insect%wing_bounding_box(6,wingID) = max(Insect%wing_bounding_box(6,wingID), zz0 + Insect%WingThickness / 2.0_pr)
                            !   !write(*,*) "x_wing(1)=",x_wing(1),"x_wing(2)=",x_wing(2),"array=",Insect%deformation_array_bbox(:,wingID)
                            !   !v_a = corrugation_a(wingID)
                            !   !v_b = corrugation_b(wingID)
                            !   !v_xmin = Insect%Wings(wingID)%corrugation_array_bbox(1,wingID)
                            !   !v_xmax = Insect%Wings(wingID)%corrugation_array_bbox(2,wingID)
                            !   !v_dx = (v_xmax-v_xmin)/(v_b-1)
                            !   !do j = 1, v_b
                            !   !       v_x = v_xmin + (j-1)*v_dx
                            !   !       do i = 1, v_a
                            !   !               !corrugation_profile(i,j,wingID) = (1.0_pr+dsin(Insect%time/0.1_rk*2_rk*pi))/20.0_pr*dsin((v_x-v_xmin)/(v_xmax-v_xmin)*2_rk*pi)
                            !   !               corrugation_profile(i,j,wingID) = dsin(Insect%time/0.01_rk*2_rk*pi)/20.0_pr*dsin((v_x-v_xmin)/(v_xmax-v_xmin)*2_rk*pi)
                            !   !       enddo
                            !   !enddo
                            !   !zz0 = interp2_nonper( x_wing(1), x_wing(2), corrugation_profile(:,:,wingID), Insect%Wings(wingID)%corrugation_array_bbox, corrugation_a(wingID), corrugation_b(wingID) )
                            ! endif
                            ! !KVN-2025<<<<<

                            zz0 = zz0 * wsign

                            ! wing thickness
                            if ( Insect%wings(wingID)%wing_thickness_distribution == "variable") then
                                ! variable wing thickness is read from an array in the wing.ini file
                                ! and interpolated linearly at the x_wing position.
                                t = interp2_nonper( x_wing(1), x_wing(2), Insect%Wings(wingID)%wing_thickness_profile, Insect%Wings(wingID)%corrugation_array_bbox )
                            else
                                ! constant thickness, read from main params.ini file
                                t = Insect%Wings(wingID)%WingThickness
                            endif
                            
                            z_tmp = step( abs(x_wing(3)-zz0), 0.5_rk*t, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int ) ! thickness
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
! Suzuki's rectangular wingshape as defined in section 3.4 of my thesis (Thomas,
! "Numerical modeling of fluid-structure interaction in bio-inspired propulsion")
! This wing has finite thickness.
!-------------------------------------------------------------------------------
subroutine draw_wing_suzuki(xx0, ddx, mask, mask_color, us,Insect,color_wing,wingID,M_g2b,M_b2w,x_pivot_b,rot_rel_wing_w)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: color_wing, wingID  !< wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
    real(kind=rk),intent(in)::M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

    integer :: ix,iy,iz
    real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
    real(kind=rk) :: R, R0, R_tmp
    real(kind=rk) :: y_tmp, x_tmp, z_tmp, y_left, y_right
    real(kind=rk) :: v_tmp(1:3), mask_tmp, theta, x_top, x_bot

    ! wing shape (determine between which x-values (x_bot, x_top) the wing is
    ! these values depend on the spanwise direction (which is y)
    x_top = 0.0667_rk
    x_bot = -0.35_rk

    y_right = 1.0_rk
    y_left = 0.1667_rk

    do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
        x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

                !-- define the various coordinate systems we are going to use
                ! disabled
                ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                x_body = matmul(M_g2b, x)
                x_wing = matmul(M_b2w, x_body-x_pivot_b)

                ! spanwise length:
                if ((x_wing(2)>=y_left-Insect%safety).and.(x_wing(2)<=y_right+Insect%safety)) then
                    ! thickness: (note left and right wing have a different orientation of the z-axis
                    ! but this does not matter since this is the same.
                    if (dabs(x_wing(3))<=0.5*Insect%Wings(wingID)%WingThickness + Insect%safety) then

                        ! in the x-direction, the actual wing shape plays.
                        if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
                            !-- smooth length
                            if ( x_wing(2) < 0.5_rk*(y_left+y_right) ) then
                                y_tmp = step(-(x_wing(2)-y_left), 0.0_rk, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            else
                                y_tmp = step( (x_wing(2)-y_left), y_right-y_left, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            endif

                            !-- smooth height
                            z_tmp = step(dabs(x_wing(3)),0.5_rk*Insect%Wings(wingID)%WingThickness, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int) ! thickness

                            !-- smooth shape
                            if (x_wing(1) < 0.0_rk) then
                                x_tmp = step(-x_wing(1),-x_bot, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            else
                                x_tmp = step( x_wing(1), x_top, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            endif

                            mask_tmp = z_tmp*y_tmp*x_tmp

                            if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0_rk)) then
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
                                us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
                            endif

                        endif
                    endif
                endif
            enddo
        enddo
    enddo
end subroutine draw_wing_suzuki

!-------------------------------------------------------------------------------

subroutine draw_wing_rectangular(xx0, ddx, mask, mask_color, us,Insect,color_wing,wingID,M_g2b,M_b2w,x_pivot_b,rot_rel_wing_w)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: color_wing, wingID  !< wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
    real(kind=rk),intent(in)::M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

    integer :: ix,iy,iz
    real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
    real(kind=rk) :: R, R0, R_tmp
    real(kind=rk) :: y_tmp, x_tmp, z_tmp, y_left, y_right
    real(kind=rk) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

    ! wing shape (determine between which x-values (x_bot, x_top) the wing is
    ! these values depend on the spanwise direction (which is y)
    x_top = 0.085_rk ! determinded from a calliphora wing, roughly estimated
    x_bot = -(0.294_rk-x_top) ! to get the same aspect ratio as in calliphora (Engels et al., RSI2020)

    y_right = 1.0_rk
    y_left = 0.0_rk

    do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
        x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

                !-- define the various coordinate systems we are going to use
                ! disabled
                ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                x_body = matmul(M_g2b,x)
                x_wing = matmul(M_b2w,x_body-x_pivot_b)

                ! spanwise length:
                if ((x_wing(2)>=y_left-Insect%safety).and.(x_wing(2)<=y_right+Insect%safety)) then
                    ! thickness: (note left and right wing have a different orientation of the z-axis
                    ! but this does not matter since this is the same.
                    if (dabs(x_wing(3))<=0.5*Insect%Wings(wingID)%WingThickness + Insect%safety) then

                        ! in the x-direction, the actual wing shape plays.
                        if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
                            !-- smooth length
                            if ( x_wing(2) < 0.5_rk*(y_left+y_right) ) then
                                y_tmp = step(-(x_wing(2)-y_left), 0.0_rk, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            else
                                y_tmp = step( (x_wing(2)-y_left), y_right-y_left, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            endif

                            !-- smooth height
                            z_tmp = step(dabs(x_wing(3)), 0.5_rk*Insect%Wings(wingID)%WingThickness, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int) ! thickness

                            !-- smooth shape
                            if (x_wing(1)<0.0_rk) then
                                x_tmp = step(-x_wing(1),-x_bot, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            else
                                x_tmp = step( x_wing(1), x_top, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            endif

                            mask_tmp = z_tmp*y_tmp*x_tmp

                            if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0_rk)) then
                                mask(ix,iy,iz) = mask_tmp
                                mask_color(ix,iy,iz) = color_wing
                                !------------------------------------------------
                                ! solid body rotation
                                ! Attention: the Matrix transpose(M) brings us back to the body
                                ! coordinate system, not to the inertial frame. this is done in
                                ! the main routine Draw_Insect
                                !------------------------------------------------

                                ! v_tmp is wing in wing reference frame
                                v_tmp(1) = rot_rel_wing_w(2)*x_wing(3)-rot_rel_wing_w(3)*x_wing(2)
                                v_tmp(2) = rot_rel_wing_w(3)*x_wing(1)-rot_rel_wing_w(1)*x_wing(3)
                                v_tmp(3) = rot_rel_wing_w(1)*x_wing(2)-rot_rel_wing_w(2)*x_wing(1)

                                ! note we set this only if it is a part of the wing
                                us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
                                ! the us velocity is now in the body system
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
end subroutine draw_wing_rectangular


subroutine draw_wing_rectangular_suzuki_butterfly(xx0, ddx, mask, mask_color, us,Insect,color_wing,wingID,M_g2b,M_b2w,x_pivot_b,rot_rel_wing_w)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: color_wing, wingID  !< wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
    real(kind=rk),intent(in)::M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)

    integer :: ix,iy,iz
    real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
    real(kind=rk) :: R, R0, R_tmp
    real(kind=rk) :: y_tmp, x_tmp, z_tmp, y_left, y_right
    real(kind=rk) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

    ! wing shape (determine between which x-values (x_bot, x_top) the wing is
    ! these values depend on the spanwise direction (which is y)
    x_top = 0.5_rk ! suzukis butterfly, freeflight configuration (a different wing than the finite-thickness classical suzuki test case)
    x_bot = -0.5_rk
    y_right = 1.0_rk
    y_left = 0.0_rk

    do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
        x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

                !-- define the various coordinate systems we are going to use
                ! disabled
                ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                x_body = matmul(M_g2b,x)
                x_wing = matmul(M_b2w,x_body-x_pivot_b)

                ! spanwise length:
                if ((x_wing(2)>=y_left-Insect%safety).and.(x_wing(2)<=y_right+Insect%safety)) then
                    ! thickness: (note left and right wing have a different orientation of the z-axis
                    ! but this does not matter since this is the same.
                    if (dabs(x_wing(3))<=0.5*Insect%Wings(wingID)%WingThickness + Insect%safety) then

                        ! in the x-direction, the actual wing shape plays.
                        if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
                            !-- smooth length
                            if ( x_wing(2) < 0.5_rk*(y_left+y_right) ) then
                                y_tmp = step(-(x_wing(2)-y_left), 0.0_rk, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            else
                                y_tmp = step( (x_wing(2)-y_left), y_right-y_left, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            endif

                            !-- smooth height
                            z_tmp = step(dabs(x_wing(3)), 0.5_rk*Insect%Wings(wingID)%WingThickness, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int) ! thickness

                            !-- smooth shape
                            if (x_wing(1)<0.0_rk) then
                                x_tmp = step(-x_wing(1),-x_bot, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            else
                                x_tmp = step( x_wing(1), x_top, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                            endif

                            mask_tmp = z_tmp*y_tmp*x_tmp

                            if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0_rk)) then
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
                                us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
end subroutine draw_wing_rectangular_suzuki_butterfly


!-------------------------------------------------------------------------------
! Draws a wing
! here, a wing is a rigid plate of constant thickness that differs from
! a rectangular plate only in the x-direction
!
! note to save a bit of computing time, we first check the easy
! conditions (thickness and spanwise length) and then the shape
! function since this saves many evaluations of the shape.
subroutine draw_wing_twoellipses(xx0, ddx, mask, mask_color, us,Insect,color_wing,wingID,M_g2b,M_b2w,x_pivot_b,rot_rel_wing_w)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2),intent(in) :: color_wing, wingID  !< wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  real(kind=rk),intent(in)::M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)


  integer :: ix,iy,iz
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=rk) :: R, R0, R_tmp,a_body
  real(kind=rk) :: y_tmp, x_tmp, z_tmp
  real(kind=rk) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

  a_body = 0.5_rk

  do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
    x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
    do iy = g, size(mask,2)-1-g
      x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
      do ix = g, size(mask,1)-1-g
        x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

        !-- define the various coordinate systems we are going to use
        ! disabled
        ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

        x_body = matmul(M_g2b,x)
        x_wing = matmul(M_b2w,x_body-x_pivot_b)

        ! spanwise length:
        if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=1.0_rk + Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%Wings(wingID)%WingThickness + Insect%safety) then
            ! wing shape (determine between which x-values (x_bot, x_top) the wing is
            ! these values depend on the spanwise direction (which is y)
            if ((1._rk - ((x_wing(2)-a_body)**2)/(a_body**2)) >= 0.0_rk) then
              x_top =  dsqrt((Insect%b_top**2)*(1._rk-((x_wing(2)-a_body)**2)/(a_body**2)))
              x_bot = -dsqrt((Insect%b_bot**2)*(1._rk-((x_wing(2)-a_body)**2)/(a_body**2)))
            else
              x_top = 0.0_rk
              x_bot = 0.0_rk
            endif

            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
              if (x_wing(2)<0.0_rk) then  ! xs is chordlength coordinate
                y_tmp = step(-x_wing(2),0.0_rk, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
              else
                y_tmp = step( x_wing(2),1.0_rk, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
              endif

              !-- smooth height
              z_tmp = step(dabs(x_wing(3)),0.5_rk*Insect%Wings(wingID)%WingThickness, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.0_rk) then
                x_tmp = step(-x_wing(1),-x_bot, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
              else
                x_tmp = step( x_wing(1), x_top, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
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
                us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_twoellipses

!-------------------------------------------------------------------------------
! Bristled wing
!-------------------------------------------------------------------------------
subroutine draw_wing_bristled(xx0, ddx, mask, mask_color, us,Insect,color_wing,wingID,M_g2b,M_b2w,x_pivot_b,rot_rel_wing_w,side)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
  real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
  integer(kind=2), intent(in) :: color_wing, wingID  !< wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  real(kind=rk),intent(in)::M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)
  character(len=1), intent(in) :: side

  integer :: ix,iy,iz,j
  real(kind=rk) :: x_body(1:3),x_wing(1:3),x(1:3),xa(1:3),xb(1:3)
  real(kind=rk) :: R,s,wsign,zz0
  real(kind=rk) :: v_tmp(1:3)

  if (side == "R") then
      wsign = +1.0_rk
  elseif (side == "L") then
      wsign = -1.0_rk
  else
      call abort(290720, "neither R nor L wing??")
  endif

  ! Draw the membranous blade using Fourier series
  call draw_blade_fourier(xx0, ddx, mask, mask_color, us,Insect,color_wing,wingID,M_g2b,M_b2w,&
       x_pivot_b,rot_rel_wing_w,side)

  ! Set the solid velocity
  s = Insect%safety
  do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2)+dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1)+dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              ! disabled
              ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
              x_body = matmul(M_g2b,x)
              x_wing = matmul(M_b2w,x_body-x_pivot_b)

              ! bounding box check: does this point lie within the bounding box? Note Insect%wing_bounding_box
              ! is set in SET_WING_BOUNDING_BOX_FOURIER
              if ( x_wing(1) >= Insect%Wings(wingID)%wing_bounding_box(1)-s &
                    .and. x_wing(1) <= Insect%Wings(wingID)%wing_bounding_box(2)+s) then
                  if ( x_wing(2) >= Insect%Wings(wingID)%wing_bounding_box(3)-s &
                        .and. x_wing(2) <= Insect%Wings(wingID)%wing_bounding_box(4)+s) then
                      if ( x_wing(3) >= Insect%Wings(wingID)%wing_bounding_box(5)-s &
                            .and. x_wing(3) <= Insect%Wings(wingID)%wing_bounding_box(6)+s) then
                          !-----------------------------------------
                          ! set new value for solid velocity us
                          !-----------------------------------------
                          ! real comparison should usually be done with a tolerance, but since we only ever set color values and do no arithmetics, this is fine
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
                              us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
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
  if (Insect%Wings(wingID)%bristles) then
      ! Loop for all bristles
      do j = 1, Insect%Wings(wingID)%n_bristles

          ! wing corrugation (i.e. deviation from a flat plate)
          if ( Insect%Wings(wingID)%corrugated ) then
              ! if the wing is corrugated, its height profile is read from ini file
              ! and interpolated at the position on the wing
              zz0 = interp2_nonper( Insect%Wings(wingID)%bristles_coords(j,1), Insect%Wings(wingID)%bristles_coords(j,2), &
              Insect%Wings(wingID)%corrugation_profile, &
              Insect%Wings(wingID)%corrugation_array_bbox )
          else
              ! no corrugation - the wing is a flat surface
              zz0 = 0.0_rk
          endif

        !   !KVN-2025>>>>>   
        !   ! NOTE: prescribed wing deformation is untested work in progress! -TE 02/2026
        !   if ( Insect%deformable(wingID) ) then
        !     ! wing deformation
        !     zz0 = interp2_nonper(Insect%Wings(wingID)%bristles_coords(,1), &
        !                 Insect%Wings(wingID)%bristles_coords(,2), &
        !                 Insect%deformation_profile(:,:,wingID), &
        !                 Insect%deformation_array_bbox(1:4,wingID), &
        !                 Insect%deformation_a(wingID), &
        !                 Insect%deformation_b(wingID))
        !   endif
        !   !KVN-2025>>>>>

          zz0 = zz0 * wsign

          ! start / end point (in wing coordinate system)
          !KVN-2025>>>>>
          if ( Insect%Wings(wingID)%bristles3D ) then
             xa = (/Insect%Wings(wingID)%bristles_coords(j,1), Insect%Wings(wingID)%bristles_coords(j,2), zz0+wsign*Insect%Wings(wingID)%bristles_coords(j,3)/)
             xb = (/Insect%Wings(wingID)%bristles_coords(j,4), Insect%Wings(wingID)%bristles_coords(j,5), zz0+wsign*Insect%Wings(wingID)%bristles_coords(j,6)/)
             R  = Insect%Wings(wingID)%bristles_coords(j,7)
          else
             xa = (/Insect%Wings(wingID)%bristles_coords(j,1), Insect%Wings(wingID)%bristles_coords(j,2), zz0/)
             xb = (/Insect%Wings(wingID)%bristles_coords(j,3), Insect%Wings(wingID)%bristles_coords(j,4), zz0/)
             R  = Insect%Wings(wingID)%bristles_coords(j,5)
          endif
          !KVN-2025<<<<<  

          ! note input to draw_bristle in in wing coordinates
          call draw_bristle(xa, xb, R, xx0, ddx, mask, mask_color, us, Insect, color_wing, M_g2b, M_b2w, x_pivot_b, rot_rel_wing_w)
      enddo
  endif

  !-----------------------------------------------------------------------------
  ! effective membrane
  !-----------------------------------------------------------------------------
  if (Insect%Wings(wingID)%bristles_simplex) then
      ! Loop for all bristles
      do j = 1, Insect%Wings(wingID)%n_bristles-1
          ! draw a membrane triangular element
          call draw_trianglular_prism(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_g2b, M_b2w, x_pivot_b, rot_rel_wing_w, &
                            Insect%Wings(wingID)%bristles_coords(j,1), Insect%Wings(wingID)%bristles_coords(j,2), &
                            Insect%Wings(wingID)%bristles_coords(j,3), Insect%Wings(wingID)%bristles_coords(j,4), &
                            Insect%Wings(wingID)%bristles_coords(j+1,1), Insect%Wings(wingID)%bristles_coords(j+1,2), wingID)
          ! draw a membrane triangular element
          call draw_trianglular_prism(xx0, ddx, mask, mask_color, us, Insect, color_wing, M_g2b, M_b2w, x_pivot_b, rot_rel_wing_w, &
                            Insect%Wings(wingID)%bristles_coords(j,3), Insect%Wings(wingID)%bristles_coords(j,4), &
                            Insect%Wings(wingID)%bristles_coords(j+1,3), Insect%Wings(wingID)%bristles_coords(j+1,4), &
                            Insect%Wings(wingID)%bristles_coords(j+1,1), Insect%Wings(wingID)%bristles_coords(j+1,2), wingID)

      enddo
  endif

end subroutine draw_wing_bristled


!-------------------------------------------------------------------------------
subroutine draw_wing_polygon(xx0, ddx, mask, mask_color, us, Insect, color_wing, wingID, M_g2b, M_b2w, x_pivot_b,rot_rel_wing_w)
    implicit none 

    type(diptera), intent(inout) :: Insect
    real(kind=rk), intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk), intent(inout) :: mask(0:,0:,0:), mask_color(0:,0:,0:)
    real(kind=rk), intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2), intent(in) :: color_wing, wingID
    real(kind=rk), intent(in) :: M_g2b(1:3,1:3), M_b2w(1:3,1:3), x_pivot_b(1:3), rot_rel_wing_w(1:3)

    integer(kind=ik) :: i, j, k, ix, iy, iz
    integer(kind=ik) :: n, iseg, idx_p1, idx_p2
    integer(kind=ik) :: Bs(1:3)
    integer(kind=ik) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer(kind=ik) :: xmin_seg, xmax_seg, ymin_seg, ymax_seg, zmin_seg, zmax_seg
    real(kind=rk) :: x, y, z
    real(kind=rk) :: x_g(1:3), x_b(1:3), x_w(1:3)
    real(kind=rk) :: corner_w(1:3), corner_b(1:3), corner_g(1:3)
    real(kind=rk) :: M_w2b(1:3,1:3), M_b2g(1:3,1:3)
    real(kind=rk) :: band_width
    real(kind=rk) :: wingbox_w(1:6), x_wingbox_w(1:2), y_wingbox_w(1:2), z_wingbox_w(1:2)
    real(kind=rk) :: segmentbox_w(1:6), x_segmentbox_w(1:2), y_segmentbox_w(1:2), z_segmentbox_w(1:2)
    real(kind=rk) :: wingbox_g_min(1:3), wingbox_g_max(1:3), segmentbox_g_min(1:3), segmentbox_g_max(1:3)
    real(kind=rk) :: p1(1:2), p2(1:2), xa(1:3), xb(1:3), R

    real(kind=rk), allocatable, save :: tmp_dist_xy(:,:,:)
    logical, allocatable, save       :: tmp_active(:,:,:)
    real(kind=rk) :: d_seg, zmin_wing, zmax_wing
    real(kind=rk) :: phi_xy, phi_z, phi, mask_value, v_tmp(1:3)
    logical :: inside_polygon


    if ((Insect%Wings(wingID)%WingShapeType) /= "polygon") call abort(01072601,"draw_wing_polygon called with non-polygon wing...")
    if (.not. allocated(Insect%Wings(wingID)%polygon_wings)) call abort(01072602, "draw_wing_polygon: polygon_wings is not allocated")
    if (Insect%Wings(wingID)%damaged) call abort(21072601, "draw_wing_polygon can currently only be used with non-damaged wings")
    if (Insect%Wings(wingID)%corrugated) call abort(21072602, "draw_wing_polygon can currently only be used with non-corrugated wings")
    if ( Insect%wings(wingID)%wing_thickness_distribution=="variable") call abort(21072603, "draw_wing_polygon can currently only be used with a constant wing thickness")

    !-----------------------------------------------------------------------------
    ! bristles (in this type of wing we need to start with them)
    !-----------------------------------------------------------------------------
    ! generic polygon wings can also have bristles: they are read from inifile.
    !
    ! NOTE: unlike the fourier wings, the polygon routines contain several layers of bounding boxes, 
    ! and this may result in incomplete bristles. Hence, we draw the bristles first, then add the wing.
    ! (TE, 30/jul/2026)
    if (Insect%Wings(wingID)%bristles) then
        ! Loop for all bristles
        do j = 1, Insect%Wings(wingID)%n_bristles
            ! start / end point (in wing coordinate system)
            xa = (/Insect%Wings(wingID)%bristles_coords(j,1), Insect%Wings(wingID)%bristles_coords(j,2), 0.0_rk/)
            xb = (/Insect%Wings(wingID)%bristles_coords(j,3), Insect%Wings(wingID)%bristles_coords(j,4), 0.0_rk/)
            R = Insect%Wings(wingID)%bristles_coords(j,5)

            ! note input to draw_bristle is in wing coordinates
            call draw_bristle(xa, xb, R, xx0, ddx, mask, mask_color, us, Insect, color_wing, M_g2b, M_b2w, &
                              x_pivot_b, rot_rel_wing_w)
        enddo
    endif

    !-----------------------------------------------------------------------------
    ! code for actual polygon wing
    !-----------------------------------------------------------------------------
    Bs(1) = size(mask,1) - 2*g
    Bs(2) = size(mask,2) - 2*g
    Bs(3) = size(mask,3) - 2*g
    n          = Insect%Wings(wingID)%n_polygon_points
    band_width = Insect%safety

   
    ! Local wing bounding box in wing coordinates extended by safety zone
    ! a bit redundant but used for readability
    wingbox_w(1) = Insect%Wings(wingID)%wing_bounding_box(1) - band_width
    wingbox_w(2) = Insect%Wings(wingID)%wing_bounding_box(2) + band_width
    wingbox_w(3) = Insect%Wings(wingID)%wing_bounding_box(3) - band_width
    wingbox_w(4) = Insect%Wings(wingID)%wing_bounding_box(4) + band_width
    wingbox_w(5) = Insect%Wings(wingID)%wing_bounding_box(5) - band_width
    wingbox_w(6) = Insect%Wings(wingID)%wing_bounding_box(6) + band_width

    ! Reverse transformations:
    ! usual point transform: global -> body -> wing
    ! here we need:          wing -> body -> global
    M_w2b = transpose(M_b2w)
    M_b2g = transpose(M_g2b)
    
    ! we need the 6 values to decribe the 8 corners of the bounding box
    x_wingbox_w = (/wingbox_w(1), wingbox_w(2)/)
    y_wingbox_w = (/wingbox_w(3), wingbox_w(4)/)
    z_wingbox_w = (/wingbox_w(5), wingbox_w(6)/)

    ! we are saving global coordinates of the bounding box 
    ! with each values wingbox_g_min(3) - x,y,z
    wingbox_g_min = 999.d9
    wingbox_g_max = -999.d9

    ! Loop over the 8 corners of the local WING bounding box
    do k = 1, 2
        do j = 1, 2
            do i = 1, 2
                ! corner
                corner_w = (/x_wingbox_w(i), y_wingbox_w(j), z_wingbox_w(k)/)
                ! wing -> body
                corner_b = x_pivot_b + matmul(M_w2b, corner_w)
                ! body -> global
                corner_g = Insect%xc_body_g + matmul(M_b2g, corner_b)

                ! Global axis-aligned box around the tilted wing box
                wingbox_g_min = min(wingbox_g_min, corner_g)
                wingbox_g_max = max(wingbox_g_max, corner_g)            
            enddo
        enddo
    enddo

    ! Convert global wing box to index range on the current block.
    xmin = floor( (wingbox_g_min(1) - xx0(1)) / ddx(1) )
    xmax = ceiling( (wingbox_g_max(1) - xx0(1)) / ddx(1) )
    ymin = floor( (wingbox_g_min(2) - xx0(2)) / ddx(2) )
    ymax = ceiling( (wingbox_g_max(2) - xx0(2)) / ddx(2) )
    zmin = floor( (wingbox_g_min(3) - xx0(3)) / ddx(3) )
    zmax = ceiling( (wingbox_g_max(3) - xx0(3)) / ddx(3) )

    ! Clip to non-ghost points of the current block (like in draw_body_superSTL)
    ! (g) and not (g+1) because of zero indexing
    xmin = max(xmin, g)
    ymin = max(ymin, g)
    zmin = max(zmin, g)

    xmax = min(xmax, size(mask,1)-1-g)
    ymax = min(ymax, size(mask,2)-1-g)
    zmax = min(zmax, size(mask,3)-1-g)

    ! If the clipped range is empty, the current block is not affected by this wing.
    if (xmax-xmin+1 <= 0) return
    if (ymax-ymin+1 <= 0) return
    if (zmax-zmin+1 <= 0) return
    
    if (.not. allocated(tmp_dist_xy)) then
        allocate(tmp_dist_xy(g:Bs(1)+g-1, &
                            g:Bs(2)+g-1, &
                            g:Bs(3)+g-1))
    endif
    if (.not. allocated(tmp_active)) then
        allocate(tmp_active(g:Bs(1)+g-1, &
                            g:Bs(2)+g-1, &
                            g:Bs(3)+g-1))
    endif

    tmp_dist_xy = band_width
    tmp_active  = .false.

    !-----------------------------------------------------------------------------
    ! Loop over all polygon segments!
    ! Each segment gets its own local 3D segment box. This box is transformed to 
    ! the global bounding box and converted to an index range.
    !-----------------------------------------------------------------------------
    do iseg = 1, n

        idx_p1 = iseg
        if (iseg < n) then
            idx_p2 = iseg + 1
        else
            idx_p2 = 1
        endif

        p1 = Insect%Wings(wingID)%polygon_wings(idx_p1, 1:2)
        p2 = Insect%Wings(wingID)%polygon_wings(idx_p2, 1:2)

        ! cheap geometry check, building a second local bounding box around segment 
        ! to check wether or not the segment is important for the specific block
        segmentbox_w(1) = min(p1(1), p2(1)) - band_width
        segmentbox_w(2) = max(p1(1), p2(1)) + band_width
        segmentbox_w(3) = min(p1(2), p2(2)) - band_width
        segmentbox_w(4) = max(p1(2), p2(2)) + band_width
        segmentbox_w(5) = Insect%Wings(wingID)%wing_bounding_box(5) - band_width
        segmentbox_w(6) = Insect%Wings(wingID)%wing_bounding_box(6) + band_width

        x_segmentbox_w = (/segmentbox_w(1), segmentbox_w(2)/)
        y_segmentbox_w = (/segmentbox_w(3), segmentbox_w(4)/)
        z_segmentbox_w = (/segmentbox_w(5), segmentbox_w(6)/)

        ! we are saving global coordinates of the bounding box with each values wingbox_g_min(3) - x,y,z
        segmentbox_g_min = 999.d9
        segmentbox_g_max = -999.d9

        ! Loop over the 8 corners of the local SEGMENT bounding box
        do k = 1, 2
            do j = 1, 2
                do i = 1, 2
                    ! corner
                    corner_w = (/x_segmentbox_w(i), y_segmentbox_w(j), z_segmentbox_w(k)/)
                    ! wing -> body
                    corner_b = x_pivot_b + matmul(M_w2b, corner_w)
                    ! body -> global
                    corner_g = Insect%xc_body_g + matmul(M_b2g, corner_b)

                    ! Global axis-aligned box around the transformed wing box.
                    segmentbox_g_min = min(segmentbox_g_min, corner_g)
                    segmentbox_g_max = max(segmentbox_g_max, corner_g)            
                enddo
            enddo
        enddo

        ! index range on the current block
        xmin_seg = floor( (segmentbox_g_min(1) - xx0(1)) / ddx(1) )
        xmax_seg = ceiling( (segmentbox_g_max(1) - xx0(1)) / ddx(1) )
        ymin_seg = floor( (segmentbox_g_min(2) - xx0(2)) / ddx(2) )
        ymax_seg = ceiling( (segmentbox_g_max(2) - xx0(2)) / ddx(2) )
        zmin_seg = floor( (segmentbox_g_min(3) - xx0(3)) / ddx(3) )
        zmax_seg = ceiling( (segmentbox_g_max(3) - xx0(3)) / ddx(3) )

        ! again clip index range if not on block
        xmin_seg = max(xmin_seg, xmin)
        ymin_seg = max(ymin_seg, ymin)
        zmin_seg = max(zmin_seg, zmin)

        xmax_seg = min(xmax_seg, xmax)
        ymax_seg = min(ymax_seg, ymax)
        zmax_seg = min(zmax_seg, zmax)

        ! contniue if segment is not on block
        if (xmin_seg > xmax_seg) cycle
        if (ymin_seg > ymax_seg) cycle
        if (zmin_seg > zmax_seg) cycle

        ! -> now we have evaluated for each point if it is important for the current 
        ! block and extracted the indizies in close range to the line segment


        !-------------------------------------------------------------------------
        ! Now only points in this small segment box are considered.
        ! For each point:
        !     global -> body -> wing
        !     optional local segment-box check
        !     compute 2D distance to current polygon segment
        !-------------------------------------------------------------------------

        do iz = zmin_seg, zmax_seg
            z = xx0(3) + dble(iz)*ddx(3)
            do iy = ymin_seg, ymax_seg
                y = xx0(2) + dble(iy)*ddx(2)
                do ix = xmin_seg, xmax_seg
                    x = xx0(1) + dble(ix)*ddx(1)

                    ! global coordinates
                    x_g = (/x,y,z/)

                    ! global -> body
                    x_b = matmul(M_g2b, x_g - Insect%xc_body_g)

                    ! body -> wing
                    x_w = matmul(M_b2w, x_b - x_pivot_b)

                    ! Extra local check.
                    ! The global bounding box is only a coarse box around a rotated box.
                    ! This check removes points that are inside the global box
                    ! but outside the actual local segment box
                    if (x_w(1) < segmentbox_w(1)) cycle
                    if (x_w(1) > segmentbox_w(2)) cycle
                    if (x_w(2) < segmentbox_w(3)) cycle
                    if (x_w(2) > segmentbox_w(4)) cycle
                    if (x_w(3) < segmentbox_w(5)) cycle
                    if (x_w(3) > segmentbox_w(6)) cycle

                    ! 2D distance in the wing plane.
                    ! Only x_w/y_w are used here.
                    d_seg = point_segment_distance_2D(x_w(1:2), p1, p2)
                    if (d_seg <= band_width) then
                        if (d_seg < tmp_dist_xy(ix,iy,iz)) then
                            tmp_dist_xy(ix,iy,iz) = d_seg
                            tmp_active(ix,iy,iz)  = .true.
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo  ! loop over polygon segments

    !-----------------------------------------------------------------------------
    ! Final SDF evaluation.
    ! - determine inside/outside in the 2D polygon
    ! - assign the sign of phi_xy
    ! - compute phi_z from the wing thickness
    ! - combine both to the 3D signed distance of the extruded polygon
    ! - convert phi to mask value
    !-----------------------------------------------------------------------------
    zmin_wing = Insect%Wings(wingID)%wing_bounding_box(5)
    zmax_wing = Insect%Wings(wingID)%wing_bounding_box(6)

    do iz = zmin, zmax
        z = xx0(3) + dble(iz)*ddx(3)
        do iy = ymin, ymax
            y = xx0(2) + dble(iy)*ddx(2)
            do ix = xmin, xmax
                x = xx0(1) + dble(ix)*ddx(1)

                ! global coordinates
                x_g = (/x,y,z/)
                ! global -> body
                x_b = matmul(M_g2b, x_g - Insect%xc_body_g)
                ! body -> wing
                x_w = matmul(M_b2w, x_b - x_pivot_b)

                ! inside the global box but outside the true local wing box.
                if (x_w(1) < wingbox_w(1)) cycle
                if (x_w(1) > wingbox_w(2)) cycle
                if (x_w(2) < wingbox_w(3)) cycle
                if (x_w(2) > wingbox_w(4)) cycle
                if (x_w(3) < wingbox_w(5)) cycle
                if (x_w(3) > wingbox_w(6)) cycle

                ! Check whether the point is inside the 2D polygon in the wing plane
                inside_polygon = point_in_polygon_2D(x_w(1:2), Insect%Wings(wingID)%polygon_wings, n)

                ! Case 1:
                ! tmp_active = true
                ! -> point was close enough to at least one polygon segment
                ! -> use computed distance
                !
                ! Case 2:
                ! tmp_active = false and point is inside polygon
                ! -> point is deep inside the polygon, away from all edges
                ! -> exact distance is not needed for truncated SDF
                ! -> use phi_xy = -band_width
                !
                ! Case 3:
                ! tmp_active = false and point is outside polygon
                ! -> point is outside and away from edges
                ! -> no mask contribution

                if (tmp_active(ix,iy,iz)) then

                    if (inside_polygon) then
                        phi_xy = -tmp_dist_xy(ix,iy,iz)
                    else
                        phi_xy =  tmp_dist_xy(ix,iy,iz)
                    endif

                else

                    if (inside_polygon) then
                        phi_xy = -band_width
                    else
                        cycle
                    endif

                endif

                ! assuming constant thickness here - incase its not maybe loop
                ! phi_z <= 0 inside the thickness interval
                ! phi_z >  0 outside the thickness interval
                phi_z = max(zmin_wing - x_w(3), x_w(3) - zmax_wing)

                ! mengenformel um echten Abstand zu bekommen
                phi = min(max(phi_xy, phi_z), 0.0_rk) + sqrt(max(phi_xy, 0.0_rk)**2 + max(phi_z, 0.0_rk)**2)
                phi = max(-band_width, min(band_width, phi))


                ! Outside the narrow band: no contribution needed
                if (phi > band_width) cycle

                !----------------------------------------------------------
                ! call step and update mask, us
                !----------------------------------------------------------
                mask_value = step(phi, 0.0_rk, Insect%L_smooth, Insect%safety,  Insect%smoothing_type_int)

                ! update only if the new mask value is higher than the one before
                if ((mask_value > mask(ix,iy,iz)) .and. (mask_value>0.0_rk)) then
                    mask(ix,iy,iz) = mask_value
                    mask_color(ix,iy,iz) = color_wing

                    v_tmp(1) = rot_rel_wing_w(2)*x_w(3)-rot_rel_wing_w(3)*x_w(2)
                    v_tmp(2) = rot_rel_wing_w(3)*x_w(1)-rot_rel_wing_w(1)*x_w(3)
                    v_tmp(3) = rot_rel_wing_w(1)*x_w(2)-rot_rel_wing_w(2)*x_w(1)

                    us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
                endif

            enddo
        enddo
    enddo

end subroutine draw_wing_polygon


function point_segment_distance_2D(x,p1,p2) result(d)

    ! x   - point in wing system (x_w, y_w)
    ! p1  - starting point of line segment
    ! p2  - end point of line segment
    ! d   - shortest distance from x to line segment p1--p2

    implicit none
    real(kind=rk), intent(in) :: x(1:2), p1(1:2), p2(1:2)
    real(kind=rk) :: d
    real(kind=rk) :: v(1:2), w(1:2), closest(1:2)
    real(kind=rk) :: vv, t

    ! p1 ---------> p2
    !       v
    v = p2 - p1
    w = x - p1

    vv = v(1)**2+v(2)**2

    ! are our segments far enough apart?
    if (vv <= 1.0e-14_rk) then
        d = sqrt(sum((x-p1)*(x-p1)))
        return
    endif

    ! Projection, restricted to finite line segment 
    ! (would be infinite line if t<0 & t>1)
    ! t = 0: projection on p1
    ! t = 1: projection on p2
    ! 0<t<1: projection between p1 and p2
    t = sum(w*v) / vv
    t = max(0.0_rk, min(1.0_rk, t))

    ! closest point on line segment
    closest = p1 + t*v

    ! finally distance
    d = sqrt(sum((x-closest)*(x-closest)))

end function point_segment_distance_2D


function point_in_polygon_2D(x, polygon, n) result(inside)
    implicit none

    real(kind=rk), intent(in) :: x(1:2), polygon(:,:)
    integer(kind=ik), intent(in) :: n
    logical :: inside
    integer(kind=ik) :: i, j
    real(kind=rk) :: xi, yi, xj, yj
    real(kind=rk) :: x_intersect

    inside = .false.
    j = n

    do i = 1, n

        xi = polygon(i,1)
        yi = polygon(i,2)

        xj = polygon(j,1)
        yj = polygon(j,2)

        ! ray-crossing test
        ! we look in positive x-direction from point x.
        ! if it crosses the polygon 2 times - outside
        ! if it crosses the polygon 1 time  - inside        
        if ( (yi > x(2)) .neqv. (yj > x(2)) ) then

            ! where does my hoizontal line cut the polygon segment?
            x_intersect = xi + (x(2)-yi) * (xj-xi) / (yj-yi)
            if (x(1) < x_intersect) then
                inside = .not. inside
            endif

        endif
        j = i

    enddo

end function point_in_polygon_2D





!-------------------------------------------------------------------------------
! evaluates the fourier series given in the ai, bi
! NOTE: the angle is [0, 2*pi] (thus the result of ATAN2 is + pi)
!-------------------------------------------------------------------------------
real(kind=rk) function Radius_Fourier( theta, Insect, wingID )
    implicit none
    integer(kind=ik) :: i,j, n_radius, j1, j2
    real(kind=rk) :: R0, theta2, dphi, area
    type(diptera),intent(inout)::Insect
    real(kind=rk), intent(in) :: theta
    integer(kind=2), intent(in) :: wingID ! wing id number
    integer(kind=ik) :: nfft_shape


    n_radius = 25000
    dphi = (2.0_rk*pi) / (dble(n_radius-1))


    !-------------------------------------------------------------------------------------
    ! evaluate the entire R(theta) once with very fine resolution, so when
    ! calling it for the second time we only need linear interpolation.
    !
    ! NOTE: this setup here is a Fourier initialization. It is also possible that the wing contour
    ! is described using a set {theta, R(theta)} with linear interpolation. In this case, the Insect%Wings(wingID)%R0_table
    ! is filled in Setup_WingShape_from_inifile, and this initialization is bypassed
    if ( .not. Insect%Wings(wingID)%wings_radius_table_ready) then
        if (.not. allocated(Insect%Wings(wingID)%R0_table) ) then
            allocate( Insect%Wings(wingID)%R0_table(1:n_radius) )
        endif

        ! initialization
        Insect%Wings(wingID)%R0_table = -9.0e9_rk

        ! number of fourier coefficients used to describe shape
        nfft_shape = size( Insect%Wings(wingID)%ai_shape, 1)

        ! loop over all thetas and compute the radius for all of them, store it
        ! in the table Insect%Wings(wingID)%R0_table
        do j = 1, n_radius
            ! zeroth mode (note unfortunate historic oddity of dividing by two..)
            R0 = Insect%Wings(wingID)%a0_shape / 2.0_rk

            theta2 = real(j-1, kind=rk) * dphi

            ! evaluate Fourier series
            do i = 1, nfft_shape
                R0 = R0 + Insect%Wings(wingID)%ai_shape(i) * cos( real(i,kind=rk)*theta2 ) &
                        + Insect%Wings(wingID)%bi_shape(i) * sin( real(i,kind=rk)*theta2 )
            enddo
            ! store the result in the table
            Insect%Wings(wingID)%R0_table(j) = R0
        enddo
        ! call this setup only once.
        ! Note: was not merged into Setup_WingShape because the hard-coded Fouier coefficients make that complicated
        Insect%Wings(wingID)%wings_radius_table_ready = .true.
        ! for debugging
        if (root) write(*,*) "Radius Fourier: pre-computation done. ", maxval(Insect%Wings(wingID)%R0_table), minval(Insect%Wings(wingID)%R0_table)
    endif

    !-------------------------------------------------------------------------------------
    ! linear interpolation, if already preprocessed the R(theta) table - that is faster than 
    ! evaluating the Fourier series.
    j1 = floor( theta / dphi ) + 1 ! note one-based indexing
    j2 = j1 + 1
    if (j2 > n_radius) j2 = 1 ! periodization
    Radius_Fourier = Insect%Wings(wingID)%R0_table(j1) + ((theta-real(j1-1,kind=rk)*dphi) / dphi) &
    * (Insect%Wings(wingID)%R0_table(j2) - Insect%Wings(wingID)%R0_table(j1))
end function


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
subroutine Setup_WingShape(Insect, wingID)
  implicit none
  real(kind=rk) :: xroot, yroot
  type(diptera),intent(inout)::Insect
  integer(kind=2), intent(in) :: wingID ! wing id number
  character(len=clong) :: wingshape_str
  integer(kind=ik) :: nfft_shape

  !-----------------------------------------
  ! hard-coded Fourier coefficients for R(theta)
  !-----------------------------------------
  select case (Insect%Wings(wingID)%WingShape)
  case ('pieris-brassicae1')
    ! butterfly (p. brassicae) wingshape, extracted from frame #633 of flight
    ! recording "160909_flight38"
    ! we have 29 points on the wing
    nfft_shape = 30
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.8262378
    Insect%Wings(wingID)%ai_shape = &
    (/0.0333388,0.0088668,0.0910448,0.0051568,-0.0252108,-0.0125628,&
    -0.0080008,-0.0054058,0.0110598,0.0090758,-0.0076998,-0.0071098,&
    0.0018768,0.0003048,0.0033128,0.0049088,-0.0021308,-0.0051658,&
    -0.0003278,0.0022188,0.0015628,0.0008548,-0.0009428,-0.0027778,&
    -0.0005728,0.0021988,0.0010718,-0.0007128,-0.0010658,-0.0008818 /)
    Insect%Wings(wingID)%bi_shape = &
    (/-0.0974478,0.0697448,0.0277868,-0.0461138,-0.0059608,0.0080508,&
    -0.0059118,0.0113048,0.0038478,-0.0157198,-0.0070308,0.0034218,&
    0.0030868,0.0040238,0.0017568,-0.0027368,-0.0050618,-0.0008578,&
    0.0033778,0.0009348,-0.0012898,-0.0011198,-0.0016878,0.0000398,&
    0.0029238,0.0009108,-0.0020888,-0.0012848,-0.0003728,0.0004738  /)
    Insect%Wings(wingID)%yc = 0.3184928
    Insect%Wings(wingID)%xc = -0.2459908
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition

  case ('drosophila')
    !********************************************
    ! Drosophila wing from Jan Gruber's png file
    !********************************************
    nfft_shape = 40
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.5140278
    Insect%Wings(wingID)%ai_shape = &
        (/0.1276258,-0.1189758,-0.0389458,0.0525938,0.0151538,-0.0247938,&
          -0.0039188,0.0104848,-0.0030638,-0.0064578,0.0042208,0.0043248,&
          -0.0026878,-0.0021458,0.0017688,0.0006398,-0.0013538,-0.0002038,&
          0.0009738,0.0002508,-0.0003548,-0.0003668,-0.0002798,0.0000568,&
          0.0003358,0.0001408,-0.0002208,0.0000028,0.0004348,0.0001218,&
          -0.0006458,-0.0003498,0.0007168,0.0003288,-0.0007078,-0.0001368,&
          0.0007828,0.0001458,-0.0007078,-0.0001358/)
    Insect%Wings(wingID)%bi_shape = &
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
    Insect%Wings(wingID)%xc =-0.1206 + xroot
    Insect%Wings(wingID)%yc = 0.3619 + yroot
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('drosophila_mutated')
    !********************************************
    ! mutated Drosophila wing from Jan Gruber's png file
    !********************************************
    nfft_shape = 70
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.4812548
    Insect%Wings(wingID)%ai_shape = &
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

    Insect%Wings(wingID)%bi_shape = &
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
    Insect%Wings(wingID)%xc =-0.1206 + xroot
    Insect%Wings(wingID)%yc = 0.3619 + yroot
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('drosophila_sandberg')
    !********************************************
    !  Drosophila wing from Ramamurti & Sandberg ( JEB 210, 881-896, 2007)
    !********************************************
    nfft_shape = 24
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.4995578
    Insect%Wings(wingID)%ai_shape = &
    (/0.0164168,-0.1621518,0.0030938,0.0601108,-0.0083988,-0.0199988,&
    0.0049048,0.0047878,-0.0005648,-0.0001108,-0.0008638,-0.0006928,&
    0.0006608,0.0001978,0.0001558,0.0006878,-0.0007498,-0.0008018,&
    0.0003878,0.0007028,0.0000408,-0.0001108,-0.0001068,-0.0003958 &
    /)
    Insect%Wings(wingID)%bi_shape = &
    (/-0.2083518,-0.0106488,0.0878308,-0.0018168,-0.0338278,0.0045768,&
    0.0113778,-0.0020678,-0.0026928,0.0002758,-0.0000838,-0.0001298,&
    0.0004118,0.0005638,-0.0001018,-0.0006918,-0.0002268,0.0005238,&
    0.0004008,-0.0001818,-0.0003038,-0.0000068,-0.0001218,0.0002008 &
    /)
    Insect%Wings(wingID)%xc =-0.0235498
    Insect%Wings(wingID)%yc = 0.1531398
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('drosophila_maeda')
    !********************************************
    !  Drosophila wing from Maeda and Liu, similar to Liu and Aono, BB2009
    !********************************************
    nfft_shape = 25
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.585432698694358
    Insect%Wings(wingID)%ai_shape = &
    (/0.113400475583443, -0.0862823485047213, -0.0346234482214816,&
    0.0237625254732323,0.00902498439287132,-0.0158926757445186,&
    -0.00549384372979449,0.0114928668063701,0.00431222381497978,&
    -0.00951270119733201,-0.00484045133879639,0.00706223174320460,&
    0.00473736389439926,-0.00449539769983697,-0.00418487169011745,&
    0.00320520052884641,0.00355631891057573,-0.00183155403463614,&
    -0.00191680264797099,0.00144768631289857,0.00135580122365068,&
    -0.000579638217642394,-0.000818378434108882,0.000132570375969864,&
    0.000683325977327827/)
    Insect%Wings(wingID)%bi_shape = &
     (/0.0939265226506824,0.0486063180327962,-0.0206591129298861,&
    -0.0136085709392758,0.0118575265347540,0.00604510770670991,&
    -0.0110263907282936,-0.00636979352727611,0.00786779216718321,&
    0.00390493804324433,-0.00797763174406198,-0.00450123591642554,&
    0.00445099872769504,0.00387237248613979,-0.00305464314668877,&
    -0.00398381251524846,0.00144450353105449,0.00257445316700965,&
    -0.00104247508055041,-0.00167946127380679,0.000577428923826108,&
    0.00114016779684690,-2.63209684213992e-05,-0.000753899380930065,&
    -0.000294894042986087/)

    !Insect%Wings(wingID)%xc = 0.0 ! original mesh
    Insect%Wings(wingID)%xc = 0.0473 ! shifted towards t.e. to 1/4 of the root chord ("+" sign here)
    !Insect%Wings(wingID)%xc = -0.0728 ! shifted towards l.e., to 0.2cmean from the l.e. (Liu and Aono BB 2009)
    !Insect%Wings(wingID)%yc = 0.7
    !Insect%Wings(wingID)%yc = 0.712 ! measured using kinematics snapshots
    Insect%Wings(wingID)%yc = 0.702 ! According to Maeda's email, Jun 21, 2014
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('drosophila_sun')
    !********************************************
    !  Drosophila virilis wing from Chen and Sun, Acta Mech Sin 2014
    !********************************************
    nfft_shape = 25
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.5427789795180327
    Insect%Wings(wingID)%ai_shape = &
    (/0.10879717599747996, -0.11445383382313232, -0.02023898255134523,&
    0.04903268079573884, 0.0012813019346402, -0.02397317767942499,&
    0.0013575396713610029, 0.0108149787395804, -0.001514114743464855,&
    -0.005364275911068656, 3.6505751634048205E-4, 0.002640180169907162,&
    -3.2673259786225535E-4, -0.0014323857426683313, 2.431115176929324E-4,&
    5.392319229992534E-4, -4.5833334881866856E-4, -1.3216432233072333E-4,&
    6.563502263270568E-4, 2.0750829321817808E-4, -4.807960800434886E-4,&
    -2.9006261005712504E-4, 2.7578746591965946E-4, 2.7519915193569374E-4,&
    -3.0570604954113513E-4/)
    Insect%Wings(wingID)%bi_shape = &
     (/-0.09385487982296374, -0.010821846776797858, 0.030052970821579587,&
    0.005312859230387492, -0.006054695188204192, -0.0015421303479118118,&
    -0.002533264559802815, -0.0014806147599133366, 0.003640199794653037,&
    0.0020416212413267134, -0.0024948946435721206, -7.83017244422372E-4,&
    0.0021574389122894035, 2.6950667683726845E-4, -0.00131044444112179,&
    6.404390762251693E-5, 2.513250728448789E-4, -4.7634735716375334E-4,&
    -1.5949800516545527E-5, 5.001276053841919E-4, 8.445613796483002E-5,&
    -5.510759077970704E-4, -3.3722093938416713E-4, 3.524656540450335E-4,&
    2.9924999100355387E-4/)
    Insect%Wings(wingID)%xc = 0.0
    Insect%Wings(wingID)%yc = 0.399446382250523
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('bumblebee')
    !********************************************
    !  Bumblebee
    !  http://www.entomology.umn.edu/museum/links/coursefiles/JPEG%20images/Hymenoptera%20web%20jpeg/Bombus-wings.jpg
    !********************************************
    nfft_shape = 25
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.594557593733011_rk
    Insect%Wings(wingID)%ai_shape = &
     (/-0.0128037920989526_rk,-0.106777418654552_rk,0.0380851289321982_rk,&
    0.0330548081081197_rk,-0.0178496286355627_rk,-0.00328588543359649_rk,&
    0.0108246924336137_rk,-0.00489302388329943_rk,-0.00708808441172961_rk,&
    0.00518244772891516_rk,0.00445979960844562_rk,-0.000108072056165527_rk,&
    0.00204437845603716_rk,0.00147176382618797_rk,-0.00229559463098105_rk,&
    -0.000514633972391526_rk,0.00134150515430486_rk,-0.000149860228261824_rk,&
    9.01456938813568d-05,0.00150639712261487_rk,0.000914624010720407_rk,&
    -0.000737650894315551_rk,-0.000843849923321745_rk,-0.000354971670482499_rk,&
    -0.000382956472432449_rk/)
    Insect%Wings(wingID)%bi_shape = &
     (/-0.0158061138788171_rk,0.0308243584184200_rk,-0.00903330410923372_rk,&
    -0.0185758334697500_rk,-0.000924452934252486_rk,-0.00242101213359519_rk,&
    -0.00204549530064489_rk,0.00291468131401423_rk,-0.000140755032337495_rk,&
    -0.00135036427128534_rk,0.00141285439042451_rk,-0.000334215276598231_rk,&
    -0.00161521722061879_rk,-0.000164055684312904_rk,-0.000256278551727569_rk,&
    -0.000740258481681094_rk,0.000847498161852221_rk,0.00157442110960973_rk,&
    -0.000559835622451578_rk,-0.000617498559228280_rk,0.00115413452523474_rk,&
    0.000322564770099778_rk,-0.000917375185844477_rk,4.44819399488798d-05,&
    0.000710028654602170_rk/)
    Insect%Wings(wingID)%xc = -0.1_rk
    Insect%Wings(wingID)%yc = 0.501549263807117_rk
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('b_ignitus')
    !********************************************
    !  Bumblebee B. ignitus
    !  Digitized from images taken at Liu Lab
    !********************************************
    nfft_shape = 25
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.536472532931637_rk
    Insect%Wings(wingID)%ai_shape = &
    (/-0.0447167394708177_rk,-0.106357727795917_rk,0.0504418160417239_rk,&
    0.0217275689429364_rk,-0.0259085955164794_rk,0.00272535910748833_rk,&
    0.00925289824790763_rk,-0.00453010629382665_rk,-0.000726647749565597_rk,&
    0.00258280099999843_rk,-0.00193033529765617_rk,-0.00121090519402499_rk,&
    0.00149872968653121_rk,0.000716207684720514_rk,-0.000205317764190544_rk,&
    0.000120507537444963_rk,-0.000381477942805165_rk,-0.000364957961985063_rk,&
    -6.70598716926467d-05,0.000166365788794039_rk,0.000332993591840758_rk,&
    -0.000225912231239784_rk,-0.000554023819155716_rk,0.000352735383706648_rk,&
    0.000650085631908143_rk/)
    Insect%Wings(wingID)%bi_shape = &
     (/-0.0580660125663764_rk,0.0271775529659247_rk,0.0178916506228727_rk,&
    -0.0196983386855655_rk,-0.00865040473524334_rk,0.0112078637630294_rk,&
    0.00505882127179290_rk,-0.00516874871678530_rk,-0.000418585234573997_rk,&
    0.00248996756589669_rk,-0.00248081765717699_rk,-0.00165307115885468_rk,&
    0.00236884835642553_rk,0.000920860396041608_rk,-0.00160449459432319_rk,&
    7.96078949775159d-05,0.000716588388745441_rk,0.000306756717543478_rk,&
    0.000310638954298390_rk,-0.000523512353114016_rk,-0.000773372382092419_rk,&
    1.97258594500968d-05,0.000261943571939630_rk,0.000262003935722642_rk,&
    0.000278542046262820_rk/)
    Insect%Wings(wingID)%xc = -0.13_rk
    Insect%Wings(wingID)%yc = 0.434820393790595_rk
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('paratuposa_flatwing')
    !********************************************
    !  Paratuposa wing simplified as a flat plate
    !  Digitized from 3D model, Moscow University entomology lab
    !********************************************
    nfft_shape = 15
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.694542662069373_rk
    Insect%Wings(wingID)%ai_shape = &
    (/-0.134736163793655_rk,0.00530251847896251_rk,-0.0345113221312334_rk,&
    0.00564308389276391_rk,-0.0151286715792430_rk,0.00702004741152472_rk,&
    0.00144649560886655_rk,-0.00185566405410384_rk,-0.00275905041561011_rk,&
    0.00217239130911607_rk,-0.00106500370428430_rk,-0.000750733476326611_rk,&
    0.00280149738648434_rk,-0.00182306390466332_rk,-0.000432849087666278_rk/)
    Insect%Wings(wingID)%bi_shape = &
    (/0.0931417992571149_rk,0.0598684008118805_rk,-0.0297533040215975_rk,&
    0.00720759394160553_rk,-0.00591704536021243_rk,-0.0103773128578859_rk,&
    0.00480389428474622_rk,-0.000497969127742664_rk,-0.000267837077569406_rk,&
    -0.00102010912896721_rk,0.00218717420988224_rk,-0.00290492852428728_rk,&
    0.000860757518054172_rk,0.00147755983849289_rk,-0.000638118479966807_rk/)
    Insect%Wings(wingID)%xc = -0.3_rk
    Insect%Wings(wingID)%yc = 0.7_rk
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('paratuposa_flatelytra')
    !********************************************
    !  Paratuposa elytra simplified as a flat plate
    !  Digitized from 3D model, Moscow University entomology lab
    !********************************************
    nfft_shape = 21
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.565146320110115_rk * 0.62_rk
    Insect%Wings(wingID)%ai_shape = &
    (/0.0977588377084158_rk,-0.134895186363991_rk,-0.0415626398480136_rk,&
    0.0462811349514129_rk,0.0116178197112501_rk,-0.0173654394721538_rk,&
    -0.00143307308642128_rk,0.00463364466230196_rk,-0.00245823642621685_rk,&
    -0.000295909979786882_rk,0.00293439703102658_rk,-0.000757888513431659_rk,&
    -0.00232308360478113_rk,0.000590866168893430_rk,0.00134477744241136_rk,&
    -0.000179797491049942_rk,-0.000579366975013497_rk,0.000298632558506859_rk,&
    0.000121897474918553_rk,-0.000129056828603867_rk,6.65231631883479d-5/) &
    * 0.62_rk
    Insect%Wings(wingID)%bi_shape = &
    (/0.0873846736908491_rk,0.0261657542571480_rk,-0.0361196299337545_rk,&
    -0.0169805639649468_rk,0.0120703518387247_rk,0.00786957912240014_rk,&
    -0.00255959633597845_rk,-0.00262639690803298_rk,-0.000425744504894171_rk,&
    -0.000897607832134308_rk,0.00110940805947636_rk,0.00119153774474213_rk,&
    -0.00114698647381269_rk,-0.000899734725370745_rk,0.000630784528183092_rk,&
    0.000278761021018954_rk,-0.000290963412817785_rk,6.98199270677237d-5,&
    0.000195892295619984_rk,-0.000190195259814270_rk,-4.82888469073568d-5/) &
    * 0.62_rk
    Insect%Wings(wingID)%xc = -0.02_rk * 0.62_rk
    Insect%Wings(wingID)%yc = 0.65_rk * 0.62_rk
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('flapper_sane')
    !********************************************
    !  Mechanical model from Sane and Dickinson, JEB 205, 2002
    !  'The aerodynamic effects...'
    !********************************************
    nfft_shape = 25
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.5379588906565078
    Insect%Wings(wingID)%ai_shape = &
     (/0.135338653455782,-0.06793162622123261,-0.0398235167675977,&
    0.006442194893963269,0.0012783260416583853,-0.007014398516674715,&
    0.0017710765408983137,0.006401601802033519,-2.970619204124993E-4,&
    -0.0038483478773981405,-6.180958756568494E-4,8.015784831786756E-4,&
    -6.957513357109226E-4,-1.4028929172227943E-4,0.0013484885717868547,&
    4.827827498543977E-4,-9.747844462919694E-4,-5.838504331939134E-4,&
    2.72834004831554E-4,2.8152492682871664E-5,-1.2802199282558645E-4,&
    4.117887216124469E-4,3.364169982438278E-4,-3.33258003686823E-4,&
    -3.5615733035757616E-4/)
    Insect%Wings(wingID)%bi_shape = &
     (/2.686408368800394E-4,0.01649582345310688,0.01288513083639708,&
    0.004711436946785864,-0.0035725088809005073,-0.00898640397179334,&
    -0.003856509905612652,0.004536524572892801,0.004849677692836578,&
    2.9194421255236984E-4,-7.512780802871473E-4,7.12685261783966E-4,&
    -1.5519932673320404E-4,-0.0012695469974603026,2.2861692091158138E-4,&
    0.0016461316319681953,5.257476721137781E-4,-7.686482830046961E-4,&
    -3.108879176661735E-4,2.2437540206568518E-4,-2.578427217327782E-4,&
    -2.5120263516966855E-4,4.1693453021778877E-4,3.9290173948150096E-4,&
    -1.9762601237675826E-4/)
    Insect%Wings(wingID)%xc = 0.0
    Insect%Wings(wingID)%yc = 0.6
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('flapper_dickinsonII')
    !********************************************
    ! Digitized from Dickinson et al 1999 Science, figure 1A, drawing
    ! of the mechanical robot
    !********************************************
    nfft_shape = 20
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.6442788
    Insect%Wings(wingID)%ai_shape = &
     (/0.0482978,-0.1208378,0.0061008,0.0356718,-0.0148328,-0.0109958,&
    0.0110268,0.0018538,-0.0061998,0.0015458,0.0025508,-0.0017538,&
    -0.0002578,0.0015018,-0.0003158,-0.0006048,0.0007168,-0.0001568,&
    -0.0005018,0.0004118/)
    Insect%Wings(wingID)%bi_shape = &
     (/-0.0521708,0.0051828,0.0369428,-0.0002868,-0.0177448,0.0023218,&
    0.0081378,-0.0036288,-0.0038168,0.0031348,0.0011858,-0.0023828,&
    -0.0001638,0.0016098,-0.0004768,-0.0007188,0.0007228,0.0002278,&
    -0.0005798,0.0001228/)
    Insect%Wings(wingID)%yc = 0.5282438
    Insect%Wings(wingID)%xc = -0.1184548
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition
  case ('robofly_dickinson')
    !********************************************
    ! Digitized from the hand drawn figure M. Dickinson sent via email, which
    ! contained the exact location of the pivot point. He also sent a CAD drawing
    ! which looks slightly different, and had no pivot point marked.
    !********************************************
    nfft_shape = 28
    allocate( Insect%Wings(wingID)%ai_shape(1:nfft_shape) )
    allocate( Insect%Wings(wingID)%bi_shape(1:nfft_shape) )

    Insect%Wings(wingID)%a0_shape = 0.5313628
    Insect%Wings(wingID)%ai_shape = &
    (/-0.0245658,-0.0842918,0.0218028,0.0105418,-0.0095288,0.0012928,&
    0.0021928,0.0000328,-0.0007648,-0.0015808,0.0013808,0.0013068,&
    -0.0010748,0.0002408,-0.0000378,-0.0010888,0.0008248,0.0004708,&
    -0.0003988,0.0002658,-0.0003178,-0.0004218,0.0002768,0.0000818,&
    0.0000318,0.0001228,-0.0001918,-0.0000558/)
    Insect%Wings(wingID)%bi_shape = &
     (/-0.0905448,0.0278058,0.0392558,-0.0125248,-0.0159598,0.0048268,&
    0.0038898,-0.0028828,0.0012618,0.0012998,-0.0019058,0.0003118,&
    0.0003198,-0.0004298,0.0006388,-0.0000648,-0.0002308,0.0002518,&
    -0.0003948,0.0000928,0.0004478,-0.0003078,-0.0000888,0.0001638,&
    -0.0002348,0.0001398,0.0001398,-0.0002358/)
    Insect%Wings(wingID)%yc = 0.4645238
    Insect%Wings(wingID)%xc = -0.0716018
    Insect%Wings(wingID)%WingShapeType = "fourier"  ! for readability only; default set in type(diptera) definition

    case default

        ! if all other options fail, we still might load coefficients from file:
        wingshape_str = Insect%Wings(wingID)%WingShape
        if (index(wingshape_str,"from_file::") /= 0) then
            !-------------------------------------------------------------------------
            ! wing shape is read from ini-file
            !-------------------------------------------------------------------------
            call Setup_WingShape_from_inifile(Insect, wingID, trim(adjustl(wingshape_str( 12:len_trim(wingshape_str) ))))

        else
            ! now we theres an error...
            write (*,*) "Insect module: trying to set up fourier descriptors for wing&
            & shape but the type Insect%WingShape is unknown! :: "// Insect%Wings(wingID)%WingShape

            call abort(554329, "Insect module: trying to set up fourier descriptors for wing&
            & shape but the type Insect%WingShape is unknown!")
        end if
    end select

  ! for many cases, it is important that Lspan and Lchord are known, but that is
  ! tedious for Fourier shapes, as the use cannot see it from the coefficients.
  ! Therefore, we compute the max / min of x / y here and store the result
  call set_wing_bounding_box( Insect, wingID )

  ! this is the old default value:
  if (maxval(Insect%Wings(wingID)%corrugation_array_bbox(:)) == 0.0_rk) then
      Insect%Wings(wingID)%corrugation_array_bbox = Insect%Wings(wingID)%wing_bounding_box(1:4)
  endif

  if (root) then
    write(*,'(30("-"))')
    write(*,'("Insect module: Setup_WingShape is done.")')
    write(*,'("Wing shape is ",A)') trim(adjustl(Insect%Wings(wingID)%WingShape))
    write(*,'(30("-"))')
  endif

end subroutine Setup_WingShape



!-------------------------------------------------------------------------------
! Initialize the wing from an ini file
!-------------------------------------------------------------------------------
! We read:
!     - Fourier coefficients for the radius wing shape
!     - Wing thickness profile (constant or variable)
!     - Wing corrugation profile (flat or corrugated)
!-------------------------------------------------------------------------------
subroutine Setup_WingShape_from_inifile( Insect, wingID, fname )
    implicit none
    type(diptera),intent(inout) :: Insect
    character(len=*), intent(in) :: fname

    type(inifile) :: ifile
    real(kind=rk), allocatable :: tmparray(:,:)
    character(len=clong) :: type_str
    integer :: a,b
    integer(kind=2), intent(in) :: wingID ! wing id number
    real(kind=rk) :: init_thickness, dphi, theta2
    real(kind=rk), dimension(:), allocatable :: theta_i, R_i
    integer(kind=ik) :: i, j, n_radius, nfft_shape
    !KVN-2025>>>>>
    integer(kind=ik) :: a1,b1,c
    !KVN-2025<<<<<

    if (root) then
        write(*,'(80("─"))')
        write(*,'("Reading wing shape for wingID=",i1," from file ",A)') wingID, fname
        write(*,'(80("─"))')
    endif

    ! instead of the hard-coded values above, read fourier coefficients for wings from
    ! an ini-file
    call read_ini_file_mpi(ifile, fname, .true.  )

    ! check if this file seem to be valid (i.e. if TYPE matches a wing that can actually
    ! be initialized from file)
    call read_param_mpi(ifile, "Wing", "type", Insect%Wings(wingID)%WingShapeType, "none")

    ! fourier:
    !       the wing shape is described in polar coordinates and the radius is encoded as fourier coefficients
    !       fourier coeffs are read as a0_wings, ai_wings, bi_wings
    !       note historic oddity that a0 is half the mean value (the 0th Fourier mode)
    !       T. Engels, D. Kolomenskiy, K. Schneider and J. Sesterhenn. FluSI: A novel parallel simulation tool for flapping insect flight using a Fourier method with volume penalization. SIAM J. Sci. Comp., 38(5), S03-S24, 2016
    ! linear:
    !       The wing contour is described in polar coordinates, just like in the "fourier" case, but the R(theta) is
    !       included as a table, not as Fourier coefficients. This is useful if the wing contains sharp edges, where the
    !       Fourier series converges badly (Gibbs ringing).
    ! polygon:
    !       the wing contour is described in cartesian coordinates of points P(:,:) given in the wing-system and can be arbitrarily shaped 

    select case(Insect%Wings(wingID)%WingShapeType)
    case ("linear")
        !-----------------------------------------------------------------------------
        ! R(theta) given as a lits of points (for linear interpolation without Fourier series)
        !-----------------------------------------------------------------------------
        ! wing mid-point (of course in wing system..)
        ! used as origin of the polar coordinates R(theta) that describe the
        ! wing contour
        call read_param_mpi(ifile, "Wing", "x0w", Insect%Wings(wingID)%xc, 0.0_rk)
        call read_param_mpi(ifile, "Wing", "y0w", Insect%Wings(wingID)%yc, 0.0_rk)


        call param_matrix_size_mpi( ifile, "Wing", "theta_i", a, b)
        ! allocate matrix
        allocate( tmparray(1:a,1:b) )
        ! allocate local work arrays
        ! the insect_module only uses R0_table in the end, one for each wing.
        allocate( theta_i(1:b) )
        allocate( R_i(1:b) )

        call param_matrix_read_mpi( ifile, "Wing", "theta_i", tmparray)
        theta_i(:) = tmparray(1,:)

        call param_matrix_read_mpi( ifile, "Wing", "R_i", tmparray)
        R_i(:) = tmparray(1,:)


        ! wing contour given as a set of {theta, R(theta)} points. Note: the line root-tip (the span) is the Y-axis,
        ! and the x coordinate runs from trailing (negative) to leading edge (positive). This is important for the
        ! definition of theta.

        ! we need equidistant data, possibly fine-sampled. For sharp-edged features, we need to bypass
        ! the Fourier series (it works not so well for such functions: Gibbs ringing...)

        ! if N<25000 we can store it in Insect%Wings(wingID)%R0_table(j,wingID) = R0
        ! upsampling is required: we cannot read so many values from ini file (max_column_width)

        ! fill the table Insect%Wings(wingID)%R0_table, this is the same as R_i but it has 25000
        ! entries and is used in the subroutine "Radius_Fourier"

        allocate( Insect%Wings(wingID)%R0_table(1:25000) )
        Insect%Wings(wingID)%R0_table(:) = 0.0_rk

        n_radius = size( Insect%Wings(wingID)%R0_table, dim=1 )

        ! loop over all thetas and compute the radius for all of them, store it
        ! in the table Insect%Wings(wingID)%R0_table
        do j = 1, n_radius
            ! in the target array (fine spacing with 25000 entries)
            dphi   = (2.0_rk*pi) / (real(n_radius-1, kind=rk))
            theta2 = real(j-1, kind=rk) * dphi

            ! In the source array (read from file), we do not even assume equidistant data.
            ! The way we do that here is not super elegant, but it is done only once anyways.
            do i = 2, size(theta_i, 1)
                if ((theta_i(i) > theta2).or.(i==size(theta_i, 1))) then
                    ! now i is the index of the first element larger than the value we look for
                    ! so we interpolate between (i-1, i)
                    ! R_j = R_i-1 + (th_j - th_i-1) / (th_i - th_i-1) * (R_i - R_i-1)
                    Insect%Wings(wingID)%R0_table(j) = R_i(i-1) + &
                    (R_i(i)-R_i(i-1)) * (theta2-theta_i(i-1)) / (theta_i(i)-theta_i(i-1))
                    ! end of inner for-loop
                    exit
                endif
            enddo
        enddo

        ! done. This flag bypasses the Fourier-series initialization in Radius_Fourier
        Insect%Wings(wingID)%wings_radius_table_ready = .true.

        ! deallocate work arrays
        deallocate(R_i, theta_i)

    case("fourier", "fourierY")
        !-----------------------------------------------------------------------------
        ! Read fourier coeffs for wing radius ("fourier") or y-coordinate of membrane (used for bristled wings)
        !-----------------------------------------------------------------------------
        call read_param_mpi( ifile, "Wing", "a0_wings", Insect%Wings(wingID)%a0_shape, 0.0_rk)

        ! NOTE: Annoyingly, the fujitsu SXF90 compiler cannot handle allocatable arrays
        ! as arguments. so we have to split the routine in one part that returns the size
        ! of the array, then let the caller allocate, then read the matrix. very tedious.
        ! fetch size of matrix
        call param_matrix_size_mpi( ifile, "Wing", "ai_wings", a, b)
        ! allocate matrix
        allocate( tmparray(1:a,1:b) )
        ! read matrix
        call param_matrix_read_mpi( ifile, "Wing", "ai_wings", tmparray)

        nfft_shape = size(tmparray,2)
        allocate( Insect%wings(wingID)%ai_shape(1:nfft_shape) )
        Insect%Wings(wingID)%ai_shape = tmparray(1,:)
        deallocate(tmparray)


        call param_matrix_size_mpi( ifile, "Wing", "bi_wings", a, b)
        ! allocate matrix
        allocate( tmparray(1:a,1:b) )
        ! read matrix
        call param_matrix_read_mpi( ifile, "Wing", "bi_wings", tmparray)

        nfft_shape = size(tmparray,2)
        allocate( Insect%wings(wingID)%bi_shape(1:nfft_shape) )
        Insect%Wings(wingID)%bi_shape = tmparray(1,:)
        deallocate(tmparray)

        ! wing mid-point (of course in wing system..)
        ! used as origin of the polar coordinates R(theta) that describe the
        ! wing contour
        call read_param_mpi(ifile,"Wing","x0w",Insect%Wings(wingID)%xc, 0.0_rk)
        call read_param_mpi(ifile,"Wing","y0w",Insect%Wings(wingID)%yc, 0.0_rk)

        if (root) then
            write(*,*) "wingID", wingID
            write(*,*) "ai_shape", Insect%Wings(wingID)%ai_shape
            write(*,*) "bi_shape", Insect%Wings(wingID)%bi_shape
        endif

    case("polygon")
        !-----------------------------------------------------------------------------
        ! Read points P(:,:), we assume: P(1,:) -> P(2,:) -> ... P(n,:) -> P(1,:)
        !-----------------------------------------------------------------------------
        ! fetch size of matrix
        call param_matrix_size_mpi( ifile, "Wing", "polygon_points", a,b)
        if (a < 3) then
            call abort(17062601, "think again :) - polygon_points must contain at least three points")
        endif
        if (b/=2) then
            call abort(17062602, "polygon_points must have exactly two columns for coordinates (x,y)")
        endif

        ! allocate matrix
        allocate( tmparray(1:a,1:b) )
        call param_matrix_read_mpi(ifile, "Wing", "polygon_points", tmparray)

        ! for now, I assume the wings are completly symmetric (in the z direction, top/down), 
        ! later maybe we can add a second variable for e.g. asymmetric wing damage (corrugation)
        if (.not. allocated(Insect%wings(wingID)%polygon_wings)) then
            allocate(Insect%wings(wingID)%polygon_wings(1:a, 1:2))
        endif

        Insect%wings(wingID)%polygon_wings(:, :) = tmparray(:, :)
        Insect%wings(wingID)%n_polygon_points = a
        
        deallocate(tmparray)

        if (root) then
            write(*,'("Polygon wing shape from ini-file: read ",i0," points from polygon_points")') a
        endif

    case default
        call abort(6652, "ini file for wing does not seem to be correct type...")

    end select

    !-----------------------------------------------------------------------------
    ! wing thickness
    !-----------------------------------------------------------------------------
    call read_param_mpi(ifile,"Wing","wing_thickness_distribution",Insect%wings(wingID)%wing_thickness_distribution, "constant")

    if ( Insect%wings(wingID)%wing_thickness_distribution == "constant") then

        if (root) write(*,*) "Wing thickness is constant along the wing"

        ! wing thickness (NOTE: overwrites settings in other params file)
        if ( (Insect%WingThickness > 0.0_rk) .and. (Insect%WingThickness < 1.0_rk) ) then
            ! Use existing value if it is reasonable
            init_thickness = Insect%WingThickness 
        else
            ! This is the default value otherwise, because we may not know dx here
            init_thickness = 0.05_rk 
        endif

        call read_param_mpi(ifile, "Wing", "wing_thickness_value", Insect%wings(wingID)%WingThickness, init_thickness)

    elseif ( Insect%wings(wingID)%wing_thickness_distribution == "variable") then

        if (root) write(*,*) "Wing thickness is variable, i.e. t = t(x,y)"

        ! read matrix from ini file, see comments on SXF90 compiler
        call param_matrix_size_mpi(ifile,"Wing","wing_thickness_profile",a,b)

        allocate( Insect%wings(wingID)%wing_thickness_profile(1:a, 1:b) )
        
        call param_matrix_read_mpi(ifile, "Wing", "wing_thickness_profile", Insect%wings(wingID)%wing_thickness_profile)

    else
        call abort(77623, " Insect wing thickness distribution is unknown (must be constant or variable)")
    endif


    !-----------------------------------------------------------------------------
    ! bristles
    !-----------------------------------------------------------------------------
    call read_param_mpi(ifile, "Wing","bristles", Insect%Wings(wingID)%bristles, .false.)
    call read_param_mpi(ifile, "Wing","bristles_simplex", Insect%Wings(wingID)%bristles_simplex, .false.)

    if (Insect%Wings(wingID)%bristles) then
        call param_matrix_size_mpi( ifile, "Wing", "bristles_coords", a, b)

        ! number of bristles on this wing
        Insect%Wings(wingID)%n_bristles = a

        allocate(Insect%Wings(wingID)%bristles_coords(1:a, 1:b))

        call param_matrix_read_mpi( ifile, "Wing", "bristles_coords", Insect%Wings(wingID)%bristles_coords(1:a, 1:b))
    endif

    !KVN-2025>>>>>
    !-----------------------------------------------------------------------------
    ! 3D-bristles
    !-----------------------------------------------------------------------------
    call read_param_mpi(ifile, "Wing", "bristles3D", Insect%Wings(wingID)%bristles3D, .false.)
    
    !-----------------------------------------------------------------------------
    ! wing deformation
    !-----------------------------------------------------------------------------
    ! NOTE: prescribed wing deformation is untested work in progress! -TE 02/2026
    ! call read_param_mpi(ifile,"Wing","deformable",Insect%deformable(wingID), .false.)
    ! if (Insect%deformable(wingID)) then
    !     if (root) write(*,*) "wing is deformable z=z(x,y,t)"
    !     call read_param_mpi(ifile,"Wing","deformation_array_bbox",Insect%deformation_array_bbox(1:4,wingID), (/0.0_rk,0.0_rk,0.0_rk,0.0_rk/))
    !     call param_matrix_size_mpi(ifile,"Wing","deformations",a1,b1)
    !     call Allocate_Arrays(Insect,"deformations",a1,b1)
    !     call param_matrix_read_mpi(ifile,"Wing","deformations",Insect%deformations(:,:,wingID))
    !     do a = 1, a1
    !       if ( Insect%deformations(a,1,wingID) /= Insect%deformations(a+1,1,wingID) ) then
    !         exit
    !       end if
    !     end do
    !     b = b1 - 1
    !     c = a1/a
    !     Insect%deformation_a(wingID) = a
    !     Insect%deformation_b(wingID) = b
    !     Insect%deformation_c(wingID) = c
    !     call Allocate_Arrays(Insect,"deformation_profile",a,b)
    ! else
    !     if (root) write(*,*) "wing is non-deformable"
    ! endif
    !KVN-2025<<<<<

    !-----------------------------------------------------------------------------
    ! wing corrugation
    !-----------------------------------------------------------------------------
    call read_param_mpi(ifile,"Wing","corrugated",Insect%Wings(wingID)%corrugated, .false.)

    if (Insect%Wings(wingID)%corrugated) then
        if (root) write(*,*) "wing is corrugated, z=z(x,y)"

        ! read matrix from ini file, see comments on SXF90 compiler
        call param_matrix_size_mpi(ifile,"Wing","corrugation_profile",a,b)
        
        allocate( Insect%Wings(wingID)%corrugation_profile( 1:a, 1:b ) )

        call param_matrix_read_mpi(ifile, "Wing", "corrugation_profile", Insect%Wings(wingID)%corrugation_profile(:,:))

    else
        if (root) write(*,*) "wing is flat (non-corrugated), z==0"
    endif

    !-----------------------------------------------------------------------------
    ! wing damage
    !-----------------------------------------------------------------------------
    call read_param_mpi(ifile, "Wing", "damaged", Insect%Wings(wingID)%damaged, .false.)

    if (Insect%Wings(wingID)%damaged) then
        if (root) write(*,*) "wing is damaged, D=D(x,y)"

        ! read matrix from ini file, see comments on SXF90 compiler
        call param_matrix_size_mpi(ifile, "Wing", "damage_mask", a, b)

        allocate( Insect%Wings(wingID)%damage_mask( 1:a, 1:b ) )

        call param_matrix_read_mpi(ifile, "Wing", "damage_mask", Insect%Wings(wingID)%damage_mask(:,:))

    else
        if (root) write(*,*) "wing is intact (non-damaged)"
    endif

    ! the 2D arrays damage_mask, corrugation_profile and wing_thickness_profile
    ! are equidistant grids, and their bounding box is all the same:
    call read_param_mpi(ifile, "Wing", "corrugation_array_bbox", Insect%Wings(wingID)%corrugation_array_bbox, (/0.0_rk,0.0_rk,0.0_rk,0.0_rk/))

end subroutine Setup_WingShape_from_inifile


!-------------------------------------------------------------------------------
! Setup extends of wing for reduction of computational time
!-------------------------------------------------------------------------------
! for many cases, it is important that Lspan and Lchord are known, but that is
! tedious for Fourier shapes, as the use cannot see it from the cooefficients.
! Therefore, we compute the max / min of x / y / z  here and store the result
! NOTE: This code is executed only once.
! 2026/06/30 EG: I added cases to include the polygon shaped wing and changed the name
! set_wing_bounding_box_fourier( Insect, wingID) -> set_wing_bounding_box( Insect, wingID)
!-------------------------------------------------------------------------------
subroutine set_wing_bounding_box( Insect, wingID)

    implicit none
    type(diptera),intent(inout) :: Insect
    real(kind=rk) :: theta, xmin,xmax, ymin, ymax, R, x, y, theta_prime, tmp
    integer(kind=2), intent(in) :: wingID ! wing id number
    integer(kind=ik) :: n

    theta = 0.0_rk
    xmin = 999.d9
    ymin = 999.d9
    xmax = -999.d9
    ymax = -999.d9

    if (root) write(*,*) "set_wing_bounding_box, wing_file_type=", trim(adjustl((Insect%Wings(wingID)%WingShapeType)))


    select case(Insect%Wings(wingID)%WingShapeType)
    case("fourier", "linear", "fourierY")
        ! construct the wing border by looping over the angle theta, look for smallest and largest x,y values
        ! note flusi uses an angle between [0, 2*pi)
        do while ( theta < 2.0_rk*pi )
            ! note this angle is [0, 2*pi)
            R = Radius_Fourier( theta, Insect, wingID )

            ! note how the usual atan2 gives angles [-pi, +pi)
            ! so here we add pi
            theta_prime = theta - pi
            x = Insect%Wings(wingID)%xc + R * cos( theta_prime )
            y = Insect%Wings(wingID)%yc + R * sin( theta_prime )

            ! NOTE: A word on theta_prime: it is the angle described with the positve x-axis
            ! and indeed rotates in positve z-direction. That means in the plane
            !
            !  ^ x_wing
            !  |
            !  |
            !  o-------> y_wing
            !
            ! It rotates CLOCKWISE, starting from the x_wing axis (zero is the x-axis, vertical)

            xmin = min( xmin, x )
            ymin = min( ymin, y )

            xmax = max( xmax, x )
            ymax = max( ymax, y )

            theta = theta + 1.0e-3_rk
        end do

        Insect%Wings(wingID)%wing_bounding_box(1:4) = (/xmin, xmax, ymin, ymax/)

    case("polygon")
        !         yw
        !         ^
        !         |
        !  ymax   |    +-----------------+
        !         |    |      P4 o       |
        !         |    |        / \      |
        !         |    |       /   o P3  |
        !         |    |  P5 o     |     |
        !         |    |     \     |     |
        !         |    |      \    o P2  |
        !         |    |       \  /      |
        !         |    |        o P1     |
        !  ymin   |    +-----------------+
        !         |
        !         o------------------------------> xw
        !            xmin               xmax     
        
        ! 2D box
        n = Insect%Wings(wingID)%n_polygon_points
        ! x = chord direction, y = span direction
        xmin = minval( Insect%Wings(wingID)%polygon_wings(1:n, 1) )
        xmax = maxval( Insect%Wings(wingID)%polygon_wings(1:n, 1) )
        ymin = minval( Insect%Wings(wingID)%polygon_wings(1:n, 2) )
        ymax = maxval( Insect%Wings(wingID)%polygon_wings(1:n, 2) )

        Insect%Wings(wingID)%wing_bounding_box(1:4) = (/xmin, xmax, ymin, ymax/)

    case default
        call abort(300626001, "ini file for wing does not seem to be correct type...")
    end select

    ! the bounding box in z-direction depends on the wing thicnkess (constant or not)
    ! and the corrugation
    if ( Insect%wings(wingID)%wing_thickness_distribution == "constant" ) then

        if ( Insect%Wings(wingID)%corrugated ) then
            Insect%Wings(wingID)%wing_bounding_box(5) = minval(Insect%Wings(wingID)%corrugation_profile) - Insect%Wings(wingID)%WingThickness / 2.0_rk
            Insect%Wings(wingID)%wing_bounding_box(6) = maxval(Insect%Wings(wingID)%corrugation_profile) + Insect%Wings(wingID)%WingThickness / 2.0_rk
        else
            ! constant thickness, no corrugation is the classical flat case:
            Insect%Wings(wingID)%wing_bounding_box(5) = -Insect%Wings(wingID)%WingThickness / 2.0_rk
            Insect%Wings(wingID)%wing_bounding_box(6) = +Insect%Wings(wingID)%WingThickness / 2.0_rk
        endif

    else

        if ( Insect%Wings(wingID)%corrugated ) then
            ! minimum of lower surface
            Insect%Wings(wingID)%wing_bounding_box(5) = minval(Insect%Wings(wingID)%corrugation_profile - Insect%Wings(wingID)%wing_thickness_profile/2.0_rk)
            ! maximum of upper surface
            Insect%Wings(wingID)%wing_bounding_box(6) = maxval(Insect%Wings(wingID)%corrugation_profile + Insect%Wings(wingID)%wing_thickness_profile/2.0_rk)
        else
            ! bounding box is +- largest thickness  simply
            Insect%Wings(wingID)%wing_bounding_box(5) = -maxval(Insect%Wings(wingID)%wing_thickness_profile / 2.0_rk)
            Insect%Wings(wingID)%wing_bounding_box(6) =  maxval(Insect%Wings(wingID)%wing_thickness_profile / 2.0_rk)
        endif

    endif

    ! Mirror the bounding box for left wings (because the z direction points in the opposite direction!)
    if ( (wingID == 1) .or. (wingID == 3) ) then
        tmp = Insect%Wings(wingID)%wing_bounding_box(6)
        Insect%Wings(wingID)%wing_bounding_box(6) = -Insect%Wings(wingID)%wing_bounding_box(5)
        Insect%Wings(wingID)%wing_bounding_box(5) = -tmp
    end if

    if (root) then
        write(*,'("Effective (=the real surface) wing lengths are:")')
        write(*,'("Lspan=",es15.8,"Lchord=",es15.8)') ymax-ymin, xmax-xmin
        write(*,'("Bounding box is:")')
        write(*,'("xwmin=",es15.8," xwmax=",es15.8)') Insect%Wings(wingID)%wing_bounding_box(1:2)
        write(*,'("ywmin=",es15.8," ywmax=",es15.8)') Insect%Wings(wingID)%wing_bounding_box(3:4)
        write(*,'("zwmin=",es15.8," zwmax=",es15.8)') Insect%Wings(wingID)%wing_bounding_box(5:6)
    endif
end subroutine set_wing_bounding_box




subroutine draw_bristle(x1w, x2w, R0, xx0, ddx, mask, mask_color, us, Insect, color_val, M_g2b, M_b2w, x_pivot_b, rot_rel_wing_w)
    implicit none

    real(kind=rk), dimension(1:3), intent(in):: x1w, x2w
    real(kind=rk),intent(in)::R0
    type(diptera),intent(inout)::Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: color_val
    real(kind=rk),intent(in) :: M_g2b(1:3,1:3), M_b2w(1:3,1:3), x_pivot_b(1:3), rot_rel_wing_w(1:3)

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
    x1 = matmul( transpose(M_b2w), x1w) + x_pivot_b
    x1 = matmul( transpose(M_g2b), x1) + Insect%xc_body_g

    x2 = matmul( transpose(M_b2w), x2w) + x_pivot_b
    x2 = matmul( transpose(M_g2b), x2) + Insect%xc_body_g

    ! Ideally, we would use this function form primitives-collection.
    ! But the bristles are slightly displaced (I don't know why, maybe the x0_indices)
    ! and this routine does not add the velocity field. As using the new function is a bigger operation
    ! I postpone it for now. After all, the old code works. -TE 30/jul/2026
    ! call draw_cylinder_rounded(mask, mask_color, xx0, ddx, g, x1, x2, R0, int(color_val, kind=ik), &
    ! Insect%smoothing_type_int, Insect%L_smooth, x0_indices=(/0,0,0/))

    !---------------------------------------------------------------------------
    ! define bristle (cylinder) coordinate system
    !---------------------------------------------------------------------------
    ! unit vector in bristle (cylinder) axis direction and bristle (cylinder) length
    e_x = x2 - x1
    clength = norm2_3d(e_x)
    e_x = e_x / clength

    ! radial unit vector
    ! use a vector perpendicular to e_x, since it is a azimuthal symmetry
    ! it does not really matter which one. however, we must be sure that the vector
    ! we use and the e_x vector are not colinear -- their cross product is the zero vector, if that is the case
    e_r = (/0.0_rk, 0.0_rk, 0.0_rk/)
    do while ( norm2_3d(e_r) <= 1.0d-12 )
        e_r = cross( (/rand_nbr(),rand_nbr(),rand_nbr()/), e_x)
    enddo
    e_r = e_r / norm2_3d(e_r)

    ! third (also radial) unit vector, simply the cross product of the others
    e_3 = cross(e_x,e_r)
    e_3 = e_3 / norm2_3d(e_3)


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
                ! disabled
                ! if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! cb is the distance to the cylinder mid-point
                cb = 0.5_rk*(x1+x2) - x_glob
                ! rb is the length of the clinder
                rb = x1 - x2

                ! this is a spherical bounding box, centered around the mid-point
                if ( sum(cb**2) < 0.25*sum(rb**2) ) then ! the 0.25 is from the 0.5 squared
                    ab = x_glob - x1
                    u = x2 - x1

                    vp = cross(ab, u)
                    R = sqrt( sum(vp**2) / sum(u**2) )

                    if (R <= R0+safety) then
                        t = step(R, R0, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                        if (t >= mask(ix,iy,iz)) then

                            mask(ix,iy,iz) = t
                            mask_color(ix,iy,iz) = color_val

                            x_body = matmul(M_g2b, x_glob - Insect%xc_body_g)
                            x_wing = matmul(M_b2w, x_body - x_pivot_b)

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
                            us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
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

                ! disabled
                ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                ! the bounding box check is incorporated in the loop bounds - no if clause!
                ! compute radius
                R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
                if ( R <= R0+Insect%safety ) then
                    t = step(R, R0, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                    if ( t >= mask(ix,iy,iz) ) then
                        ! set new value
                        mask(ix,iy,iz) = t
                        mask_color(ix,iy,iz) = color_val

                        x_body = matmul(M_g2b, x + x2 - Insect%xc_body_g)
                        x_wing = matmul(M_b2w, x_body - x_pivot_b)

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
                        us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
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

                ! disabled
                ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                ! the bounding box check is incorporated in the loop bounds - no if clause!
                ! compute radius
                R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
                if ( R <= R0+Insect%safety ) then
                    t = step(R, R0, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int)
                    if ( t >= mask(ix,iy,iz) ) then
                        ! set new value
                        mask(ix,iy,iz) = t
                        mask_color(ix,iy,iz) = color_val

                        x_body = matmul(M_g2b, x + x1 - Insect%xc_body_g)
                        x_wing = matmul(M_b2w, x_body - x_pivot_b)

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
                        us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
                    endif
                endif
            enddo
        enddo
    enddo

end subroutine


!-------------------------------------------------------------------------------
! Draw a triangular prism determined by points x1,y1 x2,y2 x3,y3
!-------------------------------------------------------------------------------
subroutine draw_trianglular_prism(xx0, ddx, mask, mask_color, us,Insect,color_wing,M_g2b,M_b2w,x_pivot_b,rot_rel_wing_w, &
                        x1,y1,x2,y2,x3,y3, wingID)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)  !< block origin and grid spacing
  real(kind=rk),intent(inout) :: mask(0:,0:,0:)  !< mask value between 0 and 1
  real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)  !< velocity field of the body
  real(kind=rk),intent(inout) :: mask_color(0:,0:,0:)  !< color, usually integer values to differentiate the parts
  integer(kind=2),intent(in) :: color_wing  !< color that will be set
  real(kind=rk),intent(in) :: M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3)  !< coordinate transformation matrices and pivot point
  real(kind=rk),intent(in) :: x1,y1,x2,y2,x3,y3  ! coordinates of the triangle in the wing system
  integer(kind=2), intent(in) :: wingID

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

  do iz = g, size(mask,3)-1-g ! note zero-based indexing in this module, which may appear odd in WABBIT (usually 1-based)
      x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
      do iy = g, size(mask,2)-1-g
          x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
          do ix = g, size(mask,1)-1-g
              x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

              !-- define the various coordinate systems we are going to use
              ! disabled
              ! if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

              x_body = matmul(M_g2b,x)
              x_wing = matmul(M_b2w,x_body-x_pivot_b)

              !-- test point
              xt = x_wing(1)
              yt = x_wing(2)

              !-- if inside bounding box
              if ( (xt>=xmin-Insect%safety) .and. (xt<=xmax+Insect%safety) .and. &
                   (yt>=ymin-Insect%safety) .and. (yt<=ymax+Insect%safety) ) then

                  if (dabs(x_wing(3))<=0.5*Insect%Wings(wingID)%WingThickness + Insect%safety) then
                      !-- determine which side
                      v1 = f_same_side_point(x1,y1,x2,y2,xc,yc,xt,yt);
                      v2 = f_same_side_point(x2,y2,x3,y3,xc,yc,xt,yt);
                      v3 = f_same_side_point(x3,y3,x1,y1,xc,yc,xt,yt);

                      !-- if inside, draw
                      if ( (v1==1) .and. (v2==1) .and. (v3==1) ) then

                           !-- smooth height
                           mask_tmp = step(dabs(x_wing(3)),0.5_rk*Insect%Wings(wingID)%WingThickness, Insect%L_smooth, Insect%safety, Insect%smoothing_type_int) ! thickness

                           if ((mask(ix,iy,iz) < mask_tmp).and.(mask_tmp>0.0_rk)) then
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
                               us(ix,iy,iz,1:3) = matmul(transpose(M_b2w), v_tmp)
                           endif
                      endif
                  endif
              endif
          enddo
      enddo
  enddo
end subroutine draw_trianglular_prism


!-------------------------------------------------------------------------------
! Given a line through two points x1,y1 and x2,y2
! determine if points x3,y3 and x4,y4 are on the same side
! Return 1 if true
! Return 0 if false
!-------------------------------------------------------------------------------
function f_same_side_point(x1,y1,x2,y2,x3,y3,x4,y4)
    implicit none

    real(kind=rk), intent(in) :: x1,y1,x2,y2,x3,y3,x4,y4
    real(kind=rk) :: rr,dx21,dy21,m1,b1,b3,b4
    integer :: f_same_side_point

    ! Initialize the result
    rr = 0

    ! Denominator in the slope of the line
    dx21 = (x2-x1)
    dy21 = (y2-y1)
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
    elseif (dabs(dx21) > 1.0d-5) then
        m1 = dx21/dy21
        b1 = x1-m1*y1
        b3 = x3-m1*y3
        b4 = x4-m1*y4
        if ( ( (b1>=b3) .and. (b1>=b4) ) .or. ( (b1<=b3) .and. (b1<=b4) ) ) then
            rr = 1
        endif
    endif
    f_same_side_point = rr
end function f_same_side_point
