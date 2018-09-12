module module_insects
  ! The following module includes required functions from FLUSI/WABBIT and hence
  ! makes this module here independent of that
  use module_insects_integration_flusi_wabbit

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  ! functions
  PUBLIC :: Draw_Insect, Update_Insect, insect_init, insect_clean, draw_fractal_tree, draw_active_grid_winglets, &
  aero_power, inert_power, read_insect_STATE_from_file, rigid_solid_init, rigid_solid_time_step, &
  BodyMotion, FlappingMotion_right, FlappingMotion_left, StrokePlane, mask_from_pointcloud, &
  body_rotation_matrix, wing_right_rotation_matrix, wing_left_rotation_matrix, write_kinematics_file
  ! type definitions
  PUBLIC :: wingkinematics, diptera

  ! we use this so only root prints write statements...
  logical :: root = .false.


  ! size (global) of domain
  real(kind=rk) :: xl, yl, zl
  ! viscosity (just for printing the Reynolds number)
  real(kind=rk) :: nu

  ! arrays for fourier coefficients are fixed size (avoiding issues with allocatable
  ! elements in derived datatypes) this is their length:
  integer, parameter :: nfft_max = 1024
  ! Maximum number of Hermite interpolation nodes (hardcoded because of sxf90 compiler requirements)
  integer, parameter :: nhrmt_max = 10000

  ! Allocatable arrays used in Insect object
  ! this will hold the surface markers and their normals used for particles:
  real(kind=rk), allocatable, dimension(:,:) :: particle_points
  ! wing thickness profile
  real(kind=rk), allocatable, dimension(:,:) :: wing_thickness_profile
  ! wing corrugation profile
  real(kind=rk), allocatable, dimension(:,:) :: corrugation_profile

  ! wing signed distance function, if the 3d-interpolation approach is used.
  ! This is useful for highly complex wings, where one generates the mask only once
  ! and then interpolates the values to the global grid. note the data allocated
  ! here is of course understood in the wing system, so the linear transformation
  ! x_g -> x_w' is used, where x_w' is not a grid-aligned value x_w
  real(kind=rk), allocatable, dimension(:,:,:), save :: mask_wing_complete
  real(kind=rk), dimension(1:3), save :: mask_wing_xl, mask_wing_x0
  integer, dimension(1:3), save :: mask_wing_nxyz
  integer, save :: mask_wing_safety=4

  !-----------------------------------------------------------------------------
  ! TYPE DEFINITIONS
  ! datatype for wing kinematics, if described by a Fourier series or kineloader
  ! For both wings such a datatype is contained in the insect.
  type wingkinematics
    ! Fourier coefficients
    real(kind=rk) :: a0_alpha=0.d0, a0_phi=0.d0, a0_theta=0.d0
    real(kind=rk), dimension(1:nfft_max) :: ai_phi=0.d0, bi_phi=0.d0, ai_theta=0.d0, &
      bi_theta=0.d0, ai_alpha=0.d0, bi_alpha=0.d0
    integer :: nfft_phi=0, nfft_alpha=0, nfft_theta=0
    ! coefficients are read only once from file (or set differently)
    logical :: initialized = .false.
    ! some details about the file, if reading from ini file
    character(len=strlen) :: infile_convention="", infile_type="", infile_units="", infile=""
    ! variables for kineloader (which uses non-periodic hermite interpolation)
    integer :: nk=0
    real(kind=rk), dimension (1:nhrmt_max) :: vec_t=0.d0, &
      vec_phi=0.d0,vec_alpha=0.d0,vec_theta=0.d0,vec_pitch=0.d0,vec_vert=0.d0,vec_horz=0.d0,  &
      vec_phi_dt=0.d0,vec_alpha_dt=0.d0,vec_theta_dt=0.d0,vec_pitch_dt=0.d0,vec_vert_dt=0.d0, &
      vec_horz_dt=0.d0
  end type


  !-----------------------------------------------------------------------------
  ! derived datatype for insect parameters (for readability)
  type diptera
    !-------------------------------------------------------------
    ! Body motion state, wing motion state and characteristic points on insect
    !-------------------------------------------------------------
    ! position of logical center, and translational velocity
    real(kind=rk), dimension(1:3) :: xc_body_g=0.d0, vc_body_g=0.d0
    ! initial or tethered position, velocity and yawpitchroll angles:
    real(kind=rk), dimension(1:3) :: x0=0.d0, v0=0.d0, yawpitchroll_0=0.d0
    ! roll pitch yaw angles and their time derivatives
    real(kind=rk) :: psi=0.d0, beta=0.d0, gamma=0.d0, psi_dt=0.d0, beta_dt=0.d0, gamma_dt=0.d0, eta0=0.d0
    ! body pitch angle, if it is constant (used in forward flight and hovering)
    real(kind=rk) :: body_pitch_const=0.d0
    ! angles of the wings (left and right)
    real(kind=rk) :: phi_r=0.d0, alpha_r=0.d0, theta_r=0.d0, phi_dt_r=0.d0, alpha_dt_r=0.d0, theta_dt_r=0.d0
    real(kind=rk) :: phi_l=0.d0, alpha_l=0.d0, theta_l=0.d0, phi_dt_l=0.d0, alpha_dt_l=0.d0, theta_dt_l=0.d0
    ! stroke plane angle
    real(kind=rk) :: eta_stroke=0.d0
    ! angular velocity vectors (wings L+R, body)
    real(kind=rk), dimension(1:3) :: rot_body_b=0.d0, rot_body_g=0.d0
    real(kind=rk), dimension(1:3) :: rot_rel_wing_l_w=0.d0, rot_rel_wing_r_w=0.d0
    real(kind=rk), dimension(1:3) :: rot_rel_wing_l_b=0.d0, rot_rel_wing_r_b=0.d0
    real(kind=rk), dimension(1:3) :: rot_rel_wing_l_g=0.d0, rot_rel_wing_r_g=0.d0
    real(kind=rk), dimension(1:3) :: rot_abs_wing_l_g=0.d0, rot_abs_wing_r_g=0.d0
    ! angular acceleration vectors (wings L+R)
    real(kind=rk), dimension(1:3) :: rot_dt_wing_l_w=0.d0, rot_dt_wing_r_w=0.d0
    real(kind=rk), dimension(1:3) :: rot_dt_wing_l_g=0.d0, rot_dt_wing_r_g=0.d0
    ! Vector from body centre to pivot points in global reference frame
    real(kind=rk), dimension(1:3) :: x_pivot_l_g=0.d0, x_pivot_r_g=0.d0
    ! vectors desribing the positoions of insect's key elements
    ! in the body coordinate system
    real(kind=rk), dimension(1:3) :: x_head=0.d0,x_eye_r=0.d0,x_eye_l=0.d0,x_pivot_l_b=0.d0,x_pivot_r_b=0.d0
    ! moments of inertia in the body reference frame
    real(kind=rk) :: Jroll_body=0.d0, Jyaw_body=0.d0, Jpitch_body=0.d0
    ! total mass of insect:
    real(kind=rk) :: mass, gravity=0.d0
    ! variables to decide whether to draw the body or not.
    character(len=strlen) :: body_moves="yes"
    logical :: body_already_drawn = .false.
    !-------------------------------------------------------------
    ! for free flight solver
    !-------------------------------------------------------------
    real(kind=rk) :: time=0.d0
    real(kind=rk), dimension(1:20) :: RHS_old=0.d0, RHS_this=0.d0
    real(kind=rk), dimension(1:20) :: STATE=0.d0
    ! STATE(1) : x-position of body
    ! STATE(2) : y-position of body
    ! STATE(3) : z-position of body
    ! STATE(4) : x-velocity of body
    ! STATE(5) : y-velocity of body
    ! STATE(6) : z-velocity of body
    ! STATE(7) : 1st component of body quaternion
    ! STATE(8) : 2nd component of body quaternion
    ! STATE(9) : 3rd component of body quaternion
    ! STATE(10) : 4th component of body quaternion
    ! STATE(11) : x-angular velocity of body (in body system)
    ! STATE(12) : y-angular velocity of body (in body system)
    ! STATE(13) : z-angular velocity of body (in body system)
    ! STATE(14) : 1st component of left wing quaternion
    ! STATE(15) : 2nd component of left wing quaternion
    ! STATE(16) : 3rd component of left wing quaternion
    ! STATE(17) : 4th component of left wing quaternion
    ! STATE(18) : x-angular velocity of left wing
    ! STATE(19) : y-angular velocity of left wing
    ! STATE(20) : z-angular velocity of left wing
    real(kind=rk), dimension(1:6) :: DoF_on_off=0.d0
    character(len=strlen) :: startup_conditioner=""
    !-------------------------------------------------------------
    ! for wing fsi solver
    !-------------------------------------------------------------
    character(len=strlen) :: wing_fsi="no"
    real(kind=rk), dimension(1:3) :: torque_muscle_l_w=0.d0, torque_muscle_r_w=0.d0
    real(kind=rk), dimension(1:3) :: torque_muscle_l_b=0.d0, torque_muscle_r_b=0.d0
    real(kind=rk), dimension(1:3) :: init_alpha_phi_theta=0.d0
    !-------------------------------------------------------------
    ! wing shape parameters
    !-------------------------------------------------------------
    ! wing shape fourier coefficients. Note notation:
    ! R = a0/2 + SUM ( ai cos(2pi*i) + bi sin(2pi*i)  )
    ! to avoid compatibility issues, the array is of fixed size, although only
    ! the first nftt_wings entries will be used
    real(kind=rk), dimension(1:nfft_max) :: ai_wings=0.d0, bi_wings=0.d0
    real(kind=rk) :: a0_wings=0.d0
    ! fill the R0(theta) array once, then only table-lookup instead of Fseries
    real(kind=rk), dimension(1:25000) :: R0_table=0.d0
    ! describes the origin of the wings system
    real(kind=rk) :: xc=0.d0,yc=0.d0
    ! number of fft coefficients for wing geometry
    integer :: nfft_wings=0
    logical :: wingsetup_done = .false.
    logical :: wings_radius_table_ready = .false.
    ! wing bounding box (xmin, xmax, ymin, ymax, zmin, zmax)
    real(kind=rk) :: wing_bounding_box(1:6) = 0.d0
    ! wing inertia
    real(kind=rk) :: Jxx=0.d0,Jyy=0.d0,Jzz=0.d0,Jxy=0.d0
    character(len=strlen) :: wing_thickness_distribution = "constant"
    character(len=strlen) :: pointcloudfile = "none"
    logical :: corrugated = .false.

    !--------------------------------------------------------------
    ! Wing kinematics
    !--------------------------------------------------------------
    ! wing kinematics Fourier coefficients
    type(wingkinematics) :: kine_wing_l, kine_wing_r
    ! the following flag makes the code write the kinematics log to either kinematics.t
    ! (regular simulation) or kinematics.dry-run.t (for a dry run). The reason for this
    ! is that during postprocessing of an existing run, the dry run would overwrite the
    ! simulation data.
    character(len=strlen) :: kinematics_file = "kinematics.t"
    ! rotation matrices for the various coordinate system for the insect
    real(kind=rk),dimension(1:3,1:3) :: M_body, M_wing_l, M_wing_r, M_body_inv

    !-------------------------------------------------------------
    ! parameters that control shape of wings, body, and motion
    !-------------------------------------------------------------
    character(len=strlen) :: WingShape="", BodyType="", BodyMotion="", HasDetails=""
    character(len=strlen) :: FlappingMotion_right="", FlappingMotion_left=""
    character(len=strlen) :: infile="", LeftWing="", RightWing=""
    ! parameters for body:
    real(kind=rk) :: L_body=0.d0, b_body=0.d0, R_head=0.d0, R_eye=0.d0
    ! parameters for wing shape:
    real(kind=rk) :: b_top=0.d0, b_bot=0.d0, L_span=0.d0, WingThickness=0.d0
    ! this is a safety distance for smoothing:
    real(kind=rk) :: safety=0.d0, smooth=0.d0
    ! parameter for hovering:
    real(kind=rk) :: distance_from_sponge=0.d0
    ! Wings and body forces (1:body, 2:left wing, 3:right wing)
    type(Integrals), dimension(1:3) :: PartIntegrals

    !-------------------------------------------------------------
    ! parameters for mask coloring
    !-------------------------------------------------------------
    integer(kind=2) :: color_body=1, color_l=2, color_r=3
  end type diptera
  !-----------------------------------------------------------------------------

contains


  !---------------------------------------
  ! note these include files also have to be specified as dependencies in the
  ! Makefile for make to check if one of them changed
  include "insect_init_clean.f90"
  include "body_geometry.f90"
  include "body_motion.f90"
  include "rigid_solid_time_stepper.f90"
  include "wings_geometry.f90"
  include "wings_motion.f90"
  include "stroke_plane.f90"
  include "kineloader.f90"
  include "pointcloud.f90"
  include "fractal_trees.f90"
  include "active_grid_winglets.f90"
  !---------------------------------------


  !-------------------------------------------------------------------------------
  ! Wrapper for dynamic array allocation within Insect structure
  ! This cannot be 'include'-ed to allow conditional compilation with cpp,
  ! as SX compiler only supports static arrays within structures
  !-------------------------------------------------------------------------------
  subroutine Allocate_Arrays ( Insect, array_name, a, b )
    implicit none

    character(len=*), intent(in) :: array_name
    integer, intent(in) :: a, b
    type(diptera), intent(inout) :: Insect

    select case ( array_name )
    case ("wing_thickness_profile")
        if (.not.allocated(wing_thickness_profile)) allocate(wing_thickness_profile(1:a,1:b))
    case ("corrugation_profile")
        if (.not.allocated(corrugation_profile)) allocate(corrugation_profile(1:a,1:b))
    case default
    endselect

  end subroutine Allocate_Arrays


    !-----------------------------------------------------------------------------
    ! Many parts of the insect mask generation are done only once per time step (i.e.
    ! per mask generation). Now, the adaptive code calls Draw_Insect several times, on each
    ! block of the grid. Draw_Insect is thus called SEVERAL times per mask generation.
    ! Therefore, we outsource the parts that need to be done only once to this routine,
    ! and call it BEFORE calling Draw_Insect. For FLUSI, this does not have any effect
    ! other than having two routines.
    !-----------------------------------------------------------------------------
    subroutine Update_Insect( time, Insect )
        implicit none

        real(kind=rk), intent(in) :: time
        type(diptera),intent(inout) :: Insect
        logical, save :: first_call = .true.

        !-----------------------------------------------------------------------------
        ! fetch current motion state
        !-----------------------------------------------------------------------------
        call BodyMotion (time, Insect)
        call FlappingMotion_right (time, Insect)
        call FlappingMotion_left (time, Insect)
        call StrokePlane (time, Insect)

        !-----------------------------------------------------------------------------
        ! define the rotation matrices to change between coordinate systems
        !-----------------------------------------------------------------------------
        call body_rotation_matrix( Insect, Insect%M_body )
        call wing_right_rotation_matrix( Insect, Insect%M_wing_r )
        call wing_left_rotation_matrix( Insect, Insect%M_wing_l )

        ! inverse of the body rotation matrices
        Insect%M_body_inv = transpose(Insect%M_body)

        ! body angular velocity vector in b/g coordinate system
        call body_angular_velocity( Insect, Insect%rot_body_b, Insect%rot_body_g, Insect%M_body )
        ! rel+abs wing angular velocities in the w/b/g coordinate system
        call wing_angular_velocities ( time, Insect, Insect%M_body )
        ! angular acceleration for wings (required for inertial power)
        call wing_angular_accel( time, Insect )

        !-----------------------------------------------------------------------------
        ! vector from body centre to left/right pivot point in global reference frame,
        ! for aerodynamic power
        !-----------------------------------------------------------------------------
        Insect%x_pivot_l_g = matmul(Insect%M_body_inv, Insect%x_pivot_l_b)
        Insect%x_pivot_r_g = matmul(Insect%M_body_inv, Insect%x_pivot_r_b)


        if (first_call) then
            ! print some important numbers, routine exectuted only once during a simulation
            call print_insect_reynolds_numbers( Insect )
            first_call = .false.
        endif

        ! save time to insect, then we can check if the update routine has been called
        ! or not (this is not necessary if Update_Insect is called, but helpful to prevent
        ! human errors)
        Insect%time = time

    end subroutine Update_Insect

  !-------------------------------------------------------------------------------
  ! Main routine for drawing insects. Loops over the entire domain, computes
  ! coordinates in various systems (global-, body-, stroke-, wing-) and calls
  ! subroutines doing the actual job of defining the mask. Note all surfaces are
  ! smoothed.
  !-------------------------------------------------------------------------------
  subroutine Draw_Insect( time, Insect, xx0, ddx, mask, mask_color, us)
      implicit none

      real(kind=rk), intent(in) :: time
      type(diptera),intent(inout) :: Insect
      real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
      real(kind=rk),intent(inout) :: mask(0:,0:,0:)
      real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
      integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

      real(kind=rk) :: t1
      real(kind=rk),dimension(1:3) :: x, x_body, v_tmp
      ! real(kind=rk),dimension(1:3,1:3) :: M_body, M_wing_l, M_wing_r, M_body_inv
      integer :: ix, iy, iz
      integer(kind=2) :: c

      if ((dabs(Insect%time-time)>1.0d-10).and.root) then
          write(*,'("error! time=",es15.8," but Insect%time=",es15.8)') time, Insect%time
          write(*,'("Did you call Update_Insect before Draw_Insect?")')
      endif

      Insect%smooth = 1.0d0*maxval(ddx)
      Insect%safety = 3.5d0*Insect%smooth

      ! delete old mask
      call delete_old_mask( time, mask, mask_color, us, Insect )

      !-----------------------------------------------------------------------------
      ! Draw individual parts of the Diptera. Separate loops are faster
      ! since the compiler can optimize them better
      !-----------------------------------------------------------------------------
      ! BODY. Now the body is special: if the insect does not move (or rotate), the
      ! body does not change in time. On the other hand, it is quite expensive to
      ! compute, since it involves a lot of points (volume), and it is a source of
      ! load balancing problems, since many cores do not draw the body at all.
      ! We thus try to draw it only once and then simply not to erase it later.
      !-----------------------------------------------------------------------------
      if (Insect%body_moves=="no" .and. avoid_drawing_static_body) then
          if (.not. Insect%body_already_drawn) then
              ! the body is at rest, but it is the first call to this routine, so
              ! draw it now.
              if (root .and. (.not. Insect%body_already_drawn) ) then
                  write(*,*) "Flag Insect%body_moves is no and we did not yet draw"
                  write(*,*) "the body once: we do that now, and skip draw_body"
                  write(*,*) "from now on. time=", time
              endif

              call draw_body( xx0, ddx, mask, mask_color, us, Insect, Insect%color_body, Insect%M_body)
              Insect%body_already_drawn = .true.
          endif
      else
          ! the body moves, draw it
          call draw_body(xx0, ddx, mask, mask_color, us, Insect, Insect%color_body, Insect%M_body)
      endif

      !-----------------------------------------------------------------------------
      ! Wings
      !-----------------------------------------------------------------------------
      if (Insect%RightWing == "yes") then
          call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%color_r, Insect%M_body, &
          Insect%M_wing_r, Insect%x_pivot_r_b, Insect%rot_rel_wing_r_w )
      endif

      if (Insect%LeftWing == "yes") then
          call draw_wing(xx0, ddx, mask, mask_color, us, Insect, Insect%color_l, Insect%M_body, &
          Insect%M_wing_l, Insect%x_pivot_l_b, Insect%rot_rel_wing_l_w )
      endif

      !-----------------------------------------------------------------------------
      ! Add solid body rotation (i.e. the velocity field that originates
      ! from the body rotation and translation). Until now, the wing velocities
      ! were the only ones set plus they are in the body reference frame
      !-----------------------------------------------------------------------------
      do iz = 0, size(mask,3)-1
          do iy = 0, size(mask,2)-1
              do ix = 0, size(mask,1)-1
                  c = mask_color(ix,iy,iz)
                  ! skip all parts that do not belong to the insect (ie they have a different color)
                  if (c==Insect%color_body .or. c==Insect%color_l .or. c==Insect%color_r ) then
                      x = (/ xx0(1)+dble(ix)*ddx(1), xx0(2)+dble(iy)*ddx(2), xx0(3)+dble(iz)*ddx(3) /)
                      x = periodize_coordinate(x - Insect%xc_body_g, (/xl,yl,zl/))
                      x_body = matmul(Insect%M_body, x)

                      ! add solid body rotation in the body-reference frame, if color
                      ! indicates that this part of the mask belongs to the insect
                      if (mask(ix,iy,iz) > 0.d0) then

                          ! translational part. we compute the rotational part in the body
                          ! reference frame, therefore, we must transform the body translation
                          ! velocity Insect%vc (which is in global coordinates) to the body frame
                          v_tmp = matmul(Insect%M_body, Insect%vc_body_g)

                          ! add solid body rotation to the translational velocity field. Note
                          ! that rot_body_b and x_body are in the body reference frame
                          v_tmp(1) = v_tmp(1)+Insect%rot_body_b(2)*x_body(3)-Insect%rot_body_b(3)*x_body(2)
                          v_tmp(2) = v_tmp(2)+Insect%rot_body_b(3)*x_body(1)-Insect%rot_body_b(1)*x_body(3)
                          v_tmp(3) = v_tmp(3)+Insect%rot_body_b(1)*x_body(2)-Insect%rot_body_b(2)*x_body(1)

                          ! the body motion is added to the wing motion, which is already in us
                          ! and they are also in the body refrence frame. However, us has to be
                          ! in the global reference frame, so M_body_inverse is applied
                          us(ix,iy,iz,1:3) = matmul( Insect%M_body_inv, us(ix,iy,iz,1:3)+v_tmp )
                      endif
                  endif
              enddo
          enddo
      enddo

      ! this is a debug test, which suceeded.
      !call check_if_us_is_derivative_of_position_wingtip(time, Insect)
  end subroutine Draw_Insect



  !-------------------------------------------------------
  ! short for the smooth step function.
  ! the smooting is defined in Insect%smooth, here we need only x, and the
  ! thickness (i.e., in the limit, steps=1 if x<t and steps=0 if x>t
  !-------------------------------------------------------
  real(kind=rk) function steps(x, t, h)
     implicit none
    real(kind=rk) :: f,x,t, h
    call smoothstep(f,x,t,h)
    steps=f
  end function


  !-------------------------------------------------------
  ! Compute angle from coefficients provided by Maeda
  !-------------------------------------------------------
  subroutine get_dangle( angles, F, a, b, shift_phase, initial_phase, dangle, dangle_dt )
    implicit none
    integer, intent(in) :: F  ! wavenumber (Dmitry, 7 Nov 2013)
    real(kind=rk), intent(in) :: angles ! 2*pi*F*time (Dmitry, 7 Nov 2013)
    real(kind=rk), intent(in) :: a
    real(kind=rk), intent(in) :: b
    real(kind=rk), intent(in) :: shift_phase
    real(kind=rk), intent(in) :: initial_phase
    real(kind=rk), intent(out) :: dangle
    real(kind=rk), intent(out) :: dangle_dt ! velocity increment (Dmitry, 7 Nov 2013)
    real(kind=rk) :: dAmp
    real(kind=rk) :: factor_amp = 1.0d0  ! Dmitry, 7 Nov 2013
    real(kind=rk) :: phase
    !!----------------------------

    !! d_amplitude
    dAmp = dsqrt(a**2 +b**2)*factor_amp

    !! phase
    if( b>0.0d0 ) then
      phase = datan(a/b)
    elseif( b<0.0d0 ) then
      phase = datan(a/b) +pi
    else !! b == 0 -> avoid division by zero
      phase = pi*0.5d0 !! sin(PI/2) = cos
    endif

    phase = phase + (shift_phase +initial_phase*2.0d0*pi)*dble(F)

    !! d_angle
    dangle = dAmp*dsin( angles +phase )

    !! velocity increment (Dmitry, 7 Nov 2013)
    dangle_dt = 2.0d0*pi*dble(F) * dAmp*dcos( angles +phase )

    return
  end subroutine get_dangle



  ! Compute aerodynamic power
  subroutine aero_power(Insect,apowtotal)
    implicit none

    integer :: color_body, color_l, color_r
    real(kind=rk), dimension(1:3) :: omrel, momrel
    real(kind=rk), intent(out) :: apowtotal
    type(diptera),intent(inout)::Insect

    ! colors for Diptera (one body, two wings)
    color_body = Insect%color_body
    color_l = Insect%color_l
    color_r = Insect%color_r

    ! body is not driven directly, therefore the power is set to zero
    Insect%PartIntegrals(color_body)%APow = 0.0d0

    !-----------
    ! left wing
    !-----------
    ! relative angular velocity, in global system
    omrel = Insect%rot_rel_wing_l_g

    ! compute moment with respect to the pivot point
    ! initialize it as the moment with respect to insect's centre point
    momrel = Insect%PartIntegrals(color_l)%Torque + &
             Insect%PartIntegrals(color_l)%Torque_unst

    ! aerodynamic power
    Insect%PartIntegrals(color_l)%APow = - sum( momrel * omrel )

    !-----------
    ! right wing
    !-----------
    ! relative angular velocity, in global system
    omrel = Insect%rot_rel_wing_r_g

    ! compute moment with respect to the pivot point
    ! initialize it as the moment with respect to insect's centre point
    momrel = Insect%PartIntegrals(color_r)%Torque + &
             Insect%PartIntegrals(color_r)%Torque_unst

    ! aerodynamic power
    Insect%PartIntegrals(color_r)%APow = - sum( momrel * omrel )

    !-----------
    ! Total aerodynamic power
    !-----------
    apowtotal = Insect%PartIntegrals(color_body)%APow + &
    Insect%PartIntegrals(color_l)%APow + Insect%PartIntegrals(color_r)%APow

  end subroutine aero_power


  !-------------------------------------------------------------------------------
  ! Compute interial power, i.e. the power the insect would have to invest
  ! when flapping its wings in vacuum.
  ! OUTPUT:
  !       ipowtotal: total inertial power
  !       Insect%PartIntegrals%IPow: (global): individual inertial power
  ! INPUT:
  !       Insect%rot_dt_wing_l_w (global): left wing angular acceleration
  !       Insect%rot_dt_wing_r_w (global): right wing angular acceleration
  !       Insect%Jxx,Jyy,Jxy,Jzz (global) Wing inertia
  ! MATHEMATICS:
  !       P_inertia = omega*( J*omega_dt + omega \cross (J*omega) )
  !                 = omega*( a + omega \cross b )
  !       The interia tensor is (it is specified in the PARAMS file)
  !           / Jxx Jxy 0   \
  !       J = | Jxy Jyy 0   |
  !           \ 0   0   Jzz /
  ! SEE ALSO
  !       Berman, Wang: Energy minimizing kinematics in hovering insect flight
  !       (JFM 582, 2007), eqn 2.22 (looks a bit different)
  !-------------------------------------------------------------------------------
  subroutine inert_power(Insect,ipowtotal)
    implicit none

    real(kind=rk), intent(out) :: ipowtotal
    real(kind=rk), dimension(1:3) :: a,b
    integer(kind=2) :: color_body, color_l,color_r
    type(diptera),intent(inout)::Insect

    ! colors for Diptera (one body, two wings)
    color_body = Insect%color_body
    color_l = Insect%color_l
    color_r = Insect%color_r

    !-- LEFT WING
    a(1) = Insect%Jxx * Insect%rot_dt_wing_l_w(1) + Insect%Jxy * Insect%rot_dt_wing_l_w(2)
    a(2) = Insect%Jxy * Insect%rot_dt_wing_l_w(1) + Insect%Jyy * Insect%rot_dt_wing_l_w(2)
    a(3) = Insect%Jzz * Insect%rot_dt_wing_l_w(3)

    b(1) = Insect%Jxx * Insect%rot_rel_wing_l_w(1) + Insect%Jxy * Insect%rot_rel_wing_l_w(2)
    b(2) = Insect%Jxy * Insect%rot_rel_wing_l_w(1) + Insect%Jyy * Insect%rot_rel_wing_l_w(2)
    b(3) = Insect%Jzz * Insect%rot_rel_wing_l_w(3)

    Insect%PartIntegrals(color_l)%IPow = &
    Insect%rot_rel_wing_l_w(1) * (a(1)+Insect%rot_rel_wing_l_w(2)*b(3)-Insect%rot_rel_wing_l_w(3)*b(2)) +&
    Insect%rot_rel_wing_l_w(2) * (a(2)+Insect%rot_rel_wing_l_w(3)*b(1)-Insect%rot_rel_wing_l_w(1)*b(3)) +&
    Insect%rot_rel_wing_l_w(3) * (a(3)+Insect%rot_rel_wing_l_w(1)*b(2)-Insect%rot_rel_wing_l_w(2)*b(1))

    !-- RIGHT WING
    a(1) = Insect%Jxx * Insect%rot_dt_wing_r_w(1) + Insect%Jxy * Insect%rot_dt_wing_r_w(2)
    a(2) = Insect%Jxy * Insect%rot_dt_wing_r_w(1) + Insect%Jyy * Insect%rot_dt_wing_r_w(2)
    a(3) = Insect%Jzz * Insect%rot_dt_wing_r_w(3)

    b(1) = Insect%Jxx * Insect%rot_rel_wing_r_w(1) + Insect%Jxy * Insect%rot_rel_wing_r_w(2)
    b(2) = Insect%Jxy * Insect%rot_rel_wing_r_w(1) + Insect%Jyy * Insect%rot_rel_wing_r_w(2)
    b(3) = Insect%Jzz * Insect%rot_rel_wing_r_w(3)

    Insect%PartIntegrals(color_r)%IPow = &
    Insect%rot_rel_wing_r_w(1) * (a(1)+Insect%rot_rel_wing_r_w(2)*b(3)-Insect%rot_rel_wing_r_w(3)*b(2)) +&
    Insect%rot_rel_wing_r_w(2) * (a(2)+Insect%rot_rel_wing_r_w(3)*b(1)-Insect%rot_rel_wing_r_w(1)*b(3)) +&
    Insect%rot_rel_wing_r_w(3) * (a(3)+Insect%rot_rel_wing_r_w(1)*b(2)-Insect%rot_rel_wing_r_w(2)*b(1))

    ipowtotal = Insect%PartIntegrals(color_r)%IPow + Insect%PartIntegrals(color_l)%IPow

  end subroutine inert_power

  !-----------------------------------------------------------------------------
  ! Body angular velocity vector
  !-----------------------------------------------------------------------------
  ! Variant (a) : free flight with quaternion solver
  !
  !    when using the quaternion based free-flight solver, the angular
  !    velocity of the body is computed dynamically, and the rotation matrix that
  !    brings us from global to body system is computed with quaternions. In body_motion
  !    the solver sets Insect%rot_body_b, so here we compute only rot_body_g
  !
  ! Variant (b) : imposed (prescribed) body dynamics
  !
  !    If the free flight solver is not active, body yaw,pitch,roll are known, so
  !    given yaw.pitch roll angles and their time derivatives, return the bodies
  !    angular velocity vector in global and body frame
  !-----------------------------------------------------------------------------
  subroutine body_angular_velocity( Insect, rot_body_b, rot_body_g, M_body )
    implicit none

    type(diptera), intent(inout) :: Insect
    real(kind=rk), intent(in) :: M_body(1:3,1:3)
    real(kind=rk), dimension(1:3), intent(out) :: rot_body_b, rot_body_g
    real(kind=rk) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt

    psi = Insect%psi
    beta = Insect%beta
    gamma = Insect%gamma
    psi_dt = Insect%psi_dt
    beta_dt = Insect%beta_dt
    gamma_dt = Insect%gamma_dt

    if ( Insect%BodyMotion == "free_flight" ) then
      ! variant (a)
      rot_body_b = Insect%rot_body_b ! copy (useless, actually, but required for interface)
      rot_body_g = matmul( transpose(M_body), rot_body_b)
    else
      ! variant (b)
      ! in global frame
      rot_body_g = (/ psi_dt*cos(beta)*cos(gamma)-beta_dt*sin(gamma) ,&
                      beta_dt*cos(gamma)+psi_dt*cos(beta)*sin(gamma) ,&
                      gamma_dt-psi_dt*sin(beta) /)
      ! in body frame
      rot_body_b = (/ psi_dt-gamma_dt*sin(beta) ,&
                      beta_dt*cos(psi)+gamma_dt*cos(beta)*sin(psi) ,&
                      gamma_dt*cos(beta)*cos(psi)-beta_dt*sin(psi) /)
    endif
  end subroutine body_angular_velocity



  !-------------------------------------------------------------------------------
  ! given the angles of each wing (and their time derivatives), compute
  ! the angular velocity vectors for both wings.
  ! output:
  !    Insect%rot_rel_wing_r_w    relative angular velocity of right wing (wing frame)
  !    Insect%rot_rel_wing_r_b    relative angular velocity of right wing (body frame)
  !    Insect%rot_rel_wing_r_g    relative angular velocity of right wing (glob frame)
  !    Insect%rot_abs_wing_r_g    absolute angular velocity of right wing (glob frame)
  !-------------------------------------------------------------------------------
  subroutine wing_angular_velocities ( time, Insect, M_body )
    implicit none

    real(kind=rk), intent(in) :: time
    real(kind=rk), intent(in) :: M_body(1:3,1:3)
    type(diptera), intent(inout) :: Insect

    real(kind=rk) :: eta_stroke
    real(kind=rk) :: phi_r, alpha_r, theta_r, phi_dt_r, alpha_dt_r, theta_dt_r
    real(kind=rk) :: phi_l, alpha_l, theta_l, phi_dt_l, alpha_dt_l, theta_dt_l
    real(kind=rk), dimension(1:3) :: rot_l_alpha, rot_l_theta, rot_l_phi, &
    rot_r_alpha, rot_r_theta, rot_r_phi
    real(kind=rk), dimension(1:3,1:3) :: M_wing_l, M_wing_r, &
    M1_tmp, M2_tmp, M1_l, M2_l, M3_l, M1_r, M2_r, M3_r, &
    M_stroke_l, M_stroke_r

    phi_r      = Insect%phi_r
    alpha_r    = Insect%alpha_r
    theta_r    = Insect%theta_r
    phi_dt_r   = Insect%phi_dt_r
    alpha_dt_r = Insect%alpha_dt_r
    theta_dt_r = Insect%theta_dt_r

    phi_l      = Insect%phi_l
    alpha_l    = Insect%alpha_l
    theta_l    = Insect%theta_l
    phi_dt_l   = Insect%phi_dt_l
    alpha_dt_l = Insect%alpha_dt_l
    theta_dt_l = Insect%theta_dt_l

    eta_stroke = Insect%eta_stroke

    !-----------------------------------------------------------------------------
    ! define the rotation matrices to change between coordinate systems
    !-----------------------------------------------------------------------------
    call Ry(M1_tmp,eta_stroke)
    M_stroke_l = M1_tmp

    call Rx(M1_tmp,pi)
    call Ry(M2_tmp,eta_stroke)
    M_stroke_r = matmul(M1_tmp,M2_tmp)

    call Ry(M1_l,alpha_l)
    call Rz(M2_l,theta_l)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3_l,phi_l)
    M_wing_l = matmul(M1_l,matmul(M2_l,matmul(M3_l,M_stroke_l)))

    ! note the coordinate system is rotated so we don't need to inverse the sign
    ! of theta, and the wings still rotate in opposite direction
    call Ry(M1_r,-alpha_r)
    call Rz(M2_r,theta_r)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3_r,-phi_r)
    M_wing_r = matmul(M1_r,matmul(M2_r,matmul(M3_r,M_stroke_r)))

    !-----------------------------------------------------------------------------
    ! angular velocity vectors (in wing system)
    !-----------------------------------------------------------------------------
    rot_l_alpha = (/ 0.0d0, alpha_dt_l, 0.0d0 /)
    rot_l_theta = (/ 0.0d0, 0.0d0, theta_dt_l /)
    rot_l_phi   = (/ phi_dt_l, 0.0d0, 0.0d0   /)
    rot_r_alpha = (/ 0.0d0, -alpha_dt_r, 0.0d0/)
    rot_r_theta = (/ 0.0d0, 0.0d0, theta_dt_r /)
    rot_r_phi   = (/ -phi_dt_r, 0.0d0, 0.0d0  /)

    ! in the wing coordinate system
    Insect%rot_rel_wing_l_w = matmul(M_wing_l,matmul(transpose(M_stroke_l),matmul(transpose(M3_l), &
    rot_l_phi+matmul(transpose(M2_l),rot_l_theta+matmul(transpose(M1_l), &
    rot_l_alpha)))))
    Insect%rot_rel_wing_r_w = matmul(M_wing_r,matmul(transpose(M_stroke_r),matmul(transpose(M3_r), &
    rot_r_phi+matmul(transpose(M2_r),rot_r_theta+matmul(transpose(M1_r), &
    rot_r_alpha)))))

    ! direct definition, equivalent to what is above.
    ! Insect%rot_rel_wing_l_w = (/phi_dt_l*cos(alpha_l)*cos(theta_l)-theta_dt_l*sin(alpha_l),&
    !   alpha_dt_l-phi_dt_l*sin(theta_l),&
    !   theta_dt_l*cos(alpha_l)+phi_dt_l*sin(alpha_l)*cos(theta_l)/)

    ! prior to the call of this routine, the routine body_angular_velocity has
    ! computed the body angular velocity (both g/b frames) so here now we can also
    ! compute global and absolute wing angular velocities.
    Insect%rot_rel_wing_l_b = matmul( transpose(M_wing_l), Insect%rot_rel_wing_l_w )
    Insect%rot_rel_wing_r_b = matmul( transpose(M_wing_r), Insect%rot_rel_wing_r_w )

    Insect%rot_rel_wing_l_g = matmul( transpose(M_body), Insect%rot_rel_wing_l_b )
    Insect%rot_rel_wing_r_g = matmul( transpose(M_body), Insect%rot_rel_wing_r_b )

    Insect%rot_abs_wing_l_g = Insect%rot_body_g + Insect%rot_rel_wing_l_g
    Insect%rot_abs_wing_r_g = Insect%rot_body_g + Insect%rot_rel_wing_r_g


    if (Insect%wing_fsi == "yes") then
      !**********************************
      !** Wing fsi model               **
      !**********************************
      ! overwrite the left wing
      Insect%rot_rel_wing_l_w = Insect%STATE(18:20)
      Insect%rot_rel_wing_l_b = matmul( transpose(M_wing_l), Insect%rot_rel_wing_l_w )
      Insect%rot_rel_wing_l_g = matmul( transpose(M_body), Insect%rot_rel_wing_l_b )
      ! the last guy is actually unused, as we have non-rotating body
      Insect%rot_abs_wing_l_g = Insect%rot_body_g + Insect%rot_rel_wing_l_g
    endif

  end subroutine wing_angular_velocities



  !-------------------------------------------------------------------------------
  ! Numerically estimate (it's a very precise estimation) the angular acceleration
  ! vectors for both wings, using one-sided finite differences (in the future)
  ! NOTE: this routine requires us to be able to evaluate both body and wing state
  !       at arbitrary times.
  !-------------------------------------------------------------------------------
  subroutine wing_angular_accel( time, Insect )
    implicit none
    real(kind=rk), intent(in) :: time
    type(diptera), intent(inout) :: Insect

    real(kind=rk) :: M_body(1:3,1:3), rot_dt_wing_g(1:3), M_wing_r(1:3,1:3), M_wing_l(1:3,1:3)
    type(diptera) :: Insect2
    real(kind=rk) :: dt,t

    dt = 1.0d-8
    Insect2 = Insect

    Insect%rot_dt_wing_l_w = 0.d0
    Insect%rot_dt_wing_r_w = 0.d0

    ! fetch motion state at time+dt
    call BodyMotion (time+dt, Insect2)
    call FlappingMotion_right(time+dt, Insect2)
    call FlappingMotion_left (time+dt, Insect2)
    call StrokePlane (time+dt, Insect2)
    call body_rotation_matrix( Insect2, M_body )
    call wing_angular_velocities ( time+dt, Insect2, M_body )
    ! this is the current state:
    call body_rotation_matrix( Insect, M_body )
    call wing_right_rotation_matrix( Insect, M_wing_r )
    call wing_left_rotation_matrix( Insect, M_wing_l )

    ! use one-sided finite differences to derive the absolute angular velocity with
    ! respect to time. NOte in older code versions, this was wrong, as we derived
    ! the ang. vel. in the wing coordinate system, which is a moving reference frame.

    ! now happily, the old results are still correct, as long as the body does not rotate
    ! see, e.g., this document https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjzl-XP6_LNAhWoC5oKHUdDCHwQFggoMAI&url=http%3A%2F%2Focw.mit.edu%2Fcourses%2Faeronautics-and-astronautics%2F16-07-dynamics-fall-2009%2Flecture-notes%2FMIT16_07F09_Lec08.pdf&usg=AFQjCNHzEB-n_NMm6K3J1eRpIaGnuKpW0Q&sig2=yEPNin3bL5DnWauNJk2hcw&bvm=bv.126993452,d.bGs&cad=rjt
    ! however, if the body moves, an additional term occurs, and this was indeed missing
    ! in previous results.
    rot_dt_wing_g = (Insect2%rot_rel_wing_l_g - Insect%rot_rel_wing_l_g) / dt
    Insect%rot_dt_wing_l_w = matmul(M_wing_l,matmul(M_body, rot_dt_wing_g))

    rot_dt_wing_g = (Insect2%rot_rel_wing_r_g - Insect%rot_rel_wing_r_g) / dt
    Insect%rot_dt_wing_r_w = matmul(M_wing_r,matmul(M_body, rot_dt_wing_g))


    ! if (root) then
    !   write(*,*) "L new code", Insect%rot_dt_wing_l_w
    !   write(*,*) "L old code", (Insect2%rot_rel_wing_l_w - Insect%rot_rel_wing_l_w)/dt
    !   write(*,*) "rot_rel_wing_l_g", Insect%rot_rel_wing_l_g
    !   write(*,*) "rot_body_g", Insect%rot_body_g
    !   write(*,*) "rot_body_b", Insect%rot_body_b
    !   write(*,*) "rot_dt_body_g", (Insect2%rot_body_g-Insect%rot_body_g)/dt
    ! endif
    !
    ! if (root) then
    !   t = 0.d0
    !   open  (17,file='test.t',status='replace')
    !   do while (t<=6.d0)
    !     call FlappingMotion_left ( t, Insect)
    !     write (17,'(7(es15.8,1x))') t,  &
    !     Insect%alpha_l, Insect%phi_l, Insect%theta_l, &
    !     Insect%alpha_dt_l, Insect%phi_dt_l, Insect%theta_dt_l
    !
    !     t = t + 1.0d-5
    !   end do
    !   close (17)
    !   call abort(7)
    ! endif

    ! ! this is the OLD CODE. It is actually wrong, since it computes the time derivative
    ! ! in the moving refrence frame of the wing, which is not the actual angular acceleration
    ! Insect%rot_dt_wing_r_w = (Insect2%rot_rel_wing_r_w - Insect%rot_rel_wing_r_w)/dt
    ! Insect%rot_dt_wing_l_w = (Insect2%rot_rel_wing_l_w - Insect%rot_rel_wing_l_w)/dt


    if (Insect%wing_fsi == "yes") then
      !**********************************
      !** Wing fsi model               **
      !**********************************
      ! overwrite the left wings acceleration. the right hand side 18,19,20 is the
      ! time derivative of the angular velocity, so the acceleration is readily available
      Insect%rot_dt_wing_l_w = Insect%RHS_THIS(18:20)
    endif

  end subroutine wing_angular_accel


  subroutine delete_old_mask( time, mask, mask_color, us, Insect )
      implicit none

      real(kind=rk), intent(in) :: time
      type(diptera),intent(in) :: Insect
      real(kind=rk),intent(inout) :: mask(0:,0:,0:)
      real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
      integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
      integer(kind=2) :: color_body, color_l, color_r
      logical, save :: cleaned_already_once = .false.

      ! colors for Diptera (one body, two wings)
      color_body = Insect%color_body
      color_l = Insect%color_l
      color_r = Insect%color_r

      !-----------------------------------------------------------------------------
      ! delete old mask
      !-----------------------------------------------------------------------------
      if (Insect%body_moves=="no" .and. avoid_drawing_static_body .and. cleaned_already_once) then
          ! the body is at rest, so we will not draw it. Delete the wings, as they move.
          where (mask_color==color_l .or. mask_color==color_r)
              mask = 0.d0
              mask_color = 0
          end where
          ! as the body rests it has no solid body velocity, which means we can safely
          ! reset the velocity everywhere (this step is actually unnessesary, but for
          ! safety we do it as well)
          us = 0.d0
      else
          ! the body of the insect moves, so we will construct the entire insect in this
          ! (and any other) call, and therefore we can safely reset the entire mask to zeros.
          mask = 0.d0
          mask_color = 0
          us = 0.d0
      endif

      cleaned_already_once = .true.

  end subroutine delete_old_mask


  !-----------------------------------------------------------------------------
  ! return the body rotation matrix
  !-----------------------------------------------------------------------------
  subroutine body_rotation_matrix( Insect, M_body )
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(out) :: M_body(1:3,1:3)
    real(kind=rk), dimension(1:3,1:3) :: M1_b, M2_b, M3_b

    if (Insect%BodyMotion=="free_flight") then
      ! entries 7,8,9,10 of the Insect%STATE vector are the body quaternion
      call rotation_matrix_from_quaternion( Insect%STATE(7:10), M_body)
    else
      ! conventional yaw, pitch, roll. Note the order of matrices is important.
      ! first we yaw, then we pitch, then we roll the insect. Note that when the
      ! free-flight solver is used, this matrix is obtained from quaternions, and
      ! not as a product of simple rotaion matrices. The latter can cause "gimbal-lock"
      call Rx(M1_b,Insect%psi)
      call Ry(M2_b,Insect%beta)
      call Rz(M3_b,Insect%gamma)
      M_body = matmul(M1_b,matmul(M2_b,M3_b))
    endif
  end subroutine body_rotation_matrix

  !-----------------------------------------------------------------------------
  ! return the rotation matrix for the right wing
  !-----------------------------------------------------------------------------
  subroutine wing_right_rotation_matrix( Insect, M_wing_r )
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(out) :: M_wing_r(1:3,1:3)
    real(kind=rk), dimension(1:3,1:3) :: M1, M2, M3, M_stroke_r


    call Rx(M1,pi)
    call Ry(M2,Insect%eta_stroke)
    M_stroke_r = matmul(M1,M2)

    ! note the coordinate system is rotated so we don't need to inverse the sign
    ! of theta, and the wings still rotate in opposite direction
    call Ry(M1,-Insect%alpha_r)
    call Rz(M2, Insect%theta_r)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3,-Insect%phi_r)
    M_wing_r = matmul(M1,matmul(M2,matmul(M3,M_stroke_r)))
  end subroutine wing_right_rotation_matrix


  !-----------------------------------------------------------------------------
  ! return the rotation matrix for the left wing
  !-----------------------------------------------------------------------------
  subroutine wing_left_rotation_matrix( Insect, M_wing_l )
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(out) :: M_wing_l(1:3,1:3)
    real(kind=rk),dimension(1:3,1:3) :: M1, M2, M3, M_stroke_l

    if ( Insect%wing_fsi /= "yes" ) then
      ! we're not using the wing fsi solver, so the wings follow a prescribed
      ! motion and we can compute the rotation matrix from the angles
      call Ry(M1,Insect%eta_stroke)
      M_stroke_l = M1

      call Ry(M1,Insect%alpha_l)
      call Rz(M2,Insect%theta_l)   ! Order changed (Dmitry, 7 Nov 2013)
      call Rx(M3,Insect%phi_l)
      M_wing_l = matmul(M1,matmul(M2,matmul(M3,M_stroke_l)))
    else
      !**********************************
      !** Wing fsi model               **
      !**********************************
      ! in the wing FSI case, a quaternion-based formulation is used to get the
      ! rotation matrix from the wing quaterion. note the wing quaternion are the
      ! entries 14,15,16,17 of the Insect%STATE vector
      ! entries 18,19,20 are the angular VELOCITY of the wing
      call rotation_matrix_from_quaternion( Insect%STATE(14:17), M_wing_l)
    endif
  end subroutine wing_left_rotation_matrix

  !-----------------------------------------------------------------------------
  ! Compute and print a couple of important numbers for insects
  !-----------------------------------------------------------------------------
  subroutine  print_insect_reynolds_numbers( Insect )
    implicit none
    type(diptera),intent(inout) :: Insect
    type(diptera) :: Insect_copy
    real(kind=rk) :: area, Re_f, Re
    real(kind=rk) :: time, dt
    real(kind=rk) :: phil_min, phil_max, phir_min, phir_max
    logical, save :: first_call = .true.

    ! the second call is just a return statement
    if ( first_call .eqv. .false.) return

    ! only root does this...
    if (root) then
      ! we need the wing area to compute the mean wing chord
      call compute_wing_surface(Insect, area)
      write(*,'(50("~"))')
      write(*,'("Wing area is A=",g15.8)') area
      write(*,'("Mean chord length is c_m=",g15.8)') area/1.d0 ! note c_m = A/R but R=1

      if (Insect%wing_fsi /= 'yes') then
        ! first we computethe stroke amplitude of the positional angle phi (for
        ! both wings). for safety, we make a copy of the insect, since the routines
        ! for the flapping motion write to this object, and we want to prevent any
        ! unwanted side effects
        Insect_copy = Insect
        time = 0.d0
        dt = 1.0d-3
        phil_min = 0.d0
        phil_max = 0.d0
        phir_min = 0.d0
        phir_max = 0.d0
        ! we use only one stroke ( the first one )
        do while (time < 1.d0)
          call FlappingMotion_left ( time, Insect_copy )
          call FlappingMotion_right ( time, Insect_copy )
          phil_min = min( phil_min, Insect_copy%phi_l )
          phil_max = max( phil_max, Insect_copy%phi_l )
          phir_min = min( phir_min, Insect_copy%phi_r )
          phir_max = max( phir_max, Insect_copy%phi_r )
          time = time + dt
        end do
        write(*,'("All following quantities are based on the first stroke 0.0 <= t <= 1.0")')
        write(*,'("Stroke amplitude is PHI_L=",g15.8)') abs(phil_max) + abs(phil_min)
        write(*,'("Stroke amplitude is PHI_R=",g15.8)') abs(phir_max) + abs(phir_min)
        write(*,'("Re_left  = 2*phi*R*f*c_m / nu =",g15.8)') 2.d0*(abs(phil_max)+abs(phil_min))*1.d0*1.d0*area/nu
        write(*,'("Re_right = 2*phi*R*f*c_m / nu =",g15.8)') 2.d0*(abs(phir_max)+abs(phir_min))*1.d0*1.d0*area/nu
        write(*,'("Re_f = R*R*f / nu =",g15.8)') 1.d0/nu
      else
        write(*,*) "In the case of wing_fsi problems, we do not know the wing"
        write(*,*) "kinematics before the computation. Therefore, we cannot tell"
        write(*,*) "the Reynolds number in advance."
        write(*,'("Re_f = R*R*f / nu =",g15.8)') 1.d0/nu
      endif
      write(*,'(50("~"))')
    endif

    first_call = .false.
  end subroutine print_insect_reynolds_numbers


  !-----------------------------------------------------------------------------
  ! In some cases, we need to reconstruct the mask or the body system in postprocessing
  ! if the free_flight solver was used, the body system cannnot be simply evaluated from closed-from
  ! expressions.
  ! In these case, we have to read the rigidsolidsolver.t file, which contains the Insect%STATE
  ! and from this the body orientation, rotation matrix, etc can be computed. So here we read this file
  ! and return Insect%STATE at the desired time (linear interpolation is used)
  !-----------------------------------------------------------------------------
  subroutine read_insect_STATE_from_file(time, Insect)
    implicit none
    real(kind=rk), intent(in) :: time
    type(diptera),intent(inout) :: Insect
    integer :: num_lines, n_header = 1, i
    character(len=maxcolumns) :: dummy
    real(kind=rk), allocatable, save :: data1(:,:)

    if ( .not. allocated(data1) ) then
      ! read rigidsolidsolver.t file
      ! skip header, count lines, read
      call count_lines_in_ascii_file_mpi('rigidsolidsolver.t', num_lines, n_header)
      ! read contents of file
      allocate( data1(1:num_lines,1:14))
      call read_array_from_ascii_file_mpi('rigidsolidsolver.t', data1 , n_header)
    endif

    ! interpolate in time
    i = 1
    do while (data1(i,1) <= time .and. i<size(data1,1)-1)
      i=i+1
    enddo
    ! we now have data1(i-1,1) <= time < data1(i,1)
    ! use linear interpolation
    Insect%STATE = 0.d0
    Insect%STATE(1:13) = data1(i-1,2:14) + (time - data1(i-1,1)) * (data1(i,2:14)-data1(i-1,2:14)) / (data1(i,1)-data1(i-1,1))

    if (root) then
        write(*,*) "The extracted Insect%STATE vector is:"
        write(*,'(21(es12.4,1x))') time, Insect%STATE(1:13)
    endif
    ! deallocate (data)
  end subroutine


  ! this routine computes  the wingteip velocity in the global system by two means:
  ! one we compute the global position vector (which we derive wrt time in postprocessing)
  ! and once the cross-products of rotation ang velocities as it is done in the actual code.
  ! we checked: both agree, also with imposed body velocity.
  subroutine check_if_us_is_derivative_of_position_wingtip(time, Insect)
    real(kind=rk), intent(in) :: time
    type(diptera), intent(inout) :: Insect

    real(kind=rk) :: M_body(1:3,1:3), M_wing_r(1:3,1:3), x_tip_w(1:3), x_tip_b(1:3), x_tip_g(1:3), &
      us_tip_g(1:3), v_tmp(1:3), v_tmp_b(1:3)
    real(kind=rk)::xd,yd,zd
    real(kind=rk)::c00,c10,c01,c11,c0,c1
    integer :: ix,iy,iz

    call body_rotation_matrix( Insect, M_body )
    call wing_right_rotation_matrix( Insect, M_wing_r )
    ! body angular velocity vector in b/g coordinate system
    call body_angular_velocity( Insect, Insect%rot_body_b, Insect%rot_body_g, M_body )
    ! rel+abs wing angular velocities in the w/b/g coordinate system
    call wing_angular_velocities ( time, Insect, M_body )

    x_tip_w = (/ 1.0d0, 1.0d0, 1.0d0 /)
    x_tip_b = matmul( transpose(M_wing_r), x_tip_w ) + Insect%x_pivot_r_b
    x_tip_g = matmul( transpose(M_body), x_tip_b ) + Insect%xc_body_g

    !-----------------------------------------------------------------------------
    ! now we extrcat how the us field is constructed
    v_tmp(1) = Insect%rot_rel_wing_r_w(2)*x_tip_w(3)-Insect%rot_rel_wing_r_w(3)*x_tip_w(2)
    v_tmp(2) = Insect%rot_rel_wing_r_w(3)*x_tip_w(1)-Insect%rot_rel_wing_r_w(1)*x_tip_w(3)
    v_tmp(3) = Insect%rot_rel_wing_r_w(1)*x_tip_w(2)-Insect%rot_rel_wing_r_w(2)*x_tip_w(1)
    v_tmp_b = matmul(transpose(M_wing_r), v_tmp) ! in body system


    ! translational part. we compute the rotational part in the body
    ! reference frame, therefore, we must transform the body translation
    ! velocity Insect%vc (which is in global coordinates) to the body frame
    v_tmp = matmul(M_body,Insect%vc_body_g)

    ! add solid body rotation to the translational velocity field. Note
    ! that rot_body_b and x_body are in the body reference frame
    v_tmp(1) = v_tmp(1) + Insect%rot_body_b(2)*x_tip_b(3)-Insect%rot_body_b(3)*x_tip_b(2)
    v_tmp(2) = v_tmp(2) + Insect%rot_body_b(3)*x_tip_b(1)-Insect%rot_body_b(1)*x_tip_b(3)
    v_tmp(3) = v_tmp(3) + Insect%rot_body_b(1)*x_tip_b(2)-Insect%rot_body_b(2)*x_tip_b(1)

    ! the body motion is added to the wing motion, which is already in us
    ! and they are also in the body refrence frame. However, us has to be
    ! in the global reference frame, so M_body_inverse is applied
    us_tip_g = matmul( transpose(M_body), v_tmp_b + v_tmp )

    open(14,file='debug_wing_us.t',status='unknown',position='append')
    write (14,'(7(es15.8,1x))') time, x_tip_g, us_tip_g
    close(14)

  end subroutine

  !-----------------------------------------------------------------------------
  ! write kinematics to disk. Note this cannot be done in the main Draw_Insect anymore
  ! because in the adaptive code, the routine is called Nblock times per time step,
  ! while FLUSI calls only once.
  !-----------------------------------------------------------------------------
  subroutine write_kinematics_file( time, Insect )
      implicit none
      real(kind=pr), intent(in) :: time
      type(diptera), intent(inout) :: Insect

      if (root) then
        open  (17, file=Insect%kinematics_file, status='unknown', position='append')
        write (17,'(26(es15.8,1x))') time, Insect%xc_body_g, Insect%psi, Insect%beta, &
        Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
        Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
        Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
        Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w
        close (17)
      endif

  end subroutine

end module module_insects
