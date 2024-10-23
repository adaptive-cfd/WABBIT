module module_insects
   ! The following module includes required functions from FLUSI/WABBIT and hence
   ! makes this module here independent of that
   use module_insects_integration_flusi_wabbit
   use module_t_files
   use module_stl_file_reader

   implicit none

   ! I usually find it helpful to use the private keyword by itself initially, which specifies
   ! that everything within the module is private unless explicitly marked public.
   PRIVATE

   ! functions
   PUBLIC :: Draw_Insect, draw_insect_body, draw_insect_wings, Update_Insect, insect_init, fractal_tree_init, &
      insect_clean, draw_fractal_tree, draw_active_grid_winglets, &
      aero_power, inert_power, read_insect_STATE_from_file, rigid_solid_init, rigid_solid_rhs, &
      BodyMotion, FlappingMotionWrap, StrokePlane, mask_from_pointcloud, &
      body_rotation_matrix, wing_rotation_matrix, write_kinematics_file
   ! type definitions
   PUBLIC :: wingkinematics, diptera

   ! we use this so only root prints write statements...
   logical :: root = .false.
   logical :: periodic_insect = .false.

   ! ghost nodes. If the insect module is used in a finite-differences code, then
   ! the data that we have often has ghost nodes, i.e. points that overlap and exist
   ! on several CPUS. On those, you normally would not create the mask (which is expensive)
   ! so we skip the first and last "g" points on the arrays used for mask creation
   integer, save :: g

   ! size (global) of domain
   real(kind=rk) :: xl, yl, zl
   ! viscosity (just for printing the Reynolds number)
   real(kind=rk) :: nu

   ! arrays for fourier coefficients are fixed size (avoiding issues with allocatable
   ! elements in derived datatypes) this is their length:
   integer, parameter :: nfft_max = 1024
   ! Maximum number of Hermite interpolation nodes (hardcoded because of sxf90 compiler requirements)
   ! JB: Array sizes resulting from this number result currently in 50% of WABBITs program size
   integer, parameter :: nhrmt_max = 10000

   ! Allocatable arrays used in Insect object
   ! this will hold the surface markers and their normals used for particles:
   real(kind=rk), allocatable, dimension(:,:) :: particle_points
   ! wing thickness profile
   real(kind=rk), allocatable, dimension(:,:,:) :: wing_thickness_profile
   ! wing thickness profile array dimensions
   integer, dimension(1:4), save :: wing_thickness_a, wing_thickness_b
   ! wing corrugation profile
   real(kind=rk), allocatable, dimension(:,:,:) :: corrugation_profile
   ! wing corrugation profile array dimensions
   integer, dimension(1:4), save :: corrugation_a, corrugation_b
   ! wing damage mask
   real(kind=rk), allocatable, dimension(:,:,:) :: damage_mask
   ! wing damage mask array dimensions
   integer, dimension(1:4), save :: damage_a, damage_b

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
   ! stuff for the fractal tree
   real(kind=rk), allocatable, save :: treedata(:,:)
   real(kind=rk), allocatable, save :: treedata_boundingbox(:,:)
   ! array for superSTL file for the body
   real(kind=rk), allocatable, save :: body_superSTL_b(:,:)
   real(kind=rk), allocatable, save :: body_superSTL_g(:,:)


   !-----------------------------------------------------------------------------
   ! TYPE DEFINITIONS
   ! datatype for wing kinematics, if described by a Fourier series or kineloader
   ! For both wings such a datatype is contained in the insect.
   type wingkinematics
      ! Fourier coefficients
      real(kind=rk) :: a0_alpha=0.0_rk, a0_phi=0.0_rk, a0_theta=0.0_rk
      real(kind=rk), dimension(1:nfft_max) :: ai_phi=0.0_rk, bi_phi=0.0_rk, ai_theta=0.0_rk, &
         bi_theta=0.0_rk, ai_alpha=0.0_rk, bi_alpha=0.0_rk
      integer :: nfft_phi=0, nfft_alpha=0, nfft_theta=0
      ! coefficients are read only once from file (or set differently)
      logical :: initialized = .false.
      ! some details about the file, if reading from ini file
      character(len=clong) :: infile_convention="", infile_type="", infile_units="", infile=""
      ! variables for kineloader (which uses non-periodic hermite interpolation)
      integer :: nk=0
      real(kind=rk), dimension (1:nhrmt_max) :: vec_t=0.0_rk, &
         vec_phi=0.0_rk,vec_alpha=0.0_rk,vec_theta=0.0_rk,vec_pitch=0.0_rk,vec_vert=0.0_rk,vec_horz=0.0_rk,  &
         vec_phi_dt=0.0_rk,vec_alpha_dt=0.0_rk,vec_theta_dt=0.0_rk,vec_pitch_dt=0.0_rk,vec_vert_dt=0.0_rk, &
         vec_horz_dt=0.0_rk
   end type


   !-----------------------------------------------------------------------------
   ! derived datatype for insect parameters (for readability)
   type diptera
      logical :: initialized = .false.
      !-------------------------------------------------------------
      ! Body motion state, wing motion state and characteristic points on insect
      !-------------------------------------------------------------
      ! position of logical center, and translational velocity
      real(kind=rk), dimension(1:3) :: xc_body_g=0.0_rk, vc_body_g=0.0_rk
      ! initial or tethered position, velocity and yawpitchroll angles:
      real(kind=rk), dimension(1:3) :: x0=0.0_rk, v0=0.0_rk, yawpitchroll_0=0.0_rk
      ! first harmonic components of the yawpitchroll angles
      real(kind=rk), dimension(1:3) :: yawpitchroll_a1=0.0_rk, yawpitchroll_b1=0.0_rk
      ! roll pitch yaw angles and their time derivatives
      real(kind=rk) :: psi=0.0_rk, beta=0.0_rk, gamma=0.0_rk, psi_dt=0.0_rk, beta_dt=0.0_rk, gamma_dt=0.0_rk, eta0=0.0_rk
      ! body pitch angle, if it is constant (used in forward flight and hovering)
      real(kind=rk) :: body_pitch_const=0.0_rk
      ! angles of the wings (left and right, second left and second right)
      real(kind=rk) :: phi_r=0.0_rk, alpha_r=0.0_rk, theta_r=0.0_rk, phi_dt_r=0.0_rk, alpha_dt_r=0.0_rk, theta_dt_r=0.0_rk
      real(kind=rk) :: phi_l=0.0_rk, alpha_l=0.0_rk, theta_l=0.0_rk, phi_dt_l=0.0_rk, alpha_dt_l=0.0_rk, theta_dt_l=0.0_rk
      real(kind=rk) :: phi_r2=0.0_rk, alpha_r2=0.0_rk, theta_r2=0.0_rk, phi_dt_r2=0.0_rk, alpha_dt_r2=0.0_rk, theta_dt_r2=0.0_rk
      real(kind=rk) :: phi_l2=0.0_rk, alpha_l2=0.0_rk, theta_l2=0.0_rk, phi_dt_l2=0.0_rk, alpha_dt_l2=0.0_rk, theta_dt_l2=0.0_rk
      ! stroke plane angle
      real(kind=rk) :: eta_stroke=0.0_rk
      ! is the body motion state described be the STATE vector? This is the case if the
      ! free-flight solver is used, and if its results are read in postprocessing or
      ! if it used used to prescribe the body motion state from a different simulation
      logical :: quaternion_solver_used = .false.
      ! angular velocity vectors (body, left and right wings, 2nd left and 2nd right wings)
      real(kind=rk), dimension(1:3) :: rot_body_b=0.0_rk, rot_body_g=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_l_w=0.0_rk, rot_rel_wing_r_w=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_l_b=0.0_rk, rot_rel_wing_r_b=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_l_g=0.0_rk, rot_rel_wing_r_g=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_abs_wing_l_g=0.0_rk, rot_abs_wing_r_g=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_l2_w=0.0_rk, rot_rel_wing_r2_w=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_l2_b=0.0_rk, rot_rel_wing_r2_b=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_l2_g=0.0_rk, rot_rel_wing_r2_g=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_abs_wing_l2_g=0.0_rk, rot_abs_wing_r2_g=0.0_rk
      ! angular acceleration vectors (left and right wings, 2nd left and 2nd right wings)
      real(kind=rk), dimension(1:3) :: rot_dt_wing_l_w=0.0_rk, rot_dt_wing_r_w=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_dt_wing_l_g=0.0_rk, rot_dt_wing_r_g=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_dt_wing_l2_w=0.0_rk, rot_dt_wing_r2_w=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_dt_wing_l2_g=0.0_rk, rot_dt_wing_r2_g=0.0_rk
      ! Vector from body centre to pivot points in global reference frame
      real(kind=rk), dimension(1:3) :: x_pivot_l_g=0.0_rk, x_pivot_r_g=0.0_rk
      real(kind=rk), dimension(1:3) :: x_pivot_l2_g=0.0_rk, x_pivot_r2_g=0.0_rk
      ! vectors desribing the positoions of insect's key elements
      ! in the body coordinate system
      real(kind=rk), dimension(1:3) :: x_head=0.0_rk,x_eye_r=0.0_rk,x_eye_l=0.0_rk,x_pivot_l_b=0.0_rk,x_pivot_r_b=0.0_rk
      real(kind=rk), dimension(1:3) :: x_pivot_l2_b=0.0_rk,x_pivot_r2_b=0.0_rk
      ! moments of inertia in the body reference frame
      real(kind=rk) :: Jroll_body=0.0_rk, Jyaw_body=0.0_rk, Jpitch_body=0.0_rk
      ! total mass of insect:
      real(kind=rk) :: mass, gravity=0.0_rk, gravity_y=0.0_rk, gravity_x=0.0_rk
      ! variables to decide whether to draw the body or not.
      character(len=clong) :: body_moves="yes"
      character(len=clong) :: BodySuperSTLfile="none.superstl"
      logical :: body_already_drawn = .false.
      ! second wing pair exists or not
      logical :: second_wing_pair
      !-------------------------------------------------------------
      ! for free flight solver
      !-------------------------------------------------------------
      real(kind=rk) :: time=0.0_rk
      real(kind=rk), allocatable :: RHS(:,:)
      real(kind=rk), dimension(1:20) :: STATE=0.0_rk
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
      ! STATE(14) : 1st component of left wing quaternion (WHY ONLY LEFT?)
      ! STATE(15) : 2nd component of left wing quaternion
      ! STATE(16) : 3rd component of left wing quaternion
      ! STATE(17) : 4th component of left wing quaternion
      ! STATE(18) : x-angular velocity of left wing
      ! STATE(19) : y-angular velocity of left wing
      ! STATE(20) : z-angular velocity of left wing
      real(kind=rk), dimension(1:6) :: DoF_on_off=0.0_rk
      character(len=clong) :: startup_conditioner=""
      !-------------------------------------------------------------
      ! for wing fsi solver
      !-------------------------------------------------------------
      real(kind=rk), dimension(1:3) :: torque_muscle_l_w=0.0_rk, torque_muscle_r_w=0.0_rk
      real(kind=rk), dimension(1:3) :: torque_muscle_l_b=0.0_rk, torque_muscle_r_b=0.0_rk
      real(kind=rk), dimension(1:3) :: init_alpha_phi_theta=0.0_rk
      !-------------------------------------------------------------
      ! wing shape parameters
      !-------------------------------------------------------------
      ! wing shape fourier coefficients. Note notation:
      ! R = a0/2 + SUM ( ai cos(2pi*i) + bi sin(2pi*i)  )
      ! to avoid compatibility issues, the array is of fixed size, although only
      ! the first nftt_wings entries will be used
      real(kind=rk), dimension(1:nfft_max,1:4) :: ai_wings=0.0_rk, bi_wings=0.0_rk
      real(kind=rk), dimension(1:4) :: a0_wings=0.0_rk
      ! fill the R0(theta) array once, then only table-lookup instead of Fseries
      ! JB: This array increases WABBIT program size quite a bit, consider making it allocatable
      real(kind=rk), dimension(1:25000,1:4) :: R0_table=0.0_rk
      ! describes the origin of the wings system
      real(kind=rk), dimension(1:4) :: xc=0.0_rk, yc=0.0_rk
      ! number of fft coefficients for wing geometry
      integer, dimension(1:4) :: nfft_wings=0
      logical, dimension(1:4) :: wingsetup_done = .false.
      logical, dimension(1:4) :: wings_radius_table_ready = .false.
      real(kind=rk), dimension(:), allocatable :: theta_i, R_i
      ! wing bounding box (xmin, xmax, ymin, ymax, zmin, zmax)
      real(kind=rk) :: wing_bounding_box(1:6,1:4) = 0.0_rk
      ! wing inertia
      real(kind=rk) :: Jxx=0.0_rk,Jyy=0.0_rk,Jzz=0.0_rk,Jxy=0.0_rk
      ! wing inertia of the second pair of wings
      real(kind=rk) :: Jxx2=0.0_rk,Jyy2=0.0_rk,Jzz2=0.0_rk,Jxy2=0.0_rk
      character(len=clong) :: wing_thickness_distribution(1:4) = "constant"
      character(len=clong) :: pointcloudfile = "none"
      character(len=clong) :: smoothing_thickness = "global", wing_file_type(1:4) = "fourier"
      logical :: corrugated(1:4) = .false.
      real(kind=rk) :: corrugation_array_bbox(1:4,1:4)
      logical :: bristles(1:4) = .false.
      logical :: bristles_simplex(1:4) = .false.
      integer :: n_bristles(1:4)
      real(kind=rk), ALLOCATABLE :: bristles_coords(:,:,:)
      ! used for rectangular part of bristled model wings (Kleemeier)
      real(kind=rk) :: B_membrane(1:4), L_membrane(1:4)
      logical :: damaged(1:4)

      !--------------------------------------------------------------
      ! Wing kinematics
      !--------------------------------------------------------------
      logical :: fractal_tree = .false.
      character(len=clong) :: fractal_tree_file = "tree_data.in"
      real(kind=rk), dimension(1:3) :: fractal_tree_x0 = (/0.0_rk, 0.0_rk, 0.0_rk/)
      real(kind=rk) :: fractal_tree_scaling = 1.0_rk

      !--------------------------------------------------------------
      ! Wing kinematics
      !--------------------------------------------------------------
      ! wing kinematics Fourier coefficients
      type(wingkinematics) :: kine_wing_l, kine_wing_r
      type(wingkinematics) :: kine_wing_l2, kine_wing_r2
      ! the following flag makes the code write the kinematics log to either kinematics.t
      ! (regular simulation) or kinematics.dry-run.t (for a dry run). The reason for this
      ! is that during postprocessing of an existing run, the dry run would overwrite the
      ! simulation data.
      character(len=clong) :: kinematics_file = "kinematics.t"
      ! rotation matrices for the various coordinate system for the insect
      real(kind=rk),dimension(1:3,1:3) :: M_g2b, M_b2w_l, M_b2w_r, M_b2g
      real(kind=rk),dimension(1:3,1:3) :: M_b2w_l2, M_b2w_r2

      !-------------------------------------------------------------
      ! parameters that control shape of wings, body, and motion
      !-------------------------------------------------------------
      character(len=clong) :: WingShape(1:4)=["","","",""] ! left, right, 2nd left, 2nd right
      character(len=clong) :: BodyType="", BodyMotion=""
      character(len=clong) :: FlappingMotion_right="", FlappingMotion_left=""
      character(len=clong) :: FlappingMotion_right2="", FlappingMotion_left2=""
      character(len=clong) :: infile="", LeftWing="", RightWing=""
      character(len=clong) :: infile2="", LeftWing2="", RightWing2=""
      ! parameters for body:
      real(kind=rk) :: L_body=0.0_rk, b_body=0.0_rk, R_head=0.0_rk, R_eye=0.0_rk
      ! parameters for wing shape:
      real(kind=rk) :: b_top=0.0_rk, b_bot=0.0_rk, L_span=0.0_rk, WingThickness=0.0_rk
      ! this is a safety distance for smoothing:
      real(kind=rk) :: safety=0.0_rk, smooth=0.0_rk, C_smooth=1.0_rk, dx_reference=0.0_rk, C_shell_thickness=5.0_rk
      ! parameter for hovering:
      real(kind=rk) :: distance_from_sponge=0.0_rk
      ! Wings and body forces (1:body, 2:left wing, 3:right wing, 4:left wing, 5:right wing)
      type(Integrals), dimension(1:5) :: PartIntegrals

      !-------------------------------------------------------------
      ! parameters for mask coloring
      !-------------------------------------------------------------
      ! available color values
      integer(kind=2) :: color_body=1, color_l=2, color_r=3, color_l2=4, color_r2=5
   end type diptera
   !-----------------------------------------------------------------------------

contains


!---------------------------------------
! note these include files also have to be specified as dependencies in the
! Makefile for make to check if one of them changed
#include "insect_init_clean.f90"
#include "body_geometry.f90"
#include "body_motion.f90"
#include "rigid_solid_time_stepper.f90"
#include "wings_geometry.f90"
#include "wings_motion.f90"
#include "stroke_plane.f90"
#include "kineloader.f90"
#include "pointcloud.f90"
#include "fractal_trees.f90"
#include "active_grid_winglets.f90"
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
      integer :: a_old, b_old
      real(kind=rk), allocatable, dimension(:,:,:) :: profile_tmp

      select case ( array_name )
       case ("wing_thickness_profile")
         if (.not.allocated(wing_thickness_profile)) then
            allocate(wing_thickness_profile(1:a,1:b,1:4))
         else
            a_old = size(wing_thickness_profile,1)
            b_old = size(wing_thickness_profile,2)
            if ( (a_old<a) .or. (b_old<b) ) then
               allocate(profile_tmp(1:a_old,1:b_old,1:4))
               profile_tmp(:,:,:) = wing_thickness_profile(:,:,:)
               deallocate(wing_thickness_profile)
               allocate(wing_thickness_profile(1:a,1:b,1:4))
               wing_thickness_profile(1:a_old,1:b_old,1:4) = profile_tmp(:,:,:)
               deallocate(profile_tmp)
            endif
         endif
       case ("corrugation_profile")
         if (.not.allocated(corrugation_profile)) then
            allocate(corrugation_profile(1:a,1:b,1:4))
         else
            a_old = size(corrugation_profile,1)
            b_old = size(corrugation_profile,2)
            if ( (a_old<a) .or. (b_old<b) ) then
               allocate(profile_tmp(1:a_old,1:b_old,1:4))
               profile_tmp(:,:,:) = corrugation_profile(:,:,:)
               deallocate(corrugation_profile)
               allocate(corrugation_profile(1:a,1:b,1:4))
               corrugation_profile(1:a_old,1:b_old,1:4) = profile_tmp(:,:,:)
               deallocate(profile_tmp)
            endif
         endif
       case ("damage_mask")
         if (.not.allocated(damage_mask)) then
            allocate(damage_mask(1:a,1:b,1:4))
         else
            a_old = size(damage_mask,1)
            b_old = size(damage_mask,2)
            if ( (a_old<a) .or. (b_old<b) ) then
               allocate(profile_tmp(1:a_old,1:b_old,1:4))
               profile_tmp(:,:,:) = damage_mask(:,:,:)
               deallocate(damage_mask)
               allocate(damage_mask(1:a,1:b,1:4))
               damage_mask(1:a_old,1:b_old,1:4) = profile_tmp(:,:,:)
               deallocate(profile_tmp)
            endif
         endif
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

      real(kind=rk), intent(in)   :: time
      type(diptera),intent(inout) :: Insect

      logical, save :: first_call = .true.
      integer(kind=ik) :: i
      integer(kind=2) :: wingID

      !-----------------------------------------------------------------------------
      ! fetch current motion state
      !-----------------------------------------------------------------------------
      call BodyMotion (time, Insect)
      wingID = 1
      call FlappingMotionWrap (time, Insect, wingID)
      wingID = 2
      call FlappingMotionWrap (time, Insect, wingID)

      if (Insect%second_wing_pair) then
         wingID = 3
         call FlappingMotionWrap (time, Insect, wingID)
         wingID = 4
         call FlappingMotionWrap (time, Insect, wingID)
      endif
      call StrokePlane (time, Insect)

      !-----------------------------------------------------------------------------
      ! define the rotation matrices to change between coordinate systems
      !-----------------------------------------------------------------------------
      call body_rotation_matrix( Insect, Insect%M_g2b )
      call wing_rotation_matrix( Insect%M_b2w_l, Insect%alpha_l, Insect%theta_l, Insect%phi_l, Insect%eta_stroke, "left" )
      call wing_rotation_matrix( Insect%M_b2w_r, Insect%alpha_r, Insect%theta_r, Insect%phi_r, Insect%eta_stroke, "right" )
      if (Insect%second_wing_pair) then
         call wing_rotation_matrix( Insect%M_b2w_l2, Insect%alpha_l2, Insect%theta_l2, Insect%phi_l2, Insect%eta_stroke, "left" )
         call wing_rotation_matrix( Insect%M_b2w_r2, Insect%alpha_r2, Insect%theta_r2, Insect%phi_r2, Insect%eta_stroke, "right" )
      endif

      ! inverse of the body rotation matrices
      Insect%M_b2g = transpose(Insect%M_g2b)

      ! body angular velocity vector in b/g coordinate system
      call body_angular_velocity( Insect, Insect%rot_body_b, Insect%rot_body_g, Insect%M_g2b )
      ! rel+abs wing angular velocities in the w/b/g coordinate system
      call wing_angular_velocities ( time, Insect, Insect%M_g2b )
      ! angular acceleration for wings (required for inertial power)
      call wing_angular_accel( time, Insect )

      !-----------------------------------------------------------------------------
      ! vector from body centre to left/right pivot point in global reference frame,
      ! for aerodynamic power
      ! NOTE: this qty is used to compute aerody. moments in FLUSI, something odd happens and
      ! we need the relative vector. In wabbit, we need the absolute vector.
      !-----------------------------------------------------------------------------
      Insect%x_pivot_l_g = matmul(Insect%M_b2g, Insect%x_pivot_l_b) + Insect%xc_body_g
      Insect%x_pivot_r_g = matmul(Insect%M_b2g, Insect%x_pivot_r_b) + Insect%xc_body_g

      if (Insect%second_wing_pair) then
         Insect%x_pivot_l2_g = matmul(Insect%M_b2g, Insect%x_pivot_l2_b) + Insect%xc_body_g
         Insect%x_pivot_r2_g = matmul(Insect%M_b2g, Insect%x_pivot_r2_b) + Insect%xc_body_g
      endif

      if (first_call) then
         ! print some important numbers, routine exectuted only once during a simulation
         !    call print_insect_reynolds_numbers( Insect ) ! ERROR: this can only be called after wing shape initialization
         first_call = .false.
      endif

      ! save time to insect, then we can check if the update routine has been called
      ! or not (this is not necessary if Update_Insect is called, but helpful to prevent
      ! human errors)
      Insect%time = time

      if (Insect%BodyType == "superSTL") then
         ! if the body is an stl file, its data (surface triangulation) is in the body coordinate
         ! system. For the mask generation, we work on the global coordinate system (Eulerian grid)
         ! and thus we require the data from the superSTL file (coordinates and a bunch of different
         ! normal vectors) in the global system.
         ! The transformation is done here (and not in draw_body_superSTL) because it needs to be
         ! done only once, and draw_body_superSTL is called for every block. The speed-up, however,
         ! is of course limited.
         do i = 1, size(body_superSTL_b, 1)
            body_superSTL_g(i, 1:3)   = matmul(Insect%M_b2g, body_superSTL_b(i, 1:3)) + Insect%xc_body_g
            body_superSTL_g(i, 4:6)   = matmul(Insect%M_b2g, body_superSTL_b(i, 4:6)) + Insect%xc_body_g
            body_superSTL_g(i, 7:9)   = matmul(Insect%M_b2g, body_superSTL_b(i, 7:9)) + Insect%xc_body_g
            body_superSTL_g(i, 10:12) = matmul(Insect%M_b2g, body_superSTL_b(i, 10:12))
            body_superSTL_g(i, 13:15) = matmul(Insect%M_b2g, body_superSTL_b(i, 13:15))
            body_superSTL_g(i, 16:18) = matmul(Insect%M_b2g, body_superSTL_b(i, 16:18))
            body_superSTL_g(i, 19:21) = matmul(Insect%M_b2g, body_superSTL_b(i, 19:21))
            body_superSTL_g(i, 22:24) = matmul(Insect%M_b2g, body_superSTL_b(i, 22:24))
            body_superSTL_g(i, 25:27) = matmul(Insect%M_b2g, body_superSTL_b(i, 25:27))
            body_superSTL_g(i, 28:30) = matmul(Insect%M_b2g, body_superSTL_b(i, 28:30))
         enddo
      endif

   end subroutine Update_Insect

   !-------------------------------------------------------------------------------
   ! Main routine for drawing insects. Draws body and wings, parameters are in "INSECT"
   !-------------------------------------------------------------------------------
   subroutine Draw_Insect( time, Insect, xx0, ddx, mask, mask_color, us )
      implicit none

      real(kind=rk), intent(in)      :: time
      type(diptera), intent(inout)   :: Insect
      real(kind=rk), intent(in)      :: xx0(1:3), ddx(1:3)
      real(kind=rk), intent(inout)   :: mask(0:,0:,0:)
      real(kind=rk), intent(inout)   :: us(0:,0:,0:,1:)
      integer(kind=2), intent(inout) :: mask_color(0:,0:,0:)

      if ((dabs(Insect%time-time)>1.0d-10).and.root) then
         write(*,'("error! time=",es15.8," but Insect%time=",es15.8)') time, Insect%time
         write(*,'("Did you call Update_Insect before Draw_Insect?")')
      endif

      ! 28/01/2019: Thomas. Discovered that this was done block based, i.e. the smoothing layer
      ! had different thickness, if some blocks happened to be at different levels (and still carry
      ! a part of the smoothing layer.) I don't know if that made sense, because the layer shrinks/expands then
      ! and because it might be discontinous. Both options are included now, default is "as before"
      ! Insect%smoothing_thickness=="local"  : smoothing_layer = c_sm * 2**-J * L/(BS-1)
      ! Insect%smoothing_thickness=="global" : smoothing_layer = c_sm * 2**-Jmax * L/(BS-1)
      if (Insect%smoothing_thickness=="local") then
         Insect%smooth = Insect%C_smooth*maxval(ddx)
         Insect%safety = 3.5_rk*Insect%smooth
      endif

      ! delete old mask
      call delete_old_mask( time, mask, mask_color, us, Insect )

      !-----------------------------------------------------------------------------
      ! BODY. Now the body is special: if the insect does not move (or rotate), the
      ! body does not change in time. On the other hand, it is quite expensive to
      ! compute, since it involves a lot of points (volume), and it is a source of
      ! load balancing problems, since many cores do not draw the body at all.
      ! We thus try to draw it only once and then simply not to erase it later.
      !-----------------------------------------------------------------------------
      if (Insect%body_moves=="no" .and. .not. grid_time_dependent) then
         if (.not. Insect%body_already_drawn) then
            ! the body is at rest, but it is the first call to this routine, so draw it now
            if (root) then
               write(*,'(80("~"))')
               write(*,*) "Flag Insect%body_moves is no and we did not yet draw"
               write(*,*) "the body once: we do that now, and skip draw_insect_body"
               write(*,*) "from now on. time=", time
               write(*,'(80("~"))')
            endif
            call draw_insect_body( time, xx0, ddx, mask, mask_color, us, Insect, delete=.false.)
            Insect%body_already_drawn = .true.
         endif
      else
         ! the body moves, draw it
         call draw_insect_body( time, xx0, ddx, mask, mask_color, us, Insect, delete=.false.)
      endif

      !-----------------------------------------------------------------------------
      ! Wings
      !-----------------------------------------------------------------------------
      call draw_insect_wings( time, xx0, ddx, mask, mask_color, us, Insect, delete=.false.)

      ! this is a debug test, which succeeded.
      !call check_if_us_is_derivative_of_position_wingtip(time, Insect)
   end subroutine Draw_Insect



   !-------------------------------------------------------
   ! short for the smooth step function.
   ! the smooting is defined in Insect%smooth, here we need only x, and the
   ! thickness (i.e., in the limit, steps=1 if x<t and steps=0 if x>t
   !-------------------------------------------------------
   real(kind=rk) function steps(x, t, h)
      implicit none
      real(kind=rk) :: x, t, h
      ! f is 1 if x<=t-h
      ! f is 0 if x>t+h
      ! f is variable (smooth) in between
      if (x<=t-h) then
         steps = 1._rk
      elseif (((t-h)<x).and.(x<(t+h))) then
         steps = 0.5_rk*(1._rk+dcos((x-t+h)*pi/(2._rk*h)) )
      else
         steps = 0.0_rk
      endif

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
      real(kind=rk) :: factor_amp = 1.0_rk  ! Dmitry, 7 Nov 2013
      real(kind=rk) :: phase
      !!----------------------------

      !! d_amplitude
      dAmp = dsqrt(a**2 +b**2)*factor_amp

      !! phase
      if( b>0.0_rk ) then
         phase = datan(a/b)
      elseif( b<0.0_rk ) then
         phase = datan(a/b) +pi
      else !! b == 0 -> avoid division by zero
         phase = pi*0.5_rk !! sin(PI/2) = cos
      endif

      phase = phase + (shift_phase +initial_phase*2.0_rk*pi)*dble(F)

      !! d_angle
      dangle = dAmp*dsin( angles +phase )

      !! velocity increment (Dmitry, 7 Nov 2013)
      dangle_dt = 2.0_rk*pi*dble(F) * dAmp*dcos( angles +phase )

      return
   end subroutine get_dangle



   ! Compute aerodynamic power
   subroutine aero_power(Insect, apowtotal)
      implicit none

      integer :: color_body, color_l, color_r, color_l2, color_r2
      real(kind=rk), dimension(1:3) :: omrel, momrel
      real(kind=rk), intent(out) :: apowtotal
      type(diptera),intent(inout)::Insect

      ! colors for Diptera (one body, two wings)
      color_body = Insect%color_body
      color_l = Insect%color_l
      color_r = Insect%color_r

      ! body is not driven directly, therefore the power is set to zero
      Insect%PartIntegrals(color_body)%APow = 0.0_rk

      !-----------
      ! left wing
      !-----------
      ! relative angular velocity, in global system
      omrel = Insect%rot_rel_wing_l_g

      ! the aerodyn moment is computed in global system and readily w.r.t. insects hinge
      ! points (in statistics_ACM.f90)
      momrel = Insect%PartIntegrals(color_l)%Torque + &
         Insect%PartIntegrals(color_l)%Torque_unst

      ! aerodynamic power
      Insect%PartIntegrals(color_l)%APow = - sum( momrel * omrel )

      !-----------
      ! right wing
      !-----------
      ! relative angular velocity, in global system
      omrel = Insect%rot_rel_wing_r_g

      ! the aerodyn moment is computed in global system and readily w.r.t. insects hinge
      ! points (in statistics_ACM.f90)
      momrel = Insect%PartIntegrals(color_r)%Torque + &
         Insect%PartIntegrals(color_r)%Torque_unst

      ! aerodynamic power
      Insect%PartIntegrals(color_r)%APow = - sum( momrel * omrel )

      !-----------
      ! Total aerodynamic power
      !-----------
      apowtotal = Insect%PartIntegrals(color_body)%APow + &
         Insect%PartIntegrals(color_l)%APow + Insect%PartIntegrals(color_r)%APow

      !-----------
      ! if second wing pair is present
      !-----------
      if (Insect%second_wing_pair) then
         ! Colors
         color_l2 = Insect%color_l2
         color_r2 = Insect%color_r2

         ! left wing
         ! relative angular velocity, in global system
         omrel = Insect%rot_rel_wing_l2_g

         momrel = Insect%PartIntegrals(color_l2)%Torque + &
            Insect%PartIntegrals(color_l2)%Torque_unst

         ! aerodynamic power
         Insect%PartIntegrals(color_l2)%APow = - sum( momrel * omrel )

         ! right wing
         ! relative angular velocity, in global system
         omrel = Insect%rot_rel_wing_r2_g

         momrel = Insect%PartIntegrals(color_r2)%Torque + &
            Insect%PartIntegrals(color_r2)%Torque_unst

         ! aerodynamic power
         Insect%PartIntegrals(color_r2)%APow = - sum( momrel * omrel )

         ! Total aerodynamic power
         apowtotal = apowtotal + Insect%PartIntegrals(color_l2)%APow + &
            Insect%PartIntegrals(color_r2)%APow
      endif

   end subroutine aero_power


   !-------------------------------------------------------------------------------
   ! Compute interial power, i.e. the power the insect would have to invest
   ! when flapping its wings in vacuum.
   !
   ! OUTPUT:
   !       ipowtotal: total inertial power
   !       iwmoment_g: inertial force moment components of the wings about the hinge in the laboratory reference frame
   !                          first index - component: x, y, z
   !                          second index - 1:body (unused), 2:left wing, 3:right wing, 4:2nd left wing, 5:2nd right wing
   !       Insect%PartIntegrals%IPow: (global): individual inertial power
   !
   ! INPUT:
   !       Insect%rot_dt_wing_l_w (in wing system): left wing angular acceleration
   !       Insect%rot_dt_wing_r_w (in wing system): right wing angular acceleration
   !       Insect%Jxx,Jyy,Jxy,Jzz (in wing system) Wing inertia
   !       Insect%Jxx2,Jyy2,Jxy2,Jzz2 (in wing system) Wing inertia of the second pair
   !
   ! MATHEMATICS:
   !       P_inertia = omega*( J*omega_dt + omega \cross (J*omega) )
   !                 = omega*( a + omega \cross b )
   !       The interia tensor is (it is specified in the PARAMS file)
   !           / Jxx Jxy 0   \
   !       J = | Jxy Jyy 0   |
   !           \ 0   0   Jzz /
   !
   ! SEE ALSO
   !       Berman, Wang: Energy minimizing kinematics in hovering insect flight
   !       (JFM 582, 2007), eqn 2.22 (looks a bit different)
   !-------------------------------------------------------------------------------
   subroutine inert_power(Insect,ipowtotal,iwmoment_g)
      implicit none

      real(kind=rk), intent(out) :: ipowtotal
      real(kind=rk), dimension(1:3,1:5), intent(out) :: iwmoment_g
      real(kind=rk), dimension(1:3) :: a,b
      real(kind=rk), dimension(1:3,1:5) :: iwmoment
      integer(kind=2) :: color_body, color_l, color_r, color_l2, color_r2
      type(diptera),intent(inout)::Insect

      iwmoment = 0
      iwmoment_g = 0

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

      iwmoment(1,color_l) = (a(1)+Insect%rot_rel_wing_l_w(2)*b(3)-Insect%rot_rel_wing_l_w(3)*b(2))
      iwmoment(2,color_l) = (a(2)+Insect%rot_rel_wing_l_w(3)*b(1)-Insect%rot_rel_wing_l_w(1)*b(3))
      iwmoment(3,color_l) = (a(3)+Insect%rot_rel_wing_l_w(1)*b(2)-Insect%rot_rel_wing_l_w(2)*b(1))

      Insect%PartIntegrals(color_l)%IPow = &
         Insect%rot_rel_wing_l_w(1) * iwmoment(1,color_l) + &
         Insect%rot_rel_wing_l_w(2) * iwmoment(2,color_l) + &
         Insect%rot_rel_wing_l_w(3) * iwmoment(3,color_l)

      !-- RIGHT WING
      a(1) = Insect%Jxx * Insect%rot_dt_wing_r_w(1) + Insect%Jxy * Insect%rot_dt_wing_r_w(2)
      a(2) = Insect%Jxy * Insect%rot_dt_wing_r_w(1) + Insect%Jyy * Insect%rot_dt_wing_r_w(2)
      a(3) = Insect%Jzz * Insect%rot_dt_wing_r_w(3)

      b(1) = Insect%Jxx * Insect%rot_rel_wing_r_w(1) + Insect%Jxy * Insect%rot_rel_wing_r_w(2)
      b(2) = Insect%Jxy * Insect%rot_rel_wing_r_w(1) + Insect%Jyy * Insect%rot_rel_wing_r_w(2)
      b(3) = Insect%Jzz * Insect%rot_rel_wing_r_w(3)

      iwmoment(1,color_r) = (a(1)+Insect%rot_rel_wing_r_w(2)*b(3)-Insect%rot_rel_wing_r_w(3)*b(2))
      iwmoment(2,color_r) = (a(2)+Insect%rot_rel_wing_r_w(3)*b(1)-Insect%rot_rel_wing_r_w(1)*b(3))
      iwmoment(3,color_r) = (a(3)+Insect%rot_rel_wing_r_w(1)*b(2)-Insect%rot_rel_wing_r_w(2)*b(1))

      Insect%PartIntegrals(color_r)%IPow = &
         Insect%rot_rel_wing_r_w(1) * iwmoment(1,color_r) + &
         Insect%rot_rel_wing_r_w(2) * iwmoment(2,color_r) + &
         Insect%rot_rel_wing_r_w(3) * iwmoment(3,color_r)

      ipowtotal = Insect%PartIntegrals(color_r)%IPow + Insect%PartIntegrals(color_l)%IPow

      !-----------
      ! if second wing pair is present
      !-----------
      if (Insect%second_wing_pair) then
         ! colors
         color_l2 = Insect%color_l2
         color_r2 = Insect%color_r2

         ! second left wing
         a(1) = Insect%Jxx2 * Insect%rot_dt_wing_l2_w(1) + Insect%Jxy2 * Insect%rot_dt_wing_l2_w(2)
         a(2) = Insect%Jxy2 * Insect%rot_dt_wing_l2_w(1) + Insect%Jyy2 * Insect%rot_dt_wing_l2_w(2)
         a(3) = Insect%Jzz2 * Insect%rot_dt_wing_l2_w(3)

         b(1) = Insect%Jxx2 * Insect%rot_rel_wing_l2_w(1) + Insect%Jxy2 * Insect%rot_rel_wing_l2_w(2)
         b(2) = Insect%Jxy2 * Insect%rot_rel_wing_l2_w(1) + Insect%Jyy2 * Insect%rot_rel_wing_l2_w(2)
         b(3) = Insect%Jzz2 * Insect%rot_rel_wing_l2_w(3)

         iwmoment(1,color_l2) = (a(1)+Insect%rot_rel_wing_l2_w(2)*b(3)-Insect%rot_rel_wing_l2_w(3)*b(2))
         iwmoment(2,color_l2) = (a(2)+Insect%rot_rel_wing_l2_w(3)*b(1)-Insect%rot_rel_wing_l2_w(1)*b(3))
         iwmoment(3,color_l2) = (a(3)+Insect%rot_rel_wing_l2_w(1)*b(2)-Insect%rot_rel_wing_l2_w(2)*b(1))

         Insect%PartIntegrals(color_l2)%IPow = &
            Insect%rot_rel_wing_l2_w(1) * iwmoment(1,color_l2) + &
            Insect%rot_rel_wing_l2_w(2) * iwmoment(2,color_l2) + &
            Insect%rot_rel_wing_l2_w(3) * iwmoment(3,color_l2)

         ! second right wing
         a(1) = Insect%Jxx2 * Insect%rot_dt_wing_r2_w(1) + Insect%Jxy2 * Insect%rot_dt_wing_r2_w(2)
         a(2) = Insect%Jxy2 * Insect%rot_dt_wing_r2_w(1) + Insect%Jyy2 * Insect%rot_dt_wing_r2_w(2)
         a(3) = Insect%Jzz2 * Insect%rot_dt_wing_r2_w(3)

         b(1) = Insect%Jxx2 * Insect%rot_rel_wing_r2_w(1) + Insect%Jxy2 * Insect%rot_rel_wing_r2_w(2)
         b(2) = Insect%Jxy2 * Insect%rot_rel_wing_r2_w(1) + Insect%Jyy2 * Insect%rot_rel_wing_r2_w(2)
         b(3) = Insect%Jzz2 * Insect%rot_rel_wing_r2_w(3)

         iwmoment(1,color_r2) = (a(1)+Insect%rot_rel_wing_r2_w(2)*b(3)-Insect%rot_rel_wing_r2_w(3)*b(2))
         iwmoment(2,color_r2) = (a(2)+Insect%rot_rel_wing_r2_w(3)*b(1)-Insect%rot_rel_wing_r2_w(1)*b(3))
         iwmoment(3,color_r2) = (a(3)+Insect%rot_rel_wing_r2_w(1)*b(2)-Insect%rot_rel_wing_r2_w(2)*b(1))

         Insect%PartIntegrals(color_r)%IPow = &
            Insect%rot_rel_wing_r2_w(1) * iwmoment(1,color_r2) + &
            Insect%rot_rel_wing_r2_w(2) * iwmoment(2,color_r2) + &
            Insect%rot_rel_wing_r2_w(3) * iwmoment(3,color_r2)

         ! total
         ipowtotal = ipowtotal + Insect%PartIntegrals(color_r2)%IPow + &
            Insect%PartIntegrals(color_l2)%IPow
      endif

      ! TODO: calculate and include rot_dt_body_b in the insect module
      !-- BODY
      !a(1) = Insect%Jroll_body * Insect%rot_dt_body_b(1)
      !a(2) = Insect%Jpitch_body * Insect%rot_dt_body_b(2)
      !a(3) = Insect%Jyaw_body * Insect%rot_dt_body_b(3)

      !b(1) = Insect%Jroll_body * Insect%rot_body_b(1)
      !b(2) = Insect%Jpitch_body * Insect%rot_body_b(2)
      !b(3) = Insect%Jyaw_body * Insect%rot_body_b(3)

      !iwmoment(1,color_b) = (a(1)+Insect%rot_body_b(2)*b(3)-Insect%rot_body_b(3)*b(2))
      !iwmoment(2,color_b) = (a(2)+Insect%rot_body_b(3)*b(1)-Insect%rot_body_b(1)*b(3))
      !iwmoment(3,color_b) = (a(3)+Insect%rot_body_b(1)*b(2)-Insect%rot_body_b(2)*b(1))

      ! transform into the laboratory reference frame
      !iwmoment_g(:,color_b) = -matmul(Insect%M_b2g,iwmoment(:,color_b))
      iwmoment_g(:,color_l) = -matmul(Insect%M_b2g,matmul(transpose(Insect%M_b2w_l),iwmoment(:,color_l)))
      iwmoment_g(:,color_r) = -matmul(Insect%M_b2g,matmul(transpose(Insect%M_b2w_r),iwmoment(:,color_r)))
      if (Insect%second_wing_pair) then
         iwmoment_g(:,color_l2) = -matmul(Insect%M_b2g,matmul(transpose(Insect%M_b2w_l2),iwmoment(:,color_l2)))
         iwmoment_g(:,color_r2) = -matmul(Insect%M_b2g,matmul(transpose(Insect%M_b2w_r2),iwmoment(:,color_r2)))
      endif

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
   subroutine body_angular_velocity( Insect, rot_body_b, rot_body_g, M_g2b )
      implicit none

      type(diptera), intent(inout) :: Insect
      real(kind=rk), intent(in) :: M_g2b(1:3,1:3)
      real(kind=rk), dimension(1:3), intent(out) :: rot_body_b, rot_body_g
      real(kind=rk) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt

      psi = Insect%psi
      beta = Insect%beta
      gamma = Insect%gamma
      psi_dt = Insect%psi_dt
      beta_dt = Insect%beta_dt
      gamma_dt = Insect%gamma_dt

      if ( Insect%quaternion_solver_used ) then
         ! variant (a)
         rot_body_b = Insect%rot_body_b ! copy (useless, actually, but required for interface)
         rot_body_g = matmul( transpose(M_g2b), rot_body_b)

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
   subroutine wing_angular_velocities ( time, Insect, M_g2b )
      implicit none

      real(kind=rk), intent(in) :: time
      real(kind=rk), intent(in) :: M_g2b(1:3,1:3)
      type(diptera), intent(inout) :: Insect

      real(kind=rk) :: eta_stroke
      real(kind=rk) :: phi_r, alpha_r, theta_r, phi_dt_r, alpha_dt_r, theta_dt_r
      real(kind=rk) :: phi_l, alpha_l, theta_l, phi_dt_l, alpha_dt_l, theta_dt_l
      real(kind=rk), dimension(1:3) :: rot_l_alpha, rot_l_theta, rot_l_phi, &
         rot_r_alpha, rot_r_theta, rot_r_phi
      real(kind=rk), dimension(1:3,1:3) :: M_b2w_l, M_b2w_r, &
         M1_tmp, M2_tmp, M1_l, M2_l, M3_l, M1_r, M2_r, M3_r, &
         M_b2s_l, M_b2s_r

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
      M_b2s_l = M1_tmp

      call Rx(M1_tmp,pi)
      call Ry(M2_tmp,eta_stroke)
      M_b2s_r = matmul(M1_tmp,M2_tmp)

      call Ry(M1_l,alpha_l)
      call Rz(M2_l,theta_l)   ! Order changed (Dmitry, 7 Nov 2013)
      call Rx(M3_l,phi_l)
      M_b2w_l = matmul(M1_l,matmul(M2_l,matmul(M3_l,M_b2s_l)))

      ! note the coordinate system is rotated so we don't need to inverse the sign
      ! of theta, and the wings still rotate in opposite direction
      call Ry(M1_r,-alpha_r)
      call Rz(M2_r,theta_r)   ! Order changed (Dmitry, 7 Nov 2013)
      call Rx(M3_r,-phi_r)
      M_b2w_r = matmul(M1_r,matmul(M2_r,matmul(M3_r,M_b2s_r)))

      !-----------------------------------------------------------------------------
      ! angular velocity vectors (in wing system)
      !-----------------------------------------------------------------------------
      rot_l_alpha = (/ 0.0_rk, alpha_dt_l, 0.0_rk /)
      rot_l_theta = (/ 0.0_rk, 0.0_rk, theta_dt_l /)
      rot_l_phi   = (/ phi_dt_l, 0.0_rk, 0.0_rk   /)
      rot_r_alpha = (/ 0.0_rk, -alpha_dt_r, 0.0_rk/)
      rot_r_theta = (/ 0.0_rk, 0.0_rk, theta_dt_r /)
      rot_r_phi   = (/ -phi_dt_r, 0.0_rk, 0.0_rk  /)

      ! in the wing coordinate system
      Insect%rot_rel_wing_l_w = matmul(M_b2w_l, matmul(transpose(M_b2s_l),matmul(transpose(M3_l), &
         rot_l_phi+matmul(transpose(M2_l),rot_l_theta+matmul(transpose(M1_l), &
         rot_l_alpha)))))

      Insect%rot_rel_wing_r_w = matmul(M_b2w_r,matmul(transpose(M_b2s_r),matmul(transpose(M3_r), &
         rot_r_phi+matmul(transpose(M2_r),rot_r_theta+matmul(transpose(M1_r), &
         rot_r_alpha)))))

      ! direct definition, equivalent to what is above.
      ! Insect%rot_rel_wing_l_w = (/phi_dt_l*cos(alpha_l)*cos(theta_l)-theta_dt_l*sin(alpha_l),&
      !   alpha_dt_l-phi_dt_l*sin(theta_l),&
      !   theta_dt_l*cos(alpha_l)+phi_dt_l*sin(alpha_l)*cos(theta_l)/)

      ! prior to the call of this routine, the routine body_angular_velocity has
      ! computed the body angular velocity (both g/b frames) so here now we can also
      ! compute global and absolute wing angular velocities.
      Insect%rot_rel_wing_l_b = matmul( transpose(M_b2w_l), Insect%rot_rel_wing_l_w )
      Insect%rot_rel_wing_r_b = matmul( transpose(M_b2w_r), Insect%rot_rel_wing_r_w )

      Insect%rot_rel_wing_l_g = matmul( transpose(M_g2b), Insect%rot_rel_wing_l_b )
      Insect%rot_rel_wing_r_g = matmul( transpose(M_g2b), Insect%rot_rel_wing_r_b )

      Insect%rot_abs_wing_l_g = Insect%rot_body_g + Insect%rot_rel_wing_l_g
      Insect%rot_abs_wing_r_g = Insect%rot_body_g + Insect%rot_rel_wing_r_g

      !-----------
      ! if second wing pair is present
      !-----------
      if (Insect%second_wing_pair) then
         ! setup the local wing parameters
         phi_r      = Insect%phi_r2
         alpha_r    = Insect%alpha_r2
         theta_r    = Insect%theta_r2
         phi_dt_r   = Insect%phi_dt_r2
         alpha_dt_r = Insect%alpha_dt_r2
         theta_dt_r = Insect%theta_dt_r2
         phi_l      = Insect%phi_l2
         alpha_l    = Insect%alpha_l2
         theta_l    = Insect%theta_l2
         phi_dt_l   = Insect%phi_dt_l2
         alpha_dt_l = Insect%alpha_dt_l2
         theta_dt_l = Insect%theta_dt_l2

         ! second wing pair rotation matrices
         call Ry(M1_l,alpha_l)
         call Rz(M2_l,theta_l)
         call Rx(M3_l,phi_l)
         M_b2w_l = matmul(M1_l,matmul(M2_l,matmul(M3_l,M_b2s_l)))
         call Ry(M1_r,-alpha_r)
         call Rz(M2_r,theta_r)
         call Rx(M3_r,-phi_r)
         M_b2w_r = matmul(M1_r,matmul(M2_r,matmul(M3_r,M_b2s_r)))

         ! angular velocity vectors (in wing system)
         rot_l_alpha = (/ 0.0_rk, alpha_dt_l, 0.0_rk /)
         rot_l_theta = (/ 0.0_rk, 0.0_rk, theta_dt_l /)
         rot_l_phi   = (/ phi_dt_l, 0.0_rk, 0.0_rk   /)
         rot_r_alpha = (/ 0.0_rk, -alpha_dt_r, 0.0_rk/)
         rot_r_theta = (/ 0.0_rk, 0.0_rk, theta_dt_r /)
         rot_r_phi   = (/ -phi_dt_r, 0.0_rk, 0.0_rk  /)

         ! in the wing coordinate system
         Insect%rot_rel_wing_l2_w = matmul(M_b2w_l,matmul(transpose(M_b2s_l),matmul(transpose(M3_l), &
            rot_l_phi+matmul(transpose(M2_l),rot_l_theta+matmul(transpose(M1_l), &
            rot_l_alpha)))))
         Insect%rot_rel_wing_r2_w = matmul(M_b2w_r,matmul(transpose(M_b2s_r),matmul(transpose(M3_r), &
            rot_r_phi+matmul(transpose(M2_r),rot_r_theta+matmul(transpose(M1_r), &
            rot_r_alpha)))))

         ! wing angular velocities
         Insect%rot_rel_wing_l2_b = matmul( transpose(M_b2w_l), Insect%rot_rel_wing_l2_w )
         Insect%rot_rel_wing_r2_b = matmul( transpose(M_b2w_r), Insect%rot_rel_wing_r2_w )
         Insect%rot_rel_wing_l2_g = matmul( transpose(M_g2b), Insect%rot_rel_wing_l2_b )
         Insect%rot_rel_wing_r2_g = matmul( transpose(M_g2b), Insect%rot_rel_wing_r2_b )
         Insect%rot_abs_wing_l2_g = Insect%rot_body_g + Insect%rot_rel_wing_l2_g
         Insect%rot_abs_wing_r2_g = Insect%rot_body_g + Insect%rot_rel_wing_r2_g
      endif
   end subroutine wing_angular_velocities



   !-------------------------------------------------------------------------------
   ! Numerically estimate (it's a very precise estimation) the angular acceleration
   ! vectors for both wings, using one-sided finite differences (in the future dt)
   ! NOTE: this routine requires us to be able to evaluate both body and wing state
   !       at arbitrary times.
   !-------------------------------------------------------------------------------
   subroutine wing_angular_accel( time, Insect )
      implicit none
      real(kind=rk), intent(in) :: time
      type(diptera), intent(inout) :: Insect

      real(kind=rk) :: M_g2b(1:3,1:3), rot_dt_wing_g(1:3), M_b2w_r(1:3,1:3), M_b2w_l(1:3,1:3)
      real(kind=rk) :: M_b2w_r2(1:3,1:3), M_b2w_l2(1:3,1:3)
      type(diptera) :: Insect2
      real(kind=rk) :: dt,t
      integer(kind=2) :: wingID

      dt = 1.0e-8_rk
      Insect2 = Insect

      Insect%rot_dt_wing_l_w = 0.0_rk
      Insect%rot_dt_wing_r_w = 0.0_rk

      ! fetch motion state at time+dt
      call BodyMotion (time+dt, Insect2)
      wingID = 1
      call FlappingMotionWrap(time+dt, Insect2, wingID)
      wingID = 2
      call FlappingMotionWrap(time+dt, Insect2, wingID)

      if (Insect%second_wing_pair) then
         wingID = 3
         call FlappingMotionWrap(time+dt, Insect2, wingID)
         wingID = 4
         call FlappingMotionWrap(time+dt, Insect2, wingID)
      endif
      call StrokePlane (time+dt, Insect2)
      call body_rotation_matrix( Insect2, M_g2b )
      call wing_angular_velocities ( time+dt, Insect2, M_g2b )

!??????????????????????????????? why not Insect2 ???!??!
!??????????????????????????????? why not Insect2 ???!??!
!??????????????????????????????? why not Insect2 ???!??!
!??????????????????????????????? why not Insect2 ???!??!
!??????????????????????????????? why not Insect2 ???!??!
!??????????????????????????????? why not Insect2 ???!??!
      ! this is the current state:
      call body_rotation_matrix( Insect, M_g2b )
      call wing_rotation_matrix( M_b2w_l, Insect%alpha_l, Insect%theta_l, Insect%phi_l, Insect%eta_stroke, "left" )
      call wing_rotation_matrix( M_b2w_r, Insect%alpha_r, Insect%theta_r, Insect%phi_r, Insect%eta_stroke, "right" )
      if (Insect%second_wing_pair) then
         call wing_rotation_matrix( M_b2w_l2, Insect%alpha_l2, Insect%theta_l2, Insect%phi_l2, Insect%eta_stroke, "left" )
         call wing_rotation_matrix( M_b2w_r2, Insect%alpha_r2, Insect%theta_r2, Insect%phi_r2, Insect%eta_stroke, "right" )
      endif

      ! use one-sided finite differences to derive the absolute angular velocity with
      ! respect to time. Note in older code versions, this was wrong, as we derived
      ! the ang. vel. in the wing coordinate system, which is a moving reference frame.

      ! now happily, the old results are still correct, as long as the body does not rotate
      ! see, e.g., this document https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjzl-XP6_LNAhWoC5oKHUdDCHwQFggoMAI&url=http%3A%2F%2Focw.mit.edu%2Fcourses%2Faeronautics-and-astronautics%2F16-07-dynamics-fall-2009%2Flecture-notes%2FMIT16_07F09_Lec08.pdf&usg=AFQjCNHzEB-n_NMm6K3J1eRpIaGnuKpW0Q&sig2=yEPNin3bL5DnWauNJk2hcw&bvm=bv.126993452,d.bGs&cad=rjt
      ! however, if the body moves, an additional term occurs, and this was indeed missing
      ! in previous results.
      rot_dt_wing_g = (Insect2%rot_rel_wing_l_g - Insect%rot_rel_wing_l_g) / dt
      Insect%rot_dt_wing_l_w = matmul(M_b2w_l,matmul(M_g2b, rot_dt_wing_g))

      rot_dt_wing_g = (Insect2%rot_rel_wing_r_g - Insect%rot_rel_wing_r_g) / dt
      Insect%rot_dt_wing_r_w = matmul(M_b2w_r,matmul(M_g2b, rot_dt_wing_g))

      ! if second wing pair is present
      if (Insect%second_wing_pair) then
         rot_dt_wing_g = (Insect2%rot_rel_wing_l2_g - Insect%rot_rel_wing_l2_g) / dt
         Insect%rot_dt_wing_l2_w = matmul(M_b2w_l2,matmul(M_g2b, rot_dt_wing_g))

         rot_dt_wing_g = (Insect2%rot_rel_wing_r2_g - Insect%rot_rel_wing_r2_g) / dt
         Insect%rot_dt_wing_r2_w = matmul(M_b2w_r2,matmul(M_g2b, rot_dt_wing_g))
      endif

   end subroutine wing_angular_accel


   subroutine delete_old_mask( time, mask, mask_color, us, Insect )
      implicit none

      real(kind=rk), intent(in) :: time
      type(diptera),intent(in) :: Insect
      real(kind=rk),intent(inout) :: mask(0:,0:,0:)
      real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
      integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
      integer(kind=2) :: color_body, color_l, color_r, color_l2, color_r2
      logical, save :: cleaned_already_once = .false.

      ! colors for Diptera (one body, two wings)
      color_body = Insect%color_body
      color_l = Insect%color_l
      color_r = Insect%color_r
      color_l2 = Insect%color_l2
      color_r2 = Insect%color_r2

      !-----------------------------------------------------------------------------
      ! delete old mask
      !-----------------------------------------------------------------------------
      if (Insect%body_moves=="no" .and. .not. grid_time_dependent .and. cleaned_already_once) then
         ! the body is at rest, so we will not draw it. Delete the wings, as they move.
         where (mask_color==color_l .or. mask_color==color_r .or. &
            mask_color==color_l2 .or. mask_color==color_r2)
            mask = 0.0_rk
            mask_color = 0
         end where
         ! as the body rests it has no solid body velocity, which means we can safely
         ! reset the velocity everywhere (this step is actually unnessesary, but for
         ! safety we do it as well)
         us = 0.0_rk
      else
         ! the body of the insect moves, so we will construct the entire insect in this
         ! (and any other) call, and therefore we can safely reset the entire mask to zeros.
         mask = 0.0_rk
         mask_color = 0
         us = 0.0_rk
      endif

      cleaned_already_once = .true.

   end subroutine delete_old_mask


   !-----------------------------------------------------------------------------
   ! return the body rotation matrix
   !-----------------------------------------------------------------------------
   subroutine body_rotation_matrix( Insect, M_g2b )
      implicit none

      type(diptera),intent(inout) :: Insect
      real(kind=rk),intent(out) :: M_g2b(1:3,1:3)
      real(kind=rk), dimension(1:3,1:3) :: M1_b, M2_b, M3_b

      if (Insect%quaternion_solver_used) then
         ! entries 7,8,9,10 of the Insect%STATE vector are the body quaternion
         call rotation_matrix_from_quaternion( Insect%STATE(7:10), M_g2b)

      else
         ! conventional yaw, pitch, roll. Note the order of matrices is important.
         ! first we yaw, then we pitch, then we roll the insect. Note that when the
         ! free-flight solver is used, this matrix is obtained from quaternions, and
         ! not as a product of simple rotaion matrices. The latter can cause "gimbal-lock"
         call Rx(M1_b,Insect%psi)
         call Ry(M2_b,Insect%beta)
         call Rz(M3_b,Insect%gamma)
         M_g2b = matmul(M1_b,matmul(M2_b,M3_b))
      endif
   end subroutine body_rotation_matrix




   !-----------------------------------------------------------------------------
   ! returns the wing rotation matrix M_b2w for a given wing
   !-----------------------------------------------------------------------------
   subroutine wing_rotation_matrix( M_b2w, alpha, theta, phi, eta_stroke, side )
      implicit none

      real(kind=rk),intent(out) :: M_b2w(1:3,1:3)
      real(kind=rk), intent(in) :: alpha, theta, phi, eta_stroke
      ! as this is called rarely okay to use string comparison
      character(len=*), intent(in) :: side ! left or right (left2 and right2 not accepted - just pass the different angles)
      real(kind=rk), dimension(1:3,1:3) :: M1, M2, M3, M_b2s


      if (side == "right") then
         call Rx(M1, pi)
         call Ry(M2, eta_stroke)
         M_b2s = matmul(M1,M2)

         ! note the coordinate system is rotated so we don't need to inverse the sign
         ! of theta, and the wings still rotate in opposite direction
         call Ry(M1, -alpha)
         call Rz(M2,  theta)
         call Rx(M3, -phi)
         M_b2w = matmul(M1,matmul(M2,matmul(M3,M_b2s)))

      elseif (side == "left") then
         call Ry(M1,  eta_stroke)
         M_b2s = M1

         call Ry(M1, alpha)
         call Rz(M2, theta)
         call Rx(M3, phi)
         M_b2w = matmul(M1,matmul(M2,matmul(M3,M_b2s)))

      else
         ! left or right (left2 and right2 not accepted - just pass the different angles)
         call abort(1071024, "insect_module: neither right nor left side ? how many sides does an insect have? seven!?")

      endif
   end subroutine wing_rotation_matrix


   !-----------------------------------------------------------------------------
   ! Compute and print a couple of important numbers for insects
   !-----------------------------------------------------------------------------
   subroutine  print_insect_reynolds_numbers( Insect )
      implicit none
      type(diptera),intent(inout) :: Insect
      !type(diptera) :: Insect_copy ! commented out due to segmentation problems
      real(kind=rk) :: area, Re_f, Re
      real(kind=rk) :: time, dt
      real(kind=rk) :: phil_min, phil_max, phir_min, phir_max
      logical, save :: first_call = .true.
      integer(kind=2) :: wingID

      ! the second call is just a return statement
      if ( first_call .eqv. .false.) return

      ! only root does this...
      if (root) then
         ! we need the wing area to compute the mean wing chord
         wingID = 1 ! use the left wing area (wingID=1)
         call compute_wing_surface(Insect, wingID, area)
         write(*,'(50("~"))')
         write(*,'("Wing area is A=",g15.8)') area
         write(*,'("Mean chord length is c_m=",g15.8)') area/1._rk ! note c_m = A/R but R=1

         ! first we compute the stroke amplitude of the positional angle phi (for
         ! both wings). for safety, we make a copy of the insect, since the routines
         ! for the flapping motion write to this object, and we want to prevent any
         ! unwanted side effects
         !Insect_copy = Insect ! save the insect state
         time = 0.0_rk
         dt = 1.0d-3
         phil_min = 0.0_rk
         phil_max = 0.0_rk
         phir_min = 0.0_rk
         phir_max = 0.0_rk
         ! we use only one stroke ( the first one )
         do while (time < 1._rk)
            wingID = 1
            call FlappingMotionWrap ( time, Insect, wingID )
!!          call FlappingMotionWrap ( time, Insect_copy, wingID )
            wingID = 2
            call FlappingMotionWrap ( time, Insect, wingID )
!!          call FlappingMotionWrap ( time, Insect_copy, wingID )
            phil_min = min( phil_min, Insect%phi_l )
            phil_max = max( phil_max, Insect%phi_l )
            phir_min = min( phir_min, Insect%phi_r )
            phir_max = max( phir_max, Insect%phi_r )
            time = time + dt
         end do
         write(*,'("All following quantities are based on the first stroke 0.0 <= t <= 1.0")')
         write(*,'("Stroke amplitude is PHI_L=",g15.8)') abs(phil_max) + abs(phil_min)
         write(*,'("Stroke amplitude is PHI_R=",g15.8)') abs(phir_max) + abs(phir_min)
         write(*,'("Re_left  = 2*phi*R*f*c_m / nu =",g15.8)') 2._rk*(abs(phil_max)+abs(phil_min))*1._rk*1._rk*area/nu
         write(*,'("Re_right = 2*phi*R*f*c_m / nu =",g15.8)') 2._rk*(abs(phir_max)+abs(phir_min))*1._rk*1._rk*area/nu
         write(*,'("Re_f = R*R*f / nu =",g15.8)') 1._rk/nu
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
   subroutine read_insect_STATE_from_file(time, Insect, fname, verbose)
      implicit none
      real(kind=rk), intent(in) :: time
      type(diptera), intent(inout) :: Insect
      character(len=*), intent(in) :: fname
      logical, intent(in) :: verbose

      integer :: num_lines, n_header = 1, i
      character(len=maxcolumns) :: dummy
      real(kind=rk), allocatable, save :: data1(:,:)

      if ( .not. allocated(data1) ) then
         if (root) write(*,*) "read_insect_STATE_from_file:", trim(adjustl(fname))
         ! read rigidsolidsolver.t file
         ! skip header, count lines, read
         call count_lines_in_ascii_file_mpi(fname, num_lines, n_header)
         ! read contents of file
         allocate( data1(1:num_lines,1:14))
         call read_array_from_ascii_file_mpi(fname, data1 , n_header)
      endif

      ! interpolate in time
      i = 1
      do while (data1(i,1) <= time .and. i<size(data1,1)-1)
         i=i+1
      enddo

      ! we now have data1(i-1,1) <= time < data1(i,1)
      ! use linear interpolation
      Insect%STATE = 0.0_rk
      Insect%STATE(1:13) = data1(i-1,2:14) + (time - data1(i-1,1)) * (data1(i,2:14)-data1(i-1,2:14)) / (data1(i,1)-data1(i-1,1))

      if (root .and. verbose) then
         write(*,*) "The extracted Insect%STATE vector is:"
         write(*,'(21(es12.4,1x))') time, Insect%STATE(1:13)
      endif
   end subroutine


   ! this routine computes  the wingteip velocity in the global system by two means:
   ! one we compute the global position vector (which we derive wrt time in postprocessing)
   ! and once the cross-products of rotation ang velocities as it is done in the actual code.
   ! we checked: both agree, also with imposed body velocity.
   subroutine check_if_us_is_derivative_of_position_wingtip(time, Insect)
      real(kind=rk), intent(in) :: time
      type(diptera), intent(inout) :: Insect

      real(kind=rk) :: M_g2b(1:3,1:3), M_b2w_r(1:3,1:3), x_tip_w(1:3), x_tip_b(1:3), x_tip_g(1:3), &
         us_tip_g(1:3), v_tmp(1:3), v_tmp_b(1:3)
      real(kind=rk)::xd,yd,zd
      real(kind=rk)::c00,c10,c01,c11,c0,c1
      integer :: ix,iy,iz

      call body_rotation_matrix( Insect, M_g2b )
      ! call wing_right_rotation_matrix( Insect, M_b2w_r )
      call wing_rotation_matrix( M_b2w_r, Insect%alpha_r, Insect%theta_r, Insect%phi_r, Insect%eta_stroke, "right" )
      ! body angular velocity vector in b/g coordinate system
      call body_angular_velocity( Insect, Insect%rot_body_b, Insect%rot_body_g, M_g2b )
      ! rel+abs wing angular velocities in the w/b/g coordinate system
      call wing_angular_velocities ( time, Insect, M_g2b )

      x_tip_w = (/ 1.0_rk, 1.0_rk, 1.0_rk /)
      x_tip_b = matmul( transpose(M_b2w_r), x_tip_w ) + Insect%x_pivot_r_b
      x_tip_g = matmul( transpose(M_g2b), x_tip_b ) + Insect%xc_body_g

      !-----------------------------------------------------------------------------
      ! now we extrcat how the us field is constructed
      v_tmp(1) = Insect%rot_rel_wing_r_w(2)*x_tip_w(3)-Insect%rot_rel_wing_r_w(3)*x_tip_w(2)
      v_tmp(2) = Insect%rot_rel_wing_r_w(3)*x_tip_w(1)-Insect%rot_rel_wing_r_w(1)*x_tip_w(3)
      v_tmp(3) = Insect%rot_rel_wing_r_w(1)*x_tip_w(2)-Insect%rot_rel_wing_r_w(2)*x_tip_w(1)
      v_tmp_b = matmul(transpose(M_b2w_r), v_tmp) ! in body system


      ! translational part. we compute the rotational part in the body
      ! reference frame, therefore, we must transform the body translation
      ! velocity Insect%vc (which is in global coordinates) to the body frame
      v_tmp = matmul(M_g2b,Insect%vc_body_g)

      ! add solid body rotation to the translational velocity field. Note
      ! that rot_body_b and x_body are in the body reference frame
      v_tmp(1) = v_tmp(1) + Insect%rot_body_b(2)*x_tip_b(3)-Insect%rot_body_b(3)*x_tip_b(2)
      v_tmp(2) = v_tmp(2) + Insect%rot_body_b(3)*x_tip_b(1)-Insect%rot_body_b(1)*x_tip_b(3)
      v_tmp(3) = v_tmp(3) + Insect%rot_body_b(1)*x_tip_b(2)-Insect%rot_body_b(2)*x_tip_b(1)

      ! the body motion is added to the wing motion, which is already in us
      ! and they are also in the body refrence frame. However, us has to be
      ! in the global reference frame, so M_b2g is applied
      us_tip_g = matmul( transpose(M_g2b), v_tmp_b + v_tmp )

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
         if (Insect%second_wing_pair) then
            write (17,'(44(es15.8,1x))') time, Insect%xc_body_g, Insect%psi, Insect%beta, &
               Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
               Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
               Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
               Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w, &
               Insect%alpha_l2, Insect%phi_l2, Insect%theta_l2, &
               Insect%alpha_r2, Insect%phi_r2, Insect%theta_r2, &
               Insect%rot_rel_wing_l2_w, Insect%rot_rel_wing_r2_w, &
               Insect%rot_dt_wing_l2_w, Insect%rot_dt_wing_r2_w
         else
            write (17,'(26(es15.8,1x))') time, Insect%xc_body_g, Insect%psi, Insect%beta, &
               Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
               Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
               Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
               Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w
         endif
         close (17)
      endif

   end subroutine

end module module_insects
