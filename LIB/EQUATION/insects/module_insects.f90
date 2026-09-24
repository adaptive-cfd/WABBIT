module module_insects
   ! The following module includes required functions from FLUSI/WABBIT and hence
   ! makes this module here independent of that
   use module_insects_integration_flusi_wabbit
   use module_t_files
   use module_stl_file_reader
   use module_geometry

   implicit none

   ! I usually find it helpful to use the private keyword by itself initially, which specifies
   ! that everything within the module is private unless explicitly marked public.
   PRIVATE

   ! functions
   PUBLIC :: Draw_Insect, draw_insect_body, draw_insect_wings, insect_geometry_indicator, Update_All_Insects, insects_array_init, initialize_insect, &
      clean_all_insects, clean_insect, init_insect_data, write_insect_data, &
      aero_power, inert_power, read_insect_STATE_from_file, rigid_solid_init, rigid_solid_rhs, &
      BodyMotion, FlappingMotionWrap, body_rotation_matrix, wing_rotation_matrix
   ! type definitions
   PUBLIC :: diptera

   ! we use this so only root prints write statements...
   logical :: root = .false.

   ! 2025-11-07: insect periodization disabled, TE. It can be tricky to compute the properly periodized mask function 
   ! (eg bounding boxes) and in WABBIT, we'd just use a large enough domain in most cases. In particular, also STL periodization
   ! is probably painful.
   ! logical :: periodic_insect = .false.

   ! ghost nodes. If the insect module is used in a finite-differences code, then
   ! the data that we have often has ghost nodes, i.e. points that overlap and exist
   ! on several CPUS. On those, you normally would not create the mask (which is expensive)
   ! so we skip the first and last "g" points on the arrays used for mask creation
   integer, save :: g

   ! size (global) of domain
   real(kind=rk) :: xl, yl, zl
   ! viscosity (just for printing the Reynolds number), this can always be global because it is not insect-specific
   real(kind=rk) :: nu



   ! **********************************************************************************************
   ! **                              DEFINITION OF A WING                                        **
   ! **********************************************************************************************
   ! we have one wing type that regroups all parameters for a wing: shape, kinematics, etc
   type wing_type
      ! is the wing used or not
      logical :: used = .false.
      
      ! is this wing on the left or right side ?
      character(len=1) :: side = "X" ! can be "R" or "L" only, to ensure initialization set to "X"

      ! coordinates of shoulder (pivot) point
      real(kind=rk), dimension(1:3) :: x_pivot_g, x_pivot_b

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~                SHAPE / GEOMETRY                    ~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character(len=clong) :: WingShape = ""
      character(len=clong) :: WingShapeType = "fourier"

      ! ---- generic parameters -----
      real(kind=rk) :: b_top=0.0_rk, b_bot=0.0_rk, L_span=0.0_rk, WingThickness=0.0_rk
      integer(kind=2) :: color

      ! ---- wing shapes with polar coordinates description ----
      ! Fourier coefficients for wing shape
      real(kind=rk), allocatable :: ai_shape(:), bi_shape(:)
      real(kind=rk) :: a0_shape
      ! fill the R0(theta) array once, then only table-lookup (linear interpolation) is used instead of Fseries
      real(kind=rk), allocatable :: R0_table(:)
      ! the center of the polar coordinates, in the wing system
      real(kind=rk) :: xc=0.0_rk, yc=0.0_rk
      logical :: wings_radius_table_ready = .false.
      ! wing bounding box (xmin, xmax, ymin, ymax, zmin, zmax)
      ! in the wing system
      real(kind=rk) :: wing_bounding_box(1:6) = 0.0_rk

      ! ---- wing shapes with polygon description (no radius) ----
      ! polygon_wings(:,1) = x_w coordinates
      ! polygon_wings(:,2) = y_w coordinates
      real(kind=rk), allocatable :: polygon_wings(:,:)
      integer(kind=ik) :: n_polygon_points = 0

      ! ---- wing damage ----
      ! used for musca domestica wing damage project
      logical :: damaged = .false.
      ! wing damage mask
      real(kind=rk), allocatable, dimension(:,:) :: damage_mask

      ! ---- wing corrugation (non-constant z values) ----
      ! Published in T. Engels, H.-N. Wehmann, F.-O. Lehmann, Three-dimensional wing structure attenuates aerodynamic efficiency in flapping fly wings. J. R. Soc. Interface 17, 20190804, link (open access), 2020
      ! wing thickness profile
      real(kind=rk), allocatable, dimension(:,:) :: wing_thickness_profile
      character(len=clong) :: wing_thickness_distribution = "constant"
      real(kind=rk) :: wing_thickness_value = 0.05
      ! wing corrugation profile
      logical :: corrugated = .false.
      real(kind=rk), allocatable, dimension(:,:) :: corrugation_profile
      ! the 2D arrays damage_mask, corrugation_profile and wing_thickness_profile
      ! are equidistant grids, and their bounding box is all the same:
      real(kind=rk) :: corrugation_array_bbox(1:4)

      ! ---- bristles ----
      logical :: bristles = .false.
      logical :: bristles3D = .false.
      logical :: bristles_simplex = .false.
      integer(kind=ik) :: n_bristles
      real(kind=rk), allocatable :: bristles_coords(:,:)

      ! !KVN-2025>>>>>
      ! ! NOTE: prescribed wing deformation is untested work in progress! -TE 02/2026
      ! real(kind=rk), allocatable, dimension(:,:) :: deformations
      ! real(kind=rk), allocatable, dimension(:,:) :: deformation_profile
      ! integer(kind=ik) :: deformation_a, deformation_b, deformation_c
      ! !KVN-2025<<<<<

      ! ! wing signed distance function, if the 3d-interpolation approach is used.
      ! ! This is useful for highly complex wings, where one generates the mask only once
      ! ! and then interpolates the values to the global grid. note the data allocated
      ! ! here is of course understood in the wing system, so the linear transformation
      ! ! x_g -> x_w' is used, where x_w' is not a grid-aligned value x_w
      ! real(kind=rk), allocatable, dimension(:,:,:) :: mask_wing_complete
      ! real(kind=rk), dimension(1:3) :: mask_wing_xl, mask_wing_x0
      ! integer, dimension(1:3) :: mask_wing_nxyz
      ! integer :: mask_wing_safety=4

      ! moment of inertia of the wing in the wing system
      real(kind=rk) :: Jxx=0.0_rk, Jyy=0.0_rk, Jzz=0.0_rk, Jxy=0.0_rk

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~                KINEMATICS                          ~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character(len=clong) :: FlappingMotion = "NOT-SET-YET"

      ! Fourier coefficients, if kinematics given by Fourier series
      real(kind=rk) :: a0_alpha=0.0_rk, a0_phi=0.0_rk, a0_theta=0.0_rk
      real(kind=rk), allocatable :: ai_phi(:), bi_phi(:), ai_theta(:), bi_theta(:), ai_alpha(:), bi_alpha(:)
      integer :: nfft_phi=0, nfft_alpha=0, nfft_theta=0
      ! coefficients are read only once from file (or set differently)
      logical :: initialized = .false.
      ! some details about the file, if reading from ini file
      character(len=clong) :: infile_convention="", infile_type="", infile_units="", infile="NOT-USED"

      ! angles of the wing (flapping, feathering and deviation)
      real(kind=rk) :: phi=0.0_rk, alpha=0.0_rk, theta=0.0_rk, phi_dt=0.0_rk, alpha_dt=0.0_rk, theta_dt=0.0_rk

      ! angular velocities
      real(kind=rk), dimension(1:3) :: rot_rel_wing_w=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_b=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_rel_wing_g=0.0_rk
      ! angular acceleration vectors (left and right wings, 2nd left and 2nd right wings)
      real(kind=rk), dimension(1:3) :: rot_dt_wing_w=0.0_rk
      real(kind=rk), dimension(1:3) :: rot_dt_wing_g=0.0_rk

      ! rotation matrices for the various coordinate system for the insect
      real(kind=rk),dimension(1:3,1:3) :: M_b2w

   end type

   ! **********************************************************************************************
   ! **                              DEFINITION OF AN INSECT                                     **
   ! **********************************************************************************************
   type diptera
      logical :: initialized = .false.

      ! most variables are hidden here, because they are wing-specific
      ! (yes, I know that diptera have only two wings.)
      ! Wings: 1=left 2=right 3=left hind 4=right hind
      type(wing_type) :: Wings(1:4)

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
      real(kind=rk) :: psi=0.0_rk, beta=0.0_rk, gamma=0.0_rk, psi_dt=0.0_rk, beta_dt=0.0_rk, gamma_dt=0.0_rk
      ! stroke plane angle
      real(kind=rk) :: eta_stroke=0.0_rk
      ! is the body motion state described be the STATE vector? This is the case if the
      ! free-flight solver is used, and if its results are read in postprocessing or
      ! if it used used to prescribe the body motion state from a different simulation
      logical :: quaternion_solver_used = .false.

      ! angular velocity vectors (body, left and right wings, 2nd left and 2nd right wings)
      real(kind=rk), dimension(1:3) :: rot_body_b=0.0_rk, rot_body_g=0.0_rk

      ! vectors desribing the positions of insect's key elements
      ! in the body coordinate system
      real(kind=rk), dimension(1:3) :: x_head=0.0_rk, x_eye_r=0.0_rk, x_eye_l=0.0_rk

      ! moments of inertia in the body reference frame
      real(kind=rk) :: Jroll_body=0.0_rk, Jyaw_body=0.0_rk, Jpitch_body=0.0_rk
      ! total mass of insect:
      real(kind=rk) :: mass, gravity=0.0_rk, gravity_y=0.0_rk, gravity_x=0.0_rk
      ! variables to decide whether to draw the body or not.
      character(len=clong) :: body_moves="yes"
      character(len=clong) :: BodySuperSTLfile="none.superstl"
      ! second wing pair exists or not
      logical :: second_wing_pair
      

      !-------------------------------------------------------------
      ! for kinematics loader 
      !-------------------------------------------------------------
      ! variables for kineloader (which uses non-periodic hermite interpolation)
      ! Describes all four wings and the body 
      integer :: nk = 0      
      ! 1   0   time
      ! 2   1   body_center_g_x
      ! 3   2   body_center_g_x_dt
      ! 4   3   body_center_g_y
      ! 5   4   body_center_g_y_dt
      ! 6   5   body_center_g_z
      ! 7   6   body_center_g_z_dt
      ! 8   7   currently unused
      ! 9   8   currently unused
      ! 10  9   currently unused
      ! 11  10  currently unused
      ! 12  11  currently unused
      ! 13  12  currently unused
      ! 14  13  psi
      ! 15  14  psi_dt
      ! 16  15  beta
      ! 17  16  beta_dt
      ! 18  17  gamma
      ! 19  18  gamma_dt
      ! 20  19  alpha_L
      ! 21  20  alpha_L_dt
      ! 22  21  phi_L
      ! 23  22  phi_L_dt
      ! 24  23  theta_L
      ! 25  24  theta_L_dt
      ! 26  25  alpha_R
      ! 27  26  alpha_R_dt
      ! 28  27  phi_R
      ! 29  28  phi_R_dt
      ! 30  29  theta_R
      ! 31  30  theta_R_dt
      ! 32  31  alpha_L2
      ! 33  32  alpha_L2_dt
      ! 34  33  phi_L2
      ! 35  34  phi_L2_dt
      ! 36  35  theta_L2
      ! 37  36  theta_L2_dt
      ! 38  37  alpha_R2
      ! 39  38  alpha_R2_dt
      ! 40  39  phi_R2
      ! 41  40  phi_R2_dt
      ! 42  41  theta_R2
      ! 43  42  theta_R2_dt
      real(kind=rk), dimension (:,:), allocatable :: data_kineloader
      logical :: kineloader_initialized = .false.
      character(len=clong) :: infile_kineloader = ""

      !-------------------------------------------------------------
      ! for free flight solver
      !-------------------------------------------------------------
      real(kind=rk) :: time=0.0_rk
      real(kind=rk), allocatable :: RHS(:,:)
      ! this is the force and moment that is applied on the insect from the fluid, it will be computed and stored during RHS computations
      real(kind=rk), dimension(1:3) :: force_g=0.0_rk, moment_g=0.0_rk
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

      ! the following flag makes the code write the kinematics log to either kinematics.t
      ! (regular simulation) or kinematics.dry-run.t (for a dry run). The reason for this
      ! is that during postprocessing of an existing run, the dry run would overwrite the
      ! simulation data.
      character(len=clong) :: kinematics_file = "kinematics.t"
      ! rotation matrices for the various coordinate system for the insect
      real(kind=rk),dimension(1:3,1:3) :: M_g2b, M_b2g

      !-------------------------------------------------------------
      ! parameters that control shape of wings, body, and motion
      !-------------------------------------------------------------
      character(len=clong) :: BodyType="", BodyMotion=""
      ! parameters for body:
      real(kind=rk) :: L_body=0.0_rk, b_body=0.0_rk, R_head=0.0_rk, R_eye=0.0_rk
      ! parameters for wing shape:
      real(kind=rk) :: b_top=0.0_rk, b_bot=0.0_rk, L_span=0.0_rk
      ! WingThickness is also a part of the diptera type (even though it is a wing property), because it
      ! sets the default for all wings. This parameter is read from the main PARAMS.INI file, while each wing
      ! may read its own thickness from the wing INI file. the value read from main is the default when reading the wing.
      real(kind=rk) :: WingThickness=0.0_rk
      ! this is a safety distance for smoothing:
      real(kind=rk) :: safety=0.0_rk, L_smooth=0.0_rk, C_smooth=1.0_rk, dx_reference=0.0_rk, C_shell_thickness=5.0_rk
      ! some more VPM parameters will be stored in the insect, that can be individual
      character(len=clong) :: smoothing_type=""
      integer :: smoothing_type_int=0  ! in point-wise loops, we have to avoid character comparisons, so we use this integer
      real(kind=rk) :: epsilon_hester=0.0_rk
      ! Wings and body forces (1:body, 2:left wing, 3:right wing, 4:left wing, 5:right wing)
      type(Integrals), dimension(1:5) :: PartIntegrals

      !-------------------------------------------------------------
      ! parameters for mask coloring
      !-------------------------------------------------------------
      ! available color values
      ! color_body - color of main body of the insect
      ! color_l - color of left wing
      ! color_r - color of right wing
      ! color_l2 - color of second left wing, if it exists
      ! color_r2 - color of second right wing, if it exists
      ! color_geometry - color of the full insect, used for computing forces / moments wrt the full insect
      integer(kind=2) :: color_body=1, color_l=2, color_r=3, color_l2=4, color_r2=5, color_geometry=0

      !-----------------------------------------------------------------------------
      ! array for superSTL file for the body
      real(kind=rk), allocatable :: body_superSTL_b(:,:)
      real(kind=rk), allocatable :: body_superSTL_g(:,:)
   end type diptera
   !-----------------------------------------------------------------------------

   ! this module contains the insects objects. All other routines have to implicitly call update routines of this module to work with them
   ! It's a bit inconvenient that this is located below the type definition, but it is necessary
   integer(kind=ik), public :: n_insects = 0
   type(diptera), allocatable, public :: Insects(:)

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
#include "kineloader.f90"
!---------------------------------------



   !> Wrapper to update all insects at once
   subroutine Update_All_Insects( time )
      implicit none

      real(kind=rk), intent(in) :: time
      integer :: i

      do i=1,n_insects
         call Update_Insect( time, Insects(i) )
      enddo

   end subroutine Update_All_Insects

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
      type(diptera), intent(inout) :: Insect

      integer(kind=ik) :: i
      integer(kind=2) :: wingID

      !-----------------------------------------------------------------------------
      ! fetch current motion state
      !-----------------------------------------------------------------------------
      call BodyMotion (time, Insect)

      ! update wing angles (alpha, phi, theta and their time derivatives for each wing)
      do wingID = 1, 4
         if (Insect%Wings(wingID)%used) then
            call FlappingMotionWrap(time, Insect, wingID)
         endif
      enddo

      !-----------------------------------------------------------------------------
      ! define the rotation matrices to change between coordinate systems
      !-----------------------------------------------------------------------------
      call body_rotation_matrix( Insect, Insect%M_g2b )
      ! inverse of the body rotation matrix
      Insect%M_b2g = transpose(Insect%M_g2b)

      ! update wing rotation matrix for all wings that are used
      do wingID = 1, 4
         if (Insect%Wings(wingID)%used) then
            call wing_rotation_matrix( Insect%Wings(wingID)%M_b2w, Insect%Wings(wingID)%alpha, &
                                       Insect%Wings(wingID)%theta, Insect%Wings(wingID)%phi, &
                                       Insect%eta_stroke, Insect%Wings(wingID)%side )
         endif
      enddo

      !-----------------------------------------------------------------------------
      ! Angular velocities & acceleration
      !-----------------------------------------------------------------------------
      ! body angular velocity vector in b/g coordinate system
      call body_angular_velocity( Insect, Insect%rot_body_b, Insect%rot_body_g, Insect%M_g2b )
      ! rel+abs wing angular velocities in the w/b/g coordinate system
      call update_all_wing_angular_velocities ( time, Insect )
      ! angular acceleration for wings (required for inertial power)
      call wing_angular_accel( time, Insect )

      !-----------------------------------------------------------------------------
      ! vector from body centre to left/right pivot point in global reference frame,
      ! for aerodynamic power
      ! NOTE: this qty is used to compute aerody. moments in FLUSI, something odd happens and
      ! we need the relative vector. In wabbit, we need the absolute vector.
      !-----------------------------------------------------------------------------
      do wingID = 1, 4
         if (Insect%Wings(wingID)%used) then
            Insect%Wings(wingID)%x_pivot_g = matmul(Insect%M_b2g, Insect%Wings(wingID)%x_pivot_b) + Insect%xc_body_g
         endif
      enddo

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
         do i = 1, size(Insect%body_superSTL_b, 1)
            Insect%body_superSTL_g(i, 1:3)   = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 1:3)) + Insect%xc_body_g
            Insect%body_superSTL_g(i, 4:6)   = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 4:6)) + Insect%xc_body_g
            Insect%body_superSTL_g(i, 7:9)   = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 7:9)) + Insect%xc_body_g
            Insect%body_superSTL_g(i, 10:12) = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 10:12))
            Insect%body_superSTL_g(i, 13:15) = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 13:15))
            Insect%body_superSTL_g(i, 16:18) = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 16:18))
            Insect%body_superSTL_g(i, 19:21) = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 19:21))
            Insect%body_superSTL_g(i, 22:24) = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 22:24))
            Insect%body_superSTL_g(i, 25:27) = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 25:27))
            Insect%body_superSTL_g(i, 28:30) = matmul(Insect%M_b2g, Insect%body_superSTL_b(i, 28:30))
         enddo
      endif

   end subroutine Update_Insect

   !-------------------------------------------------------------------------------
   ! Init kinematics log of all insects
   !-------------------------------------------------------------------------------
   subroutine Init_insect_data(overwrite)
      implicit none

      logical, intent(in) :: overwrite

      integer :: i, i_check
      logical :: do_init
      character(len=clong) :: header(1:43*n_insects+1)

      header(1) = "time"

      do i=1,n_insects
         do_init = .true.
         if (insects(i)%kinematics_file /= Insects(1)%kinematics_file) then
            write(*,'("ERROR! insects have different kinematics files: ",a," and ",a, " . We will use the first one only")') &
            Insects(i)%kinematics_file, Insects(1)%kinematics_file
         endif

         write(header(2+(i-1)*43), '(A, i0.2, A)') "insect", i, ":xc_body_g_x"
         write(header(3+(i-1)*43), '(A, i0.2, A)') "insect", i, ":xc_body_g_y"
         write(header(4+(i-1)*43), '(A, i0.2, A)') "insect", i, ":xc_body_g_z"
         write(header(5+(i-1)*43), '(A, i0.2, A)') "insect", i, ":psi (rad)"
         write(header(6+(i-1)*43), '(A, i0.2, A)') "insect", i, ":beta (rad)"
         write(header(7+(i-1)*43), '(A, i0.2, A)') "insect", i, ":gamma (rad)"
         write(header(8+(i-1)*43), '(A, i0.2, A)') "insect", i, ":eta (rad)"
         write(header(9+(i-1)*43), '(A, i0.2, A)') "insect", i, ":alpha_l (rad)"
         write(header(10+(i-1)*43), '(A, i0.2, A)') "insect", i, ":phi_l (rad)"
         write(header(11+(i-1)*43), '(A, i0.2, A)') "insect", i, ":theta_l (rad)"
         write(header(12+(i-1)*43), '(A, i0.2, A)') "insect", i, ":alpha_r (rad)"
         write(header(13+(i-1)*43), '(A, i0.2, A)') "insect", i, ":phi_r (rad)"
         write(header(14+(i-1)*43), '(A, i0.2, A)') "insect", i, ":theta_r (rad)"
         write(header(15+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_l_w_x"
         write(header(16+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_l_w_y"
         write(header(17+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_l_w_z"
         write(header(18+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_r_w_x"
         write(header(19+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_r_w_y"
         write(header(20+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_r_w_z"
         write(header(21+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_l_w_x"
         write(header(22+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_l_w_y"
         write(header(23+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_l_w_z"
         write(header(24+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_r_w_x"
         write(header(25+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_r_w_y"
         write(header(26+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_r_w_z"
         write(header(27+(i-1)*43), '(A, i0.2, A)') "insect", i, ":alpha_l2 (rad)"
         write(header(28+(i-1)*43), '(A, i0.2, A)') "insect", i, ":phi_l2 (rad)"
         write(header(29+(i-1)*43), '(A, i0.2, A)') "insect", i, ":theta_l2 (rad)"
         write(header(30+(i-1)*43), '(A, i0.2, A)') "insect", i, ":alpha_r2 (rad)"
         write(header(31+(i-1)*43), '(A, i0.2, A)') "insect", i, ":phi_r2 (rad)"
         write(header(32+(i-1)*43), '(A, i0.2, A)') "insect", i, ":theta_r2 (rad)"
         write(header(33+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_l2_w_x"
         write(header(34+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_l2_w_y"
         write(header(35+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_l2_w_z"
         write(header(36+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_r2_w_x"
         write(header(37+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_r2_w_y"
         write(header(38+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_rel_r2_w_z"
         write(header(39+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_l2_w_x"
         write(header(40+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_l2_w_y"
         write(header(41+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_l2_w_z"
         write(header(42+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_r2_w_x"
         write(header(43+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_r2_w_y"
         write(header(44+(i-1)*43), '(A, i0.2, A)') "insect", i, ":rot_dt_r2_w_z"
      enddo
      call init_t_file( Insects(1)%kinematics_file, overwrite, header=header(1:43*n_insects+1))
   end subroutine Init_insect_data

   !-------------------------------------------------------------------------------
   ! Write kinematics log of all insects - we assume that all share the same kinematics file
   !-------------------------------------------------------------------------------
   subroutine Write_insect_data( time )
      implicit none
      real(kind=rk), intent(in) :: time
      real(kind=rk), allocatable, dimension(:) :: buffer_local
      integer :: i, n_entries

      n_entries = 43

      if (allocated(buffer_local)) then
         if (size(buffer_local) < n_entries*n_insects) deallocate(buffer_local)
      endif
      if (.not. allocated(buffer_local)) allocate(buffer_local(1:n_entries*n_insects))

      do i = 1, n_insects
         buffer_local((i-1)*n_entries+1:i*n_entries) = (/&
                  Insects(i)%xc_body_g, &
                  Insects(i)%psi, &
                  Insects(i)%beta, &
                  Insects(i)%gamma, &
                  Insects(i)%eta_stroke, &
                  Insects(i)%Wings(1)%alpha, & ! left
                  Insects(i)%Wings(1)%phi, &
                  Insects(i)%Wings(1)%theta, &
                  Insects(i)%Wings(2)%alpha, & ! right
                  Insects(i)%Wings(2)%phi, &
                  Insects(i)%Wings(2)%theta, &
                  Insects(i)%Wings(1)%rot_rel_wing_w, & ! left
                  Insects(i)%Wings(2)%rot_rel_wing_w, & !right
                  Insects(i)%Wings(1)%rot_dt_wing_w, & ! left
                  Insects(i)%Wings(2)%rot_dt_wing_w, & ! right
                  Insects(i)%Wings(3)%alpha, & ! left2
                  Insects(i)%Wings(3)%phi, &
                  Insects(i)%Wings(3)%theta, &
                  Insects(i)%Wings(4)%alpha, & ! right2
                  Insects(i)%Wings(4)%phi, &
                  Insects(i)%Wings(4)%theta, &
                  Insects(i)%Wings(3)%rot_rel_wing_w, & !L2
                  Insects(i)%Wings(4)%rot_rel_wing_w, & !R2
                  Insects(i)%Wings(3)%rot_dt_wing_w, & !L2
                  Insects(i)%Wings(4)%rot_dt_wing_w & !R2
            /)
      enddo
      call append_t_file( Insects(1)%kinematics_file, (/time, buffer_local(1:n_entries*n_insects) /) )
   end subroutine Write_insect_data


   !-------------------------------------------------------------------------------
   ! Main routine for drawing insects. Draws body and wings, parameters are in "INSECT"
   !-------------------------------------------------------------------------------
   subroutine Draw_Insect( time, Insect, xx0, ddx, mask, mask_color, us)
      implicit none

      real(kind=rk), intent(in)      :: time
      type(diptera), intent(inout)   :: Insect
      real(kind=rk), intent(in)      :: xx0(1:3), ddx(1:3)
      real(kind=rk), intent(inout)   :: mask(0:,0:,0:)
      real(kind=rk), intent(inout)   :: us(0:,0:,0:,1:)
      real(kind=rk), intent(inout)   :: mask_color(0:,0:,0:)

      if ((dabs(Insect%time-time)>1.0d-10).and.root) then
         write(*,'("error! time=",es15.8," but Insect%time=",es15.8)') time, Insect%time
         write(*,'("Did you call Update_Insect before Draw_Insect?")')
      endif


      call draw_insect_body( time, xx0, ddx, mask, mask_color, us, Insect)
      call draw_insect_wings( time, xx0, ddx, mask, mask_color, us, Insect)

      ! this is a debug test, which succeeded.
      !call check_if_us_is_derivative_of_position_wingtip(time, Insect)
   end subroutine Draw_Insect

   !-------------------------------------------------------------------------------
   ! check if this block contains part of the insect
   !-------------------------------------------------------------------------------
   subroutine insect_geometry_indicator(time, insect_id, BS, g, x0, dx, geometry_in_block)
      implicit none
      real(kind=rk), intent(in) :: time
      integer(kind=ik), intent(in) :: insect_id
      integer, intent(in) :: BS(1:3), g
      real(kind=rk), intent(in) :: x0(1:3), dx(1:3)  ! positions of block from WABBIT, so situated at g+1
      logical, intent(out) :: geometry_in_block

      real(kind=rk) :: xend(1:3), block_extent(1:3), char_length, vec(1:3)

      ! compute end of blocks
      block_extent = dx * real(BS, kind=rk)
      xend = x0 + block_extent

      ! check if body center is contained - it is located around xc_body_g
      if (all( (Insects(insect_id)%xc_body_g(:) >= x0(:)) .and. all(Insects(insect_id)%xc_body_g(:) <= xend(:)) )) then
         geometry_in_block = .true.
      endif

      ! Check if wing is contained - this is tricky, as we usually only have the pivot point (which often outside the wing) and complex shapes
      ! so for now, I skip this and hope that (fingers cross), resolving the body is enough to draw the wings as well

   end subroutine insect_geometry_indicator

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
      color_body = 1  ! NOT Insect%color_body
      color_l = 2     ! NOT Insect%color_l
      color_r = 3     ! NOT Insect%color_r

      ! body is not driven directly, therefore the power is set to zero
      Insect%PartIntegrals(color_body)%APow = 0.0_rk

      !-----------
      ! left wing
      !-----------
      ! relative angular velocity, in global system
      omrel = Insect%Wings(1)%rot_rel_wing_g

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
      omrel = Insect%Wings(2)%rot_rel_wing_g

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
                  Insect%PartIntegrals(color_l)%APow + &
                  Insect%PartIntegrals(color_r)%APow

      !-----------
      ! if second wing pair is present
      !-----------
      ! Colors
      color_l2 = 4     ! NOT Insect%color_l2
      color_r2 = 5     ! NOT Insect%color_r2

      ! left wing
      ! relative angular velocity, in global system
      omrel = Insect%Wings(3)%rot_rel_wing_g

      momrel = Insect%PartIntegrals(color_l2)%Torque + &
               Insect%PartIntegrals(color_l2)%Torque_unst

      ! aerodynamic power
      Insect%PartIntegrals(color_l2)%APow = - sum( momrel * omrel )

      ! right wing
      ! relative angular velocity, in global system
      omrel = Insect%Wings(4)%rot_rel_wing_g

      momrel = Insect%PartIntegrals(color_r2)%Torque + &
               Insect%PartIntegrals(color_r2)%Torque_unst

      ! aerodynamic power
      Insect%PartIntegrals(color_r2)%APow = - sum( momrel * omrel )

      ! Total aerodynamic power
      apowtotal = apowtotal + Insect%PartIntegrals(color_l2)%APow + &
                              Insect%PartIntegrals(color_r2)%APow

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
      integer(kind=ik) :: wingID, color

      iwmoment   = 0.0_rk
      iwmoment_g = 0.0_rk
      ipowtotal  = 0.0_rk

      ! colors for Diptera (one body, four wings) (yes, I know that diptera have only two wings.)
      color_body = 1  ! NOT Insect%color_body
      color_l = 2     ! NOT Insect%color_l
      color_r = 3     ! NOT Insect%color_r
      color_l2 = 4    ! NOT Insect%color_l2
      color_r2 = 5    ! NOT Insect%color_r2

      ! computation of inertial power and moment
      do wingID = 1, 4
         color = wingID+1

         a(1) = Insect%Wings(wingID)%Jxx * Insect%Wings(wingID)%rot_dt_wing_w(1) + Insect%Wings(wingID)%Jxy * Insect%Wings(wingID)%rot_dt_wing_w(2)
         a(2) = Insect%Wings(wingID)%Jxy * Insect%Wings(wingID)%rot_dt_wing_w(1) + Insect%Wings(wingID)%Jyy * Insect%Wings(wingID)%rot_dt_wing_w(2)
         a(3) = Insect%Wings(wingID)%Jzz * Insect%Wings(wingID)%rot_dt_wing_w(3)

         b(1) = Insect%Wings(wingID)%Jxx * Insect%Wings(wingID)%rot_rel_wing_w(1) + Insect%Wings(wingID)%Jxy * Insect%Wings(wingID)%rot_rel_wing_w(2)
         b(2) = Insect%Wings(wingID)%Jxy * Insect%Wings(wingID)%rot_rel_wing_w(1) + Insect%Wings(wingID)%Jyy * Insect%Wings(wingID)%rot_rel_wing_w(2)
         b(3) = Insect%Wings(wingID)%Jzz * Insect%Wings(wingID)%rot_rel_wing_w(3)

         ! inertial moment (in wing system)
         iwmoment(1, color) = (a(1)+Insect%Wings(wingID)%rot_rel_wing_w(2)*b(3)-Insect%Wings(wingID)%rot_rel_wing_w(3)*b(2))
         iwmoment(2, color) = (a(2)+Insect%Wings(wingID)%rot_rel_wing_w(3)*b(1)-Insect%Wings(wingID)%rot_rel_wing_w(1)*b(3))
         iwmoment(3, color) = (a(3)+Insect%Wings(wingID)%rot_rel_wing_w(1)*b(2)-Insect%Wings(wingID)%rot_rel_wing_w(2)*b(1))

         ! transform inertial moment into the laboratory reference frame
         iwmoment_g(:, color) = -matmul(Insect%M_b2g, matmul(transpose(Insect%Wings(wingID)%M_b2w), iwmoment(:,color)))

         ! inertial power
         Insect%PartIntegrals(color)%IPow = Insect%Wings(wingID)%rot_rel_wing_w(1) * iwmoment(1, color) + &
                                            Insect%Wings(wingID)%rot_rel_wing_w(2) * iwmoment(2, color) + &
                                            Insect%Wings(wingID)%rot_rel_wing_w(3) * iwmoment(3, color)

         ipowtotal = ipowtotal + Insect%PartIntegrals(color)%IPow
      enddo


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


   subroutine compute_wing_angular_velocity( Insect, phi, alpha, theta, phi_dt, alpha_dt, theta_dt, eta_stroke, side, &
    rot_rel_wing_w, rot_rel_wing_b, rot_rel_wing_g)
      implicit none

      type(diptera), intent(inout) :: Insect
      real(kind=rk), intent(in) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt, eta_stroke
      character(len=1), intent(in) :: side
      real(kind=rk), dimension(1:3), intent(out) :: rot_rel_wing_w, rot_rel_wing_b, rot_rel_wing_g

      real(kind=rk), dimension(1:3,1:3) :: M_b2s, M_b2w, M_g2b
      real(kind=rk), dimension(1:3) :: rot_wing_s

      call stroke_rotation_matrix(M_b2s, eta_stroke, side)
      call wing_rotation_matrix(M_b2w, alpha, theta, phi, eta_stroke, side)
      call body_rotation_matrix(Insect, M_g2b)

      if (side == "L") then
         ! direct definition in stroke reference frame (left)
         rot_wing_s = (/ phi_dt-sin(theta)*alpha_dt, &
                         cos(phi)*cos(theta)*alpha_dt-sin(phi)*theta_dt, &
                         sin(phi)*cos(theta)*alpha_dt+cos(phi)*theta_dt /)

      elseif (side == 'R') then
         ! direct definition in stroke reference frame (right)
         rot_wing_s = (/ -phi_dt-sin(theta)*(-alpha_dt), &
                         cos(-phi)*cos(theta)*(-alpha_dt)-sin(-phi)*theta_dt, &
                         sin(-phi)*cos(theta)*(-alpha_dt)+cos(-phi)*theta_dt /)
      else
         ! left or right (left2 and right2 not accepted - just pass the different angles)
         write(*,*) "1", side
         call abort(10710242, "insect_module: neither right nor left side ? how many sides does an insect have? seven!?")

      endif

      ! bring it to body system
      rot_rel_wing_b = matmul( transpose(M_b2s), rot_wing_s )
      ! then the result to wing system
      rot_rel_wing_w = matmul( M_b2w, rot_rel_wing_b)
      ! and finally the global one
      rot_rel_wing_g = matmul( transpose(M_g2b), rot_rel_wing_b )

   end subroutine

   !-------------------------------------------------------------------------------
   ! given the angles of each wing (and their time derivatives), compute
   ! the angular velocity vectors for all wings.
   ! output:
   !    Insect%rot_rel_wing_*_w    relative angular velocity of all wings (wing frame)
   !    Insect%rot_rel_wing_*_b    relative angular velocity of all wings (body frame)
   !    Insect%rot_rel_wing_*_g    relative angular velocity of all wings (glob frame)
   !-------------------------------------------------------------------------------
   subroutine update_all_wing_angular_velocities ( time, Insect )
      implicit none

      real(kind=rk), intent(in) :: time
      type(diptera), intent(inout) :: Insect
      integer(kind=ik) :: wingID


      do wingID = 1, 4
         if (Insect%Wings(wingID)%used) then
            call compute_wing_angular_velocity( Insect, Insect%Wings(wingID)%phi, Insect%Wings(wingID)%alpha, &
                  Insect%Wings(wingID)%theta, Insect%Wings(wingID)%phi_dt, Insect%Wings(wingID)%alpha_dt, &
                  Insect%Wings(wingID)%theta_dt, Insect%eta_stroke, Insect%Wings(wingID)%side, &
                  Insect%Wings(wingID)%rot_rel_wing_w, Insect%Wings(wingID)%rot_rel_wing_b, &
                  Insect%Wings(wingID)%rot_rel_wing_g)
         endif
      enddo

   end subroutine update_all_wing_angular_velocities



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
      ! ATTENTION: This is a deep copy in Fortran, it copies all allocatables, which can be expensive
      ! Consider a new object that only copies parts that are necessary for the wing angular acceleration
      Insect2 = Insect

      ! fetch motion state at time+dt
      call BodyMotion( time+dt, Insect2 )
      call body_rotation_matrix( Insect, M_g2b ) ! (time)

      do wingID = 1, 4
         Insect%Wings(wingID)%rot_dt_wing_w = 0.0_rk

         if (Insect%Wings(wingID)%used) then
            ! get wing angles (time+dt)
            call FlappingMotionWrap(time+dt, Insect2, wingID)
            
            ! and the angular velocity (time+dt)
            call compute_wing_angular_velocity( Insect2, Insect2%Wings(wingID)%phi, Insect2%Wings(wingID)%alpha, &
               Insect2%Wings(wingID)%theta, Insect2%Wings(wingID)%phi_dt, Insect2%Wings(wingID)%alpha_dt, &
               Insect2%Wings(wingID)%theta_dt, Insect2%eta_stroke, Insect2%Wings(wingID)%side, &
               Insect2%Wings(wingID)%rot_rel_wing_w, Insect2%Wings(wingID)%rot_rel_wing_b, &
               Insect2%Wings(wingID)%rot_rel_wing_g)

            ! angular accel using one-sided finite difference
            ! use one-sided finite differences to derive the absolute angular velocity with
            ! respect to time. Note in older code versions, this was wrong, as we derived
            ! the ang. vel. in the wing coordinate system, which is a moving reference frame.
            rot_dt_wing_g = (Insect2%Wings(wingID)%rot_rel_wing_g - Insect%Wings(wingID)%rot_rel_wing_g) / dt

            ! angular accel in wing frame
            Insect%Wings(wingID)%rot_dt_wing_w = matmul( Insect%Wings(wingID)%M_b2w, matmul( M_g2b, rot_dt_wing_g))
         endif
      enddo

      ! just to be sure, we deallocate insect2 here
      call clean_insect(insect2)
   end subroutine wing_angular_accel


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
   ! returns the stroke rotation matrix M_b2s for a given side
   !-----------------------------------------------------------------------------
   subroutine stroke_rotation_matrix( M_b2s, eta_stroke, side )
      implicit none

      real(kind=rk),intent(out) :: M_b2s(1:3,1:3)
      real(kind=rk), intent(in) :: eta_stroke
      ! as this is called rarely okay to use string comparison
      character(len=1), intent(in) :: side ! left or right (left2 and right2 not accepted - just pass the different angles)
      real(kind=rk), dimension(1:3,1:3) :: M1, M2


      if (side == "R") then
         call Rx(M1, pi)
         call Ry(M2, eta_stroke)
         M_b2s = matmul(M1,M2)

      elseif (side == "L") then
         call Ry(M1, eta_stroke)
         M_b2s = M1

      else
         write(*,*) "2", side
         ! left or right (left2 and right2 not accepted - just pass the different angles)
         call abort(1071024, "insect_module: neither right nor left side ? how many sides does an insect have? seven!?")

      endif
   end subroutine stroke_rotation_matrix


   !-----------------------------------------------------------------------------
   ! returns the wing rotation matrix M_b2w for a given side
   !-----------------------------------------------------------------------------
   subroutine wing_rotation_matrix( M_b2w, alpha, theta, phi, eta_stroke, side )
      implicit none

      real(kind=rk),intent(out) :: M_b2w(1:3,1:3)
      real(kind=rk), intent(in) :: alpha, theta, phi, eta_stroke
      ! as this is called rarely okay to use string comparison
      character(len=1), intent(in) :: side ! L or R (left2 and right2 not accepted - just pass the different angles)
      real(kind=rk), dimension(1:3,1:3) :: M1, M2, M3, M_b2s

      call stroke_rotation_matrix(M_b2s, eta_stroke, side)

      if (side == "R") then
         ! note the coordinate system is rotated so we don't need to inverse the sign
         ! of theta, and the wings still rotate in opposite direction
         call Ry(M1, -alpha)
         call Rz(M2,  theta)
         call Rx(M3, -phi)
         M_b2w = matmul(M1,matmul(M2,matmul(M3,M_b2s)))

      elseif (side == "L") then
         call Ry(M1, alpha)
         call Rz(M2, theta)
         call Rx(M3, phi)
         M_b2w = matmul(M1,matmul(M2,matmul(M3,M_b2s)))

      else
         write(*,*) "3", side
         ! left or right (left2 and right2 not accepted - just pass the different angles)
         call abort(1071024, "insect_module: neither right nor left side ? how many sides does an insect have? seven!?")

      endif
   end subroutine wing_rotation_matrix


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
   ! subroutine check_if_us_is_derivative_of_position_wingtip(time, Insect)
   !    real(kind=rk), intent(in) :: time
   !    type(diptera), intent(inout) :: Insect

   !    real(kind=rk) :: M_g2b(1:3,1:3), M_b2w_r(1:3,1:3), x_tip_w(1:3), x_tip_b(1:3), x_tip_g(1:3), &
   !       us_tip_g(1:3), v_tmp(1:3), v_tmp_b(1:3)
   !    real(kind=rk)::xd,yd,zd
   !    real(kind=rk)::c00,c10,c01,c11,c0,c1
   !    integer :: ix,iy,iz

   !    call body_rotation_matrix( Insect, M_g2b )
   !    ! call wing_right_rotation_matrix( Insect, M_b2w_r )
   !    call wing_rotation_matrix( M_b2w_r, Insect%alpha_r, Insect%theta_r, Insect%phi_r, Insect%eta_stroke, "right" )
   !    ! body angular velocity vector in b/g coordinate system
   !    call body_angular_velocity( Insect, Insect%rot_body_b, Insect%rot_body_g, M_g2b )
   !    ! rel+abs wing angular velocities in the w/b/g coordinate system
   !    call update_all_wing_angular_velocities ( time, Insect )

   !    x_tip_w = (/ 1.0_rk, 1.0_rk, 1.0_rk /)
   !    x_tip_b = matmul( transpose(M_b2w_r), x_tip_w ) + Insect%Wings(2)%x_pivot_b
   !    x_tip_g = matmul( transpose(M_g2b), x_tip_b ) + Insect%xc_body_g

   !    !-----------------------------------------------------------------------------
   !    ! now we extrcat how the us field is constructed
   !    v_tmp(1) = Insect%rot_rel_wing_r_w(2)*x_tip_w(3)-Insect%rot_rel_wing_r_w(3)*x_tip_w(2)
   !    v_tmp(2) = Insect%rot_rel_wing_r_w(3)*x_tip_w(1)-Insect%rot_rel_wing_r_w(1)*x_tip_w(3)
   !    v_tmp(3) = Insect%rot_rel_wing_r_w(1)*x_tip_w(2)-Insect%rot_rel_wing_r_w(2)*x_tip_w(1)
   !    v_tmp_b = matmul(transpose(M_b2w_r), v_tmp) ! in body system


   !    ! translational part. we compute the rotational part in the body
   !    ! reference frame, therefore, we must transform the body translation
   !    ! velocity Insect%vc (which is in global coordinates) to the body frame
   !    v_tmp = matmul(M_g2b,Insect%vc_body_g)

   !    ! add solid body rotation to the translational velocity field. Note
   !    ! that rot_body_b and x_body are in the body reference frame
   !    v_tmp(1) = v_tmp(1) + Insect%rot_body_b(2)*x_tip_b(3)-Insect%rot_body_b(3)*x_tip_b(2)
   !    v_tmp(2) = v_tmp(2) + Insect%rot_body_b(3)*x_tip_b(1)-Insect%rot_body_b(1)*x_tip_b(3)
   !    v_tmp(3) = v_tmp(3) + Insect%rot_body_b(1)*x_tip_b(2)-Insect%rot_body_b(2)*x_tip_b(1)

   !    ! the body motion is added to the wing motion, which is already in us
   !    ! and they are also in the body refrence frame. However, us has to be
   !    ! in the global reference frame, so M_b2g is applied
   !    us_tip_g = matmul( transpose(M_g2b), v_tmp_b + v_tmp )

   !    open(14,file='debug_wing_us.t',status='unknown',position='append')
   !    write (14,'(7(es15.8,1x))') time, x_tip_g, us_tip_g
   !    close(14)

   ! end subroutine

end module module_insects
