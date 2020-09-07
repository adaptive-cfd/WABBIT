subroutine insect_init(time, fname_ini, Insect, resume_backup, fname_backup, box_domain, &
    viscosity, dx_reference, N_ghost_nodes, periodic)
  implicit none
  real(kind=rk), intent(in) :: time
  character(len=*), intent(in) :: fname_ini
  type(diptera),intent(inout)::Insect
  logical, intent(in) :: resume_backup
  character(len=*), intent(in) :: fname_backup
  ! why passing these parameters and not read them from the params file?? The answer is that we use this module
  ! in different codes, hence we must be sure that properties like domain size and viscosity are found in the same
  ! sections. this is not the case, so we give them to the insect module in the call here.
  real(kind=rk), intent(in) :: box_domain(1:3), viscosity
  ! as the default wing thickness is 4*dx, pass lattice spacing here. In FLUSI, this is easy
  ! but in WABBIT it requires some thought, because dx is not a constant.
  real(kind=rk), intent(in) :: dx_reference
  ! ghost nodes. If the insect module is used in a finite-differences code, then
  ! the data that we have often has ghost nodes, i.e. points that overlap and exist
  ! on several CPUS. On those, you normally would not create the mask (which is expensive)
  ! so we skip the first and last "g" points on the arrays used for mask creation
  integer, optional, intent(in) :: N_ghost_nodes
  !
  logical, optional, intent(in) :: periodic

  type(inifile) :: PARAMS
  real(kind=rk),dimension(1:3)::defaultvec
  character(len=strlen) :: DoF_string, dummystr
  integer :: j, tmp, mpirank, mpicode, ntri
  integer(kind=2) :: wingID, Nwings

  ! in this module, we use the logical ROOT to avoid the integer comparison mpirank==0
  call MPI_COMM_RANK (MPI_COMM_WORLD, mpirank, mpicode)
  if (mpirank==0) root = .true.

  ! copy parameters from the call:
  xl = box_domain(1)
  yl = box_domain(2)
  zl = box_domain(3)
  nu = viscosity
  ! header information
  if (root) then
      write(*,*) "---------------------------------------------------------------------------------------"
        write(*,*) "      .==-.                   .-==."
        write(*,*) "       \()8`-._  `.   .'  _.-'8()/"
        write(*,*) "       (88'   ::.  \./  .::   '88)"
        write(*,*) "        \_.'`-::::.(#).::::-'`._/"
        write(*,*) "          `._... .q(_)p. ..._.'        Initializing"
        write(*,*) "            ''-..-'|=|`-..-''    "
        write(*,*) "            .''' .'|=|`. `''.   Insect"
        write(*,*) "          ,':8(o)./|=|\.(o)8:`.  "
        write(*,*) "         (O :8 ::/ \_/ \:: 8: O)       Module!"
        write(*,*) "          \O `::/       \::' O/  "
        write(*,*) "           ''--'         `--''   "
      write(*,*) "---------------------------------------------------------------------------------------"
      write(*,*) "Initializing insect module!"
      write(*,*) "*.ini file is: "//trim(adjustl(fname_ini))
      write(*,'(80("<"))')
      write(*,'("Lx=",g12.4," Ly=",g12.4," Lz=",g12.4," nu=",g12.4)') xl, yl, zl, nu
      write(*,'("dx=",g12.4," nx_equidistant=",i6)') dx_reference, nint(xl/dx_reference)
  endif

  ! ghost nodes are optional..because in FLUSI, we do not have them
  if (present(N_ghost_nodes)) then
      ! g is a module global private variable.
      g = N_ghost_nodes
  else
      g = 0
  endif
  if (root) write(*,'("n_ghosts=",i2)') g

  ! is the insect periodic?
  ! attention: this functionality is not for free, so if you do not need it - disable it.
  if (present(periodic)) then
      periodic_insect = periodic
  else
      periodic_insect = .false.
  endif
  if (root) write(*,'("periodic_insect=",L1)') periodic_insect

  !-----------------------------------------------------------------------------
  ! read in parameters form ini file
  !-----------------------------------------------------------------------------

  ! read in the complete ini file, from which we initialize the insect
  call read_ini_file_mpi(PARAMS, fname_ini, verbose=.true.)

  ! determine whether the second pair of wings is present
  call read_param_mpi(PARAMS,"Insects","LeftWing2",Insect%LeftWing2,"no")
  call read_param_mpi(PARAMS,"Insects","RightWing2",Insect%RightWing2,"no")
  if ( ( Insect%LeftWing2 == "yes" ) .or. ( Insect%RightWing2 == "yes" ) ) then
    Insect%second_wing_pair = .true.
  else
    Insect%second_wing_pair = .false.
  endif

  ! read data for the first pair of wings
  call read_param_mpi(PARAMS,"Insects","WingShape",Insect%WingShape(1),"none")
  Insect%WingShape(2) = Insect%WingShape(1)
  ! The following two lines take effect if WingShapeL or WingShapeR are set
  call read_param_mpi(PARAMS,"Insects","WingShapeL",Insect%WingShape(1),Insect%WingShape(1))
  call read_param_mpi(PARAMS,"Insects","WingShapeR",Insect%WingShape(2),Insect%WingShape(2))
  ! Rectangular wing parameters
  call read_param_mpi(PARAMS,"Insects","b_top",Insect%b_top, 0.d0)
  call read_param_mpi(PARAMS,"Insects","b_bot",Insect%b_bot, 0.d0)
  call read_param_mpi(PARAMS,"Insects","L_span",Insect%L_span, 0.d0)
  ! Kinematics
  call read_param_mpi(PARAMS,"Insects","FlappingMotion_right",Insect%FlappingMotion_right,"none")
  call read_param_mpi(PARAMS,"Insects","FlappingMotion_left",Insect%FlappingMotion_left,"none")
  ! this file is used in old syntax form for both wings:
  call read_param_mpi(PARAMS,"Insects","infile",Insect%infile,"none.in")

  if ( index(Insect%FlappingMotion_right,"from_file::") /= 0 ) then
    ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
    Insect%kine_wing_r%infile = Insect%FlappingMotion_right( 12:strlen  )
    Insect%FlappingMotion_right = "from_file"

  elseif ( index(Insect%FlappingMotion_right,"kinematics_loader::") /= 0 ) then
    ! new syntax, uses the kinematics loader for non-periodic kinematics
    Insect%kine_wing_r%infile = Insect%FlappingMotion_right( 20:strlen )
    Insect%FlappingMotion_right = "kinematics_loader"

  elseif ( Insect%FlappingMotion_right == "from_file" ) then
    ! old syntax, implies symmetric periodic motion, read from *.ini file
    Insect%kine_wing_r%infile = Insect%infile

  elseif ( Insect%FlappingMotion_right == "kinematics_loader" ) then
    ! old syntax, implies symmetric non-periodic motion, read from *.dat file
    Insect%kine_wing_r%infile = Insect%infile
  endif


  if ( index(Insect%FlappingMotion_left,"from_file::") /= 0 ) then
    ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
    Insect%kine_wing_l%infile = Insect%FlappingMotion_left( 12:strlen  )
    Insect%FlappingMotion_left = "from_file"

  elseif ( index(Insect%FlappingMotion_left,"kinematics_loader::") /= 0 ) then
    ! new syntax, uses the kinematics loader for non-periodic kinematics
    Insect%kine_wing_l%infile = Insect%FlappingMotion_left( 20:strlen )
    Insect%FlappingMotion_left = "kinematics_loader"

  elseif ( Insect%FlappingMotion_left == "from_file" ) then
    ! old syntax, implies symmetric periodic motion, read from *.ini file
    Insect%kine_wing_l%infile = Insect%infile

  elseif ( Insect%FlappingMotion_left == "kinematics_loader" ) then
    ! old syntax, implies symmetric non-periodic motion, read from *.dat file
    Insect%kine_wing_l%infile = Insect%infile
  endif

  if (root) then
    write(*,*) "Left wing: "//trim(adjustl(Insect%FlappingMotion_left))
    write(*,*) "Left wing: "//trim(adjustl(Insect%kine_wing_l%infile))
    write(*,*) "Right wing: "//trim(adjustl(Insect%FlappingMotion_right))
    write(*,*) "Right wing: "//trim(adjustl(Insect%kine_wing_r%infile))
  endif

  ! these flags trigger reading the kinematics from file when the Flapping
  ! motion is first called
  Insect%kine_wing_l%initialized = .false.
  Insect%kine_wing_r%initialized = .false.

  ! read data for the second pair of wings
  if (Insect%second_wing_pair) then
    ! note that only Fourier wing shape can be different from first wings
    call read_param_mpi(PARAMS,"Insects","WingShape2",Insect%WingShape(3),"none")
    Insect%WingShape(4) = Insect%WingShape(3)
    ! The following two lines take effect if WingShape2L or WingShape2R are set
    call read_param_mpi(PARAMS,"Insects","WingShape2L",Insect%WingShape(3),Insect%WingShape(3))
    call read_param_mpi(PARAMS,"Insects","WingShape2R",Insect%WingShape(4),Insect%WingShape(4))
    ! Kinematics
    call read_param_mpi(PARAMS,"Insects","FlappingMotion_right2",Insect%FlappingMotion_right2,"none")
    call read_param_mpi(PARAMS,"Insects","FlappingMotion_left2",Insect%FlappingMotion_left2,"none")
    ! this file is used in old syntax form for both wings:
    call read_param_mpi(PARAMS,"Insects","infile2",Insect%infile2,"none.in")

    if ( index(Insect%FlappingMotion_right2,"from_file::") /= 0 ) then
      ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
      Insect%kine_wing_r2%infile = Insect%FlappingMotion_right2( 12:strlen  )
      Insect%FlappingMotion_right2 = "from_file"

    elseif ( index(Insect%FlappingMotion_right,"kinematics_loader::") /= 0 ) then
      ! new syntax, uses the kinematics loader for non-periodic kinematics
      Insect%kine_wing_r2%infile = Insect%FlappingMotion_right2( 20:strlen )
      Insect%FlappingMotion_right2 = "kinematics_loader"

    elseif ( Insect%FlappingMotion_right2 == "from_file" ) then
      ! old syntax, implies symmetric periodic motion, read from *.ini file
      Insect%kine_wing_r2%infile = Insect%infile2

    elseif ( Insect%FlappingMotion_right2 == "kinematics_loader" ) then
      ! old syntax, implies symmetric non-periodic motion, read from *.dat file
      Insect%kine_wing_r2%infile = Insect%infile2
    endif

    if ( index(Insect%FlappingMotion_left2,"from_file::") /= 0 ) then
      ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
      Insect%kine_wing_l2%infile = Insect%FlappingMotion_left2( 12:strlen  )
      Insect%FlappingMotion_left2 = "from_file"

    elseif ( index(Insect%FlappingMotion_left2,"kinematics_loader::") /= 0 ) then
      ! new syntax, uses the kinematics loader for non-periodic kinematics
      Insect%kine_wing_l2%infile = Insect%FlappingMotion_left2( 20:strlen )
      Insect%FlappingMotion_left2 = "kinematics_loader"

    elseif ( Insect%FlappingMotion_left2 == "from_file" ) then
      ! old syntax, implies symmetric periodic motion, read from *.ini file
      Insect%kine_wing_l2%infile = Insect%infile2

    elseif ( Insect%FlappingMotion_left2 == "kinematics_loader" ) then
      ! old syntax, implies symmetric non-periodic motion, read from *.dat file
      Insect%kine_wing_l2%infile = Insect%infile2
    endif

    if (root) then
      write(*,*) "Second left wing: "//trim(adjustl(Insect%FlappingMotion_left2))
      write(*,*) "Second left wing: "//trim(adjustl(Insect%kine_wing_l2%infile))
      write(*,*) "Second right wing: "//trim(adjustl(Insect%FlappingMotion_right2))
      write(*,*) "Second right wing: "//trim(adjustl(Insect%kine_wing_r2%infile))
    endif
  endif

  ! these flags trigger reading the kinematics from file when the Flapping
  ! motion is first called
  Insect%kine_wing_l2%initialized = .false.
  Insect%kine_wing_r2%initialized = .false.

  call read_param_mpi(PARAMS,"Insects","BodyType",Insect%BodyType,"ellipsoid")
  call read_param_mpi(PARAMS,"Insects","BodySuperSTLfile",Insect%BodySuperSTLfile,"none.superstl")
  call read_param_mpi(PARAMS,"Insects","HasDetails",Insect%HasDetails,"all")
  call read_param_mpi(PARAMS,"Insects","BodyMotion",Insect%BodyMotion,"tethered")
  call read_param_mpi(PARAMS,"Insects","LeftWing",Insect%LeftWing,"yes")
  call read_param_mpi(PARAMS,"Insects","RightWing",Insect%RightWing,"yes")
  call read_param_mpi(PARAMS,"Insects","b_body",Insect%b_body, 0.1d0)
  call read_param_mpi(PARAMS,"Insects","L_body",Insect%L_body, 1.d0)
  call read_param_mpi(PARAMS,"Insects","R_head",Insect%R_head, 0.1d0)
  call read_param_mpi(PARAMS,"Insects","R_eye",Insect%R_eye, 0.d1)
  call read_param_mpi(PARAMS,"Insects","mass",Insect%mass, 1.d0)
  call read_param_mpi(PARAMS,"Insects","gravity",Insect%gravity, 1.d0)
  call read_param_mpi(PARAMS,"Insects","WingThickness",Insect%WingThickness, 4.0d0*dx_reference)
  call read_param_mpi(PARAMS,"Insects","J_body_yawpitchroll",defaultvec, (/0.d0,0.d0,0.d0/))
  Insect%Jroll_body  = defaultvec(3)
  Insect%Jyaw_body   = defaultvec(1)
  Insect%Jpitch_body = defaultvec(2)
  call read_param_mpi(PARAMS,"Insects","x0",Insect%x0, (/0.5d0*xl,0.5d0*yl,0.5d0*zl/))
  call read_param_mpi(PARAMS,"Insects","v0",Insect%v0, (/0.d0, 0.d0, 0.d0/))
  call read_param_mpi(PARAMS,"Insects","yawpitchroll_0",Insect%yawpitchroll_0,&
  (/0.d0, 0.d0, 0.d0/))
  ! convert yawpitchroll to radiants
  Insect%yawpitchroll_0 = Insect%yawpitchroll_0 * (pi/180.d0)
  call read_param_mpi(PARAMS,"Insects","eta0",Insect%eta0, 0.0d0)
  Insect%eta0 = Insect%eta0*(pi/180.d0)

  call read_param_mpi(PARAMS,"Insects","pointcloudfile",Insect%pointcloudfile,"none")



  ! degrees of freedom for free flight solver. The string from ini file contains
  ! 6 characters 1 or 0 that turn on/off x,y,z,yaw,pitch,roll degrees of freedom
  ! by multiplying the respective RHS by zero, keeping the value thus constant
  call read_param_mpi(PARAMS,"Insects","DoF",DoF_string, "111111")
  do j=1,6
    read (DoF_string(j:j), '(i1)') tmp
    Insect%DoF_on_off(j) = dble(tmp)
  enddo
  if (root) write(*,'(6(f4.2,1x))') Insect%DoF_on_off


  ! wing FSI section
  call read_param_mpi(PARAMS,"Insects","wing_fsi",Insect%wing_fsi,"no")
    call read_param_mpi(PARAMS, "Insects", "init_alpha_phi_theta", Insect%init_alpha_phi_theta, (/0.d0, 0.d0, 0.d0/) )


    ! section for additional fractal tree
    call read_param_mpi(PARAMS, "Insects", "fractal_tree", Insect%fractal_tree, .false.)
    call read_param_mpi(PARAMS, "Insects", "fractal_tree_file", Insect%fractal_tree_file, "tree_data.in")
    call read_param_mpi(PARAMS, "Insects", "fractal_tree_x0", Insect%fractal_tree_x0, (/0.0d0, 0.0d0, 0.0d0/) )
    call read_param_mpi(PARAMS, "Insects", "fractal_tree_scaling", Insect%fractal_tree_scaling, 1.0_rk )


  ! wing inertia tensor (we currently assume two identical wings)
  ! this allows computing inertial power and wing FSI model
  call read_param_mpi(PARAMS,"Insects","Jxx",Insect%Jxx,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jyy",Insect%Jyy,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jzz",Insect%Jzz,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jxy",Insect%Jxy,0.d0)

  call read_param_mpi(PARAMS,"Insects","startup_conditioner",Insect%startup_conditioner,"no")

  ! 28/01/2019: Thomas. Discovered that this was done block based, i.e. the smoothing layer
  ! had different thickness, if some blocks happened to be at different levels (and still carry
  ! a part of the smoothing layer.) I don't know if that made sense, because the layer shrinks/expands then
  ! and because it might be discontinous. Both options are included now, default is "as before"
  ! Insect%smoothing_thickness=="local"  : smoothing_layer = c_sm * 2**-J * L/(BS-1)
  ! Insect%smoothing_thickness=="global" : smoothing_layer = c_sm * 2**-Jmax * L/(BS-1)
  ! NOTE: for FLUSI, this has no impact! Here, the grid is constant and equidistant.
  ! NOTE: 05/2020 Thomas, I changed the default back to local.
  call read_param_mpi(PARAMS,"Insects","smoothing_thickness",Insect%smoothing_thickness,"local")
  Insect%smooth = 1.0d0*dx_reference
  Insect%safety = 3.5d0*Insect%smooth

  ! position vector of the head
  call read_param_mpi(PARAMS,"Insects","x_head",&
  Insect%x_head, (/0.5d0*Insect%L_body,0.d0,0.d0 /) )

  ! eyes
  defaultvec = Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.0d0,+1.0d0,1.0d0/)
  call read_param_mpi(PARAMS,"Insects","x_eye_r",Insect%x_eye_r, defaultvec)

  defaultvec = Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.0d0,-1.0d0,1.0d0/)
  call read_param_mpi(PARAMS,"Insects","x_eye_l",Insect%x_eye_l, defaultvec)

  ! wing hinges (root points)
  defaultvec=(/0.d0, +Insect%b_body, 0.d0 /)
  call read_param_mpi(PARAMS,"Insects","x_pivot_l",Insect%x_pivot_l_b, defaultvec)

  defaultvec=(/0.d0, -Insect%b_body, 0.d0 /)
  call read_param_mpi(PARAMS,"Insects","x_pivot_r",Insect%x_pivot_r_b, defaultvec)

  ! read data for the second pair of wing hinges
  if (Insect%second_wing_pair) then
    defaultvec = Insect%x_pivot_l_b
    call read_param_mpi(PARAMS,"Insects","x_pivot_l2",Insect%x_pivot_l2_b, defaultvec)
    defaultvec = Insect%x_pivot_r_b
    call read_param_mpi(PARAMS,"Insects","x_pivot_r2",Insect%x_pivot_r2_b, defaultvec)
  endif

  ! default colors for body and wings
  Insect%color_body=1
  Insect%color_l=2
  Insect%color_r=3
  if (Insect%second_wing_pair) then
    Insect%color_l2=4
    Insect%color_r2=5
  endif

  ! clean ini file
  call clean_ini_file_mpi(PARAMS)

    !---------------------------------------------------------------------------
    ! initialization for superSTl body
    !---------------------------------------------------------------------------
    if (Insect%BodyType == "superSTL") then
        if (.not. allocated(xyz_nxnynz)) then
            if (root) write(*,'("INSECTS: STL: init start")')
            if (root) write(*,'("INSECTS: STL: file=",A)') Insect%BodySuperSTLfile

            call count_lines_in_ascii_file_mpi(Insect%BodySuperSTLfile, ntri, 0)

            if (root) write(*,'("INSECTS: STL: file length is ntri=", i7 )') ntri

            allocate( xyz_nxnynz(1:ntri,1:30) )

            ! No scaling or origin shift is applied: we assume you did that when generating
            ! the superSTL file. The data is thus understood in the body coordinate system.
            call read_array_from_ascii_file_mpi(Insect%BodySuperSTLfile, xyz_nxnynz, 0)

            if (root) write(*,'("INSECTS: STL: read from file...done! We are good to go.")')
        endif
    endif


  !-----------------------------------------------------------------------------
  ! other initialization
  !-----------------------------------------------------------------------------

    if (Insect%fractal_tree) then
        ! we can also simulate an insect together with a fractal tree as turbulence
        ! generators.
        call fractal_tree_init(Insect)
    endif

  ! If required, initialize rigid solid dynamics solver
  if (Insect%BodyMotion=="free_flight") then
    ! note we have to do that before init_fields as rigid_solid_init sets up
    ! the state vector without which create_mask cannot know the position and velocity
    ! of body and wings
    call rigid_solid_init( time, Insect, resume_backup, fname_backup)
  endif

  ! the update routine computes wing angles and so on, everything that is done only
  ! once per time step. Do this here as well, so we can safely call draw_insect after
  ! calling this routine.
  call Update_Insect( time, Insect )


  ! At this point, we must also initialize all wing / body data: in wabbit, it may
  ! be that the initial grid contains less blocks than we wave MPIRANKS. Therefore,
  ! not all CPUS might call the MPI_BCAST for this initialization. (This is a BUGFIX
  ! 21 Oct 2019, Yokohama, Thomas)
  Nwings = 2
  if (Insect%second_wing_pair) Nwings = 4

  ! wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
  do wingID = 1, Nwings
      ! exclude wings that are hard-coded, otherwise, call initialization routine
      if (Insect%WingShape(wingID)/="pointcloud" .and. Insect%WingShape(wingID)/="mosquito_iams" .and. &
          Insect%WingShape(wingID)/="suzuki" .and. Insect%WingShape(wingID)/="rectangular" .and. &
            Insect%WingShape(wingID)/="TwoEllipses" .and. (Insect%wing_file_type(wingID)) /= "kleemeier") then

          ! we have some pre-defined, hard-coded data, but also can read the wing shape
          ! from INI files.
          call Setup_Wing_Fourier_coefficients(Insect, wingID)
      endif

  enddo


  if (root) then
    write(*,'(80("<"))')
    write(*,*) "Insect initialization is complete."
    write(*,'(80("<"))')
  endif
    Insect%initialized = .true.

end subroutine insect_init




subroutine insect_clean(Insect)
  implicit none
  type(diptera),intent(inout)::Insect

  if (root) then
    write(*,'(80("<"))')
    write(*,*) "Finalizing insect module!"
    write(*,'(80("<"))')
  endif

  if (allocated(particle_points)) deallocate ( particle_points )
  if (allocated(wing_thickness_profile)) deallocate ( wing_thickness_profile )
  if (allocated(corrugation_profile)) deallocate ( corrugation_profile )
  if (allocated(mask_wing_complete)) deallocate(mask_wing_complete)

  call load_kine_clean( Insect%kine_wing_l )
  call load_kine_clean( Insect%kine_wing_r )
end subroutine insect_clean
