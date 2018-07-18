subroutine insect_init(time, fname_ini, Insect, resume_backup, fname_backup, box_domain, viscosity)
  implicit none
  real(kind=rk), intent(in) :: time
  character(len=*), intent(in) :: fname_ini
  type(diptera),intent(inout)::Insect
  logical, intent(in) :: resume_backup
  character(len=*), intent(in) :: fname_backup
  real(kind=rk), intent(in) :: box_domain(1:3), viscosity

  type(inifile) :: PARAMS
  real(kind=rk),dimension(1:3)::defaultvec
  character(len=strlen) :: DoF_string, dummystr
  integer :: j, tmp, mpirank, mpicode

  if (root) then
    write(*,'(80("<"))')
    write(*,*) "Initializing insect module!"
    write(*,*) "*.ini file is: "//trim(adjustl(fname_ini))
    write(*,'(80("<"))')
  endif

  ! copy parameters from the call:
  xl = box_domain(1)
  yl = box_domain(2)
  zl = box_domain(3)
  nu = viscosity

  if (root) then
      write(*,'("Lx=",g12.4," Ly=",g12.4," Lz=",g12.4," nu=",g12.4)') xl,yl,zl,nu
  endif

  call MPI_COMM_RANK (MPI_COMM_WORLD,mpirank,mpicode)
  if (mpirank==0) root = .true.

  !-----------------------------------------------------------------------------
  ! read in parameters form ini file
  !-----------------------------------------------------------------------------

  ! read in the complete ini file, from which we initialize the insect
  call read_ini_file_mpi(PARAMS, fname_ini, verbose=.true.)

  call read_param_mpi(PARAMS,"Insects","WingShape",Insect%WingShape,"none")
  call read_param_mpi(PARAMS,"Insects","b_top",Insect%b_top, 0.d0)
  call read_param_mpi(PARAMS,"Insects","b_bot",Insect%b_bot, 0.d0)
  call read_param_mpi(PARAMS,"Insects","L_span",Insect%L_span, 0.d0)
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

  call read_param_mpi(PARAMS,"Insects","BodyType",Insect%BodyType,"ellipsoid")
  call read_param_mpi(PARAMS,"Insects","HasDetails",Insect%HasDetails,"all")
  call read_param_mpi(PARAMS,"Insects","BodyMotion",Insect%BodyMotion,"yes")
  call read_param_mpi(PARAMS,"Insects","LeftWing",Insect%LeftWing,"yes")
  call read_param_mpi(PARAMS,"Insects","RightWing",Insect%RightWing,"yes")
  call read_param_mpi(PARAMS,"Insects","b_body",Insect%b_body, 0.1d0)
  call read_param_mpi(PARAMS,"Insects","L_body",Insect%L_body, 1.d0)
  call read_param_mpi(PARAMS,"Insects","R_head",Insect%R_head, 0.1d0)
  call read_param_mpi(PARAMS,"Insects","R_eye",Insect%R_eye, 0.d1)
  call read_param_mpi(PARAMS,"Insects","mass",Insect%mass, 1.d0)
  call read_param_mpi(PARAMS,"Insects","gravity",Insect%gravity, 1.d0)
  call read_param_mpi(PARAMS,"Insects","WingThickness",Insect%WingThickness, 0.05d0)
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
  call read_param_mpi(PARAMS,"Insects","init_alpha_phi_theta",&
  Insect%init_alpha_phi_theta, (/0.d0, 0.d0, 0.d0 /) )


  ! wing inertia tensor (we currently assume two identical wings)
  ! this allows computing inertial power and wing FSI model
  call read_param_mpi(PARAMS,"Insects","Jxx",Insect%Jxx,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jyy",Insect%Jyy,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jzz",Insect%Jzz,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jxy",Insect%Jxy,0.d0)

  call read_param_mpi(PARAMS,"Insects","startup_conditioner",Insect%startup_conditioner,"no")

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

  ! default colors for body and wings
  Insect%color_body=1
  Insect%color_l=2
  Insect%color_r=3

  ! clean ini file
  call clean_ini_file_mpi(PARAMS)

  !-----------------------------------------------------------------------------
  ! other initialization
  !-----------------------------------------------------------------------------

  ! If required, initialize rigid solid dynamics solver
  if (Insect%BodyMotion=="free_flight") then
    ! note we have to do that before init_fields as rigid_solid_init sets up
    ! the state vector without which create_mask cannot know the position and velocity
    ! of body and wings
    call rigid_solid_init( time, Insect, resume_backup, fname_backup)
  endif

  if (root) then
    write(*,'(80("<"))')
    write(*,*) "Insect initialization is complete."
    write(*,'(80("<"))')
  endif

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
