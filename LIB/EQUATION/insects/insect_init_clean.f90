!> We want to work with several insects, at this part we initialize the array of insects.
!! They need to be initialized individually by initialize_insect and generally are referenced by their id
subroutine insects_array_init(n_insects_in)
    implicit none
    integer, intent(in) :: n_insects_in

    n_insects = n_insects_in
    allocate(Insects(n_insects))

end subroutine insects_array_init


subroutine initialize_insect(time, fname_ini, Insect, SECTION, Insect_ID, resume_backup, fname_backup, box_domain, &
    viscosity, c_eta_in, dx_reference, smoothing_type, N_ghost_nodes, periodic, colors_default)
    implicit none
    real(kind=rk), intent(in) :: time
    character(len=*), intent(in) :: fname_ini
    logical, intent(in) :: resume_backup
    character(len=*), intent(in) :: fname_backup
    ! why passing these parameters and not read them from the params file?? The answer is that we use this module
    ! in different codes, hence we must be sure that properties like domain size and viscosity are found in the same
    ! sections. this is not the case, so we give them to the insect module in the call here.
    real(kind=rk), intent(in) :: box_domain(1:3), viscosity, c_eta_in
    ! as the default wing thickness is 4*dx, pass lattice spacing here. In FLUSI, this is easy
    ! but in WABBIT it requires some thought, because dx is not a constant.
    real(kind=rk), intent(in) :: dx_reference
    ! smoothing type for the insect mask, can be cos or hester or discontinous
    character(len=*), intent(in) :: smoothing_type
    ! section to be read in file from. [Insects] or [Insect1] for the 1st one, then [Insect2], etc
    character(len=*), intent(in) :: SECTION
    ! ID of the insect (still required, even though this routine is called on a single instance of type(diptera))
    integer(kind=ik), intent(in) :: Insect_ID
    ! the actual insect to be initialized
    type(diptera), intent(inout) :: Insect
    ! ghost nodes. If the insect module is used in a finite-differences code, then
    ! the data that we have often has ghost nodes, i.e. points that overlap and exist
    ! on several CPUS. On those, you normally would not create the mask (which is expensive)
    ! so we skip the first and last "g" points on the arrays used for mask creation
    integer, optional, intent(in) :: N_ghost_nodes
    !
    logical, optional, intent(in) :: periodic
    !> default colors for the insect geometries, in case they are not set in the ini file. This is useful to avoid color conflicts with the user defined geometries, which are colored from 0 to n_geometries-1.
    !! The default is to set the insect colors to start right after n_geometries and be unique for all parts (body, left wing, right wing, left wing 2, right wing 2, geometry / full insect) and individual insects.
    integer, optional, intent(in) :: colors_default(6)

    type(inifile) :: PARAMS
    real(kind=rk) :: defaultvec(1:3), defaultvec5(5)
    character(len=clong) :: DoF_string, dummystr
    integer :: j, tmp, mpirank, mpicode, ntri
    integer(kind=2) :: wingID
    logical :: section_exists

    ! in this module, we use the logical ROOT to avoid the integer comparison mpirank==0
    call MPI_COMM_RANK (MPI_COMM_WORLD, mpirank, mpicode)
    if (mpirank==0) root = .true.

    ! copy parameters from the call to global variables:
    xl = box_domain(1)
    yl = box_domain(2)
    zl = box_domain(3)
    nu = viscosity

    ! read in the complete ini file, from which we initialize the insect
    call read_ini_file_mpi(PARAMS, fname_ini, verbose=.true.)

    ! header information
    if (root) then
        write(*,'(A)') "---------------------------------------------------------------------------------------"
        write(*,'(A)') "      .==-.                   .-==."
        write(*,'(A)') "       \()8`-._  `.   .'  _.-'8()/"
        write(*,'(A)') "       (88'   ::.  \./  .::   '88)"
        write(*,'(A)') "        \_.'`-::::.(#).::::-'`._/"
        write(*,'(A)') "          `._... .q(_)p. ..._.'        Initializing"
        write(*,'(A)') "            ''-..-'|=|`-..-''    "
        write(*,'(A)') "            .''' .'|=|`. `''.   Insect"
        write(*,'(A)') "          ,':8(o)./|=|\.(o)8:`.  "
        write(*,'(A)') "         (O :8 ::/ \_/ \:: 8: O)       Module!"
        write(*,'(A)') "          \O `::/       \::' O/  "
        write(*,'(A)') "           ''--'         `--''   "
        write(*,'(A)') "---------------------------------------------------------------------------------------"
        write(*,'(A)') "Initializing insect module!"
        write(*,'(A)') "*.ini file is: "//trim(adjustl(fname_ini))// " and insect is in section "//trim(adjustl(SECTION))
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

    ! some VPM parameters can be insect specific
    ! We want to avoid reading them from the ini file, so they are passed along here
    ! Could be changed to be insect specific later (atleast the smoothing type, for c_eta and therefore epsilon this is more difficult)
    Insect%epsilon_hester = sqrt(nu * c_eta_in)
    Insect%smoothing_type = smoothing_type
    select case(smoothing_type)
        case("cos", "cosine")
            Insect%smoothing_type_int = STEP_METHOD_COSINE
        case ("hester")
            Insect%smoothing_type_int = STEP_METHOD_HESTER
        case("discontinuous", "dis")
            Insect%smoothing_type_int = STEP_METHOD_DISC
        case default
            call abort(260602, "Insect hit a tree - Never heard of the smoothing type "//trim(smoothing_type))
    end select

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! WINGS TO BE SIMULATED
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Which wings are to be simulated ?
    call read_param_mpi(PARAMS, SECTION, "LeftWing", Insect%Wings(1)%used, .true.)
    call read_param_mpi(PARAMS, SECTION, "RightWing", Insect%Wings(2)%used, .true.)
    call read_param_mpi(PARAMS, SECTION, "LeftWing2", Insect%Wings(3)%used, .false.)
    call read_param_mpi(PARAMS, SECTION, "RightWing2", Insect%Wings(4)%used, .false.)
    
    ! determine whether the second pair of wings is present
    ! This variable is used to control output files.
    Insect%second_wing_pair = .false.
    if ( (Insect%Wings(3)%used) .or. (Insect%Wings(4)%used) ) then
        Insect%second_wing_pair = .true.
    endif

    ! store the sides for each wing - done only once for the rest of the code
    ! Wings: 1=left 2=right 3=left hind 4=right hind
    Insect%Wings(1)%side = "L"
    Insect%Wings(2)%side = "R"
    Insect%Wings(3)%side = "L"
    Insect%Wings(4)%side = "R"

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Wing shapes
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Read wing shapes (forewings)
    call read_param_mpi(PARAMS, SECTION, "WingShapeL", Insect%Wings(1)%WingShape, 'NOT-USED')
    call read_param_mpi(PARAMS, SECTION, "WingShapeR", Insect%Wings(2)%WingShape, 'NOT-USED')

    ! Read wing shapes for second wing pair (usually hindwings or elytra)
    call read_param_mpi(PARAMS, SECTION, "WingShape2L", Insect%Wings(3)%WingShape, "NOT-USED")
    call read_param_mpi(PARAMS, SECTION, "WingShape2R", Insect%Wings(4)%WingShape, "NOT-USED")

    ! Rectangular wing parameters
    call read_param_mpi(PARAMS, SECTION, "b_top", Insect%b_top, 0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "b_bot", Insect%b_bot, 0.0_rk)

    ! Abort if deprecated parameter WingShape (was used for L/R symmetric animals) is used
    call read_param_mpi(PARAMS, SECTION, "WingShape", dummystr, "UNKNOWN")
    if (dummystr /= "UNKNOWN") then
        call abort(2307261,"Insect-module: the INI-parameter __WingShape__ is deprecated and was removed. Please set&
        & individual shapes in the Insects section of your INI file (WingShapeR, WingShapeL, WingShapeR2, WingShapeL2).&
        & Even if all wings are the same, we require a shape for each one.")
    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! WING KINEMATICS
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Wingbeat kinematics for each wing
    !-- wingID: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
    call read_param_mpi(PARAMS, SECTION, "FlappingMotion_left"  , Insect%Wings(1)%FlappingMotion, "NOT-USED")
    call read_param_mpi(PARAMS, SECTION, "FlappingMotion_right" , Insect%Wings(2)%FlappingMotion, "NOT-USED")
    call read_param_mpi(PARAMS, SECTION, "FlappingMotion_left2" , Insect%Wings(3)%FlappingMotion, "NOT-USED")
    call read_param_mpi(PARAMS, SECTION, "FlappingMotion_right2", Insect%Wings(4)%FlappingMotion, "NOT-USED")

    ! this one file contains all the kinematics (four wings and body)
    ! hence, not a wing specific qty and stored in the insect itself
    call read_param_mpi(PARAMS, SECTION, "infile_kineloader", Insect%infile_kineloader, "none")

    ! if the wing motion is read from file (the most often used case), then we extract the file to be used (from the 
    ! string "from_file::FILENAME.INI" we extract "FILENAME.INI")
    do wingID = 1, 4
        if ( index(Insect%Wings(wingID)%FlappingMotion, "from_file::") /= 0 ) then
            ! Store the file to be used in the wingkinematics 
            Insect%Wings(wingID)%infile = Insect%Wings(wingID)%FlappingMotion( 12:clong  )
            ! Set type to from_file - then the kinematics are read and initialized from (and not hardcoded coefficients)
            Insect%Wings(wingID)%FlappingMotion = "from_file"      

        elseif ( Insect%Wings(wingID)%FlappingMotion == "from_file" ) then
            ! old syntax, implies symmetric periodic motion, read from *.ini file
            call abort(2307263,"Insect-module: the INI-parameter __infile__ is deprecated and was removed. It was used&
            & for symmetric wing motion, in combination with FlappingMotion_right=from_file.&
            & This possibility has been removed, and you need to set FlappingMotion_right=from_file::XXX.ini for all&
            & wings in the simulation, even if they have the same flapping motion. Note difference from_file vs from_file::")
        endif
    enddo

    if (root) then
        !-- wingID: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
        write(*,*) "-------------- Wingbeat kinematics summary --------------"
        write(*,*) "Left wing motion:      "//trim(adjustl(Insect%Wings(1)%FlappingMotion))
        write(*,*) "Left wing infile:      "//trim(adjustl(Insect%Wings(1)%infile))
        write(*,*) "Right wing motion:     "//trim(adjustl(Insect%Wings(2)%FlappingMotion))
        write(*,*) "Right wing infile:     "//trim(adjustl(Insect%Wings(2)%infile))
        write(*,*) ""
        write(*,*) "2nd Left wing motion:  "//trim(adjustl(Insect%Wings(3)%FlappingMotion))
        write(*,*) "2nd Left wing infile:  "//trim(adjustl(Insect%Wings(3)%infile))
        write(*,*) "2nd Right wing motion: "//trim(adjustl(Insect%Wings(4)%FlappingMotion))
        write(*,*) "2nd Right wing infile: "//trim(adjustl(Insect%Wings(4)%infile))
        write(*,*) "---------------------------------------------------------"
    endif
    
    ! Wing thickness
    call read_param_mpi(PARAMS, SECTION, "WingThickness", Insect%WingThickness, 4.0_rk*dx_reference)
    ! The default for the wing thickness is from the main parameter file - that may be overwritten
    ! if the wing INI file specifies another value
    do wingID = 1, 4
        call read_param_mpi(PARAMS, SECTION, "WingThickness", Insect%Wings(wingID)%WingThickness, 4.0_rk*dx_reference)
    enddo


    call read_param_mpi(PARAMS, SECTION, "BodyType",Insect%BodyType,"ellipsoid")
    call read_param_mpi(PARAMS, SECTION, "BodyMotion",Insect%BodyMotion,"tethered")

    call read_param_mpi(PARAMS, SECTION, "mass",Insect%mass, 1._rk)
    call read_param_mpi(PARAMS, SECTION, "gravity",Insect%gravity, 0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "gravity_x",Insect%gravity_x, 0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "gravity_y",Insect%gravity_y, 0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "J_body_yawpitchroll",defaultvec, (/0.0_rk,0.0_rk,0.0_rk/))
    Insect%Jroll_body  = defaultvec(3)
    Insect%Jyaw_body   = defaultvec(1)
    Insect%Jpitch_body = defaultvec(2)
    call read_param_mpi(PARAMS, SECTION, "x0",Insect%x0, (/0.5_rk*xl,0.5_rk*yl,0.5_rk*zl/))
    call read_param_mpi(PARAMS, SECTION, "v0",Insect%v0, (/0.0_rk, 0.0_rk, 0.0_rk/))
    call read_param_mpi(PARAMS, SECTION, "yawpitchroll_0",Insect%yawpitchroll_0,(/0.0_rk, 0.0_rk, 0.0_rk/))
    call read_param_mpi(PARAMS, SECTION, "yawpitchroll_a1",Insect%yawpitchroll_a1,(/0.0_rk, 0.0_rk, 0.0_rk/))
    call read_param_mpi(PARAMS, SECTION, "yawpitchroll_b1",Insect%yawpitchroll_b1,(/0.0_rk, 0.0_rk, 0.0_rk/))
    ! convert yawpitchroll to radiants
    Insect%yawpitchroll_0 = Insect%yawpitchroll_0 * (pi/180.0_rk)
    Insect%yawpitchroll_a1 = Insect%yawpitchroll_a1 * (pi/180.0_rk)
    Insect%yawpitchroll_b1 = Insect%yawpitchroll_b1 * (pi/180.0_rk)

    ! stroke plane angle. it is constant and the same for all four wings.
    call read_param_mpi(PARAMS, SECTION, "eta0",Insect%eta_stroke, 0.0_rk)
    Insect%eta_stroke = Insect%eta_stroke*(pi/180.0_rk)

    ! degrees of freedom for free flight solver. The string from ini file contains
    ! 6 characters 1 or 0 that turn on/off x,y,z,yaw,pitch,roll degrees of freedom
    ! by multiplying the respective RHS by zero, keeping the value thus constant
    call read_param_mpi(PARAMS, SECTION, "DoF",DoF_string, "111111")
    do j=1,6
        read (DoF_string(j:j), '(i1)') tmp
        Insect%DoF_on_off(j) = dble(tmp)
    enddo
    if (root) write(*,'(6(f4.2,1x))') Insect%DoF_on_off

    ! wing inertia tensor (we currently assume two identical forewings and two identical hindwings)
    ! this allows computing inertial power and wing FSI model
    call read_param_mpi(PARAMS, SECTION, "Jxx",Insect%Wings(1)%Jxx,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jyy",Insect%Wings(1)%Jyy,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jzz",Insect%Wings(1)%Jzz,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jxy",Insect%Wings(1)%Jxy,0.0_rk)

    call read_param_mpi(PARAMS, SECTION, "Jxx",Insect%Wings(2)%Jxx,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jyy",Insect%Wings(2)%Jyy,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jzz",Insect%Wings(2)%Jzz,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jxy",Insect%Wings(2)%Jxy,0.0_rk)

    call read_param_mpi(PARAMS, SECTION, "Jxx2",Insect%Wings(3)%Jxx,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jyy2",Insect%Wings(3)%Jyy,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jzz2",Insect%Wings(3)%Jzz,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jxy2",Insect%Wings(3)%Jxy,0.0_rk)

    call read_param_mpi(PARAMS, SECTION, "Jxx2",Insect%Wings(4)%Jxx,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jyy2",Insect%Wings(4)%Jyy,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jzz2",Insect%Wings(4)%Jzz,0.0_rk)
    call read_param_mpi(PARAMS, SECTION, "Jxy2",Insect%Wings(4)%Jxy,0.0_rk)

    call read_param_mpi(PARAMS, SECTION, "startup_conditioner",Insect%startup_conditioner,"no")

    call read_param_mpi(PARAMS, SECTION, "C_smooth",Insect%C_smooth,1.0_rk)
    call read_param_mpi(PARAMS, SECTION, "BodySuperSTLfile",Insect%BodySuperSTLfile,"none.superstl")
    ! when using CT data, code computes the mask function in a shell around fluid-solid interface.
    ! The tickness of the shell is not a critical parameter, but it affects performance. Thicker shell
    ! means more points and thus more comput effort. It is given in multiples of C_smooth, that means
    ! shell_thickness = C_shell_thickness * C_smooth * dx_min. dx_min is the spacing on the finest level Jmax.
    ! Why is the shell thickness dependent on resolution? The cost to generate the mask depends on the number of
    ! triangles and the number of points in the shell. As the latter is coupled to dx and the former is constant
    ! this way the mask generation cost is constant when increasing the resolution.
    call read_param_mpi(PARAMS, SECTION, "C_shell_thickness",Insect%C_shell_thickness, 3.0_rk)

    Insect%dx_reference = dx_reference
    Insect%L_smooth = Insect%C_smooth*dx_reference
    if (Insect%smoothing_type == "hester") then
        Insect%L_smooth = Insect%epsilon_hester
        Insect%safety = max(15.0_rk*Insect%epsilon_hester, 2*dx_reference)
    else
        Insect%safety = 3.5_rk*Insect%L_smooth
    end if

    ! wing hinges (origin of rotation relative to the body)
    defaultvec=(/0.0_rk, +Insect%b_body, 0.0_rk /)
    call read_param_mpi(PARAMS, SECTION, "x_pivot_l",Insect%Wings(1)%x_pivot_b, defaultvec)

    defaultvec=(/0.0_rk, -Insect%b_body, 0.0_rk /)
    call read_param_mpi(PARAMS, SECTION, "x_pivot_r",Insect%Wings(2)%x_pivot_b, defaultvec)

    defaultvec = Insect%Wings(1)%x_pivot_b
    call read_param_mpi(PARAMS, SECTION, "x_pivot_l2",Insect%Wings(3)%x_pivot_b, defaultvec)

    defaultvec = Insect%Wings(2)%x_pivot_b
    call read_param_mpi(PARAMS, SECTION, "x_pivot_r2",Insect%Wings(4)%x_pivot_b, defaultvec)

    !---------------------------------------------------------------------------
    ! Setup of colors (used to tell different parts of the animal apart)
    !---------------------------------------------------------------------------
    ! default colors for body, left wing, right wing, left wing 2, right wing 2, geometry / full insect
    defaultvec5 = (/1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk/)
    if (present(colors_default)) defaultvec5 = real(colors_default(1:5), kind=rk)

    call read_param_mpi(PARAMS, SECTION, "colors", defaultvec5, defaultvec5)
    Insect%color_body = int(defaultvec5(1), kind=2)
    Insect%color_l    = int(defaultvec5(2), kind=2)
    Insect%color_r    = int(defaultvec5(3), kind=2)
    Insect%color_l2   = int(defaultvec5(4), kind=2)
    Insect%color_r2   = int(defaultvec5(5), kind=2)
    ! the color_geometry is used only for statistics of forces/moments wrt to the full insect
    ! physics module should provide this as input
    if (present(colors_default)) then
        Insect%color_geometry = colors_default(6)
    else
        Insect%color_geometry = 0_2
    endif

    Insect%Wings(1)%color = Insect%color_l 
    Insect%Wings(2)%color = Insect%color_r 
    Insect%Wings(3)%color = Insect%color_l2
    Insect%Wings(4)%color = Insect%color_r2

    ! clean ini file
    call clean_ini_file_mpi(PARAMS)

    !---------------------------------------------------------------------------
    ! Wing shapes setup 
    !---------------------------------------------------------------------------
    ! Processing of wing shapes (Setup from file or hardcoded shapes without a setup)
    ! wingID: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
    do wingID = 1, 4
        ! do not try to setup wings that are not used
        if (Insect%Wings(wingID)%used) then
            ! exclude wings that do not need a setup, otherwise, call setup
            if (Insect%Wings(wingID)%WingShape/="suzuki" .and. &
                Insect%Wings(wingID)%WingShape/="rectangular" .and. &
                Insect%Wings(wingID)%WingShape/="TwoEllipses" .and. &
                Insect%Wings(wingID)%WingShape/="suzuki_butterfly") then

                ! we have some pre-defined, hard-coded data, but also can read the wing shape
                ! from INI files. This is done in the routine below.
                call Setup_WingShape(Insect, wingID)
            endif
        endif
    enddo

    !---------------------------------------------------------------------------
    ! initialization for superSTl body
    !---------------------------------------------------------------------------
    if (Insect%BodyType == "superSTL") then
        if (.not. allocated(Insect%body_superSTL_b)) then
            if (root) write(*,'("INSECTS: STL: init start")')
            if (root) write(*,'("INSECTS: STL: file=",A)') Insect%BodySuperSTLfile

            call count_lines_in_ascii_file_mpi(Insect%BodySuperSTLfile, ntri, 0)

            if (root) write(*,'("INSECTS: STL: file length is ntri=", i7 )') ntri

            allocate( Insect%body_superSTL_b(1:ntri,1:30) )
            allocate( Insect%body_superSTL_g(1:ntri,1:30) )

            ! No scaling or origin shift is applied: we assume you did that when generating
            ! the superSTL file. The data is thus understood in the body coordinate system.
            call read_array_from_ascii_file_mpi(Insect%BodySuperSTLfile, Insect%body_superSTL_b, 0)

            if (root) write(*,'("INSECTS: STL: read from file...done! We are good to go.")')
        endif
    endif


    !-----------------------------------------------------------------------------
    ! other initialization
    !-----------------------------------------------------------------------------

    ! If required, initialize rigid solid dynamics solver
    if (Insect%BodyMotion=="free_flight") then
        ! note we have to do that before init_fields as rigid_solid_init sets up
        ! the state vector without which create_mask cannot know the position and velocity
        ! of body and wings
        call rigid_solid_init( time, Insect, resume_backup, Insect_ID )
    endif

    ! the update routine computes wing angles and so on, everything that is done only
    ! once per time step. Do this here as well, so we can safely call draw_insect after
    ! calling this routine.
    call Update_Insect( time, Insect )

    if (root) then
        write(*,'(80("<"))')
        write(*,*) "Insect initialization is complete."
        write(*,'(80("<"))')
    endif
    Insect%initialized = .true.

end subroutine initialize_insect


! -------------------------------------------------------------------------------
! \brief Clean up insect data structure for all insects or one specific one
subroutine clean_all_insects(insect_i)
    implicit none

    integer(kind=ik), intent(in), optional :: insect_i  ! which insect to clean?
    integer(kind=ik) :: Insects_i, insect_start, insect_end

    ! if we clean up all insects, then this mass extinction will be announced to the user
    if (root .and. .not. present(insect_i)) then
        write(*,'(80("<"))')
        write(*,*) "Finalizing insect module!"
        write(*,'(80("<"))')
    endif

    if (present(insect_i)) then
        insect_start = insect_i
        insect_end   = insect_i
    else
        insect_start = 1
        insect_end   = n_insects
    endif

    do Insects_i = insect_start, insect_end
        call clean_insect(Insects(Insects_i))
    enddo

    if (allocated(Insects)) deallocate(Insects)
   
end subroutine clean_all_insects

! -------------------------------------------------------------------------------
! \brief Clean up insect data structure of one insect
subroutine clean_insect(insect)
    implicit none

    type(diptera), intent(inout) :: insect
    integer(kind=ik) :: wingID

    do wingID = 1, 4
        if (allocated(Insect%Wings(wingID)%ai_phi)) deallocate(Insect%Wings(wingID)%ai_phi)
        if (allocated(Insect%Wings(wingID)%bi_phi)) deallocate(Insect%Wings(wingID)%bi_phi)
        if (allocated(Insect%Wings(wingID)%ai_theta)) deallocate(Insect%Wings(wingID)%ai_theta)
        if (allocated(Insect%Wings(wingID)%bi_theta)) deallocate(Insect%Wings(wingID)%bi_theta)
        if (allocated(Insect%Wings(wingID)%ai_alpha)) deallocate(Insect%Wings(wingID)%ai_alpha)
        if (allocated(Insect%Wings(wingID)%bi_alpha)) deallocate(Insect%Wings(wingID)%bi_alpha)

        if (allocated(Insect%Wings(wingID)%ai_shape)) deallocate(Insect%Wings(wingID)%ai_shape)
        if (allocated(Insect%Wings(wingID)%bi_shape)) deallocate(Insect%Wings(wingID)%bi_shape)
        if (allocated(Insect%Wings(wingID)%R0_table)) deallocate(Insect%Wings(wingID)%R0_table)
        if (allocated(Insect%Wings(wingID)%polygon_wings)) deallocate(Insect%Wings(wingID)%polygon_wings)
        if (allocated(Insect%Wings(wingID)%damage_mask)) deallocate(Insect%Wings(wingID)%damage_mask)
        if (allocated(Insect%Wings(wingID)%wing_thickness_profile)) deallocate(Insect%Wings(wingID)%wing_thickness_profile)
        if (allocated(Insect%Wings(wingID)%corrugation_profile)) deallocate(Insect%Wings(wingID)%corrugation_profile)
        if (allocated(Insect%Wings(wingID)%bristles_coords)) deallocate(Insect%Wings(wingID)%bristles_coords)
        ! if (allocated(Insect%Wings(wingID)%deformations)) deallocate(Insect%Wings(wingID)%deformations)
        ! if (allocated(Insect%Wings(wingID)%deformation_profile)) deallocate(Insect%Wings(wingID)%deformation_profile)
        ! if (allocated(Insect%Wings(wingID)%mask_wing_complete)) deallocate(Insect%Wings(wingID)%mask_wing_complete)
        ! if (allocated(Insect%Wings(wingID)%particle_points)) deallocate(Insect%Wings(wingID)%particle_points)
    enddo

    if (allocated(insect%data_kineloader)) deallocate(insect%data_kineloader)
    if (allocated(insect%RHS)) deallocate(insect%RHS)
    if (allocated(insect%body_superSTL_b)) deallocate ( insect%body_superSTL_b )
    if (allocated(insect%body_superSTL_g)) deallocate ( insect%body_superSTL_g )

end subroutine clean_insect

