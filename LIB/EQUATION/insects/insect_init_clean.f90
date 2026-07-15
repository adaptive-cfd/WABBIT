!> We want to work with several insects, at this part we initialize the array of insects.
!! They need to be initialized individually by insect_init and generally are referenced by their id
subroutine insects_array_init(n_insects_in)
    implicit none
    integer, intent(in) :: n_insects_in

    n_insects = n_insects_in
    allocate(Insects(n_insects))

end subroutine insects_array_init


subroutine insect_init(time, fname_ini, Insect_ID, resume_backup, fname_backup, box_domain, &
    viscosity, c_eta_in, dx_reference, smoothing_type, N_ghost_nodes, periodic, colors_default)
    implicit none
    real(kind=rk), intent(in) :: time
    character(len=*), intent(in) :: fname_ini
    integer, intent(in) :: Insect_ID
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
    character(len=clong) :: DoF_string, dummystr, insect_name_str
    integer :: j, tmp, mpirank, mpicode, ntri
    integer(kind=2) :: wingID, Nwings
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

    ! we differentiate different insects, legacy one is [Insects], other ones are "Insect[ID]"
    if (Insect_ID > 1) then
        write(insect_name_str, '(A,I0)') "Insect", Insect_ID
        ! check if section exists
        call param_section_exists_mpi(PARAMS, insect_name_str, section_exists)
        if (.not. section_exists) then
            call abort(260602, "Insect flew away, no section found for " // trim(insect_name_str))
        endif
    else
        ! If insect_ID==1, we either read from [Insects], which is the old format, or from the new format [Insect1]
        ! if [Insects] is found - we take that. That ensures compatibility with old INI files.
        insect_name_str = "Insects"
        ! check if insect parameters are under legacy name
        call param_section_exists_mpi(PARAMS, insect_name_str, section_exists)
        if (.not. section_exists) then
            ! if not, then we use "Insect1" as section name for the first insect, to be consistent with the other ones
            write(insect_name_str, '(A,I0)') "Insect", Insect_ID
            call param_section_exists_mpi(PARAMS, insect_name_str, section_exists)
            ! if this is also not found, then we give up
            if (.not. section_exists) then
                call abort(260602, "Insect flew away, no section found for " // trim(insect_name_str))
            endif
        endif
    endif

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
        write(*,'(A)') "*.ini file is: "//trim(adjustl(fname_ini))// " and insect is in section "//trim(adjustl(insect_name_str))
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
    Insects(Insect_ID)%epsilon_hester = sqrt(nu * c_eta_in)
    Insects(Insect_ID)%smoothing_type = smoothing_type
    select case(smoothing_type)
        case("cos", "cosine")
            Insects(Insect_ID)%smoothing_type_int = STEP_METHOD_COSINE
        case ("hester")
            Insects(Insect_ID)%smoothing_type_int = STEP_METHOD_HESTER
        case("discontinuous", "dis")
            Insects(Insect_ID)%smoothing_type_int = STEP_METHOD_DISC
        case default
            call abort(260602, "Insect hit a tree - Never heard of the smoothing type "//trim(smoothing_type))
    end select

    ! 2025-11-07: insect periodization disabled, TE. It can be tricky to compute the properly periodized mask function 
    ! (eg bounding boxes) and in WABBIT, we'd just use a large enough domain in most cases. In particular, also STL periodization
    ! is probably painful.
    ! ! is the insect periodic?
    ! ! attention: this functionality is not for free, so if you do not need it - disable it.
    ! if (present(periodic)) then
    !     periodic_insect = periodic
    ! else
    !     periodic_insect = .true.
    ! endif
    ! if (root) write(*,'("periodic_insect=",L1)') periodic_insect

    !-----------------------------------------------------------------------------
    ! read in parameters form ini file
    !-----------------------------------------------------------------------------

    ! determine whether the second pair of wings is present
    call read_param_mpi(PARAMS,insect_name_str,"LeftWing2",Insects(Insect_ID)%LeftWing2,"no")
    call read_param_mpi(PARAMS,insect_name_str,"RightWing2",Insects(Insect_ID)%RightWing2,"no")

    if ( ( Insects(Insect_ID)%LeftWing2 == "yes" ) .or. ( Insects(Insect_ID)%RightWing2 == "yes" ) ) then
        Insects(Insect_ID)%second_wing_pair = .true.
    else
        Insects(Insect_ID)%second_wing_pair = .false.
    endif

    ! read data for the first pair of wings
    call read_param_mpi(PARAMS,insect_name_str,"WingShape",Insects(Insect_ID)%WingShape(1),"none")
    Insects(Insect_ID)%WingShape(2) = Insects(Insect_ID)%WingShape(1)
    ! The following two lines take effect if WingShapeL or WingShapeR are set
    call read_param_mpi(PARAMS,insect_name_str,"WingShapeL",Insects(Insect_ID)%WingShape(1),Insects(Insect_ID)%WingShape(1))
    call read_param_mpi(PARAMS,insect_name_str,"WingShapeR",Insects(Insect_ID)%WingShape(2),Insects(Insect_ID)%WingShape(2))
    ! Rectangular wing parameters
    call read_param_mpi(PARAMS,insect_name_str,"b_top",Insects(Insect_ID)%b_top, 0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"b_bot",Insects(Insect_ID)%b_bot, 0.0_rk)
    ! Kinematics
    call read_param_mpi(PARAMS,insect_name_str,"FlappingMotion_right",Insects(Insect_ID)%FlappingMotion_right,"none")
    call read_param_mpi(PARAMS,insect_name_str,"FlappingMotion_left",Insects(Insect_ID)%FlappingMotion_left,"none")
    ! this file is used in old syntax form for both wings:
    call read_param_mpi(PARAMS,insect_name_str,"infile",Insects(Insect_ID)%infile,"none.in")

    if ( index(Insects(Insect_ID)%FlappingMotion_right,"from_file::") /= 0 ) then
        ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
        Insects(Insect_ID)%kine_wing_r%infile = Insects(Insect_ID)%FlappingMotion_right( 12:clong  )
        Insects(Insect_ID)%FlappingMotion_right = "from_file"        
    elseif ( Insects(Insect_ID)%FlappingMotion_right == "from_file" ) then
        ! old syntax, implies symmetric periodic motion, read from *.ini file
        Insects(Insect_ID)%kine_wing_r%infile = Insects(Insect_ID)%infile
    endif

    if ( index(Insects(Insect_ID)%FlappingMotion_left,"from_file::") /= 0 ) then
        ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
        Insects(Insect_ID)%kine_wing_l%infile = Insects(Insect_ID)%FlappingMotion_left( 12:clong  )
        Insects(Insect_ID)%FlappingMotion_left = "from_file"
    elseif ( Insects(Insect_ID)%FlappingMotion_left == "from_file" ) then
        ! old syntax, implies symmetric periodic motion, read from *.ini file
        Insects(Insect_ID)%kine_wing_l%infile = Insects(Insect_ID)%infile
    endif

    if (root) then
        write(*,*) "Left wing: "//trim(adjustl(Insects(Insect_ID)%FlappingMotion_left))
        write(*,*) "Left wing: "//trim(adjustl(Insects(Insect_ID)%kine_wing_l%infile))
        write(*,*) "Right wing: "//trim(adjustl(Insects(Insect_ID)%FlappingMotion_right))
        write(*,*) "Right wing: "//trim(adjustl(Insects(Insect_ID)%kine_wing_r%infile))
    endif

    ! these flags trigger reading the kinematics from file when the Flapping
    ! motion is first called
    Insects(Insect_ID)%kine_wing_l%initialized = .false.
    Insects(Insect_ID)%kine_wing_r%initialized = .false.

    ! read data for the second pair of wings
    if (Insects(Insect_ID)%second_wing_pair) then
        ! note that only Fourier wing shape can be different from first wings
        call read_param_mpi(PARAMS,insect_name_str,"WingShape2",Insects(Insect_ID)%WingShape(3),"none")
        Insects(Insect_ID)%WingShape(4) = Insects(Insect_ID)%WingShape(3)
        ! The following two lines take effect if WingShape2L or WingShape2R are set
        call read_param_mpi(PARAMS,insect_name_str,"WingShape2L",Insects(Insect_ID)%WingShape(3),Insects(Insect_ID)%WingShape(3))
        call read_param_mpi(PARAMS,insect_name_str,"WingShape2R",Insects(Insect_ID)%WingShape(4),Insects(Insect_ID)%WingShape(4))
        ! Kinematics
        call read_param_mpi(PARAMS,insect_name_str,"FlappingMotion_right2",Insects(Insect_ID)%FlappingMotion_right2,"none")
        call read_param_mpi(PARAMS,insect_name_str,"FlappingMotion_left2",Insects(Insect_ID)%FlappingMotion_left2,"none")
        ! this file is used in old syntax form for both wings:
        call read_param_mpi(PARAMS,insect_name_str,"infile2",Insects(Insect_ID)%infile2,"none.in")

        if ( index(Insects(Insect_ID)%FlappingMotion_right2,"from_file::") /= 0 ) then
            ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
            Insects(Insect_ID)%kine_wing_r2%infile = Insects(Insect_ID)%FlappingMotion_right2( 12:clong  )
            Insects(Insect_ID)%FlappingMotion_right2 = "from_file"

        elseif ( index(Insects(Insect_ID)%FlappingMotion_right2,"kinematics_loader::") /= 0 ) then
            ! new syntax, uses the kinematics loader for non-periodic kinematics
            Insects(Insect_ID)%kine_wing_r2%infile = Insects(Insect_ID)%FlappingMotion_right2( 20:clong )
            Insects(Insect_ID)%FlappingMotion_right2 = "kinematics_loader"

        elseif ( Insects(Insect_ID)%FlappingMotion_right2 == "from_file" ) then
            ! old syntax, implies symmetric periodic motion, read from *.ini file
            Insects(Insect_ID)%kine_wing_r2%infile = Insects(Insect_ID)%infile2

        elseif ( Insects(Insect_ID)%FlappingMotion_right2 == "kinematics_loader" ) then
            ! old syntax, implies symmetric non-periodic motion, read from *.dat file
            Insects(Insect_ID)%kine_wing_r2%infile = Insects(Insect_ID)%infile2
        endif

        if ( index(Insects(Insect_ID)%FlappingMotion_left2,"from_file::") /= 0 ) then
            ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
            Insects(Insect_ID)%kine_wing_l2%infile = Insects(Insect_ID)%FlappingMotion_left2( 12:clong  )
            Insects(Insect_ID)%FlappingMotion_left2 = "from_file"

        elseif ( index(Insects(Insect_ID)%FlappingMotion_left2,"kinematics_loader::") /= 0 ) then
            ! new syntax, uses the kinematics loader for non-periodic kinematics
            Insects(Insect_ID)%kine_wing_l2%infile = Insects(Insect_ID)%FlappingMotion_left2( 20:clong )
            Insects(Insect_ID)%FlappingMotion_left2 = "kinematics_loader"

        elseif ( Insects(Insect_ID)%FlappingMotion_left2 == "from_file" ) then
            ! old syntax, implies symmetric periodic motion, read from *.ini file
            Insects(Insect_ID)%kine_wing_l2%infile = Insects(Insect_ID)%infile2

        elseif ( Insects(Insect_ID)%FlappingMotion_left2 == "kinematics_loader" ) then
            ! old syntax, implies symmetric non-periodic motion, read from *.dat file
            Insects(Insect_ID)%kine_wing_l2%infile = Insects(Insect_ID)%infile2
        endif

        if (root) then
            write(*,*) "Second left wing: "//trim(adjustl(Insects(Insect_ID)%FlappingMotion_left2))
            write(*,*) "Second left wing: "//trim(adjustl(Insects(Insect_ID)%kine_wing_l2%infile))
            write(*,*) "Second right wing: "//trim(adjustl(Insects(Insect_ID)%FlappingMotion_right2))
            write(*,*) "Second right wing: "//trim(adjustl(Insects(Insect_ID)%kine_wing_r2%infile))
        endif
    endif

    ! these flags trigger reading the kinematics from file when the Flapping
    ! motion is first called
    Insects(Insect_ID)%kine_wing_l2%initialized = .false.
    Insects(Insect_ID)%kine_wing_r2%initialized = .false.

    call read_param_mpi(PARAMS,insect_name_str,"BodyType",Insects(Insect_ID)%BodyType,"ellipsoid")
    call read_param_mpi(PARAMS,insect_name_str,"BodyMotion",Insects(Insect_ID)%BodyMotion,"tethered")
    ! this one file contains all the kinematics (four wings and body)
    call read_param_mpi(PARAMS,insect_name_str,"infile_kineloader",Insects(Insect_ID)%infile_kineloader,"none")

    call read_param_mpi(PARAMS,insect_name_str,"LeftWing",Insects(Insect_ID)%LeftWing,"yes")
    call read_param_mpi(PARAMS,insect_name_str,"RightWing",Insects(Insect_ID)%RightWing,"yes")
    call read_param_mpi(PARAMS,insect_name_str,"mass",Insects(Insect_ID)%mass, 1._rk)
    call read_param_mpi(PARAMS,insect_name_str,"gravity",Insects(Insect_ID)%gravity, 0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"gravity_x",Insects(Insect_ID)%gravity_x, 0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"gravity_y",Insects(Insect_ID)%gravity_y, 0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"WingThickness",Insects(Insect_ID)%WingThickness, 4.0_rk*dx_reference)
    call read_param_mpi(PARAMS,insect_name_str,"J_body_yawpitchroll",defaultvec, (/0.0_rk,0.0_rk,0.0_rk/))
    Insects(Insect_ID)%Jroll_body  = defaultvec(3)
    Insects(Insect_ID)%Jyaw_body   = defaultvec(1)
    Insects(Insect_ID)%Jpitch_body = defaultvec(2)
    call read_param_mpi(PARAMS,insect_name_str,"x0",Insects(Insect_ID)%x0, (/0.5_rk*xl,0.5_rk*yl,0.5_rk*zl/))
    call read_param_mpi(PARAMS,insect_name_str,"v0",Insects(Insect_ID)%v0, (/0.0_rk, 0.0_rk, 0.0_rk/))
    call read_param_mpi(PARAMS,insect_name_str,"yawpitchroll_0",Insects(Insect_ID)%yawpitchroll_0,(/0.0_rk, 0.0_rk, 0.0_rk/))
    call read_param_mpi(PARAMS,insect_name_str,"yawpitchroll_a1",Insects(Insect_ID)%yawpitchroll_a1,(/0.0_rk, 0.0_rk, 0.0_rk/))
    call read_param_mpi(PARAMS,insect_name_str,"yawpitchroll_b1",Insects(Insect_ID)%yawpitchroll_b1,(/0.0_rk, 0.0_rk, 0.0_rk/))
    ! convert yawpitchroll to radiants
    Insects(Insect_ID)%yawpitchroll_0 = Insects(Insect_ID)%yawpitchroll_0 * (pi/180.0_rk)
    Insects(Insect_ID)%yawpitchroll_a1 = Insects(Insect_ID)%yawpitchroll_a1 * (pi/180.0_rk)
    Insects(Insect_ID)%yawpitchroll_b1 = Insects(Insect_ID)%yawpitchroll_b1 * (pi/180.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"eta0",Insects(Insect_ID)%eta0, 0.0_rk)
    Insects(Insect_ID)%eta0 = Insects(Insect_ID)%eta0*(pi/180.0_rk)

    call read_param_mpi(PARAMS,insect_name_str,"pointcloudfile",Insects(Insect_ID)%pointcloudfile,"none")



    ! degrees of freedom for free flight solver. The string from ini file contains
    ! 6 characters 1 or 0 that turn on/off x,y,z,yaw,pitch,roll degrees of freedom
    ! by multiplying the respective RHS by zero, keeping the value thus constant
    call read_param_mpi(PARAMS,insect_name_str,"DoF",DoF_string, "111111")
    do j=1,6
        read (DoF_string(j:j), '(i1)') tmp
        Insects(Insect_ID)%DoF_on_off(j) = dble(tmp)
    enddo
    if (root) write(*,'(6(f4.2,1x))') Insects(Insect_ID)%DoF_on_off

    ! ! section for additional fractal tree
    ! call read_param_mpi(PARAMS, "Insects", "fractal_tree", Insects(Insect_ID)%fractal_tree, .false.)
    ! call read_param_mpi(PARAMS, "Insects", "fractal_tree_file", Insects(Insect_ID)%fractal_tree_file, "tree_data.in")
    ! call read_param_mpi(PARAMS, "Insects", "fractal_tree_x0", Insects(Insect_ID)%fractal_tree_x0, (/0.0_rk, 0.0_rk, 0.0_rk/) )
    ! call read_param_mpi(PARAMS, "Insects", "fractal_tree_scaling", Insects(Insect_ID)%fractal_tree_scaling, 1.0_rk )
    ! call read_param_mpi(PARAMS, "Insects", "fractal_tree_spheres", Insects(Insect_ID)%fractal_tree_spheres, .true. )

    


    ! wing inertia tensor (we currently assume two identical forewings and two identical hindwings)
    ! this allows computing inertial power and wing FSI model
    call read_param_mpi(PARAMS,insect_name_str,"Jxx",Insects(Insect_ID)%Jxx,0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"Jyy",Insects(Insect_ID)%Jyy,0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"Jzz",Insects(Insect_ID)%Jzz,0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"Jxy",Insects(Insect_ID)%Jxy,0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"Jxx2",Insects(Insect_ID)%Jxx2,0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"Jyy2",Insects(Insect_ID)%Jyy2,0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"Jzz2",Insects(Insect_ID)%Jzz2,0.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"Jxy2",Insects(Insect_ID)%Jxy2,0.0_rk)

    call read_param_mpi(PARAMS,insect_name_str,"startup_conditioner",Insects(Insect_ID)%startup_conditioner,"no")

    ! 28/01/2019: Thomas. Discovered that this was done block based, i.e. the smoothing layer
    ! had different thickness, if some blocks happened to be at different levels (and still carry
    ! a part of the smoothing layer.) I don't know if that made sense, because the layer shrinks/expands then
    ! and because it might be discontinous. Both options are included now, default is "as before"
    ! Insect%smoothing_thickness=="local"  : smoothing_layer = c_sm * 2**-J * L/(BS-1)
    ! Insect%smoothing_thickness=="global" : smoothing_layer = c_sm * 2**-Jmax * L/(BS-1)
    ! NOTE: for FLUSI, this has no impact! Here, the grid is constant and equidistant.
    ! NOTE: 05/2020 Thomas, I changed the default back to local.
    call read_param_mpi(PARAMS,insect_name_str,"smoothing_thickness",Insects(Insect_ID)%smoothing_thickness,"local")
    call read_param_mpi(PARAMS,insect_name_str,"C_smooth",Insects(Insect_ID)%C_smooth,1.0_rk)
    call read_param_mpi(PARAMS,insect_name_str,"BodySuperSTLfile",Insects(Insect_ID)%BodySuperSTLfile,"none.superstl")
    ! when using CT data, code computes the mask function in a shell around fluid-solid interface.
    ! The tickness of the shell is not a critical parameter, but it affects performance. Thicker shell
    ! means more points and thus more comput effort. It is given in multiples of C_smooth, that means
    ! shell_thickness = C_shell_thickness * C_smooth * dx_min. dx_min is the spacing on the finest level Jmax.
    ! Why is the shell thickness dependent on resolution? The cost to generate the mask depends on the number of
    ! triangles and the number of points in the shell. As the latter is coupled to dx and the former is constant
    ! this way the mask generation cost is constant when increasing the resolution.
    call read_param_mpi(PARAMS,insect_name_str,"C_shell_thickness",Insects(Insect_ID)%C_shell_thickness, 3.0_rk)

    Insects(Insect_ID)%dx_reference = dx_reference
    Insects(Insect_ID)%smooth = Insects(Insect_ID)%C_smooth*dx_reference
    if (Insects(Insect_ID)%smoothing_type == "hester") then
        Insects(Insect_ID)%smooth = Insects(Insect_ID)%epsilon_hester
        Insects(Insect_ID)%safety = max(5.0_rk*Insects(Insect_ID)%epsilon_hester, 2*dx_reference)
    else
        Insects(Insect_ID)%safety = 3.5_rk*Insects(Insect_ID)%smooth
    end if

    ! wing hinges (root points)
    defaultvec=(/0.0_rk, +Insects(Insect_ID)%b_body, 0.0_rk /)
    call read_param_mpi(PARAMS,insect_name_str,"x_pivot_l",Insects(Insect_ID)%x_pivot_l_b, defaultvec)

    defaultvec=(/0.0_rk, -Insects(Insect_ID)%b_body, 0.0_rk /)
    call read_param_mpi(PARAMS,insect_name_str,"x_pivot_r",Insects(Insect_ID)%x_pivot_r_b, defaultvec)

    ! read data for the second pair of wing hinges
    if (Insects(Insect_ID)%second_wing_pair) then
        defaultvec = Insects(Insect_ID)%x_pivot_l_b
        call read_param_mpi(PARAMS,insect_name_str,"x_pivot_l2",Insects(Insect_ID)%x_pivot_l2_b, defaultvec)
        defaultvec = Insects(Insect_ID)%x_pivot_r_b
        call read_param_mpi(PARAMS,insect_name_str,"x_pivot_r2",Insects(Insect_ID)%x_pivot_r2_b, defaultvec)
    endif

    ! default colors for body, left wing, right wing, left wing 2, right wing 2, geometry / full insect
    defaultvec5 = (/1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk/)
    if (present(colors_default)) defaultvec5 = real(colors_default(1:5), kind=rk)

    call read_param_mpi(PARAMS,insect_name_str,"colors", defaultvec5, defaultvec5)
    Insects(Insect_ID)%color_body = int(defaultvec5(1), kind=2)
    Insects(Insect_ID)%color_l    = int(defaultvec5(2), kind=2)
    Insects(Insect_ID)%color_r    = int(defaultvec5(3), kind=2)
    Insects(Insect_ID)%color_l2   = int(defaultvec5(4), kind=2)
    Insects(Insect_ID)%color_r2   = int(defaultvec5(5), kind=2)
    ! the color_geometry is used only for statistics of forces/moments wrt to the full insect
    ! physics module should provide this as input
    if (present(colors_default)) then
        Insects(Insect_ID)%color_geometry = colors_default(6)
    else
        Insects(Insect_ID)%color_geometry = 0_2
    endif

    ! clean ini file
    call clean_ini_file_mpi(PARAMS)

    !---------------------------------------------------------------------------
    ! initialization for superSTl body
    !---------------------------------------------------------------------------
    if (Insects(Insect_ID)%BodyType == "superSTL") then
        if (.not. allocated(Insects(Insect_ID)%body_superSTL_b)) then
            if (root) write(*,'("INSECTS: STL: init start")')
            if (root) write(*,'("INSECTS: STL: file=",A)') Insects(Insect_ID)%BodySuperSTLfile

            call count_lines_in_ascii_file_mpi(Insects(Insect_ID)%BodySuperSTLfile, ntri, 0)

            if (root) write(*,'("INSECTS: STL: file length is ntri=", i7 )') ntri

            allocate( Insects(Insect_ID)%body_superSTL_b(1:ntri,1:30) )
            allocate( Insects(Insect_ID)%body_superSTL_g(1:ntri,1:30) )

            ! No scaling or origin shift is applied: we assume you did that when generating
            ! the superSTL file. The data is thus understood in the body coordinate system.
            call read_array_from_ascii_file_mpi(Insects(Insect_ID)%BodySuperSTLfile, Insects(Insect_ID)%body_superSTL_b, 0)

            if (root) write(*,'("INSECTS: STL: read from file...done! We are good to go.")')
        endif
    endif


    ! !-----------------------------------------------------------------------------
    ! ! other initialization
    ! !-----------------------------------------------------------------------------
    ! if (Insects(Insect_ID)%fractal_tree) then
    !     ! we can also simulate an insect together with a fractal tree as turbulence
    !     ! generators.
    !     call fractal_tree_init(Insects(Insect_ID))
    ! endif

    ! If required, initialize rigid solid dynamics solver
    if (Insects(Insect_ID)%BodyMotion=="free_flight") then
        ! note we have to do that before init_fields as rigid_solid_init sets up
        ! the state vector without which create_mask cannot know the position and velocity
        ! of body and wings
        call rigid_solid_init( time, Insects(Insect_ID), resume_backup, Insect_ID )
    endif

    ! the update routine computes wing angles and so on, everything that is done only
    ! once per time step. Do this here as well, so we can safely call draw_insect after
    ! calling this routine.
    call Update_Insect( time, Insects(Insect_ID) )


    ! At this point, we must also initialize all wing / body data: in wabbit, it may
    ! be that the initial grid contains less blocks than we wave MPIRANKS. Therefore,
    ! not all CPUS might call the MPI_BCAST for this initialization. (This is a BUGFIX
    ! 21 Oct 2019, Yokohama, Thomas)
    Nwings = 2
    if (Insects(Insect_ID)%second_wing_pair) Nwings = 4

    ! wing id number: 1 = left, 2 = right, 3 = 2nd left, 4 = 2nd right
    do wingID = 1, Nwings
        ! exclude wings that are hard-coded, otherwise, call initialization routine
        if (Insects(Insect_ID)%WingShape(wingID)/="pointcloud" .and. Insects(Insect_ID)%WingShape(wingID)/="mosquito_iams" .and. &
        Insects(Insect_ID)%WingShape(wingID)/="suzuki" .and. Insects(Insect_ID)%WingShape(wingID)/="rectangular" .and. &
        Insects(Insect_ID)%WingShape(wingID)/="TwoEllipses" .and. (Insects(Insect_ID)%wing_file_type(wingID)) /= "kleemeier" &
        .and. Insects(Insect_ID)%WingShape(wingID)/="suzuki_butterfly".and. Insects(Insect_ID)%WingShape(wingID)/="none") then

        ! we have some pre-defined, hard-coded data, but also can read the wing shape
        ! from INI files.
        call Setup_Wing_Fourier_coefficients(Insects(Insect_ID), wingID)
    endif

enddo


if (root) then
    write(*,'(80("<"))')
    write(*,*) "Insect initialization is complete."
    write(*,'(80("<"))')
endif
Insects(Insect_ID)%initialized = .true.

end subroutine insect_init


! -------------------------------------------------------------------------------
! \brief Clean up insect data structure for all insects or one specific one
subroutine insects_clean(insect_i)
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
        call insect_clean(Insects(Insects_i))
    enddo

    if (allocated(Insects)) deallocate(Insects)
   
end subroutine insects_clean

! -------------------------------------------------------------------------------
! \brief Clean up insect data structure of one insect
subroutine insect_clean(insect)
    implicit none

    type(diptera), intent(inout) :: insect

    call wingkinematics_clean(insect%kine_wing_l)
    call wingkinematics_clean(insect%kine_wing_r)
    call wingkinematics_clean(insect%kine_wing_l2)
    call wingkinematics_clean(insect%kine_wing_r2)

    if (allocated(insect%data_kineloader)) deallocate(insect%data_kineloader)
    if (allocated(insect%RHS)) deallocate(insect%RHS)
    if (allocated(insect%ai_wings)) deallocate(insect%ai_wings)
    if (allocated(insect%bi_wings)) deallocate(insect%bi_wings)
    if (allocated(insect%R0_table)) deallocate(insect%R0_table)
    if (allocated(insect%theta_i)) deallocate(insect%theta_i)
    if (allocated(insect%R_i)) deallocate(insect%R_i)
    if (allocated(insect%bristles_coords)) deallocate(insect%bristles_coords)
    if (allocated(insect%particle_points)) deallocate ( insect%particle_points )
    if (allocated(insect%wing_thickness_profile)) deallocate ( insect%wing_thickness_profile )
    if (allocated(insect%corrugation_profile)) deallocate ( insect%corrugation_profile )
    if (allocated(insect%damage_mask)) deallocate ( insect%damage_mask )
    if (allocated(insect%deformations)) deallocate ( insect%deformations )
    if (allocated(insect%deformation_profile)) deallocate ( insect%deformation_profile )
    if (allocated(insect%mask_wing_complete)) deallocate ( insect%mask_wing_complete )
    if (allocated(insect%body_superSTL_b)) deallocate ( insect%body_superSTL_b )
    if (allocated(insect%body_superSTL_g)) deallocate ( insect%body_superSTL_g )

end subroutine insect_clean

subroutine wingkinematics_clean(wing_kinematics)
    implicit none

    type(wingkinematics), intent(inout) :: wing_kinematics

    if (allocated(wing_kinematics%ai_phi)) deallocate(wing_kinematics%ai_phi)
    if (allocated(wing_kinematics%bi_phi)) deallocate(wing_kinematics%bi_phi)
    if (allocated(wing_kinematics%ai_theta)) deallocate(wing_kinematics%ai_theta)
    if (allocated(wing_kinematics%bi_theta)) deallocate(wing_kinematics%bi_theta)
    if (allocated(wing_kinematics%ai_alpha)) deallocate(wing_kinematics%ai_alpha)
    if (allocated(wing_kinematics%bi_alpha)) deallocate(wing_kinematics%bi_alpha)

end subroutine wingkinematics_clean