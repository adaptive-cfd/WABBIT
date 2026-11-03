subroutine init_beams( fname_ini, dx_min_in, time, backup_resuming )
    use module_ini_files_parser_mpi

    !---------------------------------------------------
    ! initializes an array of beams. the initial state is always
    ! straight lines, possible oriented with different angles, at rest.
    !---------------------------------------------------
    implicit none
    character(len=*), intent(in) :: fname_ini
    real(kind=rk), intent(in) :: dx_min_in, time
    logical, intent(in) :: backup_resuming

    integer :: n, i, a
    real (kind=rk) :: alpha, alpha_t, alpha_tt
    ! LeadingEdge: x, y, vx, vy, ax, ay (Array)
    real (kind=rk), dimension(1:6) :: LeadingEdge
    character(len=1)  :: beamstr
    character(len=16) :: frmt
    integer(kind=ik) :: mpicode
    character(len=clong) :: fname
    type(inifile) :: PARAMS

    dx_min = dx_min_in

    ! in this module, we use the logical ROOT to avoid the integer comparison mpirank==0
    call MPI_COMM_RANK (MPI_COMM_WORLD, mpirank, mpicode)
    if (mpirank==0) root = .true.

    if (mpirank==0) then
        write (*,'("Initializing ",i1," beams")')  nBeams
    endif

    call lapack_unit_test()

    !-------------------------------------------
    ! allocate beam storage for each beam
    !-------------------------------------------
    do i = 1, nBeams
        !----------------------------------------------
        ! read in beam parameters from INI file
        !----------------------------------------------
        ! read in the complete ini file, from which we initialize the insect
        call read_ini_file_mpi(PARAMS, fname_ini, verbose=.true.)

        !-- beam resolution
        call read_param_mpi(PARAMS,"SolidModel","ns",ns, 32)

        allocate ( beams(i)%x(0:ns-1) )
        allocate ( beams(i)%y(0:ns-1) )
        allocate ( beams(i)%vx(0:ns-1) )
        allocate ( beams(i)%vy(0:ns-1) )
        allocate ( beams(i)%ax(0:ns-1) )
        allocate ( beams(i)%ay(0:ns-1) )
        ! note both theta and T include their ghosts [actually: theta(-1:ns+1), T(-1:ns-1)]
        allocate ( beams(i)%theta(-1:ns+1) )
        allocate ( beams(i)%theta_dot(0:ns-1) )
        allocate ( beams(i)%T(-1:ns-1) )
        ! external forces:
        allocate ( beams(i)%pressure_old(0:ns-1) )
        allocate ( beams(i)%pressure_new(0:ns-1) )
        allocate ( beams(i)%tau_old(0:ns-1) )
        allocate ( beams(i)%tau_new(0:ns-1) )
        ! values of theta and theta_dot at times t^(n) and t^(n-1)
        allocate ( beams(i)%theta_old(0:ns-1) )
        allocate ( beams(i)%theta_oldold(0:ns-1) )
        allocate ( beams(i)%theta_dot_old(0:ns-1) )
        allocate ( beams(i)%theta_dot_oldold(0:ns-1) )
        ! material properties and their derivatives:
        allocate ( beams(i)%mu(0:ns-1) )
        allocate ( beams(i)%mu_s(0:ns-1) )
        allocate ( beams(i)%mu_star(0:ns-1) )
        allocate ( beams(i)%zeta(0:ns-1) )
        allocate ( beams(i)%zeta_s(0:ns-1) )
        allocate ( beams(i)%zeta_ss(0:ns-1) )
        allocate ( beams(i)%zeta_sss(0:ns-1) )
        ! grid s and width in rigid direction:
        allocate ( beams(i)%s(0:ns-1) )
        allocate ( beams(i)%L_rigid(0:ns-1) )

        !-- density / stiffness / gravity
        call read_param_mpi(PARAMS,"SolidModel","mue",mue0,1.0d0)
        call read_param_mpi(PARAMS,"SolidModel","eta",eta0,1.0d0)
        call read_param_mpi(PARAMS,"SolidModel","f",frequ,1.0d0)
        call read_param_mpi(PARAMS,"SolidModel","angle",AngleBeam,0.0d0)
        call read_param_mpi(PARAMS,"SolidModel","gravity",grav,0.0d0)
        !-- beam thickness (from midline to surface, so top to bottom thickness is twice that value)
        call read_param_mpi(PARAMS,"SolidModel","t_beam", t_beam, 0.05d0)
        !-- damping
        call read_param_mpi(PARAMS,"SolidModel","sigma",sigma,0.0d0)
        !-- timing
        call read_param_mpi(PARAMS,"SolidModel","T_release",T_release,0.0d0)
        call read_param_mpi(PARAMS,"SolidModel","tau",tau,0.0d0)
        call read_param_mpi(PARAMS,"SolidModel","L_span",L_span,1.0d0)
        !-- time marching method for the solid
        call read_param_mpi(PARAMS,"SolidModel","TimeMethodSolid",TimeMethodSolid,"BDF2")

        call read_param_mpi(PARAMS,"SolidModel","imposed_motion_leadingedge", imposed_motion_leadingedge,"fixed_middle")
        call read_param_mpi(PARAMS,"SolidModel","plate_shape",plate_shape,"rectangular")

        call read_param_mpi(PARAMS,"SolidModel", "x0_origin", beams(i)%x0_origin, (/0.0_rk, 0.0_rk, 0.0_rk/) )

        call read_param_mpi(PARAMS,"VPM", "C_smooth", C_smooth, 1.5_rk )

        write (*,'("Beam width (2*t_beam) is ",es12.4)') 2.0_rk*t_beam
        write (*,'("Beam width (2*t_beam) covers ",(es12.4)," points")') 2.0_rk*t_beam/dx_min

        ! clean ini file
        call clean_ini_file_mpi(PARAMS)

        !-- grid spacing
        ds = 1.0_rk/dble(ns-1)


        !---------------------------------------------
        ! define adjustable parameters for each beam
        ! this is position and motion protocoll
        !--------------------------------------------
        beams(i)%AngleBeam = AngleBeam
        beams(i)%phase = 0.0_rk

        !--------------------------------------
        !-- initialize beam
        !--------------------------------------
        ! fetch leading edge position
        call mouvement(0.0_rk, alpha, alpha_t, alpha_tt, LeadingEdge, beams(i) )

        ! initialize as zero
        beams(i)%x = 0.0_rk
        beams(i)%y = 0.0_rk
        beams(i)%vx = 0.0_rk
        beams(i)%vy = 0.0_rk
        beams(i)%ax = 0.0_rk
        beams(i)%ay = 0.0_rk
        beams(i)%T  = 0.0_rk
        beams(i)%theta = 0.0_rk
        beams(i)%theta_dot = 0.0_rk
        beams(i)%theta_old = 0.0_rk
        beams(i)%theta_dot_old = 0.0_rk
        beams(i)%theta_oldold = 0.0_rk
        beams(i)%theta_dot_oldold = 0.0_rk
        beams(i)%pressure_old = 0.0_rk
        beams(i)%pressure_new = 0.0_rk
        beams(i)%tau_old = 0.0_rk
        beams(i)%tau_new = 0.0_rk
        beams(i)%Inertial_Force = 0.0_rk
        ! this is used for unst correction computation:
        beams(i)%Force = 0.0_rk
        beams(i)%Force_unst = 0.0_rk
        beams(i)%Force_press = 0.0_rk
        beams(i)%Inertial_Force = 0.0_rk
        beams(i)%E_kinetic = 0.0_rk
        beams(i)%E_elastic = 0.0_rk
        beams(i)%E_pot = 0.0_rk
        ! To use CN2 time stepper on the very first instant, because BDF2 is not available
        beams(i)%StartupStep = .true.
        beams(i)%dt_old = 0.0_rk

        ! remember that what follows is described in relative coordinates
        do n = 0, ns-1
            beams(i)%x(n) = LeadingEdge(1) + dble(n)*ds
            beams(i)%y(n) = LeadingEdge(2)
            beams(i)%vx(n)= 0.0_rk
            beams(i)%vy(n)= 0.0_rk
            beams(i)%ax(n)= 0.0_rk
            beams(i)%ay(n)= 0.0_rk
            ! the grid:
            beams(i)%s(n) = real(n, kind=rk)*ds
            ! width in rigid (usually span) direction:
            beams(i)%L_rigid(n) = z_top(beams(i)%s(n)) + z_bottom(beams(i)%s(n))
        enddo

        !---------------------------------------------------------------------------
        ! set up material and derivatives
        !---------------------------------------------------------------------------
        beams(i)%zeta(0:ns-1) = eta0 * beams(i)%L_rigid(0:ns-1)
        beams(i)%mu(0:ns-1)   = mue0 * beams(i)%L_rigid(0:ns-1)
        call Differentiate1D (beams(i)%zeta, beams(i)%zeta_s, ns, ds, 1)
        call Differentiate1D (beams(i)%zeta, beams(i)%zeta_ss, ns, ds, 2)
        call Differentiate1D (beams(i)%zeta, beams(i)%zeta_sss, ns, ds, 3)
        call Differentiate1D (beams(i)%mu, beams(i)%mu_s, ns, ds, 1)
        beams(i)%mu_star(0:ns-1) = beams(i)%mu_s(0:ns-1) / beams(i)%mu(0:ns-1)


        if (root) then
            write(*,'(80("-"))')
            write(*,'("setting up solid material with mu=",es12.4," and eta=",es12.4)') mue0, eta0
            write(frmt,'("(",i3.3,"(es12.4,1x))")') ns
            write(*,*) "Length in rigid direction:"
            write(*,frmt) beams(i)%L_rigid(0:ns-1)
            write(*,*) "density coefficient mu:"
            write(*,frmt) beams(i)%mu(0:ns-1)
            write(*,*) "stiffness coefficient eta:"
            write(*,frmt) beams(i)%zeta(0:ns-1)
            write(*,'(80("-"))')
        endif

        ! to take static angle into account
        call integrate_position (0.0_rk, beams(i))
    enddo


    !-------------------------------------------
    ! If we resume a backup, read from file (all ranks do that)
    !-------------------------------------------
    if ( backup_resuming ) then
        write( fname,'(a, "_", a, ".dat")') "solid", trim(adjustl(timestr(time)))

        call read_solid_backup( trim(adjustl(fname)) )
    endif

    call integrate_position (0.0_rk, beams(1))


end subroutine init_beams
