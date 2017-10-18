!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name ini_file_to_params.f90
!> \version 0.5
!> \author msr
!
!> \brief distribute blocks at start => create light data array
!
!>
!! input:    - filename \n
!! output:   - filled parameter struct \n
!!
!!
!! = log ======================================================================================
!! \n
!! 25/01/17    - create \n
!! 29/01/17    - add filter parameter \n
!! 30/01/17    - add automatic memory management
!
! ********************************************************************************************

subroutine ini_file_to_params( params, filename )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)               :: params
    !> inifile name
    character(len=*), intent(in)                   :: filename

    ! process rank
    integer(kind=ik)                                :: rank
    ! number of processes
    integer(kind=ik)                                :: number_procs
    ! inifile structure
    type(inifile)                                   :: FILE
    ! maximum memory available on all cpus
    real(kind=rk)                                   :: maxmem
    ! power used for dimensionality (d=2 or d=3)
    integer(kind=ik)                                :: d,i
    ! string read from command line call
    character(len=80)                               :: memstring

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

!---------------------------------------------------------------------------------------------
! main body

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)


    !***************************************************************************
    ! read BLOCK parameters
    ! read number_block_nodes
    call read_param_mpi(FILE, 'Blocks', 'number_block_nodes', params%number_block_nodes, 1 )
    ! read number_ghost_nodes
    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%number_ghost_nodes, 1 )
    ! read number_blocks
    call read_param_mpi(FILE, 'Blocks', 'number_blocks', params%number_blocks, 1 )

    ! read number_data_fields
    call read_param_mpi(FILE, 'Blocks', 'number_data_fields', params%number_data_fields, 1 )
    ! set number of fields in heavy work data
    ! every datafield has 5 additional fields: old, k1, k2, k3, k4
    params%number_fields = params%number_data_fields*5
    ! read threshold value
    call read_param_mpi(FILE, 'Blocks', 'eps', params%eps, 1e-3_rk )
    ! read treelevel bounds
    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params%max_treelevel, 5 )
    call read_param_mpi(FILE, 'Blocks', 'min_treelevel', params%min_treelevel, 1 )
    ! read switch to turn on|off mesh refinement
    call read_param_mpi(FILE, 'Blocks', 'adapt_mesh', params%adapt_mesh, .true. )
    ! block distribution
    call read_param_mpi(FILE, 'Blocks', 'block_dist', params%block_distribution, "---" )
    ! use non-uniform mesh correction
    call read_param_mpi(FILE, 'Blocks', 'non_uniform_mesh_correction', params%non_uniform_mesh_correction, .true. )

    !***************************************************************************
    ! read TIME parameters
    !
    ! read time_max
    call read_param_mpi(FILE, 'Time', 'time_max', params%time_max, 1.0_rk )
    ! read output write frequency
    call read_param_mpi(FILE, 'Time', 'write_freq', params%write_freq, 25 )
    ! read method to calculate time step
    call read_param_mpi(FILE, 'Time', 'time_step_calc', params%time_step_calc, "CFL_cond" )
    ! read value of fixed time step
    call read_param_mpi(FILE, 'Time', 'dt', params%dt, 1e-6_rk )
    ! read CFL number
    call read_param_mpi(FILE, 'Time', 'CFL', params%CFL, 0.5_rk )
    ! read butcher tableau (set default value to RK4)
    call read_param_mpi(FILE, 'Time', 'butcher_tableau', params%butcher_tableau, &
    reshape((/ 0.0_rk, 0.5_rk, 0.5_rk, 1.0_rk, 0.0_rk, &
    0.0_rk, 0.5_rk, 0.0_rk, 0.0_rk, 1.0_rk/6.0_rk, &
    0.0_rk, 0.0_rk, 0.5_rk, 0.0_rk, 1.0_rk/3.0_rk,&
    0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk/3.0_rk,&
    0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/6.0_rk /), (/ 5,5 /)))

    !***************************************************************************
    ! read PHYSICS parameters
    !
    ! first: read physics type
    call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )

    select case(params%physics_type)

        case('2D_convection_diffusion')

            ! domain size
            call read_param_mpi(FILE, 'Physics', 'Lx', params%Lx, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Ly', params%Ly, 1.0_rk )
            ! set third dimension to zero
            params%Lz = 0.0_rk

            ! allocate memory in params structure (need 2*data_fields for velocity
            ! and 1*data_fields for diffusion coefficient)
            allocate( params%physics%u0( 2*params%number_data_fields ) )
            allocate( params%physics%nu( params%number_data_fields ) )

            ! reset values, use as default values
            params%physics%u0 = 0.0_rk
            params%physics%nu = 0.0_rk

            ! read velocity
            call read_param_mpi(FILE, 'Physics', 'u0', params%physics%u0, params%physics%u0 )
            ! read diffusion
            call read_param_mpi(FILE, 'Physics', 'nu', params%physics%nu, params%physics%nu )

            ! read variable names
            ! allocate names list
            allocate( params%physics%names( params%number_data_fields ) )

            params%physics%names = "---"
            ! read file
            call read_param_mpi(FILE, 'Physics', 'names', params%physics%names, params%physics%names )

        case('2D_navier_stokes')

            ! error case: try to solve navier stokes equation with less or more than 4 datafields
!            if ( params%number_data_fields /= 4) then
!                write(*,'(80("_"))')
!                write(*,'("ERROR: try to solve navier stokes equation with", i3, " datafield(s)")') params%number_data_fields
!                stop
!            end if

            ! domain size
            call read_param_mpi(FILE, 'Physics', 'Lx', params%Lx, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Ly', params%Ly, 1.0_rk )
            ! set third dimension to zero
            params%Lz = 0.0_rk

            ! physics parameter
            ! read adiabatic coefficient
            call read_param_mpi(FILE, 'Physics', 'gamma_', params%physics_ns%gamma_, 0.0_rk )
            ! read specific gas constant
            call read_param_mpi(FILE, 'Physics', 'Rs', params%physics_ns%Rs, 0.0_rk )
            ! calculate isochoric heat capacity
            params%physics_ns%Cv = params%physics_ns%Rs/(params%physics_ns%gamma_-1.0_rk)
            ! calculate isobaric heat capacity
            params%physics_ns%Cp = params%physics_ns%Rs*params%physics_ns%gamma_
            ! read prandtl number
            call read_param_mpi(FILE, 'Physics', 'Pr', params%physics_ns%Pr, 0.0_rk )
            ! read dynamic viscosity
            call read_param_mpi(FILE, 'Physics', 'mu0', params%physics_ns%mu0, 0.0_rk )
            ! read switch to turn on|off dissipation
            call read_param_mpi(FILE, 'Blocks', 'dissipation', params%physics_ns%dissipation, .true. )

            ! read variable names
            ! allocate names list
            allocate( params%physics_ns%names( params%number_data_fields ) )
            params%physics_ns%names = "---"
            ! read file
            call read_param_mpi(FILE, 'Physics', 'names_ns', params%physics_ns%names, params%physics_ns%names )

        case('3D_convection_diffusion')

            ! domain size
            call read_param_mpi(FILE, 'Physics', 'Lx', params%Lx, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Ly', params%Ly, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Lz', params%Lz, 1.0_rk )

            ! allocate memory in params structure (need 3*data_fields for velocity
            ! and 1*data_fields for diffusion coefficient)
            allocate( params%physics%u0( 3*params%number_data_fields ) )
            allocate( params%physics%nu( params%number_data_fields ) )

            ! reset values, use as default values
            params%physics%u0 = 0.0_rk
            params%physics%nu = 0.0_rk

            ! read velocity
            call read_param_mpi(FILE, 'Physics', 'u0', params%physics%u0, params%physics%u0 )
            ! read diffusion
            call read_param_mpi(FILE, 'Physics', 'nu', params%physics%nu, params%physics%nu )

            ! read variable names
            ! allocate names list
            allocate( params%physics%names( params%number_data_fields ) )

            params%physics%names = "---"
            ! read file
            call read_param_mpi(FILE, 'Physics', 'names', params%physics%names, params%physics%names )

        case('3D_navier_stokes')

            ! error case: try to solve navier stokes equation with less or more than 5 datafields
            if ( params%number_data_fields /= 5) then
                write(*,'(80("_"))')
                write(*,'("ERROR: try to solve navier stokes equation with", i3, " datafield(s)")') params%number_data_fields
                stop
            end if

            ! domain size
            call read_param_mpi(FILE, 'Physics', 'Lx', params%Lx, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Ly', params%Ly, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Lz', params%Lz, 1.0_rk )

            ! physics parameter
            ! read adiabatic coefficient
            call read_param_mpi(FILE, 'Physics', 'gamma_', params%physics_ns%gamma_, 0.0_rk )
            ! read specific gas constant
            call read_param_mpi(FILE, 'Physics', 'Rs', params%physics_ns%Rs, 0.0_rk )
            ! calculate isochoric heat capacity
            params%physics_ns%Cv = params%physics_ns%Rs/(params%physics_ns%gamma_-1.0_rk)
            ! calculate isobaric heat capacity
            params%physics_ns%Cp = params%physics_ns%Rs*params%physics_ns%gamma_
            ! read prandtl number
            call read_param_mpi(FILE, 'Physics', 'Pr', params%physics_ns%Pr, 0.0_rk )
            ! read dynamic viscosity
            call read_param_mpi(FILE, 'Physics', 'mu0', params%physics_ns%mu0, 0.0_rk )
            ! read switch to turn on|off dissipation
            call read_param_mpi(FILE, 'Blocks', 'dissipation', params%physics_ns%dissipation , .true. )

            ! read variable names
            ! allocate names list
            allocate( params%physics_ns%names( params%number_data_fields ) )
            params%physics_ns%names = "---"
            ! read file
            call read_param_mpi(FILE, 'Physics', 'names_ns', params%physics_ns%names, params%physics_ns%names )

        case('2D_advection')
            ! domain size
            call read_param_mpi(FILE, 'Physics', 'Lx', params%Lx, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Ly', params%Ly, 1.0_rk )
            ! set third dimension to zero
            params%Lz = 0.0_rk

            ! use convection velocity for time step calculation
            allocate( params%physics%u0( 2*params%number_data_fields ) )


            ! reset values, use as default values
            params%physics%u0 = 0.0_rk
            ! read velocity
            call read_param_mpi(FILE, 'Physics', 'u0', params%physics%u0, params%physics%u0 )

            ! read variable names
            ! allocate names list
            ! use convection-diffusion physics struct !
            allocate( params%physics%names( params%number_data_fields ) )


            params%physics%names = "---"
            ! read file
            call read_param_mpi(FILE, 'Physics', 'names', params%physics%names, params%physics%names )

        case('2D_acm')
            ! domain size
            call read_param_mpi(FILE, 'Physics', 'Lx', params%Lx, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Ly', params%Ly, 1.0_rk )
            ! set third dimension to zero
            params%Lz = 0.0_rk
            if ( params%number_data_fields /= 3) then
                write(*,'(80("_"))')
                write(*,'("ERROR: try to solve acm with", i3, " datafield(s)")') params%number_data_fields
                stop
            end if
            allocate( params%physics_acm%names( params%number_data_fields ) )
            params%physics_acm%names = "---"
            ! read file
            call read_param_mpi(FILE, 'Physics', 'names_acm', params%physics_acm%names, params%physics_acm%names )
            ! speed of sound for acm
            call read_param_mpi(FILE, 'Physics', 'c_0', params%physics_acm%c_0, 10.0_rk)

            ! inverse of Re
            call read_param_mpi(FILE, 'Physics', 'nu', params%physics_acm%nu, 1e-1_rk)
            
            ! gamma_p
            call read_param_mpi(FILE, 'Physics', 'gamma_p', params%physics_acm%gamma_p, 1.0_rk)

        case('3D_acm')
            ! domain size
            call read_param_mpi(FILE, 'Physics', 'Lx', params%Lx, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Ly', params%Ly, 1.0_rk )
            call read_param_mpi(FILE, 'Physics', 'Lz', params%Lz, 1.0_rk )
            if ( params%number_data_fields /= 4) then
                write(*,'(80("_"))')
                write(*,'("ERROR: try to solve acm with", i3, " datafield(s)")') params%number_data_fields
                stop
            end if
            allocate( params%physics_acm%names( params%number_data_fields ) )
            params%physics_acm%names = "---"
            ! read file
            call read_param_mpi(FILE, 'Physics', 'names_acm', params%physics_acm%names, params%physics_acm%names )
            ! speed of sound for acm
            call read_param_mpi(FILE, 'Physics', 'c_0', params%physics_acm%c_0, 10.0_rk)

            ! inverse of Re
            call read_param_mpi(FILE, 'Physics', 'nu', params%physics_acm%nu, 1e-1_rk)
            
            ! gamma_p
            call read_param_mpi(FILE, 'Physics', 'gamma_p', params%physics_acm%gamma_p, 1.0_rk)

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    !**************************************************************************
    ! read INITIAL CONDITION parameters

    ! initial condition
    call read_param_mpi(FILE, 'Physics', 'initial_cond', params%initial_cond, "---" )

    ! width of initial condition (e.g. Gauss-Blob, depends on Lx and Ly)
    call read_param_mpi(FILE, 'Physics', 'inicond_width', params%inicond_width, 1e-2_rk)
   
    if (params%initial_cond == 'read_from_files') then
        ! read variable names
        allocate( params%input_files( params%number_data_fields ) )

        params%input_files = "---"
        call read_param_mpi(FILE, 'Physics', 'input_files', params%input_files, params%input_files)
    end if

    !***************************************************************************
    ! read VOLUME PENALIZATION METHOD parameters

    ! penalization flag
    call read_param_mpi(FILE, 'VPM', 'penalization', params%penalization, .false.)
    ! penalization factor
    call read_param_mpi(FILE, 'VPM', 'eps_penal', params%eps_penal, 1e-4_rk)
    ! smooth mask for penalization term
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', params%smooth_mask, .true.)

    !***************************************************************************
    ! read DISCRETIZATION parameters
    !
    ! discretization order
    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params%order_discretization, "---" )
    ! order of predictor for refinement
    call read_param_mpi(FILE, 'Discretization', 'order_predictor', params%order_predictor, "---" )
    ! boundary condition
    call read_param_mpi(FILE, 'Discretization', 'boundary_cond', params%boundary_cond, "---" )

    ! filter type
    call read_param_mpi(FILE, 'Discretization', 'filter_type', params%filter_type, "no-filter" )
    ! filter frequency
    call read_param_mpi(FILE, 'Discretization', 'filter_freq', params%filter_freq, -1 )
    ! bogey shock detector threshold
    call read_param_mpi(FILE, 'Discretization', 'r_th', params%r_th, 1e-3_rk )


    !***************************************************************************
    ! read MPI parameters
    !
    ! data exchange method
    call read_param_mpi(FILE, 'MPI', 'mpi_data_exchange', params%mpi_data_exchange, "---" )

    !***************************************************************************
    ! read DEBUG parameters
    !
    ! debug flag
    call read_param_mpi(FILE, 'Debug', 'debug', params%debug, .true. )
    ! unit test time_stepper flag
    call read_param_mpi(FILE, 'Debug', 'test_time_stepper', params%test_time_stepper, .false.)
    ! unit test spatial flag
    call read_param_mpi(FILE, 'Debug', 'test_spatial', params%test_spatial, .false.)
    ! unit test ghost nodes synchronization flag
    call read_param_mpi(FILE, 'Debug', 'test_ghost_nodes_synch', params%test_ghost_nodes_synch, .false.)
    ! unit test wavelet compression flag
    call read_param_mpi(FILE, 'Debug', 'test_wavelet_comp', params%test_wavelet_comp, .false.)
    ! unit test treecode flag
    call read_param_mpi(FILE, 'Debug', 'test_treecode', params%test_treecode, .false.)

    !---------------------------------------------------------------------------
    ! Automatic memory management. If specified --memory=0.3GB in the call line,
    ! wabbit will automatically select the number of blocks per rank to be allocated
    ! to consume this amount of memory. helpful on local machines.
    !---------------------------------------------------------------------------
    ! loop over all arguments until you find the string "--memory" in them (or not)
    ! this ensures compatibility with future versions, as the argument may be anywhere in the call.
    do i = 1, command_argument_count()
        call get_command_argument(i, memstring)
        ! is memory limitation used?
        if ( index(memstring,"--memory=")==1 ) then
            ! read memory from command line and convert to bytes
            read(memstring(10:len_trim(memstring)-2),* ) maxmem
            ! how much memory is reserved, in bytes?
            maxmem = maxmem * 1000.0 * 1000.0 * 1000.0 ! in bytes
            if ( params%threeD_case ) then
                d = 3
            else
                d = 2
            endif
            params%number_blocks = nint( maxmem /( 8.0 * params%number_procs*(6*params%number_data_fields+1)*(params%number_block_nodes+2*params%number_ghost_nodes)**d )  )
            if (params%rank==0) write(*,'(80("-"))')
            if (params%rank==0) write(*,'("INIT: automatic selection of blocks per rank is active!")')

            if (params%rank==0) write(*,'("INIT: we allocated ",i6," blocks per rank (total: ",i7," blocks) ")') params%number_blocks, params%number_blocks*params%number_procs
            if (params%rank==0) write(*,'("INIT: consuming total memory of",f12.4,"GB")') maxmem/1000.0/1000.0/1000.0

        endif
    end do


    ! clean up
    call clean_ini_file_mpi(FILE)


    ! check ghost nodes number
    if ( (params%number_ghost_nodes < 4) .and. (params%order_predictor == 'multiresolution_4th') ) then
        call error_msg("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%number_ghost_nodes < 2) .and. (params%order_predictor == 'multiresolution_2nd') ) then
        call error_msg("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%number_ghost_nodes < 2) .and. (params%order_discretization == 'FD_4th_central_optimized') ) then
        call error_msg("ERROR: need more ghost nodes for given derivative order")
    end if

end subroutine ini_file_to_params
