! ********************************************************************************************
!> Module for 2D/3D NSPP (Navier-Stokes Pressure Poisson) physics
! ********************************************************************************************
module module_nspp

  use mpi
  use module_insects
  use module_t_files

  use module_globals
  ! ini file parser module, used to read parameters. note: in principle, you can also
  ! just use any reader you feel comfortable with, as long as you can read the parameters
  ! from a file.
  use module_ini_files_parser_mpi
  use module_operators, only : compute_vorticity, compute_divergence, compute_derivative, compute_dissipation, compute_gradient, compute_laplacian, subtract_scalar_gradient
  use module_helpers
  use module_timing
  use module_geometry

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: READ_PARAMETERS_NSPP, PREPARE_SAVE_DATA_NSPP, RHS_NSPP, GET_DT_BLOCK_NSPP, &
  INICOND_NSPP, BOUNDCOND_NSPP, FIELD_NAMES_NSPP, STATISTICS_NSPP, TIME_STATISTICS_NSPP, create_mask_2D_NSPP, &
  create_mask_3D_NSPP, geometry_indicator_nspp, PREPARE_THRESHOLDFIELD_NSPP, &
  INITIALIZE_ASCII_FILES_NSPP
  !**********************************************************************************************

  ! for 2d wing section optimization
  type :: wingsection
      logical :: initialized = .false.
      real(kind=rk), allocatable :: ai_x0(:), bi_x0(:)
      real(kind=rk), allocatable :: ai_y0(:), bi_y0(:)
      real(kind=rk), allocatable :: ai_alpha(:), bi_alpha(:)
      real(kind=rk) :: a0_x0, section_thickness
      real(kind=rk) :: a0_y0
      real(kind=rk) :: a0_alpha
      real(kind=rk) :: time
      integer(kind=ik) :: nfft_x0, nfft_y0, nfft_alpha
      character(len=cshort) :: kinematics_type
  end type

  type(wingsection) :: wingsections(2)

! how many different parts of the mask can be distinguished at max
  integer(kind=ik) :: ncolors=8

  ! user defined data structure for time independent parameters, settings, constants
  ! and the like. only visible here.
  type :: type_params_nspp
    real(kind=rk) :: CFL, T_end, CFL_eta, CFL_nu=0.094
    real(kind=rk) :: C_eta, beta, C_eta_temp, C_eta_start, C_eta_ring, penalization_startup_tau, penalization_startup_time
    integer(kind=ik) :: penalization_startup_colors
    logical :: use_free_flight_solver = .false., soft_penalization_startup=.false.
    ! this is the force and moment that is applied on the insect from the fluid, it will be computed during RHS computations
    ! as it is computed by the physics module, it was designed to be a member of the physics module at first
    ! ToDo: the insects should contain the information themselves
    real(kind=rk), allocatable :: force_insect_g(:, :), moment_insect_g(:, :)
    ! nu
    real(kind=rk) :: nu
    real(kind=rk) :: dx_min = -1.0_rk
    ! Forces for the different colors
    ! These are computed and used only in statistics output
    real(kind=rk), allocatable :: force_color(:,:), moment_color(:,:)
    logical :: penalization, compute_flow=.true.
    ! sponge term:
    logical :: use_sponge = .false.
    real(kind=rk) :: C_sponge, L_sponge, p_sponge=20.0, C_smooth=1.5_rk
    character(len=cshort) :: eps_norm
    logical :: symmetry_BC(1:3) = .false., periodic_BC(1:3) = .true.

    ! linear forcing
    logical :: HIT_linear_forcing = .false., skew_symmetry=.false.
    real(kind=rk) :: HIT_energy = 1.0_rk
    real(kind=rk) :: HIT_gain = 100.0_rk

    ! channel flow
    logical :: use_channel_forcing = .false.
    real(kind=rk), allocatable :: mask_volume(:)
    real(kind=rk) :: meanflow_channel(1:3) = 0.0_rk

    logical :: use_passive_scalar = .false.
    integer(kind=ik) :: N_scalars = 0
    real(kind=rk), allocatable :: schmidt_numbers(:), x0source(:), y0source(:), &
    z0source(:), scalar_Ceta(:), widthsource(:)
    character(len=cshort), allocatable :: scalar_inicond(:), scalar_source_type(:)
    ! when computing passive scalars, we require derivatives of the mask function, which
    ! is not too difficult on paper. however, in wabbit, ghost node syncing is not a physics
    ! module task so the NSPP module cannot do it. Note it has to be done only if scalars are used.
    ! Its a waste of resources otherwise. Hence, we have the flag to set masks on ghost nodes as well
    ! to set the mask on all points of a block (incl ghost nodes)
    logical :: set_mask_on_ghost_nodes = .false.
    logical :: absorbing_sponge = .true.

    logical :: time_statistics = .false.
    integer(kind=ik) :: N_time_statistics = 0
    character(len=cshort), allocatable :: time_statistics_names(:)
    real(kind=rk), allocatable :: time_statistics_mean(:), time_statistics_maxabs(:)

    logical :: read_from_files = .false.

    integer(kind=ik) :: dim, N_fields_saved
    real(kind=rk), dimension(3) :: domain_size=0.0_rk
    character(len=clong) :: inicond="", discretization=""

    ! VPM section
    real(kind=rk) :: x_cntr(1:3), u_cntr(1:3), R_cyl, length, thickness, u_mean_set(1:3), freq, h_channel
    integer(kind=ik) :: n_geometries = 1
    character(len=clong) :: geometry_legacy="", geometry_string=""
    character(len=clong), allocatable :: geometries(:), geometry_files(:)
    integer(kind=ik), allocatable :: geometry_colors(:)
    character(len=cshort) :: sponge_type=""
    character(len=cshort) :: coarsening_indicator=""
    character(len=cshort) :: wingsection_inifiles(1:2)
    character(len=cshort) :: scalar_BC_type="neumann"
    character(len=cshort), allocatable :: names(:)
    logical :: inicond_vorticity_formulation = .false.
    ! the mean flow, as required for some forcing terms. it is computed in the RHS
    real(kind=rk) :: mean_flow(1:3)=0.0_rk, mean_p=0.0_rk, umax=0.0_rk, umag=0.0_rk
    real(kind=rk) :: start_time = 0.0_rk
    ! integral quantities like kinetic energy and enstrophy
    real(kind=rk) :: e_kin=0.0_rk, enstrophy=0.0_rk, max_vort=0.0_rk, helicity=0.0_rk, u_residual(1:3)=0.0_rk, &
    sponge_volume=0.0_rk, dissipation=0.0_rk, scalar_removal=0.0_rk, urms(1:3), div_max, div_min, penal_power(1:3)
    ! we need to know which mpirank prints output..
    integer(kind=ik) :: mpirank, mpisize
    !
    integer(kind=ik) :: Jmax, n_ghosts = -9999999
    integer(kind=ik), dimension(3) :: Bs = -1

    ! stuff for lamballais cylinder
    real(kind=rk) :: R0, R1, R2
    character(len=clong) :: file_usx, file_usy, file_usp
    character(len=cshort), allocatable :: smoothing_type(:)
    integer(kind=ik), allocatable :: smoothing_type_int(:)
    real(kind=rk), allocatable :: smoothing_width(:), smoothing_safety(:)
    real(kind=rk), allocatable :: u_lamballais(:,:,:)

    logical :: initialized = .false.
  end type type_params_nspp

  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.
  ! HACK: made them public for FSI time stepper (18 Feb 2021)
  type(type_params_nspp), public, save :: params_nspp

contains

#include "rhs_NSPP.f90"
#include "create_mask.f90"
#include "inicond_NSPP.f90"
#include "boundcond_NSPP.f90"
#include "sponge.f90"
#include "save_data_NSPP.f90"
#include "statistics_NSPP.f90"
#include "time_statistics_NSPP.f90"
#include "2D_wingsection.f90"

  !-----------------------------------------------------------------------------
  ! main level wrapper routine to read parameters in the physics module. It reads
  ! from the same ini file as wabbit, and it reads all it has to know. note in physics modules
  ! the parameter struct for wabbit is not available.
  subroutine READ_PARAMETERS_NSPP( filename, N_mask_components, g )
    implicit none

    character(len=*), intent(in) :: filename
    integer(kind=ik) :: mpicode, nx_max, n_entries
    real(kind=rk) :: dx_min, dt_min_vpm, dt_min_nu
    character(len=cshort) :: Bs_str, Bs_conc
    character(len=400) :: input_files
    character(len=12) :: timestamp
    character(:), allocatable :: Bs_short
    real(kind=rk), dimension(3) :: ddx
    integer(kind=ik), intent(out) :: N_mask_components
    integer(kind=ik), intent(in) :: g
    integer(kind=ik) :: num_lines
    real(kind=rk), allocatable :: buffer_array(:,:)

    type(inifile) :: FILE
    integer :: Neqn, i, insect_id

    N_mask_components = 0
    ! WABBIT decides how many ghost nodes we have (because the versions >=2024 determine G 
    ! automatically depending on the wavelet). Just store the number.
    params_nspp%n_ghosts = g

    ! we still need to know about mpirank and mpisize, occasionally
    call MPI_COMM_SIZE (WABBIT_COMM, params_nspp%mpisize, mpicode)
    call MPI_COMM_RANK (WABBIT_COMM, params_nspp%mpirank, mpicode)

    if (params_nspp%mpirank==0) then
        write(*,'(80("~"))')
        write(*,'(A)') " _   _  _____ _____  ____  "
        write(*,'(A)') "| \ | |/ ____|  __ \|  _ \ "
        write(*,'(A)') "|  \| | (___ | |__) | |_) |"
        write(*,'(A)') "| . ` |\___ \|  ___/|  __/ "
        write(*,'(A)') "| |\  |____) | |    | |    "
        write(*,'(A)') "|_| \_|_____/|_|    |_|    "
        write(*,'(80("~"))')
        write(*,'(A)') "Initializing Navier-Stokes Pressure-Poisson module!"
        write(*,'(80("~"))')
    endif

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    call read_param_mpi(FILE, 'Domain', 'dim', params_nspp%dim, 2 )
    call read_param_mpi(FILE, 'Domain', 'domain_size', params_nspp%domain_size(1:params_nspp%dim), (/ 1.0_rk, 1.0_rk, 1.0_rk /) )
    params_nspp%periodic_BC = .true.
    call read_param_mpi(FILE, 'Domain', 'periodic_BC', params_nspp%periodic_BC, params_nspp%periodic_BC )

    params_nspp%symmetry_BC = .not. params_nspp%periodic_BC
    call read_param_mpi(FILE, 'Domain', 'symmetry_BC', params_nspp%symmetry_BC, params_nspp%symmetry_BC )

    ! --- saving ----
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_nspp%N_fields_saved, 3 )
    if (allocated(params_nspp%names)) deallocate(params_nspp%names)
    allocate( params_nspp%names(1:params_nspp%N_fields_saved) )
    call read_param_mpi(FILE, 'Saving', 'field_names', params_nspp%names, (/"ux","uy","p "/) )


    ! viscosity
    call read_param_mpi(FILE, 'NSPP', 'nu', params_nspp%nu, 1e-1_rk)
    call read_param_mpi(FILE, 'NSPP', 'u_mean_set', params_nspp%u_mean_set, (/1.0_rk, 0.0_rk, 0.0_rk/) )
    call read_param_mpi(FILE, 'NSPP', 'beta', params_nspp%beta, 0.05_rk )
    call read_param_mpi(FILE, 'NSPP', 'compute_flow', params_nspp%compute_flow, .true. )
    call read_param_mpi(FILE, 'NSPP', 'HIT_linear_forcing', params_nspp%HIT_linear_forcing, .false. )
    call read_param_mpi(FILE, 'NSPP', 'HIT_energy', params_nspp%HIT_energy, 1.0_rk )
    call read_param_mpi(FILE, 'NSPP', 'HIT_gain', params_nspp%HIT_gain, 100.0_rk )
    call read_param_mpi(FILE, 'NSPP', 'skew_symmetry', params_nspp%skew_symmetry, .false. )

    ! channel flow
    call read_param_mpi(FILE, 'NSPP', 'use_channel_forcing', params_nspp%use_channel_forcing, .false. )

    ! initial condition
    call read_param_mpi(FILE, 'NSPP', 'inicond', params_nspp%inicond, "meanflow")
    call read_param_mpi(FILE, 'Physics', 'inicond_vorticity_formulation', params_nspp%inicond_vorticity_formulation, .false.)
    ! the free flight FSI solver needs to know if it resumes a backup or not
    call read_param_mpi(FILE, 'Physics', 'read_from_files', params_nspp%read_from_files, .false.)
    ! free flight also requires the time at which we resume (the structure of wabbit main does no allow to pass it to this routine...)
    if (params_nspp%read_from_files) then
        ! read in all files as one string (so no check for file existence), then hack-extract the timestamp, which is used for insect_init
        call read_param_mpi(FILE, 'Physics', 'input_files', input_files, "")
        timestamp = input_files( scan(input_files,'_', back=.true.)+1:scan(input_files,'.h5', back=.true.)-3)
        read(timestamp,*) params_nspp%start_time
        ! note this requires to have timestamp in the filename (so we cannot rename files...)
        params_nspp%start_time = params_nspp%start_time * 1.0e-6
    endif


    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_nspp%discretization, "FD_4th_central_optimized")

    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params_nspp%coarsening_indicator, "threshold-state-vector")
    call read_param_mpi(FILE, 'Blocks', 'eps_norm', params_nspp%eps_norm, "Linfty")

    ! penalization:
    call read_param_mpi(FILE, 'VPM', 'penalization', params_nspp%penalization, .true.)
    ! store the same value in two variables
    call read_param_mpi(FILE, 'VPM', 'C_eta', params_nspp%C_eta, 1.0_rk)
    params_nspp%C_eta_temp = params_nspp%C_eta  ! C_eta is sometimes modified, so here we store the actual constant value
    call read_param_mpi(FILE, 'VPM', 'soft_penalization_startup', params_nspp%soft_penalization_startup, .false.)  ! if soft_penalization is used at all
    if (params_nspp%soft_penalization_startup) then
        call read_param_mpi(FILE, 'VPM', 'penalization_startup_tau', params_nspp%penalization_startup_tau, 0.20_rk)  ! how long the penalization startup is happening
        call read_param_mpi(FILE, 'VPM', 'penalization_startup_time', params_nspp%penalization_startup_time, 0.0_rk)  ! when the penalization startup starts
        call read_param_mpi(FILE, 'VPM', 'C_eta_start', params_nspp%C_eta_start, 1.0_rk)  ! the c_eta with which we start with soft-approach
        call read_param_mpi(FILE, 'VPM', 'penalization_startup_colors', params_nspp%penalization_startup_colors, 1)  ! all colors >= this get the soft startup approach applied
    else
        params_nspp%penalization_startup_tau = 0.0_rk
        params_nspp%penalization_startup_time = 0.0_rk
        params_nspp%C_eta_start = params_nspp%C_eta_temp
        params_nspp%penalization_startup_colors = 1
    endif
    call read_param_mpi(FILE, 'VPM', 'C_eta_ring', params_nspp%C_eta_ring, params_nspp%C_eta)  ! C_eta for ring in lamballais test-case

    ! geometry read_in: first we check if only one geometry should be set. This is in parity with the old ini-files
    ! if not, we read in an array of geometries and can sert in several geometry files or strings
    call read_param_mpi(FILE, 'VPM', 'geometry', params_nspp%geometry_legacy, "")
    if (params_nspp%geometry_legacy /= "") then
        params_nspp%n_geometries = 1
        allocate(params_nspp%geometries(1))
        allocate(params_nspp%geometry_files(1))
        allocate(params_nspp%geometry_colors(1))
        params_nspp%geometries(1) = params_nspp%geometry_legacy
        params_nspp%geometry_colors(:) = 1
    else
        call read_param_mpi(FILE, 'VPM', 'n_geometries', params_nspp%n_geometries, 1)
        allocate(params_nspp%geometries(params_nspp%n_geometries))
        allocate(params_nspp%geometry_files(params_nspp%n_geometries))
        allocate(params_nspp%geometry_colors(params_nspp%n_geometries))
        params_nspp%geometries(:) = "none"
        call read_param_mpi(FILE, 'VPM', 'geometries', params_nspp%geometries, defaultvalue=params_nspp%geometries )
    endif
    params_nspp%geometry_files = ""
    params_nspp%geometry_string = ""
    call read_param_mpi(FILE, 'VPM', 'geometry_files', params_nspp%geometry_files, defaultvalue=params_nspp%geometry_files )
    call read_param_mpi(FILE, 'VPM', 'geometry_string', params_nspp%geometry_string, defaultvalue=params_nspp%geometry_string )
    ! make colors to be steadily increasing from 0
    do i=1,params_nspp%n_geometries
        params_nspp%geometry_colors(i) = i
    enddo
    call read_param_mpi(FILE, 'VPM', 'geometry_colors', params_nspp%geometry_colors, defaultvalue=params_nspp%geometry_colors )
    
    ! fixed geometry parameters. They are read only once and are the same for all geometries
    call read_param_mpi(FILE, 'VPM', 'x_cntr', params_nspp%x_cntr, (/0.5*params_nspp%domain_size(1), 0.5*params_nspp%domain_size(2), 0.5*params_nspp%domain_size(3)/)  )
    call read_param_mpi(FILE, 'VPM', 'R_cyl', params_nspp%R_cyl, 0.5_rk )
    call read_param_mpi(FILE, 'VPM', 'length', params_nspp%length, 1.0_rk )
    call read_param_mpi(FILE, 'VPM', 'thickness', params_nspp%thickness, 1.0_rk )
    call read_param_mpi(FILE, 'VPM', 'freq', params_nspp%freq, 1.0_rk )
    call read_param_mpi(FILE, 'VPM', 'C_smooth', params_nspp%C_smooth, 1.5_rk )

    ! stuff for lamballais cylinder
    call read_param_mpi(FILE, 'VPM', 'R0', params_nspp%R0, 0.50_rk)
    call read_param_mpi(FILE, 'VPM', 'R1', params_nspp%R1, 1.75_rk)
    call read_param_mpi(FILE, 'VPM', 'R2', params_nspp%R2, 2.25_rk)
    call read_param_mpi(FILE, 'VPM', 'file_usx', params_nspp%file_usx, 'none')
    call read_param_mpi(FILE, 'VPM', 'file_usy', params_nspp%file_usy, 'none')
    call read_param_mpi(FILE, 'VPM', 'file_usp', params_nspp%file_usp, 'none')
    allocate(params_nspp%smoothing_type(1:params_nspp%n_geometries), params_nspp%smoothing_type_int(1:params_nspp%n_geometries), params_nspp%smoothing_width(1:params_nspp%n_geometries), params_nspp%smoothing_safety(1:params_nspp%n_geometries) )
    params_nspp%smoothing_type(:) = "cos"
    call read_param_mpi(FILE, 'VPM', 'smoothing_type', params_nspp%smoothing_type, params_nspp%smoothing_type )
    do i=1,params_nspp%n_geometries
        select case(params_nspp%smoothing_type(i))
            case("cos", "cosine")
                params_nspp%smoothing_type_int(i) = STEP_METHOD_COSINE
            case("hester")
                params_nspp%smoothing_type_int(i) = STEP_METHOD_HESTER
            case("dis", "disc", "discontinuous")
                params_nspp%smoothing_type_int(i) = STEP_METHOD_DISC
            case default
                call abort(62371118, "unknown smoothing type: "//params_nspp%smoothing_type(i))
        end select
    enddo

    ! stuff for channel flow
    call read_param_mpi(FILE, 'nspp-new', 'use_channel_forcing', params_nspp%use_channel_forcing, .false. )
    call read_param_mpi(FILE, 'nspp-new', 'u_mean_set', params_nspp%u_mean_set, (/1.0_rk, 0.0_rk, 0.0_rk/) )
    call read_param_mpi(FILE, 'VPM', 'h_channel', params_nspp%h_channel, 0.25_rk)

    do i=1,params_nspp%n_geometries
        if (params_nspp%geometries(i) == "2D-wingsection" .or. params_nspp%geometries(i) == "two-moving-cylinders") then
            call read_param_mpi(FILE, 'VPM', 'wingsection_inifiles', params_nspp%wingsection_inifiles, (/"", ""/))
        endif
    enddo

    call read_param_mpi(FILE, 'Sponge', 'use_sponge', params_nspp%use_sponge, .false. )
    call read_param_mpi(FILE, 'Sponge', 'L_sponge', params_nspp%L_sponge, 0.0_rk )
    call read_param_mpi(FILE, 'Sponge', 'C_sponge', params_nspp%C_sponge, 1.0e-2_rk )
    call read_param_mpi(FILE, 'Sponge', 'sponge_type', params_nspp%sponge_type, "rect" )
    call read_param_mpi(FILE, 'Sponge', 'p_sponge', params_nspp%p_sponge, 20.0_rk )

    call read_param_mpi(FILE, 'Time', 'CFL', params_nspp%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'CFL_eta', params_nspp%CFL_eta, 0.99_rk   )

    ! default value for CFL_nu (diffusion time step restriction) depends on the scheme and the dimension.
    ! The defaults are valid only for RK4 (this is the factor 2.79). Note this condition is relevant only for
    ! small reynolds numbers, as the diffusion is not important for the time step at higher Re.
    if (params_nspp%discretization(4:4) == "2") then
        call read_param_mpi(FILE, 'Time', 'CFL_nu', params_nspp%CFL_nu, 0.95_rk*2.79_rk/(4.000_rk*dble(params_nspp%dim)) )

    elseif (params_nspp%discretization(4:4) == "4") then
        call read_param_mpi(FILE, 'Time', 'CFL_nu', params_nspp%CFL_nu, 0.95_rk*2.79_rk/(5.333_rk*dble(params_nspp%dim)) )

    elseif (params_nspp%discretization(4:4) == "6") then
        call read_param_mpi(FILE, 'Time', 'CFL_nu', params_nspp%CFL_nu, 0.95_rk*2.79_rk/(6.0444_rk*dble(params_nspp%dim)) )

    else
        call abort(62371118, "unknown NSPP discretization:"//params_nspp%discretization )
    endif

    call read_param_mpi(FILE, 'Time', 'time_max', params_nspp%T_end, 1.0_rk   )

    call read_param_mpi(FILE, 'FreeFlightSolver', 'use_free_flight_solver', params_nspp%use_free_flight_solver, .false.   )

    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params_nspp%Jmax, 1   )

    ! passive scalars
    call read_param_mpi(FILE, 'NSPP', 'use_passive_scalar', params_nspp%use_passive_scalar, .false.)
    if (params_nspp%use_passive_scalar) then
        call read_param_mpi(FILE, 'ConvectionDiffusion', 'N_scalars', params_nspp%N_scalars, 1)

        allocate( params_nspp%schmidt_numbers(1:params_nspp%N_scalars) )
        allocate( params_nspp%x0source(1:params_nspp%N_scalars) )
        allocate( params_nspp%y0source(1:params_nspp%N_scalars) )
        allocate( params_nspp%z0source(1:params_nspp%N_scalars) )
        allocate( params_nspp%widthsource(1:params_nspp%N_scalars) )
        allocate( params_nspp%scalar_Ceta(1:params_nspp%N_scalars) )

        allocate( params_nspp%scalar_inicond(1:params_nspp%N_scalars) )
        allocate( params_nspp%scalar_source_type(1:params_nspp%N_scalars) )

        params_nspp%scalar_inicond = "dummy"
        params_nspp%scalar_source_type = "dummy"

        call read_param_mpi( FILE, 'ConvectionDiffusion', 'Sc', params_nspp%schmidt_numbers )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'x0source', params_nspp%x0source )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'y0source', params_nspp%y0source )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'z0source', params_nspp%z0source )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'widthsource', params_nspp%widthsource )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'C_eta', params_nspp%scalar_Ceta)

        call read_param_mpi( FILE, 'ConvectionDiffusion', 'inicond', params_nspp%scalar_inicond, params_nspp%scalar_inicond )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'source', params_nspp%scalar_source_type, params_nspp%scalar_source_type )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'scalar_BC_type', params_nspp%scalar_BC_type, "neumann" )

        if (params_nspp%use_sponge) then
            call read_param_mpi( FILE, 'ConvectionDiffusion', 'absorbing_sponge', params_nspp%absorbing_sponge, .true. )
        else
            params_nspp%absorbing_sponge = .false.
        endif

        ! when computing passive scalars, we require derivatives of the mask function, which
        ! is not too difficult on paper. however, in wabbit, ghost node syncing is not a physics
        ! module task so the NSPP module cannot do it. Note it has to be done only if scalars are used.
        ! Its a waste of resources otherwise. Hence, we have the flag to set masks on ghost nodes as well
        ! to set the mask on all points of a block (incl ghost nodes)
        params_nspp%set_mask_on_ghost_nodes = .true.
    else
        params_nspp%set_mask_on_ghost_nodes = .false.
        params_nspp%N_scalars = 0
    endif

    ! time statistics (averaging or similar)
    call read_param_mpi(FILE, 'Time-Statistics', 'time_statistics', params_nspp%time_statistics, .false.)
    if (params_nspp%time_statistics) then
        call read_param_mpi(FILE, 'Time-Statistics', 'N_time_statistics', params_nspp%N_time_statistics, 1)
        allocate( params_nspp%time_statistics_names(1:params_nspp%N_time_statistics) )
        allocate( params_nspp%time_statistics_mean(1:params_nspp%N_time_statistics) )
        allocate( params_nspp%time_statistics_maxabs(1:params_nspp%N_time_statistics) )
        call read_param_mpi(FILE, 'Time-Statistics', 'time_statistics_names', params_nspp%time_statistics_names, (/"none"/))
    endif

    ! set defaults
    if (params_nspp%dim==3) then
      params_nspp%Bs=(/17,17,17/)
    else
      params_nspp%Bs=(/17,17,1/)
    endif
    params_nspp%Bs = read_bs(FILE,'Blocks', 'number_block_nodes', params_nspp%Bs, params_nspp%dim, params_nspp%mpirank)

    call clean_ini_file_mpi( FILE )

    if (Neqn /= params_nspp%dim + 1 + params_nspp%N_scalars) then
        ! call abort(220819, "the state vector length is not appropriate. number_equation must be DIM+1+N_scalars")
    endif

    ! uniqueGrid modification
    ddx(1:params_nspp%dim) = 2.0_rk**(-params_nspp%Jmax) * (params_nspp%domain_size(1:params_nspp%dim) / real(params_nspp%Bs(1:params_nspp%dim), kind=rk))

    dx_min = minval( ddx(1:params_nspp%dim) )
    ! uniqueGrid modification
    nx_max = maxval( (params_nspp%Bs) * 2**(params_nspp%Jmax) )

    ! compute some settings for VPM
    do i=1,params_nspp%n_geometries
        select case(params_nspp%smoothing_type(i))
            case("cos", "cosine")
                params_nspp%smoothing_width(i) = dx_min * params_nspp%C_smooth
                params_nspp%smoothing_safety(i) = 1.0_rk * params_nspp%smoothing_width(i)
            case ("hester")
                params_nspp%smoothing_width(i) = sqrt(params_nspp%nu * params_nspp%C_eta)
                params_nspp%smoothing_safety(i) = max(5.0_rk * params_nspp%smoothing_width(i), 2*maxval( ddx(1:params_nspp%dim) ))
            case("discontinuous", "dis")
                params_nspp%smoothing_width(i) = dx_min * params_nspp%C_smooth
                params_nspp%smoothing_safety(i) = 3.0_rk * params_nspp%smoothing_width(i)
            case default
                call abort(260602, "Never heard of the smoothing type "//trim(params_nspp%smoothing_type(i)))
        end select
    enddo

    ! print some time-step related information
    if (params_nspp%nu > 0.0_rk) then
        dt_min_nu = params_nspp%CFL_nu*dx_min**2/params_nspp%nu
    else
        dt_min_nu = 0.0_rk
    endif
    if (params_nspp%penalization .and. params_nspp%C_eta > 0.0_rk) then
        dt_min_vpm = params_nspp%CFL_eta*params_nspp%C_eta
    else
        dt_min_vpm = 0.0_rk
    endif

    ! nice to have this elsewhere in the NSPP module:
    params_nspp%dx_min = dx_min

    ! at most, we need 6 components: mask, usx, usy, usz, color, sponge
    ! in 2d, less arrays could be used, but its easier to just go ahead and use all of them.
    N_mask_components = 6

    if (params_nspp%mpirank==0) then
      write(*,'(80("<"))')
      write(*,'(A)') "Some information:"
      write(*,'("   CFL=",g12.4, "CFL_eta=",g12.4, "CFL_nu=",g12.4)') params_nspp%CFL, params_nspp%CFL_eta, params_nspp%CFL_nu
      write(*,'("   dx_min=",g12.4)') dx_min
      write(*,'("   dt(CFL_nu,nu,dx_min)=",g12.4)') dt_min_nu
      write(*,'("   dt(CFL_vpm,C_eta,dx_min)=",g12.4)') dt_min_vpm
      write(*,'("   if all blocks were at Jmax, the resolution would be nx=",i5)') nx_max
      if (params_nspp%penalization) then
          write(*,'("   C_eta=",es12.4," K_eta=",es12.4)') params_nspp%C_eta, sqrt(params_nspp%C_eta*params_nspp%nu)/dx_min
      endif
      write(*,'("N_mask_components=",i1, " N_scalars=",i1, " N_time_statistics=",i1)') N_mask_components, params_nspp%N_scalars, params_nspp%N_time_statistics
      write(*,'(80("<"))')
    endif

    ! before we init the insects, we have to count how many there are
    do i=1,params_nspp%n_geometries
        if (strings_are_similar(params_nspp%geometries(i), "insect") .or. strings_are_similar(params_nspp%geometries(i), "active-grid") .or. &
            strings_are_similar(params_nspp%geometries(i), "cylinder-free") .or. strings_are_similar(params_nspp%geometries(i), "sphere-free") .or. &
            strings_are_similar(params_nspp%geometries(i), "plate-free")) then
            n_insects = n_insects + 1
        endif
    enddo
    ! now we initialize the insects array
    call insects_array_init(n_insects)
    allocate( params_nspp%force_insect_g(1:3, 0:n_insects))
    allocate( params_nspp%moment_insect_g(1:3, 0:n_insects))

    ! Loop over all geometries and do geometry-specific initialization
    ncolors = 1
    do i=1,params_nspp%n_geometries

        ncolors = max( ncolors, params_nspp%geometry_colors(i) )

        ! if used, setup insect. Note active grid is part of the insects: they require the same init module
        !
        ! NOTE: there are several testing geometries used to test the free-flight solver: cylinder-free, sphere-free and plate-free (2D)
        if (strings_are_similar(params_nspp%geometries(i), "insect") .or. strings_are_similar(params_nspp%geometries(i), "active-grid") .or. &
            strings_are_similar(params_nspp%geometries(i), "cylinder-free") .or. strings_are_similar(params_nspp%geometries(i), "sphere-free") .or. &
            strings_are_similar(params_nspp%geometries(i), "plate-free")) then
            ! when computing passive scalars, we require derivatives of the mask function, which
            ! is not too difficult on paper. however, in wabbit, ghost node syncing is not a physics
            ! module task so the NSPP module cannot do it. Note it has to be done only if scalars are used.
            ! Its a waste of resources otherwise. Hence, we have the flag to set masks on ghost nodes as well
            ! to set the mask on all points of a block (incl ghost nodes)
            call get_insect_id( i, insect_id )
            if (params_nspp%set_mask_on_ghost_nodes) then
                call insect_init( params_nspp%start_time, filename, insect_id, params_nspp%read_from_files, "", params_nspp%domain_size, &
                params_nspp%nu, params_nspp%C_eta, dx_min, params_nspp%smoothing_type(i), N_ghost_nodes=0, colors_default=(/params_nspp%n_geometries+1+5*(insect_id-1),params_nspp%n_geometries+2+5*(insect_id-1),params_nspp%n_geometries+3+5*(insect_id-1),params_nspp%n_geometries+4+5*(insect_id-1),params_nspp%n_geometries+5+5*(insect_id-1), params_nspp%geometry_colors(i)/))
            else
                call insect_init( params_nspp%start_time, filename, insect_id, params_nspp%read_from_files, "", params_nspp%domain_size, &
                params_nspp%nu, params_nspp%C_eta, dx_min, params_nspp%smoothing_type(i), N_ghost_nodes=g, colors_default=(/params_nspp%n_geometries+1+5*(insect_id-1),params_nspp%n_geometries+2+5*(insect_id-1),params_nspp%n_geometries+3+5*(insect_id-1),params_nspp%n_geometries+4+5*(insect_id-1),params_nspp%n_geometries+5+5*(insect_id-1), params_nspp%geometry_colors(i)/))
            endif

            ! compute maximum color
            ncolors = max( ncolors, int(maxval((/ insects(insect_id)%color_body, insects(insect_id)%color_l, insects(insect_id)%color_r, insects(insect_id)%color_l2, insects(insect_id)%color_r2/)), kind=ik) )

            ! check if this insect is in free-flight mode but we did not activate free flight for physics
            if (.not. params_nspp%use_free_flight_solver .and. strings_are_similar(insects(insect_id)%BodyMotion, "free_flight")) then
                call abort(62371118, "You seem to use an insect in free-flight mode, but you did not activate the free-flight solver in the physics settings. Please check your ini-file.")
            endif
        endif

        if (params_nspp%geometries(i)=="2D-wingsection" .or. params_nspp%geometries(i)=="two-moving-cylinders") then
            call init_wingsection_from_file(params_nspp%wingsection_inifiles(1), wingsections(1), 0.0_rk)
            call init_wingsection_from_file(params_nspp%wingsection_inifiles(2), wingsections(2), 0.0_rk)
        endif


        ! read lamballais reference fields, see
        ! Gautier, R., Biau, D., Lamballais, E.: A reference solution of the flow over a circular cylinder at Re = 40 , Computers & Fluids 75, 103–111, 2013
        if ((params_nspp%geometries(i) == "lamballais").or.(params_nspp%geometries(i)=="lamballais-local")) then
            if (params_nspp%dim /= 2) call abort(1409241, "lamballais is a 2D test case")

            ! read us field
            call count_lines_in_ascii_file_mpi(params_nspp%file_usx, num_lines, 0)
            ! avoid maxcolumns restriction (read in a single long column and reshape)
            allocate(buffer_array(1:num_lines,1:1))
            call read_array_from_ascii_file_mpi(params_nspp%file_usx, buffer_array, 0)

            write(*,*) "lamballais num_lines", num_lines, nx_max


            if (num_lines /= nx_max**2) then
                call abort(2410011, "Lamballais: you seem to read the wrong field (size mismatch?!)")
            endif

            allocate(params_nspp%u_lamballais(0:nx_max-1, 0:nx_max-1, 1:3))
            ! reshape
            params_nspp%u_lamballais(:,:,1) = reshape(buffer_array, (/nx_max, nx_max/))

            call read_array_from_ascii_file_mpi(params_nspp%file_usy, buffer_array, 0)
            params_nspp%u_lamballais(:,:,2) = reshape(buffer_array, (/nx_max, nx_max/))

            call read_array_from_ascii_file_mpi(params_nspp%file_usp, buffer_array, 0)
            params_nspp%u_lamballais(:,:,3) = reshape(buffer_array, (/nx_max, nx_max/))
        endif
    enddo

    ! now initialze force arrays for colors at last, because we know how many colors we have
    allocate( params_nspp%force_color(1:3, 1:ncolors), params_nspp%moment_color(1:3, 1:ncolors), params_nspp%mask_volume(1:ncolors) )

    params_nspp%initialized = .true.
  end subroutine READ_PARAMETERS_NSPP


  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_NSPP( time, iteration, u, Bs, g, x0, dx, dt )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent(in) :: time
    integer(kind=ik), intent(in) :: iteration

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! the dt for this block is returned to the caller:
    real(kind=rk), intent(out) :: dt
    ! temporary array. note this is just one block and hence not important for overall memory consumption
    real(kind=rk) :: u_eigen, kappa, uu
    integer :: iscalar, dim, ix, iy ,iz

    if (.not. params_nspp%initialized) write(*,*) "WARNING: GET_DT_BLOCK_NSPP called but NSPP not initialized"

    dim = params_nspp%dim

    ! NSPP: only convective velocity matters (no acoustic waves, pressure is elliptic)
    ! In contrast to nspp, there are no fast modes with velocity u
    ! The pressure Poisson equation is solved implicitly and has no time-step restriction
    ! Compute maximum velocity magnitude directly
    if (params_nspp%dim == 2) then
        u_eigen = maxval(sqrt(u(g+1:Bs(1)+g,g+1:Bs(2)+g,1,1)**2 + u(g+1:Bs(1)+g,g+1:Bs(2)+g,1,2)**2))
    else
        u_eigen = maxval(sqrt(u(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,1)**2 + u(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,2)**2 + &
                        u(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,3)**2))
    endif

    ! CFL condition based on convective velocity only
    ! if the characteristic velocity is very small, avoid division by zero
    if ( u_eigen >= 1.0e-6_rk ) then
        dt = params_nspp%CFL * minval(dx(1:dim)) / u_eigen
    else
        dt = 1.0e-2_rk
    endif

    ! explicit diffusion (NOTE: factor 0.5 is valid only for RK4, other time steppers have more
    ! severe restrictions)
    if (params_nspp%nu>1.0e-13_rk) then
        dt = min(dt,  params_nspp%CFL_nu * minval(dx(1:dim))**2 / params_nspp%nu)
    endif

    ! diffusivity of scalars, if they are in use.(NOTE: factor 0.5 is valid only for RK4, other time steppers have more
    ! severe restrictions)
    if (params_nspp%use_passive_scalar)  then
        do iscalar = 1, params_nspp%N_scalars
            kappa = params_nspp%nu * params_nspp%schmidt_numbers(iscalar)
            if (kappa>1.0e-13_rk) dt = min(dt, params_nspp%CFL_nu * minval(dx(1:dim))**2 / kappa)
        enddo
    endif

    ! penalization
    if (params_nspp%penalization) dt = min( dt, params_nspp%CFL_eta*params_nspp%C_eta )

    ! sponge
    if (params_nspp%use_sponge) dt = min( dt, params_nspp%CFL_eta*params_nspp%C_sponge )

  end subroutine GET_DT_BLOCK_NSPP



  subroutine continue_periodic(x,L)
      !> position x
      real(kind=rk), intent(inout)     :: x
      !> domain length
      real(kind=rk), intent(in)     :: L

      real(kind=rk)                  :: min_dx

      if ( x>L ) then
        x=x-L
      elseif( x<0 ) then
        ! note it is actually x=L-abs(x) but since x is negative its
        x=L+x
      endif

      min_dx = 2.0_rk**(-params_nspp%Jmax) * min(params_nspp%domain_size(1)/real(params_nspp%Bs(1)-1, kind=rk) &
                                                  , params_nspp%domain_size(2)/real(params_nspp%Bs(2)-1, kind=rk))
      ! u(x=0) should be set equal to u(x=L)
      if ( abs(x-L)<min_dx*0.5_rk ) then
          x = 0.0_rk
      end if

  end subroutine continue_periodic

  !-----------------------------------------------------------------------------
  ! Adaptation is dependent on the different physics application.
  ! Every physics module can choose its own coarsening indicator.
  !-----------------------------------------------------------------------------
  subroutine PREPARE_THRESHOLDFIELD_NSPP( u, g, x0, dx, threshold_field, N_thresholding_components )
      implicit none

      ! block data, containg the state vector. In general a 4D field (3 dims+components)
      ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
      real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

      ! as you are allowed to compute the RHS only in the interior of the field
      ! you also need to know where 'interior' starts: so we pass the number of ghost points
      integer, intent(in) :: g

      ! for each block, you'll need to know where it lies in physical space. The first
      ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
      real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

      ! output. Note assumed-shape arrays
      real(kind=rk), intent(inout) :: threshold_field(1:,1:,1:,1:)

      ! WABBIT needs to know how many components are involved in thresholding,
      ! so return this number as well.
      integer(kind=ik), intent(out):: N_thresholding_components

      ! this function should not be called for the regular statevector thresholding
      call abort(27022017,"The NSPP module supports only coarsening_indicator=threshold-state-vector; !")
  end subroutine


  ! the statistics are written to ascii files (usually *.t files) with the help
  ! of module_t_files. In any case, the files have to be intialized: ideally, they
  ! are equipped with a header and resetted on the very first call, and they must not be deleted
  ! if the simuation is resumed from a backup. We therefore provide this function so that all physics
  ! modules can initialize those files.
  subroutine INITIALIZE_ASCII_FILES_NSPP( time, overwrite )
      implicit none

      ! it may happen that some source terms have an explicit time-dependency
      ! therefore the general call has to pass time
      real(kind=rk), intent (in) :: time
      logical, intent(in) :: overwrite
      logical :: is_insect, has_two_wings
      integer :: i_insect, i_color
      character(len=cshort) :: headers(1:100)  ! we can use this to create headers

      is_insect = .false.
      if (any(strings_are_similar(params_nspp%geometries(:), "insect")) .or. any(strings_are_similar(params_nspp%geometries(:), "active-grid")) .or. &
        any(strings_are_similar(params_nspp%geometries(:), "cylinder-free")) .or. any(strings_are_similar(params_nspp%geometries(:), "sphere-free")) .or. &
        any(strings_are_similar(params_nspp%geometries(:), "plate-free"))) is_insect = .true.

      has_two_wings = .false.
      do i_insect = 1, n_insects
        has_two_wings = has_two_wings .or. (insects(i_insect)%second_wing_pair)
      enddo

      headers(1) = "time"

      call init_t_file('meanflow.t', overwrite)
      headers(2) = "e_kin"
      call init_t_file('e_kin.t', overwrite, headers(1:2))
      headers(2) = "enstrophy"
      headers(3) = "max(omega)"
      call init_t_file('enstrophy.t', overwrite, headers(1:3))
      if (params_nspp%dim == 3) then
        headers(2) = "helicity"
        call init_t_file('helicity.t', overwrite, headers(1:2))
      endif
      headers(2) = "nu u laplace u"
      call init_t_file('dissipation.t', overwrite, headers(1:2))
      headers(2) = "max(div)"
      headers(3) = "min(div)"
      call init_t_file('div.t', overwrite,  headers(1:3))
      headers(2) = "max(|u|)=um"
      call init_t_file('umag.t', overwrite, headers(1:2))
      if (params_nspp%HIT_linear_forcing) then
        call init_t_file('turbulent_statistics.t', overwrite, (/"           time", "    dissipation", "         energy", "          u_RMS", &
      "    kolm_length", "      kolm_time", "  kolm_velocity", "   taylor_micro", "reynolds_taylor"/))
      endif
      call init_t_file('CFL.t', overwrite, (/&
      "           time", &
      "            CFL", &
      "         CFL_nu", &
      "        CFL_eta"/) )
      if (params_nspp%time_statistics) then
        call init_t_file('time_statistics_mean.t', overwrite)
        call init_t_file('time_statistics_maxabs.t', overwrite)
      endif
      if (params_nspp%penalization .or. params_nspp%use_sponge) then
        call init_t_file('forces.t', overwrite)

        ! dynamic initialization of force array so that it makes sense
        headers(1) = "           time"
        do i_color = 1, ncolors
            write(headers((i_color-1)*3 + 2),"(A,i0.3,A)") "c", i_color, ":force_X"
            write(headers((i_color-1)*3 + 3),"(A,i0.3,A)") "c", i_color, ":force_Y"
            write(headers((i_color-1)*3 + 4),"(A,i0.3,A)") "c", i_color, ":force_Z"
        enddo
        call init_t_file('forces_color.t', overwrite, headers(1:3*ncolors+1) )
        headers(1) = "           time"
        do i_color = 1, ncolors
            write(headers((i_color-1)*3 + 2),"(A,i0.3,A)") "c", i_color, ":moment_X"
            write(headers((i_color-1)*3 + 3),"(A,i0.3,A)") "c", i_color, ":moment_Y"
            write(headers((i_color-1)*3 + 4),"(A,i0.3,A)") "c", i_color, ":moment_Z"
        enddo
        call init_t_file('moments_color.t', overwrite, headers(1:3*ncolors+1) )

        if (is_insect) then
            call init_t_file('moments.t', overwrite)
            
            ! headers for aero power file
            headers(1) = "           time"
            do i_insect = 1, n_insects
                write(headers((i_insect-1)*2 + 2),"(A,i0.2,A)") "I", i_insect, ":apow"
                write(headers((i_insect-1)*2 + 3),"(A,i0.2,A)") "I", i_insect, ":ipow"
            enddo
            call init_t_file('aero_power.t', overwrite, headers(1:2*n_insects+1) )

            call init_t_file('forces_body.t', overwrite)
            call init_t_file('moments_body.t', overwrite)
            call init_t_file('forces_leftwing.t', overwrite)
            call init_t_file('moments_leftwing.t', overwrite)
            call init_t_file('forces_rightwing.t', overwrite)
            call init_t_file('moments_rightwing.t', overwrite)
            
            ! headers for state vector file
            headers(1) = "           time"
            do i_insect = 1, n_insects
                write(headers((i_insect-1)*23 + 2),"(A,i0.2,A)") "I", i_insect, ":x-pos"
                write(headers((i_insect-1)*23 + 3),"(A,i0.2,A)") "I", i_insect, ":y-pos"
                write(headers((i_insect-1)*23 + 4),"(A,i0.2,A)") "I", i_insect, ":z-pos"
                write(headers((i_insect-1)*23 + 5),"(A,i0.2,A)") "I", i_insect, ":x-vel"
                write(headers((i_insect-1)*23 + 6),"(A,i0.2,A)") "I", i_insect, ":y-vel"
                write(headers((i_insect-1)*23 + 7),"(A,i0.2,A)") "I", i_insect, ":z-vel"
                write(headers((i_insect-1)*23 + 8),"(A,i0.2,A)") "I", i_insect, ":q1-body"
                write(headers((i_insect-1)*23 + 9),"(A,i0.2,A)") "I", i_insect, ":q2-body"
                write(headers((i_insect-1)*23 + 10),"(A,i0.2,A)") "I", i_insect, ":q3-body"
                write(headers((i_insect-1)*23 + 11),"(A,i0.2,A)") "I", i_insect, ":q4-body"
                write(headers((i_insect-1)*23 + 12),"(A,i0.2,A)") "I", i_insect, ":w-x-body"
                write(headers((i_insect-1)*23 + 13),"(A,i0.2,A)") "I", i_insect, ":w-y-body"
                write(headers((i_insect-1)*23 + 14),"(A,i0.2,A)") "I", i_insect, ":w-z-body"
                write(headers((i_insect-1)*23 + 15),"(A,i0.2,A)") "I", i_insect, ":q1-l"
                write(headers((i_insect-1)*23 + 16),"(A,i0.2,A)") "I", i_insect, ":q2-l"
                write(headers((i_insect-1)*23 + 17),"(A,i0.2,A)") "I", i_insect, ":q3-l"
                write(headers((i_insect-1)*23 + 18),"(A,i0.2,A)") "I", i_insect, ":q4-l"
                write(headers((i_insect-1)*23 + 19),"(A,i0.2,A)") "I", i_insect, ":w-x-l"
                write(headers((i_insect-1)*23 + 20),"(A,i0.2,A)") "I", i_insect, ":w-y-l"
                write(headers((i_insect-1)*23 + 21),"(A,i0.2,A)") "I", i_insect, ":w-z-l"
                write(headers((i_insect-1)*23 + 22),"(A,i0.2,A)") "I", i_insect, ":force-g-x"
                write(headers((i_insect-1)*23 + 23),"(A,i0.2,A)") "I", i_insect, ":force-g-y"
                write(headers((i_insect-1)*23 + 24),"(A,i0.2,A)") "I", i_insect, ":force-g-z"
            enddo
            call init_t_file('insect_state_vector.t', overwrite, headers(1:23*n_insects+1) )

            if (has_two_wings) then
                call init_t_file('forces_leftwing2.t', overwrite)
                call init_t_file('moments_leftwing2.t', overwrite)
                call init_t_file('forces_rightwing2.t', overwrite)
                call init_t_file('moments_rightwing2.t', overwrite)
            endif

            call init_insect_data(overwrite)
        endif
        do i_color = 1, ncolors
            write(headers(i_color+1),"(A,i0.3,A)") "c", i_color, ":mask_volume"
        enddo
        headers(ncolors+2) = "sponge_volume"
        call init_t_file('mask_volume.t', overwrite, headers(1:ncolors+2) )
        call init_t_file('u_residual.t', overwrite)
        call init_t_file('forces_rk.t', overwrite)
        call init_t_file('penal_power.t', overwrite, (/&
        "           time", &
        "  E_dot_f_solid"/))
    endif


  end subroutine INITIALIZE_ASCII_FILES_NSPP


  !> In geometry string, we might have [primitives-collection, insect, primitives-collection, insect].
  !! Now, when looping over all geometries, we want to now that geometry 4 corresponds to the second insect.
  !! This routine does exactly that.
  subroutine get_insect_id(i_geom, insect_id)
    implicit none
    integer, intent(in) :: i_geom  !< index of the geometry in the geometries array
    integer, intent(inout) :: insect_id  !< the insect_id corresponding to the geometry

    integer :: i, count_insects
    count_insects = 0
    do i = 1, i_geom
        if (strings_are_similar(params_nspp%geometries(i), "insect") .or. strings_are_similar(params_nspp%geometries(i), "active-grid") .or. &
            strings_are_similar(params_nspp%geometries(i), "cylinder-free") .or. strings_are_similar(params_nspp%geometries(i), "sphere-free") .or. &
            strings_are_similar(params_nspp%geometries(i), "plate-free")) then
            count_insects = count_insects + 1
        endif
    enddo
    insect_id = count_insects
  end subroutine

end module module_nspp
