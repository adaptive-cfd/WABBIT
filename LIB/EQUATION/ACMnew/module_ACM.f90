! ********************************************************************************************
!> Module for 2D/3D acm physics
! ********************************************************************************************
module module_acm

  use mpi
  use module_insects
  use module_t_files

  use module_globals
  ! ini file parser module, used to read parameters. note: in principle, you can also
  ! just use any reader you feel comfortable with, as long as you can read the parameters
  ! from a file.
  use module_ini_files_parser_mpi
  use module_operators, only : compute_vorticity, divergence
  use module_params, only : read_bs
  use module_helpers, only : startup_conditioner, smoothstep, random_data, fseries_eval
  use module_timing

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: READ_PARAMETERS_ACM, PREPARE_SAVE_DATA_ACM, RHS_ACM, GET_DT_BLOCK_ACM, &
  INICOND_ACM, BOUNDCOND_ACM, FIELD_NAMES_ACM, STATISTICS_ACM, FILTER_ACM, create_mask_2D_ACM, &
  create_mask_3D_ACM, PREPARE_THRESHOLDFIELD_ACM, &
  INITIALIZE_ASCII_FILES_ACM, WRITE_INSECT_DATA, Update_Insect_wrapper
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

  ! user defined data structure for time independent parameters, settings, constants
  ! and the like. only visible here.
  type :: type_params_acm
    real(kind=rk) :: CFL, T_end, CFL_eta, CFL_nu=0.094
    real(kind=rk) :: c_0
    real(kind=rk) :: C_eta, beta
    logical :: use_free_flight_solver = .false.
    real(kind=rk),dimension(1:3) :: force_insect_g=0.0_rk, moment_insect_g=0.0_rk
    ! nu
    real(kind=rk) :: nu, nu_p=0.0_rk
    real(kind=rk) :: dx_min = -1.0_rk
    real(kind=rk) :: x_cntr(1:3), u_cntr(1:3), R_cyl, length, thickness, u_mean_set(1:3),  &
                     urms(1:3), div_max, div_min, freq, u_vert=0.0_rk, z_vert, penal_power
    ! forces for the different colors
    real(kind=rk) :: force_color(1:3,0:6), moment_color(1:3,0:6)
    ! gamma_p
    real(kind=rk) :: gamma_p
    ! want to add forcing?
    logical :: penalization, smooth_mask=.True., compute_flow=.true.
    ! sponge term:
    logical :: use_sponge = .false.
    real(kind=rk) :: C_sponge, L_sponge, p_sponge=20.0, C_smooth=1.5_rk
    character(len=cshort) :: eps_norm
    logical :: symmetry_BC(1:3) = .false., periodic_BC(1:3) = .true.

    ! linear forcing
    logical :: HIT_linear_forcing = .false.
    real(kind=rk) :: HIT_energy = 1.0_rk
    real(kind=rk) :: HIT_gain = 100.0_rk

    logical :: use_passive_scalar = .false.
    integer(kind=ik) :: N_scalars = 0, nsave_stats = 999999
    real(kind=rk), allocatable :: schmidt_numbers(:), x0source(:), y0source(:), &
    z0source(:), scalar_Ceta(:), widthsource(:)
    character(len=cshort), allocatable :: scalar_inicond(:), scalar_source_type(:)
    ! when computing passive scalars, we require derivatives of the mask function, which
    ! is not too difficult on paper. however, in wabbit, ghost node syncing is not a physics
    ! module task so the ACM module cannot do it. Note it has to be done only if scalars are used.
    ! Its a waste of resources otherwise. Hence, we have the flag to set masks on ghost nodes as well
    ! to set the mask on all points of a block (incl ghost nodes)
    logical :: set_mask_on_ghost_nodes = .false.
    logical :: absorbing_sponge = .true.

    logical :: read_from_files = .false.

    integer(kind=ik) :: dim, N_fields_saved
    real(kind=rk), dimension(3) :: domain_size=0.0_rk
    character(len=cshort) :: inicond="", discretization="", filter_type="", geometry="cylinder", order_predictor=""
    character(len=cshort) :: sponge_type=""
    character(len=cshort) :: nonlinear_formulation="convective" ! or "divergence"
    character(len=cshort) :: coarsening_indicator=""
    character(len=cshort) :: wingsection_inifiles(1:2)
    character(len=cshort) :: scalar_BC_type="neumann"
    character(len=cshort), allocatable :: names(:)
    ! the mean flow, as required for some forcing terms. it is computed in the RHS
    real(kind=rk) :: mean_flow(1:3), mean_p, umax, umag
    real(kind=rk) :: start_time = 0.0_rk
    ! kinetic energy and enstrophy (both integrals)
    real(kind=rk) :: e_kin, enstrophy, mask_volume, u_residual(1:3), sponge_volume, dissipation, scalar_removal=0.0_rk
    ! we need to know which mpirank prints output..
    integer(kind=ik) :: mpirank, mpisize
    !
    integer(kind=ik) :: Jmax, n_ghosts = -9999999
    integer(kind=ik), dimension(3) :: Bs = -1

    ! stuff for lamballais cylinder
    real(kind=rk) :: R0, R1, R2
    character(len=clong) :: file_usx, file_usy, file_usp
    character(len=cshort) :: smoothing_type
    real(kind=rk), allocatable :: u_lamballais(:,:,:)

    logical :: initialized = .false.
  end type type_params_acm

  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.
  ! HACK: made them public for FSI time stepper (18 Feb 2021)
  type(type_params_acm), public, save :: params_acm

  ! all parameters for insects go here:
  ! HACK: made them public for FSI time stepper (18 Feb 2021)
  type(diptera), public, save :: insect

contains

#include "rhs_ACM.f90"
#include "create_mask.f90"
#include "inicond_ACM.f90"
#include "boundcond_ACM.f90"
#include "sponge.f90"
#include "save_data_ACM.f90"
#include "statistics_ACM.f90"
#include "filter_ACM.f90"
#include "2D_wingsection.f90"

! this routine is public, even though it is non-standard for all physics modules.
! it is used in "dry-run" mode, which we use to create insect mask functions without
! solving the fluid equations. This is an incredibly useful mode, and it needs to write
! the kinematics to *.t file. until a permanent solution is found, this is a HACK
subroutine WRITE_INSECT_DATA(time)
    implicit none
    real(kind=rk), intent(in) :: time

    if (.not. params_acm%initialized) write(*,*) "WARNING: WRITE_INSECT_DATA called but ACM not initialized"

    if (Insect%second_wing_pair) then
        call append_t_file( 'kinematics.t', (/time, Insect%xc_body_g, Insect%psi, Insect%beta, &
        Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
        Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
        Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
        Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w, &
        Insect%alpha_l2, Insect%phi_l2, Insect%theta_l2, &
        Insect%alpha_r2, Insect%phi_r2, Insect%theta_r2, &
        Insect%rot_rel_wing_l2_w, Insect%rot_rel_wing_r2_w, &
        Insect%rot_dt_wing_l2_w, Insect%rot_dt_wing_r2_w/) )
    else
        call append_t_file( 'kinematics.t', (/time, Insect%xc_body_g, Insect%psi, Insect%beta, &
        Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
        Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
        Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
        Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w/) )
    endif

end subroutine

  !-----------------------------------------------------------------------------
  ! main level wrapper routine to read parameters in the physics module. It reads
  ! from the same ini file as wabbit, and it reads all it has to know. note in physics modules
  ! the parameter struct for wabbit is not available.
  subroutine READ_PARAMETERS_ACM( filename, N_mask_components, g )
    implicit none

    character(len=*), intent(in) :: filename
    integer(kind=ik) :: mpicode, nx_max, n_entries
    real(kind=rk) :: dx_min, dt_min
    character(len=cshort) :: Bs_str, Bs_conc
    character(len=16834) :: input_files
    character(len=12) :: timestamp
    character(:), allocatable :: Bs_short
    real(kind=rk), dimension(3) :: ddx
    integer(kind=ik), intent(out) :: N_mask_components
    integer(kind=ik), intent(in) :: g
    integer(kind=ik) :: num_lines
    real(kind=rk), allocatable :: buffer_array(:,:)

    type(inifile) :: FILE
    integer :: Neqn, i

    N_mask_components = 0
    ! WABBIT decides how many ghost nodes we have (because the versions >=2024 determine G 
    ! automatically depending on the wavelet). Just store the number.
    params_acm%n_ghosts = g

    ! we still need to know about mpirank and mpisize, occasionally
    call MPI_COMM_SIZE (WABBIT_COMM, params_acm%mpisize, mpicode)
    call MPI_COMM_RANK (WABBIT_COMM, params_acm%mpirank, mpicode)

    if (params_acm%mpirank==0) then
        write(*,'(80("~"))')
        write(*,'(A)') "    _    ____ __  __       _       _ _   "
        write(*,'(A)') "   / \  / ___|  \/  |     (_)_ __ (_) |_ "
        write(*,'(A)') "  / _ \| |   | |\/| |_____| | '_ \| | __|"
        write(*,'(A)') " / ___ \ |___| |  | |_____| | | | | | |_ "
        write(*,'(A)') "/_/   \_\____|_|  |_|     |_|_| |_|_|\__|"
        write(*,'(80("~"))')
        write(*,'(A)') "Initializing artificial compressibility module!"
        write(*,'(80("~"))')
    endif

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    call read_param_mpi(FILE, 'Domain', 'dim', params_acm%dim, 2 )
    call read_param_mpi(FILE, 'Domain', 'domain_size', params_acm%domain_size(1:params_acm%dim), (/ 1.0_rk, 1.0_rk, 1.0_rk /) )
    params_acm%periodic_BC = .true.
    call read_param_mpi(FILE, 'Domain', 'periodic_BC', params_acm%periodic_BC, params_acm%periodic_BC )

    params_acm%symmetry_BC = .not. params_acm%periodic_BC
    call read_param_mpi(FILE, 'Domain', 'symmetry_BC', params_acm%symmetry_BC, params_acm%symmetry_BC )

    ! --- saving ----
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_acm%N_fields_saved, 3 )
    if (allocated(params_acm%names)) deallocate(params_acm%names)
    allocate( params_acm%names(1:params_acm%N_fields_saved) )
    call read_param_mpi(FILE, 'Saving', 'field_names', params_acm%names, (/"ux","uy","p "/) )


    ! speed of sound for acm
    call read_param_mpi(FILE, 'ACM-new', 'c_0', params_acm%c_0, 10.0_rk)
    ! viscosity
    call read_param_mpi(FILE, 'ACM-new', 'nu', params_acm%nu, 1e-1_rk)
    call read_param_mpi(FILE, 'ACM-new', 'nu_p', params_acm%nu_p, 0.0_rk)
    ! gamma_p
    call read_param_mpi(FILE, 'ACM-new', 'gamma_p', params_acm%gamma_p, 1.0_rk)
    call read_param_mpi(FILE, 'ACM-new', 'u_mean_set', params_acm%u_mean_set, (/1.0_rk, 0.0_rk, 0.0_rk/) )
    call read_param_mpi(FILE, 'ACM-new', 'beta', params_acm%beta, 0.05_rk )
    call read_param_mpi(FILE, 'ACM-new', 'compute_flow', params_acm%compute_flow, .true. )
    call read_param_mpi(FILE, 'ACM-new', 'HIT_linear_forcing', params_acm%HIT_linear_forcing, .false. )
    call read_param_mpi(FILE, 'ACM-new', 'HIT_energy', params_acm%HIT_energy, 1.0_rk )
    call read_param_mpi(FILE, 'ACM-new', 'HIT_gain', params_acm%HIT_gain, 100.0_rk )
    call read_param_mpi(FILE, 'ACM-new', 'nonlinear_formulation', params_acm%nonlinear_formulation, "convective" )


    ! initial condition
    call read_param_mpi(FILE, 'ACM-new', 'inicond', params_acm%inicond, "meanflow")
    ! the free flight FSI solver needs to know if it resumes a backup or not
    call read_param_mpi(FILE, 'Physics', 'read_from_files', params_acm%read_from_files, .false.)
    ! free flight also requires the time at which we resume (the structure of wabbit main does no allow to pass it to this routine...)
    if (params_acm%read_from_files) then
        call read_param_mpi(FILE, 'Physics', 'input_files', input_files, "")
        timestamp = input_files( index(input_files,'_')+1:index(input_files,'.h5') )
        read(timestamp,*) params_acm%start_time
        ! note this requires to have timestamp in the filename (so we cannot rename files...)
        params_acm%start_time = params_acm%start_time * 1.0e-6
    endif


    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_acm%discretization, "FD_4th_central_optimized")
    call read_param_mpi(FILE, 'Discretization', 'filter_type', params_acm%filter_type, "no_filter")
    call read_param_mpi(FILE, 'Discretization', 'order_predictor', params_acm%order_predictor, "multiresolution_4th")

    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params_acm%coarsening_indicator, "threshold-state-vector")
    call read_param_mpi(FILE, 'Blocks', 'eps_norm', params_acm%eps_norm, "Linfty")

    ! penalization:
    call read_param_mpi(FILE, 'VPM', 'penalization', params_acm%penalization, .true.)
    call read_param_mpi(FILE, 'VPM', 'C_eta', params_acm%C_eta, 1.0_rk)
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', params_acm%smooth_mask, .true.)
    call read_param_mpi(FILE, 'VPM', 'geometry', params_acm%geometry, "cylinder")
    
    call read_param_mpi(FILE, 'VPM', 'x_cntr', params_acm%x_cntr, (/0.5*params_acm%domain_size(1), 0.5*params_acm%domain_size(2), 0.5*params_acm%domain_size(3)/)  )
    call read_param_mpi(FILE, 'VPM', 'R_cyl', params_acm%R_cyl, 0.5_rk )
    call read_param_mpi(FILE, 'VPM', 'length', params_acm%length, 1.0_rk )
    call read_param_mpi(FILE, 'VPM', 'thickness', params_acm%thickness, 1.0_rk )
    call read_param_mpi(FILE, 'VPM', 'freq', params_acm%freq, 1.0_rk )
    call read_param_mpi(FILE, 'VPM', 'C_smooth', params_acm%C_smooth, 1.5_rk )

    ! stuff for lamballais cylinder
    call read_param_mpi(FILE, 'VPM', 'R0', params_acm%R0, 0.50_rk)
    call read_param_mpi(FILE, 'VPM', 'R1', params_acm%R1, 1.75_rk)
    call read_param_mpi(FILE, 'VPM', 'R2', params_acm%R2, 2.25_rk)
    call read_param_mpi(FILE, 'VPM', 'file_usx', params_acm%file_usx, 'none')
    call read_param_mpi(FILE, 'VPM', 'file_usy', params_acm%file_usy, 'none')
    call read_param_mpi(FILE, 'VPM', 'file_usp', params_acm%file_usp, 'none')
    call read_param_mpi(FILE, 'VPM', 'smoothing_type', params_acm%smoothing_type, 'hester')

    if (params_acm%geometry=="2D-wingsection" .or. params_acm%geometry=="two-moving-cylinders") then
        call read_param_mpi(FILE, 'VPM', 'wingsection_inifiles', params_acm%wingsection_inifiles, (/"", ""/))
    endif

    call read_param_mpi(FILE, 'Sponge', 'use_sponge', params_acm%use_sponge, .false. )
    call read_param_mpi(FILE, 'Sponge', 'L_sponge', params_acm%L_sponge, 0.0_rk )
    call read_param_mpi(FILE, 'Sponge', 'C_sponge', params_acm%C_sponge, 1.0e-2_rk )
    call read_param_mpi(FILE, 'Sponge', 'sponge_type', params_acm%sponge_type, "rect" )
    call read_param_mpi(FILE, 'Sponge', 'p_sponge', params_acm%p_sponge, 20.0_rk )

    call read_param_mpi(FILE, 'Time', 'CFL', params_acm%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'CFL_eta', params_acm%CFL_eta, 0.99_rk   )
    call read_param_mpi(FILE, 'Time', 'CFL_nu', params_acm%CFL_nu, 0.99_rk*2.79_rk/(dble(params_acm%dim)*pi**2) )
    call read_param_mpi(FILE, 'Time', 'time_max', params_acm%T_end, 1.0_rk   )
    call read_param_mpi(FILE, 'Statistics', 'nsave_stats', params_acm%nsave_stats, 999999   )
    call read_param_mpi(FILE, 'FreeFlightSolver', 'use_free_flight_solver', params_acm%use_free_flight_solver, .false.   )

    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params_acm%Jmax, 1   )


    call read_param_mpi(FILE, 'ACM-new', 'use_passive_scalar', params_acm%use_passive_scalar, .false.)
    if (params_acm%use_passive_scalar) then
        call read_param_mpi(FILE, 'ConvectionDiffusion', 'N_scalars', params_acm%N_scalars, 1)

        allocate( params_acm%schmidt_numbers(1:params_acm%N_scalars) )
        allocate( params_acm%x0source(1:params_acm%N_scalars) )
        allocate( params_acm%y0source(1:params_acm%N_scalars) )
        allocate( params_acm%z0source(1:params_acm%N_scalars) )
        allocate( params_acm%widthsource(1:params_acm%N_scalars) )
        allocate( params_acm%scalar_Ceta(1:params_acm%N_scalars) )

        allocate( params_acm%scalar_inicond(1:params_acm%N_scalars) )
        allocate( params_acm%scalar_source_type(1:params_acm%N_scalars) )

        params_acm%scalar_inicond = "dummy"
        params_acm%scalar_source_type = "dummy"

        call read_param_mpi( FILE, 'ConvectionDiffusion', 'Sc', params_acm%schmidt_numbers )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'x0source', params_acm%x0source )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'y0source', params_acm%y0source )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'z0source', params_acm%z0source )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'widthsource', params_acm%widthsource )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'C_eta', params_acm%scalar_Ceta)

        call read_param_mpi( FILE, 'ConvectionDiffusion', 'inicond', params_acm%scalar_inicond, params_acm%scalar_inicond )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'source', params_acm%scalar_source_type, params_acm%scalar_source_type )
        call read_param_mpi( FILE, 'ConvectionDiffusion', 'scalar_BC_type', params_acm%scalar_BC_type, "neumann" )

        if (params_acm%use_sponge) then
            call read_param_mpi( FILE, 'ConvectionDiffusion', 'absorbing_sponge', params_acm%absorbing_sponge, .true. )
        else
            params_acm%absorbing_sponge = .false.
        endif

        ! when computing passive scalars, we require derivatives of the mask function, which
        ! is not too difficult on paper. however, in wabbit, ghost node syncing is not a physics
        ! module task so the ACM module cannot do it. Note it has to be done only if scalars are used.
        ! Its a waste of resources otherwise. Hence, we have the flag to set masks on ghost nodes as well
        ! to set the mask on all points of a block (incl ghost nodes)
        params_acm%set_mask_on_ghost_nodes = .true.
    else
        params_acm%set_mask_on_ghost_nodes = .false.
        params_acm%N_scalars = 0
    endif

    ! set defaults
    if (params_acm%dim==3) then
      params_acm%Bs=(/17,17,17/)
    else
      params_acm%Bs=(/17,17,1/)
    endif
    params_acm%Bs = read_bs(FILE,'Blocks', 'number_block_nodes', params_acm%Bs, params_acm%dim)

    call clean_ini_file_mpi( FILE )

    if (Neqn /= params_acm%dim + 1 + params_acm%N_scalars) then
        ! call abort(220819, "the state vector length is not appropriate. number_equation must be DIM+1+N_scalars")
    endif

    ! uniqueGrid modification
    ddx(1:params_acm%dim) = 2.0_rk**(-params_acm%Jmax) * (params_acm%domain_size(1:params_acm%dim) / real(params_acm%Bs(1:params_acm%dim), kind=rk))

    dx_min = minval( ddx(1:params_acm%dim) )
    ! uniqueGrid modification
    nx_max = maxval( (params_acm%Bs) * 2**(params_acm%Jmax) )

    if (params_acm%c_0 > 0.0_rk) then
        dt_min = params_acm%CFL*dx_min/params_acm%c_0
    else
        dt_min = 0.0_rk
    endif
    ! nice to have this elsewhere in the ACM module:
    params_acm%dx_min = dx_min

    ! at most, we need 6 components: mask, usx, usy, usz, color, sponge
    ! in 2d, less arrays could be used, but its easier to just go ahead and use all of them.
    N_mask_components = 6

    if (params_acm%mpirank==0) then
      write(*,'(80("<"))')
      write(*,*) "Some information:"
      write(*,'("c0=",g12.4," C_eta=",g12.4," CFL=",g12.4)') params_acm%c_0, params_acm%C_eta, params_acm%CFL
      write(*,'("dx_min=",g12.4," dt(CFL,c0,dx_min)=",g12.4)') dx_min, dt_min
      write(*,'("if all blocks were at Jmax, the resolution would be nx=",i5)') nx_max
      if (params_acm%penalization) then
          write(*,'("C_eta=",es12.4," K_eta=",es12.4)') params_acm%C_eta, sqrt(params_acm%C_eta*params_acm%nu)/dx_min
      endif
      write(*,'("N_mask_components=",i1)') N_mask_components
      write(*,'("N_scalars=",i2)') params_acm%N_scalars
      write(*,'(80("<"))')
    endif

    ! if used, setup insect. Note fractal tree and active grid are part of the insects: they require the same init module
    !
    ! NOTE: there are several testing geometries used to test the free-flight solver: cylinder-free, sphere-free and plate-free (2D)
    if (params_acm%geometry == "Insect".or.params_acm%geometry=="fractal_tree".or.params_acm%geometry=="active_grid" &
    .or. params_acm%geometry=="cylinder-free".or.params_acm%geometry=="sphere-free".or.params_acm%geometry=="plate-free") then
        ! when computing passive scalars, we require derivatives of the mask function, which
        ! is not too difficult on paper. however, in wabbit, ghost node syncing is not a physics
        ! module task so the ACM module cannot do it. Note it has to be done only if scalars are used.
        ! Its a waste of resources otherwise. Hence, we have the flag to set masks on ghost nodes as well
        ! to set the mask on all points of a block (incl ghost nodes)
        if (params_acm%set_mask_on_ghost_nodes) then
            call insect_init( params_acm%start_time, filename, insect, params_acm%read_from_files, "", params_acm%domain_size, &
            params_acm%nu, dx_min, N_ghost_nodes=0)
        else
            call insect_init( params_acm%start_time, filename, insect, params_acm%read_from_files, "", params_acm%domain_size, &
            params_acm%nu, dx_min, N_ghost_nodes=g)
        endif
    endif

    if (params_acm%geometry=="fractal_tree") then
        call fractal_tree_init( Insect )
    endif

    if (params_acm%geometry=="2D-wingsection" .or. params_acm%geometry=="two-moving-cylinders") then
        call init_wingsection_from_file(params_acm%wingsection_inifiles(1), wingsections(1), 0.0_rk)
        call init_wingsection_from_file(params_acm%wingsection_inifiles(2), wingsections(2), 0.0_rk)
    endif


    ! read lamballais reference fields, see
    ! Gautier, R., Biau, D., Lamballais, E.: A reference solution of the flow over a circular cylinder at Re = 40 , Computers & Fluids 75, 103–111, 2013 
    if (params_acm%geometry == "lamballais") then
        if (params_acm%dim /= 2) call abort(1409241, "lamballais is a 2D test case")
        
        ! read us field
        call count_lines_in_ascii_file_mpi(params_acm%file_usx, num_lines, 0)
        ! avoid maxcolumns restriction (read in a single long column and reshape)
        allocate(buffer_array(1:num_lines,1:1))
        call read_array_from_ascii_file_mpi(params_acm%file_usx, buffer_array, 0)

        write(*,*) "lamballais num_lines", num_lines, nx_max


        if (num_lines /= nx_max**2) then
            call abort(2410011, "Lamballais: you seem to read the wrong field (size mismatch?!)")
        endif

        allocate(params_acm%u_lamballais(0:nx_max-1, 0:nx_max-1, 1:3))
        ! reshape
        params_acm%u_lamballais(:,:,1) = reshape(buffer_array, (/nx_max, nx_max/))

        call read_array_from_ascii_file_mpi(params_acm%file_usy, buffer_array, 0)
        params_acm%u_lamballais(:,:,2) = reshape(buffer_array, (/nx_max, nx_max/))

        call read_array_from_ascii_file_mpi(params_acm%file_usp, buffer_array, 0)
        params_acm%u_lamballais(:,:,3) = reshape(buffer_array, (/nx_max, nx_max/))
    endif

    params_acm%initialized = .true.
  end subroutine READ_PARAMETERS_ACM


  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_ACM( time, iteration, u, Bs, g, x0, dx, dt )
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
    real(kind=rk) :: u_eigen, kappa, uu, u_mag
    integer :: iscalar, dim, ix, iy ,iz

    if (.not. params_acm%initialized) write(*,*) "WARNING: GET_DT_BLOCK_ACM called but ACM not initialized"

    dim = params_acm%dim

    ! compute square of velocity magnitude
    if (params_acm%dim == 2) then
        u_mag = 0.0_rk
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                uu = u(ix,iy,1,1)*u(ix,iy,1,1)+u(ix,iy,1,2)*u(ix,iy,1,2)
                u_mag = max( u_mag, uu)
            enddo
        enddo

    else
        u_mag = 0.0_rk
        do iz=g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    uu = u(ix,iy,iz,1)*u(ix,iy,iz,1)+u(ix,iy,iz,2)*u(ix,iy,iz,2)+u(ix,iy,iz,3)*u(ix,iy,iz,3)
                    u_mag = max( u_mag, uu)
                enddo
            enddo
        enddo
    endif

    ! the velocity of the fast modes is u +- W and W= sqrt(c0^2 + u^2)
    u_eigen = sqrt(u_mag) + sqrt(params_acm%c_0**2 + u_mag )

    ! ususal CFL condition
    ! if the characteristic velocity is very small, avoid division by zero
    if ( u_eigen >= 1.0e-6_rk ) then
        dt = params_acm%CFL * minval(dx(1:dim)) / u_eigen
    else
        dt = 1.0e-2_rk
    endif

    ! explicit diffusion (NOTE: factor 0.5 is valid only for RK4, other time steppers have more
    ! severe restrictions)
    if (params_acm%nu>1.0e-13_rk) then
        dt = min(dt,  params_acm%CFL_nu * minval(dx(1:dim))**2 / params_acm%nu)
    endif

    ! diffusivity of scalars, if they are in use.(NOTE: factor 0.5 is valid only for RK4, other time steppers have more
    ! severe restrictions)
    if (params_acm%use_passive_scalar)  then
        do iscalar = 1, params_acm%N_scalars
            kappa = params_acm%nu * params_acm%schmidt_numbers(iscalar)
            dt = min(dt, params_acm%CFL_nu * minval(dx(1:dim))**2 / kappa)
        enddo
    endif

    ! just for completeness...this condition should never be active (gamma ~ 1)
    if (params_acm%gamma_p>0) dt = min( dt, params_acm%CFL_eta*params_acm%gamma_p )

    ! penalization
    if (params_acm%penalization) dt = min( dt, params_acm%CFL_eta*params_acm%C_eta )

    ! sponge
    if (params_acm%use_sponge) dt = min( dt, params_acm%CFL_eta*params_acm%C_sponge )

  end subroutine GET_DT_BLOCK_ACM



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

      min_dx = 2.0_rk**(-params_acm%Jmax) * min(params_acm%domain_size(1)/real(params_acm%Bs(1)-1, kind=rk) &
                                                  , params_acm%domain_size(2)/real(params_acm%Bs(2)-1, kind=rk))
      ! u(x=0) should be set equal to u(x=L)
      if ( abs(x-L)<min_dx*0.5_rk ) then
          x = 0.0_rk
      end if

  end subroutine continue_periodic

  !-----------------------------------------------------------------------------
  ! Adaptation is dependent on the different physics application.
  ! Every physics module can choose its own coarsening indicator.
  !-----------------------------------------------------------------------------
  subroutine PREPARE_THRESHOLDFIELD_ACM( u, g, x0, dx, threshold_field, N_thresholding_components )
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
      call abort(27022017,"The ACM module supports only coarsening_indicator=threshold-state-vector; !")
  end subroutine


  ! the statistics are written to ascii files (usually *.t files) with the help
  ! of module_t_files. In any case, the files have to be intialized: ideally, they
  ! are equipped with a header and resetted on the very first call, and they must not be deleted
  ! if the simuation is resumed from a backup. We therefore provide this function so that all physics
  ! modules can initialize those files.
  subroutine INITIALIZE_ASCII_FILES_ACM( time, overwrite )
      implicit none

      ! it may happen that some source terms have an explicit time-dependency
      ! therefore the general call has to pass time
      real(kind=rk), intent (in) :: time
      logical, intent(in) :: overwrite
      logical :: is_insect, use_color

      is_insect = .false.
      use_color = .false.
      if (params_acm%geometry == "Insect") is_insect = .true.
      if (params_acm%geometry == "two-moving-cylinders") use_color = .true.



      call init_t_file('meanflow.t', overwrite)
      call init_t_file('forces.t', overwrite)
      call init_t_file('e_kin.t', overwrite, (/"           time", "          e_kin"/))
      call init_t_file('turbulent_statistics.t', overwrite, (/"           time", "    dissipation", "         energy", "          u_RMS", &
      "    kolm_length", "      kolm_time", "  kolm_velocity", "   taylor_micro", "reynolds_taylor"/))
      call init_t_file('enstrophy.t', overwrite)
      call init_t_file('div.t', overwrite)
      call init_t_file('umag.t', overwrite)
      ! write(44,'(5(A15,1x))') "%          time","u_max","c0","MachNumber","u_eigen"
      call init_t_file('CFL.t', overwrite, (/&
      "           time", &
      "            CFL", &
      "         CFL_nu", &
      "        CFL_eta"/) )
      if (is_insect) then
        call init_t_file('moments.t', overwrite)
        call init_t_file('aero_power.t', overwrite)
        call init_t_file('forces_body.t', overwrite)
        call init_t_file('moments_body.t', overwrite)
        call init_t_file('forces_leftwing.t', overwrite)
        call init_t_file('moments_leftwing.t', overwrite)
        call init_t_file('forces_rightwing.t', overwrite)
        call init_t_file('moments_rightwing.t', overwrite)
        call init_t_file('insect_state_vector.t', overwrite)
      endif
      if (use_color) then
        call init_t_file('forces_1.t', overwrite)
        call init_t_file('forces_2.t', overwrite)
      endif
      call init_t_file('mask_volume.t', overwrite)
      call init_t_file('u_residual.t', overwrite)
      call init_t_file('kinematics.t', overwrite)
      call init_t_file('forces_rk.t', overwrite)
      call init_t_file('penal_power.t', overwrite, (/&
      "           time", &
      "  E_dot_f_solid"/))

      if (Insect%second_wing_pair) then
          call init_t_file('forces_leftwing2.t', overwrite)
          call init_t_file('moments_leftwing2.t', overwrite)
          call init_t_file('forces_rightwing2.t', overwrite)
          call init_t_file('moments_rightwing2.t', overwrite)
          call init_t_file('kinematics.t', overwrite, (/&
          "           time", &
          "    xc_body_g_x", &
          "    xc_body_g_y", &
          "    xc_body_g_z", &
          "      psi (rad)", &
          "     beta (rad)", &
          "    gamma (rad)", &
          "      eta (rad)", &
          "  alpha_l (rad)", &
          "    phi_l (rad)", &
          "  theta_l (rad)", &
          "  alpha_r (rad)", &
          "    phi_r (rad)", &
          "  theta_r (rad)", &
          "  rot_rel_l_w_x", &
          "  rot_rel_l_w_y", &
          "  rot_rel_l_w_z", &
          "  rot_rel_r_w_x", &
          "  rot_rel_r_w_y", &
          "  rot_rel_r_w_z", &
          "   rot_dt_l_w_x", &
          "   rot_dt_l_w_y", &
          "   rot_dt_l_w_z", &
          "   rot_dt_r_w_x", &
          "   rot_dt_r_w_y", &
          "   rot_dt_r_w_z", &
          " alpha_l2 (rad)", &
          "   phi_l2 (rad)", &
          " theta_l2 (rad)", &
          " alpha_r2 (rad)", &
          "   phi_r2 (rad)", &
          " theta_r2 (rad)", &
          " rot_rel_l2_w_x", &
          " rot_rel_l2_w_y", &
          " rot_rel_l2_w_z", &
          " rot_rel_r2_w_x", &
          " rot_rel_r2_w_y", &
          " rot_rel_r2_w_z", &
          "  rot_dt_l2_w_x", &
          "  rot_dt_l2_w_y", &
          "  rot_dt_l2_w_z", &
          "  rot_dt_r2_w_x", &
          "  rot_dt_r2_w_y", &
          "  rot_dt_r2_w_z"/) )
      else
          call init_t_file('kinematics.t', overwrite, (/&
          "           time", &
          "    xc_body_g_x", &
          "    xc_body_g_y", &
          "    xc_body_g_z", &
          "      psi (rad)", &
          "     beta (rad)", &
          "    gamma (rad)", &
          "      eta (rad)", &
          "  alpha_l (rad)", &
          "    phi_l (rad)", &
          "  theta_l (rad)", &
          "  alpha_r (rad)", &
          "    phi_r (rad)", &
          "  theta_r (rad)", &
          "  rot_rel_l_w_x", &
          "  rot_rel_l_w_y", &
          "  rot_rel_l_w_z", &
          "  rot_rel_r_w_x", &
          "  rot_rel_r_w_y", &
          "  rot_rel_r_w_z", &
          "   rot_dt_l_w_x", &
          "   rot_dt_l_w_y", &
          "   rot_dt_l_w_z", &
          "   rot_dt_r_w_x", &
          "   rot_dt_r_w_y", &
          "   rot_dt_r_w_z"/) )
      endif

  end subroutine INITIALIZE_ASCII_FILES_ACM

  ! this tiny wrapper avoids us to use "module_insects" as public module in some places.
  subroutine Update_Insect_wrapper(time)
      implicit none
      real(kind=rk), intent(in) :: time
      call Update_Insect(time, Insect)
  end subroutine

end module module_acm
