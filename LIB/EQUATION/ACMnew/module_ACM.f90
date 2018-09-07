!> \dir
!> \brief Implementation of 3d/2d acm physics

! ********************************************************************************************
!> Module for 2D/3D acm physics
! ********************************************************************************************
!> \details
!> \version 0.5
!> \author engels
!! \date pls add creation date
!!
! ********************************************************************************************

module module_acm

  !---------------------------------------------------------------------------------------------
  ! modules

  use mpi
  use module_insects

  use module_precision
  ! ini file parser module, used to read parameters. note: in principle, you can also
  ! just use any reader you feel comfortable with, as long as you can read the parameters
  ! from a file.
  use module_ini_files_parser_mpi
  use module_operators, only : compute_vorticity, divergence
  use module_helpers, only : startup_conditioner, smoothstep

  !---------------------------------------------------------------------------------------------
  ! variables

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: READ_PARAMETERS_ACM, PREPARE_SAVE_DATA_ACM, RHS_ACM, GET_DT_BLOCK_ACM, &
  INICOND_ACM, FIELD_NAMES_ACM, STATISTICS_ACM
  !**********************************************************************************************

  ! user defined data structure for time independent parameters, settings, constants
  ! and the like. only visible here.
  type :: type_params
    real(kind=rk) :: CFL, T_end
    real(kind=rk) :: c_0
    real(kind=rk) :: C_eta, beta
    ! nu
    real(kind=rk) :: nu
    real(kind=rk) :: x_cntr(1:3), u_cntr(1:3), R_cyl, u_mean_set(1:3), force(1:3)
    ! gamma_p
    real(kind=rk) :: gamma_p
    ! want to add forcing?
    logical :: forcing, penalization, smooth_mask=.True.
    ! the mean pressure has no meaning in incompressible fluids, but sometimes it can
    ! be nice to ensure the mean is zero, e.g., for comparison wit other codes. if set to true
    ! wabbit removes the mean pressure at every time step.
    logical :: p_mean_zero
    ! sponge term:
    logical :: use_sponge=.false.
    real(kind=rk) :: C_sponge, L_sponge

    integer(kind=ik) :: dim, N_fields_saved
    real(kind=rk), dimension(3)      :: domain_size=0.0_rk
    character(len=80) :: inicond, discretization, geometry="cylinder"
    character(len=80), allocatable :: names(:), forcing_type(:)
    ! the mean flow, as required for some forcing terms. it is computed in the RHS
    real(kind=rk) :: mean_flow(1:3), mean_p
    ! the error compared to an analytical solution (e.g. taylor-green)
    real(kind=rk) :: error(1:6)
    ! kinetic energy and enstrophy (both integrals)
    real(kind=rk) :: e_kin, enstrophy
    ! we need to know which mpirank prints output..
    integer(kind=ik) :: mpirank, mpisize
    !
    integer(kind=ik) :: Jmax, Bs

  end type type_params

  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.
  type(type_params), save :: params_acm

  ! all parameters for insects go here:
  type(diptera), save :: insect

  !---------------------------------------------------------------------------------------------
  ! variables initialization

  !---------------------------------------------------------------------------------------------
  ! main body

contains

  include "rhs.f90"
  include "create_mask.f90"
  include "inicond_ACM.f90"
  include "sponge.f90"
  include "save_data_ACM.f90"
  include "statistics_ACM.f90"

  !-----------------------------------------------------------------------------
  ! main level wrapper routine to read parameters in the physics module. It reads
  ! from the same ini file as wabbit, and it reads all it has to know. note in physics modules
  ! the parameter struct for wabbit is not available.
  subroutine READ_PARAMETERS_ACM( filename )
    implicit none

    character(len=*), intent(in) :: filename
    integer(kind=ik) :: mpicode, nx_max
    real(kind=rk) :: dx_min, dt_min

    ! inifile structure
    type(inifile) :: FILE


    ! we still need to know about mpirank and mpisize, occasionally
    call MPI_COMM_SIZE (WABBIT_COMM, params_acm%mpisize, mpicode)
    call MPI_COMM_RANK (WABBIT_COMM, params_acm%mpirank, mpicode)

    if (params_acm%mpirank==0) then
      write(*,'(80("<"))')
      write(*,*) "Initializing artificial compressibility module!"
      write(*,'(80("<"))')
    endif

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    call read_param_mpi(FILE, 'Dimensionality', 'dim', params_acm%dim, 2 )
    call read_param_mpi(FILE, 'Domain', 'domain_size', params_acm%domain_size(1:params_acm%dim), (/ 1.0_rk, 1.0_rk, 1.0_rk /) )
    ! call read_param_mpi(FILE, 'DomainSize', 'Lx', params_acm%domain_size(1), 1.0_rk )
    ! call read_param_mpi(FILE, 'DomainSize', 'Ly', params_acm%domain_size(2), 1.0_rk )
    ! call read_param_mpi(FILE, 'DomainSize', 'Lz', params_acm%domain_size(3), 0.0_rk )

    ! --- saving ----
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_acm%N_fields_saved, 3 )
    allocate( params_acm%names(1:params_acm%N_fields_saved) )
    call read_param_mpi(FILE, 'Saving', 'field_names', params_acm%names, (/"ux","uy","p "/) )


    ! speed of sound for acm
    call read_param_mpi(FILE, 'ACM-new', 'c_0', params_acm%c_0, 10.0_rk)
    ! viscosity
    call read_param_mpi(FILE, 'ACM-new', 'nu', params_acm%nu, 1e-1_rk)
    ! gamma_p
    call read_param_mpi(FILE, 'ACM-new', 'gamma_p', params_acm%gamma_p, 1.0_rk)
    ! want to add a forcing term?
    call read_param_mpi(FILE, 'ACM-new', 'forcing', params_acm%forcing, .false.)
    allocate( params_acm%forcing_type(1:3) )
    call read_param_mpi(FILE, 'ACM-new', 'forcing_type', params_acm%forcing_type, (/"accelerate","none      ","none      "/) )
    call read_param_mpi(FILE, 'ACM-new', 'u_mean_set', params_acm%u_mean_set, (/1.0_rk, 0.0_rk, 0.0_rk/) )
    call read_param_mpi(FILE, 'ACM-new', 'p_mean_zero', params_acm%p_mean_zero, .false. )
    call read_param_mpi(FILE, 'ACM-new', 'beta', params_acm%beta, 0.05_rk )


    ! initial condition
    call read_param_mpi(FILE, 'ACM-new', 'inicond', params_acm%inicond, "meanflow")

    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_acm %discretization, "FD_2nd_central")
    ! penalization:
    call read_param_mpi(FILE, 'VPM', 'penalization', params_acm%penalization, .true.)
    call read_param_mpi(FILE, 'VPM', 'C_eta', params_acm%C_eta, 1.0e-3_rk)
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', params_acm%smooth_mask, .true.)
    call read_param_mpi(FILE, 'VPM', 'geometry', params_acm%geometry, "cylinder")
    call read_param_mpi(FILE, 'VPM', 'x_cntr', params_acm%x_cntr, (/0.5*params_acm%domain_size(1), 0.5*params_acm%domain_size(2), 0.5*params_acm%domain_size(3)/)  )
    call read_param_mpi(FILE, 'VPM', 'R_cyl', params_acm%R_cyl, 0.5_rk )

    call read_param_mpi(FILE, 'Sponge', 'use_sponge', params_acm%use_sponge, .false. )
    call read_param_mpi(FILE, 'Sponge', 'L_sponge', params_acm%L_sponge, 0.0_rk )
    call read_param_mpi(FILE, 'Sponge', 'C_sponge', params_acm%C_sponge, 1.0e-2_rk )

    call read_param_mpi(FILE, 'Time', 'CFL', params_acm%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'time_max', params_acm%T_end, 1.0_rk   )


    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params_acm%Jmax, 1   )
    call read_param_mpi(FILE, 'Blocks', 'number_block_nodes', params_acm%Bs, 1   )


    call clean_ini_file_mpi( FILE )


    if (params_acm%mpirank==0) then
      write(*,'(80("<"))')
      write(*,*) "Some information:"
      write(*,'("c0=",g12.4," C_eta=",g12.4," CFL=",g12.4)') params_acm%c_0, params_acm%C_eta, params_acm%CFL

      dx_min = 2.0_rk**(-params_acm%Jmax) * params_acm%domain_size(1) / real(params_acm%Bs-1, kind=rk)
      nx_max = (params_acm%Bs-1) * 2**(params_acm%Jmax)
      dt_min = params_acm%CFL*dx_min/params_acm%c_0

      write(*,'("dx_min=",g12.4," dt(CFL,c0,dx_min)=",g12.4)') dx_min, dt_min
      write(*,'("if all blocks were at Jmax, the resolution would be nx=",i5)') nx_max

      if (params_acm%penalization) then
          write(*,'("C_eta=",g12.4," K_eta=",g12.4)') params_acm%C_eta, sqrt(params_acm%C_eta*params_acm%nu)/dx_min
      endif

      write(*,'(80("<"))')
    endif

    ! if used, setup insect
    if (params_acm%geometry == "Insect") then
        call insect_init( 0.0_rk, filename, insect, .false.,"" , (/params_acm%domain_size(1), params_acm%domain_size(2), params_acm%domain_size(3)/), params_acm%nu)
    endif
  end subroutine READ_PARAMETERS_ACM


  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_ACM( time, u, Bs, g, x0, dx, dt )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: Bs, g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! the dt for this block is returned to the caller:
    real(kind=rk), intent(out) :: dt
    ! temporary array. note this is just one block and hence not important for overall memory consumption
    real(kind=rk), allocatable, save :: u_mag(:,:,:)
    real(kind=rk) :: u_eigen

    if (.not.allocated(u_mag)) allocate(u_mag(1:size(u,1), 1:size(u,2), 1:size(u,3)))

    ! compute square of velocity magnitude
    if (params_acm%dim == 2) then
        u_mag = u(:,:,:,1)*u(:,:,:,1) + u(:,:,:,2)*u(:,:,:,2)
    else
        u_mag = u(:,:,:,1)*u(:,:,:,1) + u(:,:,:,2)*u(:,:,:,2) + u(:,:,:,3)*u(:,:,:,3)
    endif
    ! the velocity of the fast modes is u +- W and W= sqrt(c0^2 + u^2)
    u_eigen = sqrt(maxval(u_mag)) + sqrt(params_acm%c_0**2 + maxval(u_mag) )

    ! ususal CFL condition
    dt = params_acm%CFL * minval(dx(1:params_acm%dim)) / u_eigen

    ! explicit diffusion (NOTE: factor 0.5 is valid only for RK4, other time steppers have more
    ! severe restrictions)
    if (params_acm%nu>1.0e-13_rk) dt = min(dt, 0.5_rk * minval(dx(1:params_acm%dim))**2 / params_acm%nu)

    ! just for completeness...this condition should never be active (gamma ~ 1)
    if (params_acm%gamma_p>0) dt = min( dt, 0.99_rk*params_acm%gamma_p )

    ! penalization
    if (params_acm%penalization) dt = min( dt, 0.99_rk*params_acm%C_eta )

    ! sponge
    if (params_acm%use_sponge) dt = min( dt, 0.99_rk*params_acm%C_sponge )

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

      min_dx = 2.0_rk**(-params_acm%Jmax) * min(params_acm%domain_size(1),params_acm%domain_size(2))&
                        / real(params_acm%Bs-1, kind=rk)
      ! u(x=0) should be set equal to u(x=L)
      if ( abs(x-L)<min_dx*0.5_rk ) then
        x = 0.0_rk
      end if

  end subroutine continue_periodic

end module module_acm
