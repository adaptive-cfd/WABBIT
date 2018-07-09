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

module module_acm_new

  !---------------------------------------------------------------------------------------------
  ! modules

  use module_precision
  ! ini file parser module, used to read parameters. note: in principle, you can also
  ! just use any reader you feel comfortable with, as long as you can read the parameters
  ! from a file.
  use module_ini_files_parser_mpi
  use module_operators, only : compute_vorticity, divergence
  use mpi
  !---------------------------------------------------------------------------------------------
  ! variables

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: READ_PARAMETERS_ACM, PREPARE_SAVE_DATA_ACM, RHS_ACM, GET_DT_BLOCK_ACM, INICOND_ACM, FIELD_NAMES_ACM, STATISTICS_ACM
  !**********************************************************************************************

  ! user defined data structure for time independent parameters, settings, constants
  ! and the like. only visible here.
  type :: type_params
    real(kind=rk) :: CFL, T_end
    real(kind=rk) :: c_0
    real(kind=rk) :: C_eta
    ! nu
    real(kind=rk) :: nu, Lx, Ly, Lz
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



  !---------------------------------------------------------------------------------------------
  ! variables initialization

  !---------------------------------------------------------------------------------------------
  ! main body

contains

  include "rhs.f90"
  include "create_mask_new.f90"
  ! include "dt.f90"
  include "iniconds.f90"
  include "sponge_new.f90"

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

    call read_param_mpi(FILE, 'DomainSize', 'Lx', params_acm%Lx, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', params_acm%Ly, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Lz', params_acm%Lz, 0.0_rk )

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


    ! initial condition
    call read_param_mpi(FILE, 'ACM-new', 'inicond', params_acm%inicond, "meanflow")

    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_acm %discretization, "FD_2nd_central")
    ! penalization:
    call read_param_mpi(FILE, 'VPM', 'penalization', params_acm%penalization, .true.)
    call read_param_mpi(FILE, 'VPM', 'C_eta', params_acm%C_eta, 1.0e-3_rk)
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', params_acm%smooth_mask, .true.)
    call read_param_mpi(FILE, 'VPM', 'geometry', params_acm%geometry, "cylinder")
    call read_param_mpi(FILE, 'VPM', 'x_cntr', params_acm%x_cntr, (/0.5*params_acm%Lx, 0.5*params_acm%Ly, 0.5*params_acm%Lz/)  )
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
      dx_min = 2.0_rk**(-params_acm%Jmax) * params_acm%Lx / real(params_acm%Bs-1, kind=rk)
      nx_max = (params_acm%Bs-1) * 2**(params_acm%Jmax)
      dt_min = params_acm%CFL*dx_min/params_acm%c_0
      write(*,'("dx_min=",g12.4," dt(CFL,c0,dx_min)=",g12.4)') dx_min, dt_min
      write(*,'("if all blocks were at Jmax, the resolution would be nx=",i5)') nx_max
      write(*,'(80("<"))')
    endif
  end subroutine READ_PARAMETERS_ACM


  !-----------------------------------------------------------------------------
  ! save data. Since you might want to save derived data, such as the vorticity,
  ! the divergence etc., which are not in your state vector, this routine has to
  ! copy and compute what you want to save to the work array.
  !
  ! In the main code, save_fields than saves the first N_fields_saved components of the
  ! work array to file.
  !
  ! NOTE that as we have way more work arrays than actual state variables (typically
  ! for a RK4 that would be >= 4*dim), you can compute a lot of stuff, if you want to.
  !-----------------------------------------------------------------------------
  subroutine PREPARE_SAVE_DATA_ACM( time, u, g, x0, dx, work )
    implicit none
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! output in work array.
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

    ! local variables
    integer(kind=ik)  :: neqn, nwork, Bs, k
    character(len=80) :: name

    ! number of state variables
    neqn = size(u,4)
    ! number of available work array slots
    nwork = size(work,4)

    Bs = size(u,1)-2*g

    ! copy state vector
    work(:,:,:,1:size(u,4)) = u(:,:,:,:)

    do k = neqn, size(params_acm%names,1)
        name = params_acm%names(k)
        select case(name(1:3))
            case('vor')
                ! vorticity
                call compute_vorticity(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), &
                    dx, Bs, g, params_acm%discretization, work(:,:,:,k:k+3))
            case('div')
                ! div(u)
                call divergence(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, &
                    g, params_acm%discretization,work(:,:,:,k))
            case('mas')
                ! mask
                call create_mask_2D_NEW(work(:,:,1,k), x0, dx, Bs, g )
            case('spo')
                ! mask for sponge
                call sponge_2D_NEW(work(:,:,1,k), x0, dx, Bs, g )
        end select
    end do

  end subroutine


  !-----------------------------------------------------------------------------
  ! when savig to disk, WABBIT would like to know how you named you variables.
  ! e.g. u(:,:,:,1) is called "ux"
  !
  ! the main routine save_fields has to know how you label the stuff you want to
  ! store from the work array, and this routine returns those strings
  !-----------------------------------------------------------------------------
  subroutine FIELD_NAMES_ACM( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=80), intent(out) :: name

    if (allocated(params_acm%names)) then
      name = params_acm%names(N)
    else
      call abort(5554,'Something ricked')
    endif

  end subroutine FIELD_NAMES_ACM


  !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid any an MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine RHS_ACM( time, u, g, x0, dx, rhs, stage )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

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
    real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in) :: stage

    ! local variables
    integer(kind=ik) :: Bs, mpierr
    real(kind=rk) :: tmp(1:3), tmp2

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    select case(stage)
    case ("init_stage")
      !-------------------------------------------------------------------------
      ! 1st stage: init_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.
      ! performs initializations in the RHS module, such as resetting integrals

      params_acm%mean_flow = 0.0_rk
      params_acm%mean_p = 0.0_rk

    case ("integral_stage")
      !-------------------------------------------------------------------------
      ! 2nd stage: init_stage.
      !-------------------------------------------------------------------------
      ! For some RHS, the eqn depend not only on local, block based qtys, such as
      ! the state vector, but also on the entire grid, for example to compute a
      ! global forcing term (e.g. in FSI the forces on bodies). As the physics
      ! modules cannot see the grid, (they only see blocks), in order to encapsulate
      ! them nicer, two RHS stages have to be defined: integral / local stage.
      !
      ! called for each block.

      if (maxval(abs(u))>1.0e5) then
        call abort(6661,"ACM fail: very very large values in state vector.")
      endif

      if (params_acm%dim == 2) then
        params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1))*dx(1)*dx(2)
        params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 2))*dx(1)*dx(2)
        params_acm%mean_p = params_acm%mean_p + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 3))*dx(1)*dx(2)
      else
        params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 1))*dx(1)*dx(2)*dx(3)
        params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 2))*dx(1)*dx(2)*dx(3)
        params_acm%mean_flow(3) = params_acm%mean_flow(3) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 3))*dx(1)*dx(2)*dx(3)
        params_acm%mean_p = params_acm%mean_p + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 4))*dx(1)*dx(2)*dx(3)
      endif ! NOTE: MPI_SUM is perfomed in the post_stage.

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.

      tmp = params_acm%mean_flow
      call MPI_ALLREDUCE(tmp, params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      tmp2 = params_acm%mean_p
      call MPI_ALLREDUCE(tmp2, params_acm%mean_p, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

      if (params_acm%dim == 2) then
        params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly)
        params_acm%mean_p = params_acm%mean_p / (params_acm%Lx*params_acm%Ly)
      else
        params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
        params_acm%mean_p = params_acm%mean_p / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
      endif

    case ("local_stage")
      !-------------------------------------------------------------------------
      ! 4th stage: local evaluation of RHS on all blocks
      !-------------------------------------------------------------------------
      ! the second stage then is what you would usually do: evaluate local differential
      ! operators etc.
      !
      ! called for each block.

      if (size(u,4) == 3) then
        ! this is a 2d case (ux,uy,p)
        call RHS_2D_acm_new(g, Bs, dx(1:2), x0(1:2), u(:,:,1,:), params_acm%discretization, &
        (/1.0_rk,0.0_rk,0.0_rk/), time, rhs(:,:,1,:))

      elseif (size(u,4) == 4) then
        ! this is a 3d case (ux,uy,uz,p)
        call abort(888,"module_ACM-new.f90: 3d not yet implemented")
      endif

    case default
      call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")
    end select


  end subroutine RHS_ACM


  !-----------------------------------------------------------------------------
  ! main level wrapper to compute statistics (such as mean flow, global energy,
  ! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
  ! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
  ! on individual blocks. This requires one to use the same staging concept as for the RHS.
  !-----------------------------------------------------------------------------
  subroutine STATISTICS_ACM( time, u, g, x0, dx, stage, work )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! work data, for mask, vorticity etc. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in) :: stage

    ! local variables
    integer(kind=ik) :: Bs, mpierr, ix, iy
    real(kind=rk) :: tmp(1:6)
    real(kind=rk) :: x, y
    real(kind=rk) :: eps_inv

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    select case(stage)
    case ("init_stage")
      !-------------------------------------------------------------------------
      ! 1st stage: init_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, NOT for each block.
      ! performs initializations in the RHS module, such as resetting integrals
      params_acm%mean_flow = 0.0_rk
      if (params_acm%forcing_type(1) .eq. "taylor_green") params_acm%error = 0.0_rk
      params_acm%force = 0.0_rk
      params_acm%e_kin = 0.0_rk
      params_acm%enstrophy = 0.0_rk

    case ("integral_stage")
      !-------------------------------------------------------------------------
      ! 2nd stage: integral_stage.
      !-------------------------------------------------------------------------
      ! This stage contains all operations which are running on the blocks
      !
      ! called for each block.

      if (maxval(abs(u))>1.0e5) then
        call abort(6661,"ACM fail: very very large values in state vector.")
      endif

      !-------------------------------------------------------------------------
      ! compute mean flow for output in statistics
      if (params_acm%dim == 2) then
        params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1))*dx(1)*dx(2)
        params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 2))*dx(1)*dx(2)
      else
        params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 1))*dx(1)*dx(2)*dx(3)
        params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 2))*dx(1)*dx(2)*dx(3)
        params_acm%mean_flow(3) = params_acm%mean_flow(3) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 3))*dx(1)*dx(2)*dx(3)
      endif ! NOTE: MPI_SUM is perfomed in the post_stage.

      !-------------------------------------------------------------------------
      ! if the forcing is taylor-green, then we know the exact solution in time. Therefore
      ! we compute the error w.r.t. this solution heres
      if (params_acm%forcing_type(1) .eq. "taylor_green") then
        do iy = g+1,Bs+g
          do ix = g+1, Bs+g
              x = x0(1) + dble(ix-g-1)*dx(1)
              y = x0(2) + dble(iy-g-1)*dx(2)
              tmp(1) = params_acm%u_mean_set(1) + dsin(x-params_acm%u_mean_set(1)*time)*&
                  dcos(y-params_acm%u_mean_set(2)*time)*dcos(time)
              tmp(2) = params_acm%u_mean_set(2) - dcos(x-params_acm%u_mean_set(1)*time)*&
                  dsin(y-params_acm%u_mean_set(2)*time)*dcos(time)
              tmp(3) = 0.25_rk*(dcos(2.0_rk*(x-params_acm%u_mean_set(1)*time)) +&
                  dcos(2.0_rk*(y-params_acm%u_mean_set(2)*time)))*dcos(time)**2
              params_acm%error(1:3) = params_acm%error(1:3) + abs(u(ix,iy,1,:)-&
                  tmp(1:3))
              params_acm%error(4:6) = params_acm%error(4:6) + sqrt(tmp(1:3)**2)
          end do
        end do
        params_acm%error = params_acm%error*dx(1)*dx(2)
      end if

      !-------------------------------------------------------------------------
      ! compute fluid force on penalized obstacle. The force can be computed by
      ! volume integration (which is much easier than surface integration), see
      ! Angot et al. 1999
      call create_mask_2D_NEW(work(:,:,1,1), x0, dx, Bs, g)
      eps_inv = 1.0_rk / params_acm%C_eta

      if (params_acm%dim == 2) then
        params_acm%force(1) = params_acm%force(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1)*work(g+1:Bs+g-1, g+1:Bs+g-1,1,1)*eps_inv)*dx(1)*dx(2)
        params_acm%force(2) = params_acm%force(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 2)*work(g+1:Bs+g-1, g+1:Bs+g-1,1,1)*eps_inv)*dx(1)*dx(2)
        params_acm%force(3) = 0.d0
      else
        call abort(6661,"ACM 3D not implemented.")
      endif

      !-------------------------------------------------------------------------
      ! compute kinetic energy in the whole domain (including penalized regions)
      if (params_acm%dim == 2) then
          params_acm%e_kin = params_acm%e_kin + 0.5_rk*sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1:2)**2)*dx(1)*dx(2)
      else
          params_acm%e_kin = params_acm%e_kin + 0.5_rk*sum(u(g+1:Bs+g-1,g+1:Bs+g-1,g+1:Bs+g-1, 1:3)**2)*dx(1)*dx(2)*dx(3)
      end if

      !-------------------------------------------------------------------------
      ! compute enstrophy in the whole domain (including penalized regions)
      call compute_vorticity(u(:,:,:,1), u(:,:,:,2), work(:,:,:,2), dx, Bs, g, params_acm%discretization, work(:,:,:,:))
      if (params_acm%dim ==2) then
          params_acm%enstrophy = params_acm%enstrophy + sum(work(g+1:Bs+g-1,g+1:Bs+g-1,1,1)**2)*dx(1)*dx(2)
      else
          call abort(6661,"ACM 3D not implemented.")
      end if

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, NOT for each block.


      !-------------------------------------------------------------------------
      ! mean flow
      tmp(1:3) = params_acm%mean_flow
      call MPI_ALLREDUCE(tmp(1:3), params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      if (params_acm%dim == 2) then
        params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly)
      else
        params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
      endif

      !-------------------------------------------------------------------------
      ! force
      tmp(1:3) = params_acm%force
      call MPI_ALLREDUCE(tmp(1:3), params_acm%force, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

      !-------------------------------------------------------------------------
      ! kinetic energy
      tmp(1) = params_acm%e_kin
      call MPI_ALLREDUCE(tmp(1), params_acm%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

      !-------------------------------------------------------------------------
      ! kinetic enstrophy
      tmp(1)= params_acm%enstrophy
      call MPI_ALLREDUCE(tmp(1), params_acm%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

      !-------------------------------------------------------------------------
      ! write statistics to ascii files.
      if (params_acm%mpirank == 0) then
        ! write mean flow to disk...
        open(14,file='meanflow.t',status='unknown',position='append')
        write (14,'(4(es15.8,1x))') time, params_acm%mean_flow
        close(14)

        ! write forces to disk...
        open(14,file='forces.t',status='unknown',position='append')
        write (14,'(4(es15.8,1x))') time, params_acm%force
        close(14)

        ! write kinetic energy to disk...
        open(14,file='e_kin.t',status='unknown',position='append')
        write (14,'(2(es15.8,1x))') time, params_acm%e_kin
        close(14)

        ! write enstrophy to disk...
        open(14,file='enstrophy.t',status='unknown',position='append')
        write (14,'(2(es15.8,1x))') time, params_acm%enstrophy
        close(14)
      end if

      if (params_acm%forcing_type(1) .eq. "taylor_green") then
        tmp = params_acm%error
        call MPI_REDUCE(tmp, params_acm%error, 6, MPI_DOUBLE_PRECISION, MPI_SUM, 0, WABBIT_COMM,mpierr)
        !params_acm%error(1:3) = params_acm%error(1:3)/params_acm%error(4:6)
        params_acm%error(1:3) = params_acm%error(1:3)/(params_acm%Lx*params_acm%Ly)

        if (params_acm%mpirank == 0) then
          ! write error to disk...
          open(15,file='error_taylor_green.t',status='unknown',position='append')
          write (15,'(4(es15.8,1x))') time, params_acm%error(1:3)
          close(15)
        end if

      end if
    case default
      call abort(7772,"the STATISTICS wrapper requests a stage this physics module cannot handle.")
    end select


  end subroutine STATISTICS_ACM


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

    dt = params_acm%CFL * dx(1) / params_acm%c_0

    if (params_acm%penalization) dt = min( dt, params_acm%C_eta )

  end subroutine GET_DT_BLOCK_ACM


  !-----------------------------------------------------------------------------
  ! main level wrapper for setting the initial condition on a block
  !-----------------------------------------------------------------------------
  subroutine INICOND_ACM( time, u, g, x0, dx, work, adapting )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! work data, for mask, vorticity etc. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! if we are still adapting the initial condition, we may use penalization for refinement.
    ! if the initial grid is adapted we set our initial condition without penalization (impulsive start).
    logical, intent(in) :: adapting

    real(kind=rk)    :: x,y
    integer(kind=ik) :: Bs, ix, iy

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    select case (params_acm%inicond)
    case("meanflow")
      u = 0.0_rk
      u(:,:,:,1) = params_acm%u_mean_set(1)
      u(:,:,:,2) = params_acm%u_mean_set(2)
      if (params_acm%dim == 3) then
        u(:,:,:,3) = params_acm%u_mean_set(3)
      endif
    case("taylor_green")
      do iy= 1,Bs+2*g
        do ix= 1, Bs+2*g
          x = x0(1) + dble(ix-g-1)*dx(1)
          y = x0(2) + dble(iy-g-1)*dx(2)
          call continue_periodic(x,params_acm%Lx)
          call continue_periodic(y,params_acm%Ly)
          u(ix,iy,1,1) = params_acm%u_mean_set(1) + dsin(x)*dcos(y)
          u(ix,iy,1,2) = params_acm%u_mean_set(2) - dcos(x)*dsin(y)
          u(ix,iy,1,3) = 0.25_rk*(dcos(2.0_rk*x) + dcos(2.0_rk*y))
        end do
      end do
    case default
      write(*,*) "errorrroororor"
    end select
    ! if we use volume penalization, the mask is first used for refinement of the grid.
    ! In a second stage, the initial condition without penalization is then applied to the refined grid.
    if (adapting .and. params_acm%penalization) then
        call create_mask_2D_NEW(work(:,:,1,1), x0, dx, Bs, g )
        u(:,:,:,1) = (1.0_rk-work(:,:,:,1))*u(:,:,:,1)
        u(:,:,:,2) = (1.0_rk-work(:,:,:,1))*u(:,:,:,2)
        if (params_acm%dim == 3) then
            u(:,:,:,3) = (1.0_rk-work(:,:,:,1))*u(:,:,:,3)
        end if
    end if

  end subroutine INICOND_ACM

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

      min_dx = 2.0_rk**(-params_acm%Jmax) * min(params_acm%Lx,params_acm%Ly)&
                        / real(params_acm%Bs-1, kind=rk)
      ! u(x=0) should be set equal to u(x=L)
      if ( abs(x-L)<min_dx*0.5_rk ) then
        x = 0.0_rk
      end if

  end subroutine continue_periodic

end module module_acm_new
