!> \file
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_acm.f90
!> \version 0.5
!> \author engels
!!
!! \brief module for 2D/3D acm physics
!!
!!
!! = log ======================================================================================
!! \n
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
  use module_operators, only : compute_vorticity
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
    real(kind=rk) :: x_cntr(1:3), u_cntr(1:3), R_cyl, u_mean_set(1:3)
    ! gamma_p
    real(kind=rk) :: gamma_p
    ! want to add forcing?
    logical :: forcing, penalization,smooth_mask=.True., sponge_layer
    ! alpha for sponge
    real(kind=rk) :: alpha
    integer(kind=ik) :: dim, N_fields_saved
    character(len=80) :: inicond, discretization, geometry="cylinder"
    character(len=80), allocatable :: names(:), forcing_type(:)
    ! the mean flow, as required for some forcing terms. it is computed in the RHS
    real(kind=rk) :: mean_flow(1:3)
    ! we need to know which mpirank prints output..
    integer(kind=ik) :: mpirank, mpisize
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
    integer(kind=ik) :: mpicode

    ! inifile structure
    type(inifile) :: FILE

    ! we still need to know about mpirank and mpisize, occasionally
    call MPI_COMM_SIZE (MPI_COMM_WORLD, params_acm%mpisize, mpicode)
    call MPI_COMM_RANK (MPI_COMM_WORLD, params_acm%mpirank, mpicode)

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

    call read_param_mpi(FILE, 'Time', 'CFL', params_acm%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'time_max', params_acm%T_end, 1.0_rk   )
    ! sponge layer:
    call read_param_mpi(FILE, 'Physics', 'sponge_layer', params_acm%sponge_layer, .false.)
    call read_param_mpi(FILE, 'Physics', 'alpha', params_acm%alpha, 100.0_rk)

    call clean_ini_file_mpi( FILE )
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
    integer(kind=ik) :: neqn, nwork, Bs

    ! number of state variables
    neqn = size(u,4)
    ! number of available work array slots
    nwork = size(work,4)

    Bs = size(u,1)-2*g

    ! copy state vector
    work(:,:,:,1:size(u,4)) = u(:,:,:,:)

    ! vorticity
    call compute_vorticity(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization,&
    work(:,:,:,4:6))

    ! mask
    call create_mask_2D_NEW(work(:,:,1,5), x0, dx, Bs, g )

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
    real(kind=rk) :: tmp(1:3)

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

      if (params_acm%forcing) then
        if (params_acm%dim == 2) then
          params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g, g+1:Bs+g, 1, 1))*dx(1)*dx(2)
          params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g, g+1:Bs+g, 1, 2))*dx(1)*dx(2)
        else
          params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 1))*dx(1)*dx(2)*dx(3)
          params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 2))*dx(1)*dx(2)*dx(3)
          params_acm%mean_flow(3) = params_acm%mean_flow(3) + sum(u(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 3))*dx(1)*dx(2)*dx(3)
        endif ! NOTE: MPI_SUM is perfomed in the post_stage.
      endif

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.

      if (params_acm%forcing) then
        tmp = params_acm%mean_flow
        call MPI_ALLREDUCE(tmp, params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
        if (params_acm%dim == 2) then
          params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly)
        else
          params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
        endif
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
  subroutine STATISTICS_ACM( time, u, g, x0, dx, rhs, stage )
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
    real(kind=rk) :: tmp(1:3)

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

      if (params_acm%dim == 2) then
        params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g, g+1:Bs+g, 1, 1))*dx(1)*dx(2)
        params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g, g+1:Bs+g, 1, 2))*dx(1)*dx(2)
      else
        params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 1))*dx(1)*dx(2)*dx(3)
        params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 2))*dx(1)*dx(2)*dx(3)
        params_acm%mean_flow(3) = params_acm%mean_flow(3) + sum(u(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 3))*dx(1)*dx(2)*dx(3)
      endif ! NOTE: MPI_SUM is perfomed in the post_stage.

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.

      tmp = params_acm%mean_flow
      call MPI_ALLREDUCE(tmp, params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
      if (params_acm%dim == 2) then
        params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly)
      else
        params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
      endif

      if (params_acm%mpirank == 0) then
        ! write mean flow to disk...
        open(14,file='meanflow.t',status='unknown',position='append')
        write (14,'(4(es15.8,1x))') time, params_acm%mean_flow
        close(14)
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
  subroutine INICOND_ACM( time, u, g, x0, dx )
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


    select case (params_acm%inicond)
    case("meanflow")
      u = 0.0_rk
      u(:,:,:,1) = 1.0_rk
    case default
      write(*,*) "errorrroororor"
    end select

  end subroutine INICOND_ACM



end module module_acm_new
