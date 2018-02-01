!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_navier_stokes_new.f90
!> \version 0.5
!> \author Pkrah
!!
!! \brief module for 2D/3D navier_stokes
!!
!!
!! = log ======================================================================================
!! \n
!! + creation 23.1.2018
!!
! ********************************************************************************************

module module_navier_stokes_new

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
  PUBLIC :: READ_PARAMETERS_NSTOKES, PREPARE_SAVE_DATA_NSTOKES, RHS_NSTOKES, GET_DT_BLOCK_NSTOKES, INICOND_NSTOKES, FIELD_NAMES_NStokes
  !**********************************************************************************************

  ! user defined data structure for time independent parameters, settings, constants
  ! and the like. only visible here.
  type :: type_params_ns
        ! Courant-Friedrichs-Lewy 
        real(kind=rk) :: CFL, T_end
        ! spatial domain
        real(kind=rk) :: Lx, Ly, Lz
        ! number data fields
        integer(kind=ik)                            :: number_data_fields
        ! dimension
        integer(kind=ik)                            :: dim, N_fields_saved
        ! adiabatic coefficient
        real(kind=rk)                               :: gamma_
        ! specific gas constant
        real(kind=rk)                               :: Rs
        ! isochoric heat capacity
        real(kind=rk)                               :: Cv
        ! isobaric heat capacity
        real(kind=rk)                               :: Cp
        ! prandtl number
        real(kind=rk)                               :: Pr
        ! dynamic viscosity
        real(kind=rk)                               :: mu0
        ! dissipation switch
        logical                                     :: dissipation
        ! variable names
        character(len=80), allocatable              :: names(:)
        ! discretization
        character(len=80)                           :: discretization

  end type type_params_ns

  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.
  type(type_params_ns), save :: params_ns
  ! statevector
  integer(kind=ik),save      :: rhoF,UxF,UyF,UzF,pF
      


  !---------------------------------------------------------------------------------------------
  ! variables initialization

  !---------------------------------------------------------------------------------------------
  ! main body

contains

  !include "rhs_navier_stokes.f90"
  include "RHS_3D_navier_stokes.f90"
  include "RHS_2D_navier_stokes.f90"

  !-----------------------------------------------------------------------------
  ! main level wrapper routine to read parameters in the physics module. It reads
  ! from the same ini file as wabbit, and it reads all it has to know. note in physics modules
  ! the parameter struct for wabbit is not available.
  subroutine READ_PARAMETERS_NStokes( filename )
    implicit none

    character(len=*), intent(in) :: filename

    ! inifile structure
    type(inifile) :: FILE
    integer(kind=ik)            :: dF
      

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)


     ! read number_data_fields
    call read_param_mpi(FILE, 'Blocks', 'number_data_fields', params_ns%number_data_fields, 1 )
    !error case: try to solve navier stokes equation with less or more than 5 datafields
    if ( params_ns%number_data_fields /= 5) then
        write(*,'(80("_"))')
        write(*,'("ERROR: try to solve navier stokes equation with", i3, " datafield(s)")') params_ns%number_data_fields
        stop
    end if
    !===================================================================
    ! physics parameter
    !===================================================================
    ! dimension
    call read_param_mpi(FILE, 'Dimensionality', 'dim', params_ns%dim, 2 )
    ! spatial domain size
    call read_param_mpi(FILE, 'DomainSize', 'Lx', params_ns%Lx, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', params_ns%Ly, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Lz', params_ns%Lz, 0.0_rk )
    ! read adiabatic coefficient
    call read_param_mpi(FILE, 'Physics', 'gamma_', params_ns%gamma_, 0.0_rk )
    ! read specific gas constant
    call read_param_mpi(FILE, 'Physics', 'Rs', params_ns%Rs, 0.0_rk )
    ! calculate isochoric heat capacity
    params_ns%Cv = params_ns%Rs/(params_ns%gamma_-1.0_rk)
    ! calculate isobaric heat capacity
    params_ns%Cp = params_ns%Rs*params_ns%gamma_
    ! read prandtl number
    call read_param_mpi(FILE, 'Physics', 'Pr', params_ns%Pr, 0.0_rk )
    ! read dynamic viscosity
    call read_param_mpi(FILE, 'Physics', 'mu0', params_ns%mu0, 0.0_rk )
    ! read switch to turn on|off dissipation
    call read_param_mpi(FILE, 'Physics', 'dissipation', params_ns%dissipation, .true. )

  ! read variable names
    ! allocate names list
    allocate( params_ns%names( params_ns%number_data_fields ) )
    params_ns%names = "---"
    ! read file
    call read_param_mpi(FILE, 'Physics', 'names_ns', params_ns%names, params_ns%names )

    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_ns%discretization, "FD_2nd_central")

    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_ns%N_fields_saved, 1 )
    allocate( params_ns%names(1:params_ns%N_fields_saved))
          ! allocate names list
    allocate( params_ns%names( params_ns%number_data_fields ) )
    params_ns%names = "---"
                ! read file
    call read_param_mpi(FILE, 'Physics', 'names_ns', params_ns%names, params_ns%names )

    call read_param_mpi(FILE, 'Time', 'CFL', params_ns%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'time_max', params_ns%T_end, 1.0_rk   )

    ! set global parameters pF,rohF, UxF etc
    UzF=-1
    do dF = 1, params_ns%number_data_fields
                if ( params_ns%names(dF) == "p" ) pF = dF
                if ( params_ns%names(dF) == "rho" ) rhoF = dF
                if ( params_ns%names(dF) == "Ux" ) UxF = dF
                if ( params_ns%names(dF) == "Uy" ) UyF = dF
                if ( params_ns%names(dF) == "Uz" ) UzF = dF
    end do


    call clean_ini_file_mpi( FILE )
  end subroutine READ_PARAMETERS_NStokes


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
  subroutine PREPARE_SAVE_DATA_NStokes( time, u, g, x0, dx, work )
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
    integer(kind=ik) :: neqn, nwork, Bs, dF

    Bs = size(u,1)-2*g

    ! copy state vector
    !work(:,:,:,1:size(u,4)) = u(:,:,:,:)
    if (size(u,3)==1) then
          work(:,:,:,1) = u(:,:,:,UxF)/u(:,:,:,rhoF)**2 !u
          work(:,:,:,2) = u(:,:,:,UyF)/u(:,:,:,rhoF)**2 !v
    else
        work(:,:,:,1) = u(:,:,:,UxF)/u(:,:,:,rhoF)**2 !u
        work(:,:,:,2) = u(:,:,:,UyF)/u(:,:,:,rhoF)**2 !v
        work(:,:,:,3) = u(:,:,:,UzF)/u(:,:,:,rhoF)**2 !w
    endif

    ! vorticity
    call compute_vorticity(work(:,:,:,1), work(:,:,:,2), work(:,:,:,3), dx, Bs, g, params_ns%discretization,work(:,:,:,4:6))

    ! mask
    !call create_mask_2D_NEW(work(:,:,1,5), x0, dx, Bs, g )

  end subroutine


  !-----------------------------------------------------------------------------
  ! when savig to disk, WABBIT would like to know how you named you variables.
  ! e.g. u(:,:,:,1) is called "ux"
  !
  ! the main routine save_fields has to know how you label the stuff you want to
  ! store from the work array, and this routine returns those strings
  !-----------------------------------------------------------------------------
  subroutine FIELD_NAMES_NStokes( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=80), intent(out) :: name

    if (allocated(params_ns%names)) then
      name = params_ns%names(N)
    else
      call abort(5554,'Something ricked')
    endif

  end subroutine FIELD_NAMES_NStokes


  !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid any an MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine RHS_NStokes( time, u, g, x0, dx, rhs, stage )
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

    ! output. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actuall RHS evaluation.
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

      return

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

      return

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.

      return

    case ("local_stage")
      !-------------------------------------------------------------------------
      ! 4th stage: local evaluation of RHS on all blocks
      !-------------------------------------------------------------------------
      ! the second stage then is what you would usually do: evaluate local differential
      ! operators etc.
      !
      ! called for each block.
      if (size(u,3)==1) then
        !call RHS_2D_navier_stokes(g, Bs, dx(1),dx(2), u(:,:,1,:), rhs)
      else  
        call RHS_3D_navier_stokes(g, Bs, dx(1),dx(2),dx(3), u, rhs)
      endif


    case default
      call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")
    end select


  end subroutine RHS_NStokes

  !-----------------------------------------------------------------------------
  ! subroutine statistics_NStokes()
  !   implicit none
  ! end subroutine


  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_NStokes( time, u, Bs, g, x0, dx, dt )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g, bs

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! the dt for this block is returned to the caller:
    real(kind=rk), intent(out) :: dt

    ! loop variables
    integer(kind=ik)                    :: dF, N_dF, UxF, UyF, UzF
    integer(kind=ik) :: i, ix, iy
    real(kind=rk) :: x,y,unorm,deltax

    ! get smallest spatial seperation 
    deltax=minval(dx)
    dt = 0.0_rk
    if (UzF==-1) then
        unorm = maxval( u(:,:,:,UxF)*u(:,:,:,UxF) + u(:,:,:,UyF)*u(:,:,:,UyF))  
    else
        unorm = maxval( u(:,:,:,UxF)*u(:,:,:,UxF) + u(:,:,:,UyF)*u(:,:,:,UyF)+u(:,:,:,UzF)*u(:,:,:,UzF) )
    endif
    ! max velocity in one block
    if (unorm < 1e-12_rk) unorm = 9e9_rk
    dt = min(dt, params_ns%CFL * deltax / unorm)
  end subroutine GET_DT_BLOCK_NStokes


  !-----------------------------------------------------------------------------
  ! main level wrapper for setting the initial condition on a block
  !-----------------------------------------------------------------------------
  subroutine INICOND_NStokes( time, u, g, x0, dx )
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

    integer(kind=ik) :: ix, iy, Bs,i
    real(kind=rk) :: x,y,c0x,c0y

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    u = 0.0_rk

    ! do i = 1, params_ns%N_scalars
    !   c0x = params_NStokes%x0(i)
    !   c0y = params_NStokes%y0(i)

    !   select case (params_NStokes%inicond(i))
    !   case("blob")
    !     ! create gauss pulse
    !     do ix = g+1,Bs+g
    !       do iy = g+1,Bs+g
    !         ! compute x,y coordinates from spacing and origin
    !         x = dble(ix-(g+1)) * dx(1) + x0(1) - c0x
    !         y = dble(iy-(g+1)) * dx(2) + x0(2) - c0y

    !         if (x<-params_NStokes%Lx/2.0) x = x + params_NStokes%Lx
    !         if (x>params_NStokes%Lx/2.0) x = x - params_NStokes%Lx

    !         if (y<-params_NStokes%Ly/2.0) y = y + params_NStokes%Ly
    !         if (y>params_NStokes%Ly/2.0) y = y - params_NStokes%Ly

    !         ! set actual inicond gauss blob
    !         u(ix,iy,:,i) = dexp( -( (x)**2 + (y)**2 ) / params_NStokes%blob_width(i) )
    !       end do
    !     end do
    !   case default
    !     write(*,*) "errorrroororor"
    !   end select
    !enddo


  end subroutine INICOND_NStokes



end module module_navier_stokes_new
 