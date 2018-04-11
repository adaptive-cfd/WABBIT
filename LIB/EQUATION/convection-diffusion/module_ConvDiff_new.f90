 !> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_convdiff.f90
!> \version 0.5
!> \author engels
!!
!! \brief module for 2D/3D convdiff physics
!!
!!
!! = log ======================================================================================
!! \n
!!
! ********************************************************************************************

module module_convdiff_new

  !---------------------------------------------------------------------------------------------
  ! modules

  use module_precision
  ! ini file parser module, used to read parameters. note: in principle, you can also
  ! just use any reader you feel comfortable with, as long as you can read the parameters
  ! from a file.
  use module_ini_files_parser_mpi
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
  PUBLIC :: READ_PARAMETERS_convdiff, PREPARE_SAVE_DATA_convdiff, RHS_convdiff, GET_DT_BLOCK_convdiff, INICOND_convdiff, FIELD_NAMES_convdiff
  !**********************************************************************************************

  ! user defined data structure for time independent parameters, settings, constants
  ! and the like. only visible here.
  type :: type_paramsb
    real(kind=rk) :: CFL, T_end
    real(kind=rk) :: Lx, Ly, Lz
    real(kind=rk), allocatable, dimension(:) :: nu, u0x,u0y,u0z,blob_width,x0,y0,z0
    integer(kind=ik) :: dim, N_scalars, N_fields_saved
    character(len=80), allocatable :: names(:), inicond(:), velocity(:)
    character(len=80) :: discretization
  end type type_paramsb

  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.
  type(type_paramsb), save :: params_convdiff



  !---------------------------------------------------------------------------------------------
  ! variables initialization

  !---------------------------------------------------------------------------------------------
  ! main body

contains

  include "rhs_convdiff.f90"

  !-----------------------------------------------------------------------------
  ! main level wrapper routine to read parameters in the physics module. It reads
  ! from the same ini file as wabbit, and it reads all it has to know. note in physics modules
  ! the parameter struct for wabbit is not available.
  subroutine READ_PARAMETERS_convdiff( filename )
    implicit none

    character(len=*), intent(in) :: filename

    ! inifile structure
    type(inifile) :: FILE

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'N_scalars', params_convdiff%N_scalars, 1  )
    allocate( params_convdiff%nu(1:params_convdiff%N_scalars))
    allocate( params_convdiff%u0x(1:params_convdiff%N_scalars))
    allocate( params_convdiff%u0y(1:params_convdiff%N_scalars))
    allocate( params_convdiff%u0z(1:params_convdiff%N_scalars))

    allocate( params_convdiff%x0(1:params_convdiff%N_scalars))
    allocate( params_convdiff%y0(1:params_convdiff%N_scalars))
    allocate( params_convdiff%z0(1:params_convdiff%N_scalars))

    allocate( params_convdiff%inicond(1:params_convdiff%N_scalars))
    allocate( params_convdiff%velocity(1:params_convdiff%N_scalars))

    allocate( params_convdiff%blob_width(1:params_convdiff%N_scalars))

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'nu', params_convdiff%nu )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'u0x', params_convdiff%u0x )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'u0y', params_convdiff%u0y )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'u0z', params_convdiff%u0z )

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'x0', params_convdiff%x0 )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'y0', params_convdiff%y0 )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'z0', params_convdiff%z0 )

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'blob_width', params_convdiff%blob_width )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'inicond', params_convdiff%inicond, (/'gauss_blob'/) )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'velocity', params_convdiff%velocity, (/'constant'/) )

    call read_param_mpi(FILE, 'Dimensionality', 'dim', params_convdiff%dim, 2 )
    call read_param_mpi(FILE, 'DomainSize', 'Lx', params_convdiff%Lx, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', params_convdiff%Ly, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Lz', params_convdiff%Lz, 0.0_rk )

    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_convdiff%discretization, "FD_2nd_central")

    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_convdiff%N_fields_saved, 1 )
    allocate( params_convdiff%names(1:params_convdiff%N_fields_saved))
    call read_param_mpi(FILE, 'Saving', 'field_names', params_convdiff%names, (/"phi1","phi2","phi3"/) )


    call read_param_mpi(FILE, 'Time', 'CFL', params_convdiff%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'time_max', params_convdiff%T_end, 1.0_rk   )

    call clean_ini_file_mpi( FILE )

    if ( params_convdiff%dim == 2) params_convdiff%u0z=0.0_rk
  end subroutine READ_PARAMETERS_convdiff


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
  subroutine PREPARE_SAVE_DATA_convdiff( time, u, g, x0, dx, work )
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

    Bs = size(u,1)-2*g

    ! copy state vector
    work(:,:,:,1:size(u,4)) = u(:,:,:,:)

    call create_velocity_field_2d( time, g, Bs, dx, x0, work(:,:,1,2:3), 1 )

  end subroutine


  !-----------------------------------------------------------------------------
  ! when savig to disk, WABBIT would like to know how you named you variables.
  ! e.g. u(:,:,:,1) is called "ux"
  !
  ! the main routine save_fields has to know how you label the stuff you want to
  ! store from the work array, and this routine returns those strings
  !-----------------------------------------------------------------------------
  subroutine FIELD_NAMES_convdiff( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=80), intent(out) :: name

    if (allocated(params_convdiff%names)) then
      name = params_convdiff%names(N)
    else
      call abort(5554,'Something ricked')
    endif

  end subroutine FIELD_NAMES_convdiff


  !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid any an MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine RHS_convdiff( time, u, g, x0, dx, rhs, stage )
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

      call RHS_convdiff_new(time, g, Bs, dx, x0, u, rhs)


    case default
      call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")
    end select


  end subroutine RHS_convdiff

  !-----------------------------------------------------------------------------
  ! subroutine statistics_convdiff()
  !   implicit none
  ! end subroutine


  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_convdiff( time, u, Bs, g, x0, dx, dt )
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

    real(kind=rk) :: u0(1:Bs+2*g, 1:Bs+2*g, 1:2)
    integer(kind=ik) :: i, ix, iy
    real(kind=rk) :: x,y,unorm

    dt = 9.9e9_rk


    do i = 1, params_convdiff%N_scalars
      call create_velocity_field_2d( time, g, Bs, dx, x0, u0, i )

      unorm = maxval( u0(:,:,1)*u0(:,:,1) + u0(:,:,2)*u0(:,:,2) )

      if ( unorm < 1.0e-5_rk) then
        ! if the value of u is very small, which may happen if it is time dependent
        ! we choose some fixed value in order not to miss the instant when u becomes
        ! large again.
        dt = min( 1.0e-3_rk, params_convdiff%CFL * dx(1) / sqrt(unorm) )
      else
        dt = min(params_convdiff%CFL * dx(1) / sqrt(unorm), dt)
      endif

    enddo


  end subroutine GET_DT_BLOCK_convdiff


  !-----------------------------------------------------------------------------
  ! main level wrapper for setting the initial condition on a block
  !-----------------------------------------------------------------------------
  subroutine INICOND_convdiff( time, u, g, x0, dx )
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

    integer(kind=ik) :: ix, iy, iz, Bs,i
    real(kind=rk) :: x,y,c0x,c0y,z,c0z

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    u = 0.0_rk

    do i = 1, params_convdiff%N_scalars
      c0x = params_convdiff%x0(i)
      c0y = params_convdiff%y0(i)
      c0z = params_convdiff%z0(i)

      select case (params_convdiff%inicond(i))
      case("blob")
          if (params_convdiff%dim==2) then
            ! create gauss pulse
            do ix = g+1,Bs+g
              do iy = g+1,Bs+g
                ! compute x,y coordinates from spacing and origin
                x = dble(ix-(g+1)) * dx(1) + x0(1) - c0x
                y = dble(iy-(g+1)) * dx(2) + x0(2) - c0y

                if (x<-params_convdiff%Lx/2.0) x = x + params_convdiff%Lx
                if (x>params_convdiff%Lx/2.0) x = x - params_convdiff%Lx

                if (y<-params_convdiff%Ly/2.0) y = y + params_convdiff%Ly
                if (y>params_convdiff%Ly/2.0) y = y - params_convdiff%Ly

                ! set actual inicond gauss blob
                u(ix,iy,:,i) = dexp( -( (x)**2 + (y)**2 ) / params_convdiff%blob_width(i) )
              end do
            end do
        else
            ! create gauss pulse
            do ix = g+1,Bs+g
              do iy = g+1,Bs+g
                  do iz = g+1,Bs+g
                    ! compute x,y coordinates from spacing and origin
                    x = dble(ix-(g+1)) * dx(1) + x0(1) - c0x
                    y = dble(iy-(g+1)) * dx(2) + x0(2) - c0y
                    z = dble(iz-(g+1)) * dx(3) + x0(3) - c0z

                    if (x<-params_convdiff%Lx/2.0) x = x + params_convdiff%Lx
                    if (x>params_convdiff%Lx/2.0) x = x - params_convdiff%Lx

                    if (y<-params_convdiff%Ly/2.0) y = y + params_convdiff%Ly
                    if (y>params_convdiff%Ly/2.0) y = y - params_convdiff%Ly

                    if (z<-params_convdiff%Lz/2.0) z = z + params_convdiff%Lz
                    if (z>params_convdiff%Lz/2.0) z = z - params_convdiff%Lz

                    ! set actual inicond gauss blob
                    u(ix,iy,iz,i) = dexp( -( (x)**2 + (y)**2 + (z)**2 ) / params_convdiff%blob_width(i) )
                end do
              end do
            end do
        end if
      case default
        call abort(72637,"Error. Inital conditon for conv-diff is unkown: "//trim(adjustl(params_convdiff%inicond(i))))
      end select


    enddo


  end subroutine INICOND_convdiff



end module module_convdiff_new
