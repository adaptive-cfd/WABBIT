!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_ConvDiff.f90
!> \version 0.5
!> \author sm
!!
!! \brief module for 2D/3D acm physics
!!
!!
!! = log ======================================================================================
!! \n
!!
! ********************************************************************************************

module module_ConvDiff_new

!---------------------------------------------------------------------------------------------
! modules

    use module_precision
    ! ini file parser module, used to read parameters. note: in principle, you can also
    ! just use any reader you feel comfortable with, as long as you can read the parameters
    ! from a file.
    use module_ini_files_parser_mpi
!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    !**********************************************************************************************
    PUBLIC :: READ_PARAMETERS_ConvDiff, PREPARE_SAVE_DATA_ConvDiff, RHS_ConvDiff, GET_DT_BLOCK_ConvDiff, INICOND_ConvDiff, FIELD_NAMES_ConvDiff
    !**********************************************************************************************

    ! user defined data structure for time independent parameters, settings, constants
    ! and the like. only visible here.
    type :: type_params
      ! how many scalar fields do you want to solve?
      integer(kind=ik) :: N_scalars=1
      ! note you need to specify one value per scalar field:
      real(kind=rk), allocatable :: nu(:)
      real(kind=rk), allocatable :: u0x(:)
      real(kind=rk), allocatable :: u0y(:)
      real(kind=rk), allocatable :: u0z(:)
      character(len=80), allocatable :: inicond(:)
      real(kind=rk), allocatable :: blob_width(:)
      ! variable names
      character(len=80), allocatable :: names(:)
    end type type_params

    ! parameters for this module. they should not be seen outside this physics module
    ! in the rest of the code. WABBIT does not need to know them.
    type(type_params), save :: params

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

  !-----------------------------------------------------------------------------
  ! main level wrapper routine to read parameters in the physics module. It is
  ! called in
  subroutine READ_PARAMETERS_ConvDiff( file )
    implicit none
    character(len=*), intent(in) :: file

    ! local variables
    ! inifile structure
    type(inifile) :: FILE
    ! abbreviation:
    integer(kind=ik) :: N

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'N_scalars', N, 1 )
    params%N_scalars = N
    allocate( params%nu(1:N), params%u0x(1:N), params%u0y(1:N), params%u0z(1:N) )
    allocate( params%blob_width(1:N), params%inicond(1:N), params%names(1:N) )

    call clean_ini_file_mpi( FILE )
  end subroutine READ_PARAMETERS_ConvDiff

  !-----------------------------------------------------------------------------

  subroutine SAVE_DATA_ConvDiff()
    implicit none
  end subroutine

  !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid any an MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine RHS_ConvDiff( time, u, g, x0, dx, rhs )
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
    ! output
    real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

    rhs = 0.0

  end subroutine RHS_ConvDiff

  !-----------------------------------------------------------------------------
  ! subroutine statistics_ConvDiff()
  !   implicit none
  ! end subroutine


  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_ConvDiff( time, u, g, x0, dx, dt )
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
    ! the dt for this block is returned to the caller:
    real(kind=rk), intent(out) :: dt

    dt = 0.0
  end subroutine GET_DT_BLOCK_ConvDiff

  !-----------------------------------------------------------------------------
  ! main level wrapper for setting the initial condition on a block
  !-----------------------------------------------------------------------------
  subroutine INICOND_ConvDiff( time, u, g, x0, dx )
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


  end subroutine INICOND_ConvDiff

  !-----------------------------------------------------------------------------
  ! when savig to disk, WABBIT would like to know how you named you variables.
  ! e.g. u(:,:,:,1) is called "ux"
  !-----------------------------------------------------------------------------
  subroutine FIELD_NAMES_ConvDiff( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=80), intent(out) :: name

    name = params%names(N)

  end subroutine FIELD_NAMES_ConvDiff

end module module_ConvDiff_new
