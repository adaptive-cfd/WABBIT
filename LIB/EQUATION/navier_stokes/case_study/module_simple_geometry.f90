
!-----------------------------------------------------------------
!> Implementation of simple penalized geometries (cylinder,triangle,rhombus)
!> \details
!> \version 23.2.2018
!> \author P.Krah
!-----------------------------------------------------------------

module module_simple_geometry

  use module_navier_stokes_params
  use module_precision
  use module_ini_files_parser_mpi
  use module_ns_penalization

  implicit none

  !**********************************************************************************************
  ! make everything private if not explicitly marked public
  PRIVATE
  !**********************************************************************************************
  ! only this functions are visible outside this module
  PUBLIC :: init_geometry,add_geometry,draw_geometry
  !**********************************************************************************************

  character(len=80),                 save        :: MASK_GEOMETRY
    !------------------------------------------------
  !> \file
  !> \details
  !> Available mask geometrys
  !! ------------------------
  !!
  !!    | geometry    | geometrical properties                             |
  !!    |-------------|----------------------------------------------------|
  !!    | \b cylinder |   \c x_cntr(1:2),\c radius                         |
  !!    | \b triangle | \c x_cntr, \c length, \c angle                     |
  !!    | \b rhombus  | is a double triangle                               |
  !!     -------------------------------------------------------------------
  type :: type_cylinder
      real(kind=rk), dimension(2) :: x_cntr
      real(kind=rk)               :: radius
  end type type_cylinder


  type :: type_triangle
      real(kind=rk), dimension(3) :: x_cntr
      real(kind=rk)               :: length
      real(kind=rk)               :: angle
      logical                     :: rhombus
  end type type_triangle

  type :: type_shock_params
      real(kind=rk) :: rho_left,u_left,p_left
      real(kind=rk) :: rho_right,u_right,p_right
  end type type_shock_params



  !------------------------------------------------
  type(type_cylinder) , save :: cyl
  type(type_triangle) , save :: triangle
  type(type_shock_params) , save :: shock_params
  !------------------------------------------------

contains

include "vortex_street.f90"
include "triangle.f90"



!> \brief reads parameters for mask function from file
subroutine init_geometry( params,FILE )
    use module_navier_stokes_params
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params

    ! add sponge parameters for the in and outflow of the domain
    call init_simple_sponge(params,FILE)
    MASK_GEOMETRY=params%geometry
    ! geometry of the solid obstacle
    select case(mask_geometry)
    case ('triangle','rhombus')
      call init_triangle(params,FILE)
    case ('vortex_street','cylinder')
      call init_vortex_street(FILE)
    case default
      call abort(8546501," ERROR: geometry for VPM is unknown"//mask_geometry)
    end select

end subroutine init_geometry






!==========================================================================
!> This function adds a penalization term
!> to navier stokes equations
subroutine add_geometry(penalization, x0, dx, Bs, g, phi)
      implicit none
      !-------------------------------------------------------
      !> grid parameter
      integer(kind=ik), intent(in)    :: g, Bs
      !> rhs
      real(kind=rk), intent(inout)    :: penalization(:,:,:,:)
      !> state variables
      real(kind=rk), intent(in)       :: phi(:,:,:,:)
      !> spacing and origin of block
      real(kind=rk), intent(in)       :: x0(3), dx(3)
      !--------------------------------------------------------

      select case( MASK_GEOMETRY )
      case('vortex_street','cylinder')
        call add_cylinder(penalization(:,:,1,:), x0, dx, Bs, g, phi(:,:,1,:) )
      case('rhombus','triangle')
        call add_triangle(penalization(:,:,1,:),x0,dx,Bs,g,phi(:,:,1,:))
      case default
        call abort(120401,"ERROR: geometry for VPM is unknown"//mask_geometry)
      end select


end subroutine add_geometry
!==========================================================================




!==========================================================================
subroutine draw_geometry(x0, dx, Bs, g, mask)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g        !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
    real(kind=rk), intent(inout), allocatable  :: mask(:,:,:)    !< mask function
    ! -----------------------------------------------------------------

    select case(MASK_GEOMETRY)
    case('vortex_street','cylinder')
      call draw_cylinder(mask(:,:,1), x0, dx, Bs, g )
    case('rhombus','triangle')
      call draw_triangle(mask(:,:,1),x0,dx,Bs,g)
    case default
      call abort(120601,"ERROR: geometry is unknown"//mask_geometry)
    end select

end subroutine draw_geometry
!==========================================================================





end module module_simple_geometry
