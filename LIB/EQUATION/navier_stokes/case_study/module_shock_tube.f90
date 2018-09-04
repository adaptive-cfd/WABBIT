!----------------------------------------------------------------
!> This module implements different shock tests to validate the
!> euler equations of NS and to test the filters in a multiresolution
!> framework.
!> \details
!> \author PKrah
!> \date 10.08.2018 - creation of the module
!----------------------------------------------------------------

module module_shock_tube

    use module_navier_stokes_params
    use module_precision
    use module_ini_files_parser_mpi
    use module_ns_penalization

    implicit none

    !**********************************************************************************************
    ! make everything private if not explicitly marked public
    PRIVATE
    character(len=80),                 save        :: SHOCK_TYPE
    !**********************************************************************************************
    ! only this functions are visible outside this module
    PUBLIC :: init_shock_tube,add_shock_tube,draw_simple_shock
    !**********************************************************************************************

  !> \file
  !> \details
  !> Available shock tubes
  !! ------------------------
  !!
  !!    | type                | short description                                  |
  !!    |---------------------|----------------------------------------------------|
  !!    |  sod's shock tube |   see publication  <a href="https://doi.org/10.1016/0021-9991(78)90023-2">Sod,G. A (1978)</a>                          |
  !!    |  moving shock     |  shock moves in positiv x direction using Rankine-Hugoniot conditions |
  !!     --------------------------------------------------------------------------------


  type :: type_shock_params
      real(kind=rk) :: rho_left,u_left,p_left
      real(kind=rk) :: rho_right,u_right,p_right
  end type type_shock_params



  !------------------------------------------------
  type(type_shock_params) , save :: shock_params
  !------------------------------------------------

contains


  include "sod_shock_tube.f90"
  include "simple_shock.f90"



!> \brief reads parameters for mask function from file
subroutine init_shock_tube( params,FILE )
    use module_navier_stokes_params
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params

    call init_simple_sponge(params,FILE)
    SHOCK_TYPE=params%geometry
    select case(SHOCK_TYPE)
    case ('moving-shock')
      call init_simple_shock(params,FILE)
    case('sod_shock_tube')
      ! nothing to do<
    case default
      call abort(8546501,"[module_ns_penalization.f90] ERROR: geometry for VPM is unknown"//SHOCK_TYPE)

    end select

end subroutine init_shock_tube





!==========================================================================
!> This function adds a penalization term
!> to navier stokes equations
subroutine add_shock_tube(penalization, x0, dx, Bs, g, phi)
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

      select case( SHOCK_TYPE )
      case('sod_shock_tube')
        call add_sod_shock_tube(penalization(:,:,1,:), x0, dx, Bs, g, phi(:,:,1,:) )
      case('moving-shock')
        call add_simple_shock(penalization(:,:,1,:),x0,dx,Bs,g,phi(:,:,1,:))
      case default
        call abort(12071,"ERROR: geometry for VPM is unknown"//SHOCK_TYPE)
      end select


end subroutine add_shock_tube
!==========================================================================








end module module_shock_tube
