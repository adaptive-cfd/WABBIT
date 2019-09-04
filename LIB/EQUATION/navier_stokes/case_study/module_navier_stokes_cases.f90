
!----------------------------------------------------------------
!> Implementation of all the case studies in Navier Stokes
!> \details
!> \author PKrah
!> \date 10.08.2018 - creation of the module
!----------------------------------------------------------------
module module_navier_stokes_cases

    use module_navier_stokes_params
    use module_funnel
    use module_ns_penalization
    use module_simple_geometry
    use module_shock
    use module_pipe_flow


    implicit none
    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE
    !***************************************************************
    ! routines visible outside this module
    !***************************************************************
    PUBLIC :: read_case_parameters, set_inicond_case, get_mask, &
              compute_mask_and_ref2D, compute_mask_and_ref3D
    !***************************************************************
contains

  ! include 'sod_shock_tube.f90'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> In this subroutine you should read in the necessary parameters of your study.
!> Furthermore you can choose a geometry for the volume penalization.
 subroutine read_case_parameters( FILE )
  implicit none
  !--------------------------------------------------------
  type(inifile),intent(inout)   :: FILE         !< filepointer
  !--------------------------------------------------------

  select case ( params_ns%case )
  case ('funnel')
    call read_params_funnel(params_ns,FILE)
  case('simple_geometry')
    call read_params_geometry(params_ns,FILE)
  case('shock_tube')
    call read_params_shock_tube(params_ns,FILE)
  case('pipe_flow')
    call read_params_pipe_flow(params_ns,FILE)
  case('no')

  case default
    call abort(1992132,'Computer says NOOOOOOOOOOOOOO! ~Little Britain~')
  end select

end subroutine read_case_parameters
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !> \brief Set the initial condition of a specific case
 !> If sucessfull set_inidcond returns true, if not then set_inicond is false.
 !> In this case we assume that the case is using a generic initial condition provided in the
 !> inicond function
 !> \details
  function set_inicond_case(params,x0, dx, Bs, g, phi) result(set_inicond)
      implicit none
      ! -----------------------------------------------------------------
      integer(kind=ik), intent(in)  :: g          !< grid parameter
      integer(kind=ik), dimension(3), intent(in) :: Bs
      real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
      real(kind=rk), intent(inout)  :: phi(:,:,:,:)    !< Statevector for t=0
      type(type_params_ns),intent(inout)   :: params    !< NStokes Params structure
      logical :: set_inicond
      ! -----------------------------------------------------------------

      set_inicond=.true.
      select case(params%case)
      case('funnel')
        call set_inicond_funnel(x0, dx, Bs, g, phi )
        return
      case('shock_tube')
        call set_inicond_shock_tube(x0, dx, Bs, g, phi )
        return
      case default
        set_inicond=.false.
      end select

  end function set_inicond_case
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 !==========================================================================
 !> \brief This function computes a penalization term for different geometries
 !! in Navier Stokes physics
 subroutine compute_mask_and_ref2D(params, Bs, g, x0, dx, phi, mask, phi_ref)
     implicit none
     !-------------------------------------------------------
     integer(kind=ik), intent(in)  :: g          !< grid parameter
     integer(kind=ik), dimension(3), intent(in) :: Bs
     real(kind=rk), intent(in)       :: phi(:,:,:)     !< primary state variables
     real(kind=rk), intent(inout)    :: mask(:,:,:)    !< rhs
     real(kind=rk), intent(inout)    :: phi_ref(:,:,:) !< reference state variables
     real(kind=rk), intent(in)       :: x0(2), dx(2)   !< spacing and origin of block
     type(type_params_ns),intent(inout)   :: params    !< NStokes Params structure
     !--------------------------------------------------------
     ! compute mask and reference statevector for volume penalization
     ! term for the different case studies
     select case( params%CASE )
     case('shock_tube')
       call shock_tube_penalization2D(Bs, g, x0, dx, mask, phi_ref)
     case('simple_geometry')
       call geometry_penalization2D(Bs, g, x0, dx, phi(:,:,rhoF), mask, phi_ref)
     case('funnel')
       call funnel_penalization2D(Bs, g, x0, dx, phi, mask, phi_ref)
     case('pipe_flow')
       call pipe_flow_penalization2D(Bs, g, x0, dx, mask, phi_ref)
     case('no')
       return
     case default
       call abort(1201,"ERROR: Come down to earth and tell me what is that:"//params%CASE)
     end select

 end subroutine compute_mask_and_ref2D
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 !==========================================================================
 !> \brief This function computes a penalization term for different geometries
 !! in Navier Stokes physics
 subroutine compute_mask_and_ref3D(params, Bs, g, x0, dx, phi, mask, phi_ref)
     implicit none
     !-------------------------------------------------------
     integer(kind=ik), intent(in)  :: g          !< grid parameter
     integer(kind=ik), dimension(3), intent(in) :: Bs
     real(kind=rk), intent(in)       :: phi(:,:,:,:)     !< primary state variables
     real(kind=rk), intent(inout)    :: mask(:,:,:,:)    !< rhs
     real(kind=rk), intent(inout)    :: phi_ref(:,:,:,:) !< reference state variables
     real(kind=rk), intent(in)       :: x0(3), dx(3)   !< spacing and origin of block
     type(type_params_ns),intent(inout)   :: params    !< NStokes Params structure
     !--------------------------------------------------------

     ! compute mask and reference statevector for volume penalization
     ! term for the different case studies
     select case( params%CASE )
     case('shock_tube')
       call shock_tube_penalization3D(Bs, g, x0, dx, mask, phi_ref)
     case('funnel')
       call funnel_penalization3D(Bs, g, x0, dx, phi, mask, phi_ref)
     case('no')
       return
     case default
       call abort(1201014,"ERROR: Do not know what that is:"//params%CASE)
     end select

 end subroutine compute_mask_and_ref3D
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !> Computes the mask field.
 !> \details Many routines depend on the geometry and need to know the
 !>          mask function. For example for calculation of the force on an obstacle or the
 !>          density within the obstacle, as well as outside the obstacle.
  subroutine get_mask(params,x0, dx, Bs, g, mask, mask_is_colored)
      implicit none
      ! -----------------------------------------------------------------
      integer(kind=ik), intent(in)  :: g          !< grid parameter
      integer(kind=ik), dimension(3), intent(in) :: Bs
      real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
      real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
      type(type_params_ns),intent(inout)   :: params    !< NStokes Params structure
      logical, optional, intent(in) :: mask_is_colored
      logical, save :: is_colored =.false.
      ! -----------------------------------------------------------------
      if( present(mask_is_colored)) is_colored=mask_is_colored

      select case(params%CASE)
      case('pipe_flow')
        call draw_pipe_sponges(mask, x0, dx, Bs, g )
      case('simple_geometry')
        call draw_geometry(x0, dx, Bs, g, mask)
      case('funnel')
        call draw_funnel(x0, dx, Bs, g, mask,is_colored)
      case('shock_tube')
        call draw_simple_shock(mask(:,:,:), x0, dx, Bs, g )
      case('no')
        return
      case default
        call abort(120601,"ERROR: say whaaaaaaaaaaaaaat, don't know that case:"//params%CASE)
      end select

  end subroutine get_mask
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





end module module_navier_stokes_cases
