
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
    use module_shock_tube


    implicit none
    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE
    !***************************************************************
    ! routines visible outside this module
    !***************************************************************
    PUBLIC :: set_case_parameters, set_inicond_case, add_constraints, get_mask
    !***************************************************************
contains

  ! include 'sod_shock_tube.f90'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> In this subroutine you should read in the necessary parameters of your study.
!> Furthermore you can choose a geometry for the volume penalization.
 subroutine set_case_parameters( params , FILE )
  implicit none
  !--------------------------------------------------------
  type(inifile)       ,intent(inout)   :: FILE         !< filepointer
  type(type_params_ns),intent(inout)   :: params    !< NStokes Params structure
  !--------------------------------------------------------

  select case ( params%case )
  case ('funnel')
    call read_params_funnel(params,FILE)
  case('simple_geometry')
    call read_params_geometry(params,FILE)
  case('shock_tube')
    call read_params_shock_tube(params,FILE)
  case('no')

  case default
    call abort(1992132,'Computer says NOOOOOOOOOOOOOO! ~Little Britain~')
  end select

end subroutine set_case_parameters
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
      integer(kind=ik), intent(in)  :: Bs, g        !< grid parameter
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




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !> Use this subroutine to add case specific constraints like
 !>             - boundary conditions
 !>             - volume penalization
 !> to the right hand side
 !> \details
 !> This routine is called in the local stage (i.e. for every block)
 !> after the right hand side was computed

 !==========================================================================
 !> \brief This function computes a penalization term for different geometries
 !! in Navier Stokes physics
 subroutine add_constraints(params, rhs ,Bs , g, x0, dx, phi)
     implicit none
     !-------------------------------------------------------
     integer(kind=ik), intent(in)    :: g, Bs          !< grid parameter
     real(kind=rk), intent(inout)    :: rhs(:,:,:,:)   !< rhs
     real(kind=rk), intent(in)       :: phi(:,:,:,:)   !< state variables
     real(kind=rk), intent(in)       :: x0(3), dx(3)   !< spacing and origin of block
     type(type_params_ns),intent(inout)   :: params    !< NStokes Params structure
     !--------------------------------------------------------
      real(kind=rk), allocatable, save :: penalization(:,:,:,:)

      ! allocates the penalization array if not allready allocated
      call allocate_statevector_ns(penalization, Bs, g)

     ! 1. compute volume penalization term for the different case studies
     select case( params%CASE )
     case('shock_tube')
       call add_shock_tube(penalization, x0, dx, Bs, g, phi )
     case('simple_geometry')
       call add_geometry2D(penalization(:,:,1,:),x0,dx,Bs,g,phi(:,:,1,:))
     case('funnel')
       call add_funnel(penalization, x0, dx, Bs, g, phi)
     case default
       call abort(1201,"ERROR: Come down to earth and tell me what is that:"//params%CASE)
     end select

     ! 2. add penalty to the right hand side
     call add_penalization_term(rhs, penalization, phi)

     ! 3. add boundary conditions (set_boundary conditions) here

 end subroutine add_constraints
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !> Computes the mask field.
 !> \details Many routines depend on the geometry and need to know the
 !>          mask function. For example for calculation of the force on an obstacle or the
 !>          density within the obstacle, as well as outside the obstacle.
  subroutine get_mask(params,x0, dx, Bs, g, mask,mask_is_colored)
      implicit none
      ! -----------------------------------------------------------------
      integer(kind=ik), intent(in)  :: Bs, g        !< grid parameter
      real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
      real(kind=rk), intent(inout), allocatable  :: mask(:,:,:)    !< mask function
      type(type_params_ns),intent(inout)   :: params    !< NStokes Params structure
      logical, optional, intent(in) :: mask_is_colored
      logical, save :: is_colored =.false.
      ! -----------------------------------------------------------------
      if( present(mask_is_colored)) is_colored=mask_is_colored

      select case(params%CASE)
      case('simple_geometry')
        call draw_geometry(x0, dx, Bs, g, mask)
      case('funnel')
        call draw_funnel(x0, dx, Bs, g, mask,is_colored)
      case('shock_tube')
        call draw_simple_shock(mask(:,:,1), x0, dx, Bs, g,'boundary' )
      case default
        call abort(120601,"ERROR: say whaaaaaaaaaaaaaat, don't know that case:"//params%CASE)
      end select

  end subroutine get_mask
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





end module module_navier_stokes_cases
