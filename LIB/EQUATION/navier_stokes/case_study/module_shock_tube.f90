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
    !**********************************************************************************************
    ! only this functions are visible outside this module
    PUBLIC :: read_params_shock_tube,add_shock_tube,draw_simple_shock, set_inicond_shock_tube, &
              set_inicond_moving_shock,moving_shockVals, shockVals
    !**********************************************************************************************

  !> \file
  !> \details
  !> Available shock tubes
  !! ------------------------
  !!
  !!    | type                | short description                                  |
  !!    |---------------------|----------------------------------------------------|
  !!    |  sod's shock tube |   see publication  <a href="https://doi.org/10.1016/0021-9991(78)90023-2">Sod,G. A (1978)</a>                          |
  !!    |  moving shock     |  shock moves in positiv x direction  (Rankine-Hugoniot conditions) |
  !!    |  standing shock   |  shockfront is not moving   (Rankine-Hugoniot conditions) |
  !!     --------------------------------------------------------------------------------


  type :: type_shock_params
      character(len=80) :: name
      real(kind=rk)     :: rho_left,u_left,p_left
      real(kind=rk)     :: rho_right,u_right,p_right
      real(kind=rk)     :: machnumber
  end type type_shock_params

  !---------------------------------------
  type(type_shock_params) , save :: shock
  !---------------------------------------

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> \brief reft and right shock values for 1D shock moving with mach to the right
  !> \detail This function converts with the Rankine-Hugoniot Conditions
  !>  values \f$\rho_L,p_L,Ma\f$ to the values of the right of the shock
  !>  \f$\rho_R,u_R,p_R\f$ and \f$u_L\f$ .
  !> See: formula 3.51-3.56 in Riemann Solvers and Numerical Methods for Fluid Dynamics
  !> author F.Toro
  subroutine moving_shockVals(rhoL,uL,pL,rhoR,uR,pR,gamma,mach)
      implicit none
      !> one side of the shock (density, pressure)
      real(kind=rk), intent(in)      ::rhoL,pL
      !> shock speed
      real(kind=rk), intent(in)      :: mach
      !> speed on
      real(kind=rk), intent(inout)      :: uL
      !> other side of the shock (density, velocity, pressure)
      real(kind=rk), intent(out)      ::rhoR,uR,pR
      !> heat capacity ratio
      real(kind=rk), intent(in)      ::gamma

      real(kind=rk)                ::c_R


       uR    =   0
       rhoR  =   ((gamma-1)*mach**2+2)/((gamma+1)*mach**2)*rhoL
       pR    = (gamma+1)/(2*gamma*mach**2-gamma+1)*pL
       c_R   = sqrt(gamma*pR/rhoR)
       uL    = (1-rhoR/rhoL)*mach*c_R
  end subroutine moving_shockVals
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> \brief This function calculates from \f$\rho_1,u_1,p_1\f$
  !> values \f$\rho_2,u_2,p_2\f$ on the ohter side
  !> of the shock
  subroutine shockVals(rho1,u1,p1,rho2,u2,p2,gamma)
      implicit none
      !> one side of the shock (density, velocity, pressure)
      real(kind=rk), intent(in)      ::rho1,u1,p1
      !> other side of the shock (density, velocity, pressure)
      real(kind=rk), intent(out)      ::rho2,u2,p2
      !> heat capacity ratio
      real(kind=rk), intent(in)      ::gamma
      real(kind=rk)                ::cstar_sq

      cstar_sq = 2*(gamma-1)/(gamma+1)*( p1/rho1*(gamma/(gamma-1))+u1**2/2 ) ;
      !sqrt(cstar_sq)
      u2 = cstar_sq /u1;
      rho2 = (rho1*u1)/u2;
      p2= (p1+ rho1*u1**2 )-rho2*u2**2;
  end subroutine shockVals
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \brief reads parameters for mask function from file
subroutine read_params_shock_tube( params,FILE )
    use module_navier_stokes_params
    implicit none
    !-------------------------------------------------------------------
    type(inifile) ,intent(inout)       :: FILE    !< pointer to inifile
    type(type_params_ns),intent(inout)  :: params !< params structure of navier stokes
    !-------------------------------------------------------------------

    call read_param_mpi(FILE, 'Shock_Tube', 'name', shock%name, 'no_shock' )

    select case(shock%name)
    case ('moving-shock')
      call read_param_mpi(FILE, 'Shock_Tube', 'pressure_left', shock%p_left,1.0_rk )
      call read_param_mpi(FILE, 'Shock_Tube', 'velocity_left', shock%u_left, 3.0_rk )
      call read_param_mpi(FILE, 'Shock_Tube', 'density_left', shock%rho_left, 2.0_rk )
      call read_param_mpi(FILE, 'Shock_Tube', 'shock_wave_speed', shock%machnumber, 2.0_rk )
      ! speed of sound
      params%c0 =sqrt(params%gamma_ * shock%p_left / shock%rho_left)
      ! shock wave speed as mach number
      shock%machnumber    =shock%machnumber / params%c0
      ! return the shock wave speed to Navier Stokes module (only for logfiles additional info)
      params%machnumber=shock%machnumber
      ! from the initial values we calculate the values of the other side of the shock
      call moving_shockVals(shock%rho_left ,shock%u_left ,shock%p_left , &
                            shock%rho_right,shock%u_right,shock%p_right, &
                            params%gamma_,shock%machnumber)
    case ('standing-shock')
        call read_param_mpi(FILE, 'Shock_Tube', 'pressure_left', shock%p_left,1.0_rk )
        call read_param_mpi(FILE, 'Shock_Tube', 'velocity_left', shock%u_left, 3.0_rk )
        call read_param_mpi(FILE, 'Shock_Tube', 'density_left', shock%rho_left, 2.0_rk )
        ! from the initial values we calculate the values of the other side of the shock
        call shockVals( shock%rho_left,shock%u_left,shock%p_left, &
                        shock%rho_right,shock%u_right,shock%p_right,params%gamma_)
    case('sod_shock_tube')
      ! Sods test case: shock tube
      ! ---------------------------
      ! following values are imposed and smoothed with tangens:
      ! ---------------------------------
      !   rhoL    | rhoR         | rhoL
      !   uL      | uR           | uL
      !   pL      | pR           | pL
      ! 0------------------------------Lx
      !          Lx/2
      ! Test case for shock capturing filter
      ! The initial condition is devided into
      ! Left part x<= Lx/2
      shock%rho_left=1
      shock%p_left  =1
      shock%u_left  =0
      !
      ! Rigth part x> Lx/2
      shock%rho_right=0.125
      shock%p_right  =0.1
      shock%u_right  =0
    case default
      call abort(8546501,"ERROR: Mr. Spock: Insufficient facts always invite danger. "//shock%name)

    end select

end subroutine read_params_shock_tube
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> \brief Set the initial condition of a specific case
  !> \details
   subroutine set_inicond_shock_tube(x0, dx, Bs, g, u)
       implicit none
       ! -----------------------------------------------------------------
       integer(kind=ik), intent(in)  :: Bs, g        !< grid parameter
       real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
       real(kind=rk), intent(inout)  :: u(:,:,:,:)    !< Statevector for t=0
       ! -----------------------------------------------------------------
       real(kind=rk),allocatable:: mask(:,:,:)
       real(kind=rk)            :: x,p_init, rho_init,u_init(3),T_init
       integer(kind=ik)         :: ix


       p_init    =params_ns%initial_pressure
       rho_init  =params_ns%initial_density
       u_init    =params_ns%initial_velocity
       T_init    =params_ns%initial_temp

       select case(shock%name)
       case ('moving-shock','standing-shock')

       case('sod_shock_tube')
         ! Sods test case: shock tube
         ! ---------------------------
         !
         ! Test case for shock capturing filter
         ! The initial condition is devided into
         ! Left part x<= Lx/2
         !
         ! rho=1
         ! p  =1
         ! u  =0
         !
         ! Rigth part x> Lx/2
         !
         ! rho=0.125
         ! p  =0.1
         ! u  =0
         do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)
           call continue_periodic(x,params_ns%domain_size(1))
           if (x <= params_ns%domain_size(1)*0.5_rk) then
             ! left domain half
             u( ix, :, :, rhoF) = sqrt(shock%rho_left)
             u( ix, :, :, pF)   = shock%p_left
             u( :, :, :, UxF)   = shock%u_left
           else
             ! right domain half
             u( ix, :, :, rhoF) = sqrt(shock%rho_right)
             u( ix, :, :, pF)   = shock%p_right
             u( :, :, :, UxF)   = shock%u_right
           endif
         end do
         ! velocity set to 0
         u( :, :, :, UyF) = 0.0_rk
       case default
         call abort(8546501,"ERROR: Mr. Spock: Insufficient facts always invite danger. "//shock%name)
       end select

    end subroutine set_inicond_shock_tube
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





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

      call add_simple_shock(penalization(:,:,1,:),x0,dx,Bs,g,phi(:,:,1,:))

end subroutine add_shock_tube
!==========================================================================




!==========================================================================
!> \brief Compute mask function of sod shock tube
!> \detail
!>      +For boundary_type=='left'
!>             -> mask will be generatet for 0<x<0.1*Lx
!>      +For boundary_type=='rigth'
!>             -> mask will be generatet for 0.9Lx<x<Lx
!>      +For boundary_type else
!>             -> mask will be generatet for both sides
subroutine draw_simple_shock(mask, x0, dx, Bs, g, boundary_type )


    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx
    !> boundary_type={'left','rigth'}
     character(len=*),          intent(in)                    :: boundary_type
    ! auxiliary variables
    real(kind=rk)                                             :: x, h, x_boundary(2),boundary_width
    ! loop variables
    integer(kind=ik)                                          :: ix

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk
    ! left and right boundary of shock tube
    x_boundary(1)   =0.1_rk*params_ns%domain_size(1)
    x_boundary(2)   =params_ns%domain_size(1)*0.9_rk
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)

              mask(ix,:) = smoothstep(x-x_boundary(1),h) &
                         + smoothstep(x_boundary(2)-x,h)

       end do

end subroutine draw_simple_shock
!==========================================================================



!==========================================================================
!> \brief Adds penalization terms for a simple shock_tube
!> \detail
!>        A mask will be generatet at 0<x<0.1*Lx and 0.9Lx<x<Lx and
!>        the left and right values of the shock tube will be assigned to
!>        the penalized domains.
!>        Be aware that this function uses the C_sp_inv parameter globaly defined
!>        when importing the module_penalization.
subroutine add_simple_shock(penalization, x0, dx, Bs, g ,phi)

    implicit none
    ! grid
    integer(kind=ik), intent(in)                     :: Bs, g
    !> penalization term including mask
    real(kind=rk), dimension(:,:,:), intent(inout)   :: penalization
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)          :: x0, dx
    !> statevector
    real(kind=rk), dimension(:,:,:), intent(in)      :: phi
      !> statevector
    real(kind=rk)                                    :: mask

    ! auxiliary variables
    real(kind=rk)                                    :: x, y, h,domain_size(3)
    ! preasure,density velocities
    real(kind=rk)                                    :: p(Bs+2*g,Bs+2*g),rho(Bs+2*g,Bs+2*g), &
                                                        u(Bs+2*g,Bs+2*g),v(Bs+2*g,Bs+2*g)
    ! loop variables
    integer(kind=ik)                                 :: ix, iy,n
    ! left and right boundary
    real(kind=rk)                                    :: u_ref,u_R,u_L,rho_R,rho_L,p_L,p_R, &
                                                        x_L,x_R,rho_ref,p_ref

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(penalization,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask         = 0.0_rk
    penalization = 0.0_rk

    ! parameter for smoothing function (width)
    h       = 1.5_rk*max(dx(1), dx(2))

    rho         = phi(:,:,1)**2
    u           = phi(:,:,2)/phi(:,:,1)
    v           = phi(:,:,3)/phi(:,:,1)
    p           = phi(:,:,4)
    domain_size = params_ns%domain_size

    ! left boundary
    x_L      =0.1_rk*domain_size(1)
    p_L      = shock%p_left
    rho_L    =  shock%rho_left
    u_L      =  shock%u_left

    ! right boundary
    x_R      =domain_size(1)*0.9_rk
    p_R      = shock%p_right
    rho_R    =  shock%rho_right
    u_R      =  shock%u_right

    ! parameter for smoothing function (width)

    do ix=1, Bs+2*g
      x = dble(ix-(g+1)) * dx(1) + x0(1)
      if (x<domain_size(1)*0.5_rk) then
        ! values of the left sponge domain
        mask      = smoothstep(x-x_L,h)
        rho_ref   = rho_L
        p_ref     = p_L
        u_ref   = u_L
      else
        mask      = smoothstep(x_R-x,h)
        ! Because WABBIT provides only periodic BC, we have to ensure a smooth transition
        ! from the sponge values on the left side of the domain (rho_L, p_L ,u_L)
        ! to the values on the right side (rho_R, p_R, u_R),
        ! For this reason we use the transition function for assigning the values
        ! inside the right sponge domain.
        ! The transition is located in the right sponge (x>x_R) and starts
        ! from x0=0.925_rk*domain_size(1) to x0+widht=(0.925_rk+0.05)*domain_size(1)
        rho_ref   = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),rho_R,rho_L)
        p_ref     = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),p_R  ,p_L    )
        u_ref     = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),u_R  ,u_L    )
      endif

      ! density
      penalization(ix,:,1)= penalization(ix,:,1) + C_sp_inv*mask * ( rho(ix,:) - rho_ref )
      ! x-velocity
      penalization(ix,:,2)= penalization(ix,:,2) + C_sp_inv*mask * ( rho(ix,:)*u(ix,:)- rho_ref*u_ref )
      ! y-velocity
      penalization(ix,:,3)= penalization(ix,:,3) + C_sp_inv*mask * ( rho(ix,:)*v(ix,:) )
      ! preasure
      penalization(ix,:,4)= penalization(ix,:,4) + C_sp_inv*mask *( p(ix,:) - p_ref )
    end do



end subroutine add_simple_shock



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \brief
!>    + Sets the initial condition of a moving shock wave along a one axis, using the Rankine-Hugoniot
!>      conditions.
!>    + If shock\_speed is not passed to the routine, the shock wave is a standing-shock
!>    + example: call set_inicond_moving_shock(x0(2), dx(2), Bs, g, phi(:,1,:,:), (/1, 1, 0/), (/2, 3, 0/), 10)
!>               for a shock wave along the y axis
!> \details
subroutine set_inicond_moving_shock(x0, dx, Bs, g, phi, phi_left, phi_right)
     implicit none
     ! -----------------------------------------------------------------
     integer(kind=ik), intent(in)       :: Bs, g        !< grid parameter
     real(kind=rk), intent(in)          :: x0(3), dx(3) !< coordinates of block and block spacinf
     real(kind=rk), intent(inout)       :: phi(:,:)    !< Statevector in skew-symetric form : phi([ix|iy|iz],dF)
     !> Statevector values for the left and right side of the shock
     !>       + phi_[left/right](1) density at the left or right
     !>       + phi_[left/right](2) velocity at the left or right
     !>       + phi_[left/right](3) pressure at the left or right
     real(kind=rk), intent(in)          :: phi_left(3),phi_right(3)
     ! -----------------------------------------------------------------
     real(kind=rk)            :: width, max_R, x, x0_inicond, b, radius
     real(kind=rk)            :: rho_left,u_left,p_left, rho_right, u_right, p_right
     integer(kind=ik)         :: ix

     ! density at the left and right of the shock
     rho_left = phi_left(1)
     rho_right= phi_right(1)
     ! velocity at hte left and right of the shock
     u_left = phi_left(2)
     u_right= phi_right(2)
     ! pressure at hte left and right of the shock
     p_left = phi_left(3)
     p_right= phi_right(3)

       ! check for usefull inital values
       if ( rho_left<0 .or. rho_right<0 .or.&
            p_left  <0 .or. p_right  <0) then
         call abort(3572,"ERROR initial values are insufficient")
       end if
       ! following values are imposed and smoothed with tangens:
       ! ------------------------------------------
       !   rhoL    | rhoR                  | rhoL
       !   uL      | uR                    | uL
       !   pL      | pR                    | pL
       ! 0-----------------------------------------Lx
       !           x0_inicond             x0_inicond+width
       width       = params_ns%domain_size(1)*(1-params_ns%inicond_width-0.1)
       x0_inicond  = params_ns%inicond_width*params_ns%domain_size(1)
       max_R       = width*0.5_rk
       do ix=g+1, Bs+g
          x = dble(ix-(g+1)) * dx(1) + x0(1)
          call continue_periodic(x,params_ns%domain_size(1))
          ! left region
          radius = abs(x-x0_inicond-width*0.5_rk)
          b      = 0.5_rk*(1-tanh((radius-(max_R-10*dx(1)))*2*PI/(10*dx(1)) ))
          phi(ix, rhoF) = dsqrt(rho_left)-b*(dsqrt(rho_left)-dsqrt(rho_right))
          phi(ix, UxF)  = phi(ix,rhoF)*(u_left-b*(u_left-u_right))
          phi(ix, UyF)  = 0.0_rk
          phi(ix, pF)   = p_left-b*(p_left - p_right)
       end do



  end subroutine set_inicond_moving_shock
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




end module module_shock_tube
