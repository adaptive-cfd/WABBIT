!----------------------------------------------------------------
!> This module implements different shock tests to validate the
!> euler equations of NS and to test the filters in a multiresolution
!> framework.
!> \details
!> \author PKrah
!> \date 10.08.2018 - creation of the module
!----------------------------------------------------------------

module module_shock

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
    PUBLIC :: read_params_shock_tube,draw_simple_shock, set_inicond_shock_tube, &
              set_shock_in_direction,moving_shockVals, standing_shockVals, &
              shock_tube_penalization2D, shock_tube_penalization3D
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
      character(len=80) :: name                      !< name of the shock (sod_shock_tube,moving-shock,etc.)
      real(kind=rk)     :: rho_left,u_left,p_left    !< left values (behind the shock front)
      real(kind=rk)     :: rho_right,u_right,p_right !< right values of the shock front
      real(kind=rk)     :: machnumber     !< machnumber of a moving shock
      integer(kind=ik)   :: normalvector=-1   !< normalvector of the shock front
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
  !> of the shock to generate a standing shock front
  subroutine standing_shockVals(rho1,u1,p1,rho2,u2,p2,gamma)
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
  end subroutine standing_shockVals
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
        call standing_shockVals( shock%rho_left,shock%u_left,shock%p_left, &
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
      ! read in the normal vector of the shock wavefront
      ! 1: direction e_x
      ! 2: direction e_y
      ! 3: direction e_z only in 3D
      call read_param_mpi(FILE, 'Shock_Tube', 'shock_front_normal', shock%normalvector,1_ik )
      ! for this case we want to have fixed sponge layer size
      params_ns%L_sponge=0.1*params_ns%domain_size(shock%normalvector)
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
       integer(kind=ik), intent(in) :: g
       integer(kind=ik), dimension(3), intent(in) :: Bs !< grid parameter
       real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
       real(kind=rk), intent(inout)  :: u(:,:,:,:)    !< Statevector for t=0
       ! -----------------------------------------------------------------
       real(kind=rk),allocatable:: mask(:,:,:)
       real(kind=rk)            :: x,p_init, rho_init,u_init(3),T_init
       real(kind=rk)            :: u_left(size(u,4)),u_right(size(u,4))
       integer(kind=ik)         :: ix,UshockF


       select case(shock%name)
       case ('moving-shock','standing-shock')
         ! not implemented yet
       case('sod_shock_tube')
         ! Sods test case: shock tube
         ! ---------------------------
         ! Left part x<= Lx/2
         ! rho=1
         ! p  =1
         ! u  =0
         !
         u_left(:)      =0
         u_left(rhoF)   =sqrt(shock%rho_left)
         u_left(pF)     =shock%p_left
         !
         ! Rigth part x> Lx/2
         ! rho=0.125
         ! p  =0.1
         ! u  =0
         u_right(:)=0
         u_right(rhoF)   =sqrt(shock%rho_right)
         u_right(pF)     =shock%p_right
         call set_shock_in_direction(x0, dx, Bs, g, u, u_left, u_right,&
                                      0.5_rk*params_ns%domain_size(1), shock%normalvector)
       case default
         call abort(8546501,"ERROR: Mr. Spock: Insufficient facts always invite danger. "//shock%name)
       end select

    end subroutine set_inicond_shock_tube
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!==========================================================================
!> \brief Compute mask function of sod shock tube
subroutine draw_simple_shock(mask, x0, dx, Bs, g )
    implicit none
    !-----------------------------------------------------------
    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(inout)   :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(3), intent(in)          :: x0, dx
    !-----------------------------------------------------------

    if (size(mask,1) /= Bs(1)+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    if ( params_ns%dim==2 ) then
      call sponge_2D(mask(:,:,1), x0(1:2), dx(1:2), Bs, g, shock%normalvector)
    else
      call wall_3D(mask, x0, dx, Bs, g, shock%normalvector)
    end if

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
subroutine shock_tube_penalization3D(Bs, g, x0, dx, mask, phi_ref)

    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs  !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3)   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:,:) !< reference values of penalized volume
    real(kind=rk), intent(inout)  :: mask(:,:,:,:)    !< mask function
    ! -----------------------------------------------------------------
    ! coordinate systems
    real(kind=rk)       :: X(3),domain_size(3)
    integer(kind=ik)    :: ix, iy, iz, alpha,UshockF,neq
    ! left and right boundary
    real(kind=rk)       :: u_ref,u_R,u_L,rho_R,rho_L,p_L,p_R,width,rho_ref,p_ref

    if (size(mask,3) /= Bs(3)+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! index x_alpha of the coordinate vector which is the normal of the shock front
    alpha=shock%normalvector

    ! draw a shock on the boundaries of the x_alhpa coordinate
    ! all compontents are penalized
    do neq = 1, params_ns%n_eqn
      call sponge_3D(mask(:,:,:,neq), x0, dx, Bs, g, alpha)
      !we allready add the penalization strength to the mask
      mask(:,:,:,neq)=C_sp_inv*mask(:,:,:,neq)
    end do


    if ( alpha==1 ) then
      UshockF=UxF
    elseif (alpha==2) then
      UshockF=UyF
    else
      UshockF=UzF
    end if

    domain_size = params_ns%domain_size
    width    = params_ns%L_sponge* 0.25_rk ! width of the transition area
    !                    shockfront
    ! -----------------------------------------
    ! |sponge  \ rhoL    ||     rhoR  \sponge |
    ! |sponge  \ uL      ||-->  uR    \sponge |
    ! |sponge  \ pL      ||     pR    \sponge |
    ! 0-------------------_--------------------L_alpha
    !          L_alpha/2

    ! left side of the shock
    p_L      = shock%p_left
    rho_L    =  shock%rho_left
    u_L      =  shock%u_left

    ! right side of the shock
    p_R      = shock%p_right
    rho_R    =  shock%rho_right
    u_R      =  shock%u_right

    ! parameter for smoothing function (width)

    do iz = 1, Bs(3)+2*g
        ! z-coordinate
        X(3) = dble(iz-(g+1)) * dx(3) + x0(3)

        do iy = 1, Bs(2)+2*g
            ! y-coordinate
            X(2) = dble(iy-(g+1)) * dx(2) + x0(2)

            do ix = 1, Bs(1)+2*g
                ! x-coordinate
                X(1) = dble(ix-(g+1)) * dx(1) + x0(1)

                !default is 0
                phi_ref(ix,iy,iz,:)=0.0_rk
                ! at this stage we deside which coordinate x_alpha (alpha=1,2,3) the shock front travels along
                if (X(alpha)<domain_size(alpha)*0.5_rk) then
                  phi_ref(ix,iy,iz,rhoF)   = rho_L
                  phi_ref(ix,iy,iz,pF)     = p_L
                  phi_ref(ix,iy,iz,UshockF)= u_L
                else
                  ! Because WABBIT provides only periodic BC, we have to ensure a smooth transition
                  ! from the sponge values on the left side of the domain (rho_L, p_L ,u_L)
                  ! to the values on the right side (rho_R, p_R, u_R),
                  ! For this reason we use the transition function for assigning the values
                  ! inside the right sponge domain.
                  ! The transition is located in the right sponge (x>x_R) and starts
                  ! from xstart to xstart+widht
                  phi_ref(ix,iy,iz,rhoF)   = transition(X(alpha),domain_size(alpha)-width,width,rho_R,rho_L)
                  phi_ref(ix,iy,iz,pF)     = transition(X(alpha),domain_size(alpha)-width,width,p_R  ,p_L    )
                  phi_ref(ix,iy,iz,UshockF)= transition(X(alpha),domain_size(alpha)-width,width,u_R  ,u_L    )
                endif
            end do
        end do
    end do

end subroutine shock_tube_penalization3D



!==========================================================================
!> \brief Adds penalization terms for a simple shock_tube
!> \detail
!>        A mask will be generatet at 0<x<0.1*Lx and 0.9Lx<x<Lx and
!>        the left and right values of the shock tube will be assigned to
!>        the penalized domains.
!>        Be aware that this function uses the C_sp_inv parameter globaly defined
!>        when importing the module_penalization.
subroutine shock_tube_penalization2D(Bs, g, x0, dx, mask, phi_ref)

    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs!< grid parameter
    real(kind=rk), intent(in)     :: x0(2), dx(2)   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:) !< reference values of penalized volume
    real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
    ! -----------------------------------------------------------------
    ! auxiliary variables
    real(kind=rk)       :: X(2),domain_size(2)
    integer(kind=ik)    :: ix, iy, alpha, UshockF,neq
    ! left and right boundary
    real(kind=rk)       :: u_R,u_L,rho_R,rho_L,p_L,p_R,width

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs(1)+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! index x_alpha of the coordinate vector which is the normal of the shock front
    alpha=shock%normalvector

    ! all compontents are penalized
    do neq = 1, params_ns%n_eqn
      call sponge_2D(mask(:,:,neq), x0, dx, Bs, g, alpha)
      !we allready add the penalization strength to the mask
      mask(:,:,neq)=C_sp_inv*mask(:,:,neq)
    end do

    if ( alpha==1 ) then
      UshockF=UxF
    else
      UshockF=UyF
    end if

    domain_size = params_ns%domain_size(1:2)
    width       = params_ns%L_sponge *0.5_rk
    !                    shockfront
    ! -----------------------------------------
    ! |sponge  \ rhoL    ||     rhoR  \sponge |
    ! |sponge  \ uL      ||-->  uR    \sponge |
    ! |sponge  \ pL      ||     pR    \sponge |
    ! 0-------------------_--------------------L_alpha
    !          L_alpha/2

    ! left side of the shock
    p_L      = shock%p_left
    rho_L    =  shock%rho_left
    u_L      =  shock%u_left

    ! right side of the shock
    p_R      = shock%p_right
    rho_R    =  shock%rho_right
    u_R      =  shock%u_right

    ! parameter for smoothing function (width)

        do iy = 1, Bs(2)+2*g
            ! y-coordinate
            X(2) = dble(iy-(g+1)) * dx(2) + x0(2)

            do ix = 1, Bs(1)+2*g
                ! x-coordinate
                X(1) = dble(ix-(g+1)) * dx(1) + x0(1)

                !default is 0
                phi_ref(ix,iy,:)=0.0_rk
                ! at this stage we deside which coordinate x_alpha (alpha=1,2,3) the shock front travels along
                if (X(alpha)<domain_size(alpha)*0.5_rk) then
                  phi_ref(ix,iy,rhoF)   = rho_L
                  phi_ref(ix,iy,pF)     = p_L
                  phi_ref(ix,iy,UshockF)= u_L
                else
                  ! Because WABBIT provides only periodic BC, we have to ensure a smooth transition
                  ! from the sponge values on the left side of the domain (rho_L, p_L ,u_L)
                  ! to the values on the right side (rho_R, p_R, u_R),
                  ! For this reason we use the transition function for assigning the values
                  ! inside the right sponge domain.
                  ! The transition is located in the right sponge (x>x_R) and starts
                  ! from x0=0.925_rk*domain_size(1) to x0+widht=(0.925_rk+0.05)*domain_size(1)
                  phi_ref(ix,iy,rhoF)   = transition(X(alpha),domain_size(alpha)-width,width,rho_R,rho_L)
                  phi_ref(ix,iy,pF)     = transition(X(alpha),domain_size(alpha)-width,width,p_R  ,p_L    )
                  phi_ref(ix,iy,UshockF)= transition(X(alpha),domain_size(alpha)-width,width,u_R  ,u_L    )
                endif
            end do
        end do

end subroutine shock_tube_penalization2D

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> \brief
  !>    + Sets the initial condition of a shock with shock values phi_left and phi_right
  !>      along the axis alpha,
  !>    + example: set_shock_in_direction(x0, dx, Bs, g, phi, phi_left, phi_right, 0.4*domain_size(3), 3)
  !>               for a shock wave along the z axis
  subroutine set_shock_in_direction(x0, dx, Bs, g, phi, phi_left, phi_right, x0_shock, alpha)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs !< grid parameter
    real(kind=rk), intent(in)          :: x0(1:3), dx(1:3) !< block coordinate x0 and spacing dx
    !> direction \f$\vec{e}_\alpha\f$ along the shock front is moving (perpendicular to the initial shock front)
    !>          * alpha=1 : shock is put in x direction
    !>          * alpha=2 : shock is put in y direction
    !>          * alpha=3 : shock is put in z direction
    integer(kind=ik),intent(in)        :: alpha
    ! shockfront
    !       |
    !       |
    !       |----------------> x0(alpha),dx(alpha)
    ! left  |  right
    !       |
    !
    real(kind=rk), intent(in)          :: x0_shock    !< x0_shock is the position of the shock front
    real(kind=rk), intent(inout)       :: phi(:,:,:,:)    !< Statevector in skew-symetric form : phi(ix,iy,dF)
    !> Statevector values for the left and right side of the shock
    !>       + phi_[left/right](rhoF) density at the left or right
    !>       + phi_[left/right](UxF) velocity at the left or right
    !>       + phi_[left/right](pF) pressure at the left or right
    real(kind=rk), intent(in)          :: phi_left(size(phi,2)),phi_right(size(phi,2))
    ! -----------------------------------------------------------------
    integer(kind=ik)    :: ix,iy,iz,nF ! loop variables
    real(kind=rk)       :: X(3),b,h
    ! check if the direction is allowed in the 2D case
    if ( alpha > params_ns%dim ) call abort(2410181,"Oh no, I am sorry to tell you that the direction is not available!")

  h        = 3*dx(alpha) !smoothing width
  ! loop over all fields
  if (params_ns%dim==3) then
    do nF = 1, params_ns%n_eqn
      do iz = 1, Bs(3)+2*g
          ! z-coordinate
          X(3) = dble(iz-(g+1)) * dx(3) + x0(3)

          do iy = 1, Bs(2)+2*g
              ! y-coordinate
              X(2) = dble(iy-(g+1)) * dx(2) + x0(2)

              do ix = 1, Bs(1)+2*g
                  ! x-coordinate
                  X(1) = dble(ix-(g+1)) * dx(1) + x0(1)
                  !default is 0
                  phi(ix,iy,iz,nF)=0.0_rk
                  ! at this stage we deside which coordinate x_alpha (alpha=1,2,3) the shock front travels along
                  b      = smoothstep(x0_shock -X(alpha),h)
                  phi(ix,iy,iz,nF) = phi_left(nF)-b*(phi_left(nF)-phi_right(nF))
              end do !ix
          end do !iy
      end do !iz
    end do! field
  else
    do nF = 1, params_ns%n_eqn
          do iy = 1, Bs(2)+2*g
              ! y-coordinate
              X(2) = dble(iy-(g+1)) * dx(2) + x0(2)

              do ix = 1, Bs(1)+2*g
                  ! x-coordinate
                  X(1) = dble(ix-(g+1)) * dx(1) + x0(1)
                  !default is 0
                  phi(ix,iy,1,nF)=0.0_rk
                  ! at this stage we deside which coordinate x_alpha (alpha=1,2,3) the shock front travels along
                  b      = smoothstep(x0_shock -X(alpha),h)
                  phi(ix,iy,1,nF) = phi_left(nF)-b*(phi_left(nF)-phi_right(nF))
              end do !ix
          end do !iy
    end do! field
  end if


end subroutine set_shock_in_direction
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




end module module_shock
