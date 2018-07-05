
!-----------------------------------------------------------------
!> \file
!> \brief
!! Module for volume penaliza \f$\chi(x,t)\f$
!> \details
!> This module implements the mask function on each Block
!!          \f[
!!                           \chi(x,t)\quad \forall\;  x\in\mathcal{–}^l
!!            \f]
!!                        for Volume penalization
!> To increase peformance there are certain tasks:
!> \todo
!>       + include update flag: block_updated
!>                - if true the mask needs to be computed again
!>                - if false the grid has not changed on the proc rank and the mask function stays
!!                  the same
!!       + maybe it is better to have a global mask on the finest grid level and coarsen it for
!!         lower levels on the specific Blocks
!!
!> \version 23.2.2018
!> \author P.Krah
!-----------------------------------------------------------------

module module_ns_penalization

  !---------------------------------------------------------------------------------------------
  ! modules
  use module_navier_stokes_params
  use module_precision
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
  PUBLIC :: shockVals,init_penalization,add_constraints,get_mask,mean_quantity,integrate_over_pump_area
  !**********************************************************************************************

!  real(kind=rk),    allocatable,     save        :: mask(:,:,:)
  character(len=80),                 save        :: mask_geometry!273.15_rk
  logical           ,                save        :: smooth_mask, use_sponge
  real(kind=rk)                    , save        :: C_eta_inv,C_sp_inv,L_sponge
  real(kind=rk),                     save        :: domain_size(3)
  ! radius of domain (Ly/2)
  real(kind=rk),                     save        :: R_domain
  real(kind=rk),                     save        :: Rs,gamma_
  !------------------------------------------------
  !> \file
  !> \details
  !> Available mask geometrys
  !! ------------------------
  !!
  !!    | geometry    | geometrical properties                             |
  !!    |-------------|----------------------------------------------------|
  !!    | \b cylinder |   \c x_cntr(1:2),\c radius                              |
  !!    | \b funnel   |   \c offset(1:2),\c Nr_plates, \c dx_plates,\c diameter_slope  \c dout, \c dmin \c dmax,\c Nr_plates, \c dx_plates  |
  !!    | \b fplate   |   \c x_0(1:2),\c R_in, \c R_out , \c width         |
  !!     ------------------------------------------
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



  type :: type_funnel_plate
    real(kind=rk) :: x0(1:2)
    real(kind=rk) :: width
    real(kind=rk) :: r_in
    real(kind=rk) :: r_out
  end type type_funnel_plate


  type :: type_funnel
      real(kind=rk)       ::outer_diameter         ! outer diameter
      real(kind=rk)       ::max_inner_diameter     ! maximal inner diameter
      real(kind=rk)       ::min_inner_diameter    =-1.0_rk ! minimal inner diameter
      integer(kind=ik)    ::nr_plates             =0.0_ik ! Number of plates
      real(kind=rk)       ::plates_distance       =-1.0_rk ! distance between origin of plates
      real(kind=rk)       ::plates_thickness      =-1.0_rk !
      real(kind=rk)       ::first_plate_thickness =-1.0_rk
      real(kind=rk)       ::temperatur            =-1.0_rk ! temperatur of plates

      real(kind=rk)       ::length                =-1.0_rk ! total length of funnel
      real(kind=rk)       ::slope                 =-1.0_rk ! slope of funnel
      real(kind=rk)       ::offset(2)             =-1.0_rk ! offset of funnel in x and y

      ! parameters of flow inlet outlet
      real(kind=rk)       ::pump_diameter  =-1.0_rk
      real(kind=rk)       ::pump_x_center  =-1.0_rk
      real(kind=rk)       ::jet_radius     =-1.0_rk        ! cappilary inner Radius
      real(kind=rk)       ::r_out_cappilary=-1.0_rk         ! cappilary outer Radus
      real(kind=rk)       ::wall_thickness =-1.0_rk           !

      real(kind=rk)       ::inlet_velocity(2)       !
      real(kind=rk)       ::inlet_density       !
      real(kind=rk)       ::inlet_pressure       !
      real(kind=rk)       ::outlet_pressure       !
      real(kind=rk)       ::pump_speed       !
      real(kind=rk)       ::pump_density      !
      real(kind=rk)       ::pump_pressure     !


      type(type_funnel_plate), allocatable:: plate(:)
  end type type_funnel


  !------------------------------------------------
  type(type_funnel)   , save :: funnel
  type(type_cylinder) , save :: cyl
  type(type_triangle) , save :: triangle
  type(type_shock_params) , save :: shock_params
  !------------------------------------------------

contains

include "funnel.f90"
include "vortex_street.f90"
include "sod_shock_tube.f90"
include "simple_shock.f90"
include "triangle.f90"



!> \brief reads parameters for mask function from file
subroutine init_penalization( params,FILE )
    use module_navier_stokes_params
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params

     if (params%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: penalization and geometries!"
      write(*,'(" -----------------------------------")')
    endif
    ! =============================================================================
    ! parameters needed for ns_physics module
    ! -----------------------------------------------------------------------------
    call read_param_mpi(FILE, 'VPM', 'penalization', params%penalization, .false.)


    if (params%penalization) then
      call read_param_mpi(FILE, 'VPM', 'geometry', params%geometry, "cylinder")
      call read_param_mpi(FILE, 'Sponge', 'C_sponge',  params%C_sp, 0.01_rk )
      call read_param_mpi(FILE, 'VPM', 'C_eta', params%C_eta, 0.01_rk )
    else
      return
    endif

    ! =============================================================================
    ! parameters needed for penalization only
    ! -----------------------------------------------------------------------------
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', smooth_mask, .true.)
    call read_param_mpi(FILE, 'VPM', 'geometry', mask_geometry, "cylinder")
    call read_param_mpi(FILE, 'DomainSize', 'Lx', domain_size(1), 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', domain_size(2), 1.0_rk )
    call read_param_mpi(FILE,'Sponge', 'C_sponge', C_sp_inv, 0.01_rk )
        ! read adiabatic coefficient
    call read_param_mpi(FILE, 'Navier_Stokes', 'gamma_', gamma_, 0.0_rk )
    ! read specific gas constant
    call read_param_mpi(FILE, 'Navier_Stokes', 'Rs', Rs, 0.0_rk )

    C_eta_inv=1.0_rk/params%C_eta
    C_sp_inv =1.0_rk/C_sp_inv
    R_domain =domain_size(2)*0.5_rk

    select case(mask_geometry)
    case ('triangle','rhombus')
      call init_simple_sponge(params,FILE)
      call init_triangle(params,FILE)
    case ('moving-shock')
      call init_simple_shock(params,FILE)
    case('sod_shock_tube')
      ! nothing to do
    case ('vortex_street','cylinder')
      call init_simple_sponge(params,FILE)
      call init_vortex_street(FILE)
    case ('funnel')
      call init_funnel(FILE)
    case default
      call abort(8546501,"[module_ns_penalization.f90] ERROR: geometry for VPM is unknown"//mask_geometry)

    end select

end subroutine init_penalization





!> \brief reads parameters for sponge initializing a sipmle sponge form file
!> \details The simple sponge starts at x=0 and ends at x=L_sponge
subroutine init_simple_sponge( params,FILE )
    use module_navier_stokes_params
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params

     if (params%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: init simple sponge!"
      write(*,'(" -----------------------------------")')
    endif

    call read_param_mpi(FILE, 'Sponge', 'use_sponge', use_sponge, .false. )

    if (use_sponge) then
      call read_param_mpi(FILE, 'Sponge', 'L_sponge', L_sponge, 0.1_rk )
      if ( L_sponge<1e-10 ) then
        L_sponge=0.1*params%Lx
      end if
      call read_param_mpi(FILE, 'Sponge', 'C_sponge', C_sp_inv, 1.0e-2_rk )
      C_sp_inv=1.0_rk/C_sp_inv
    else
      return
    endif

end subroutine init_simple_sponge






















subroutine get_mask(mask, x0, dx, Bs, g )


    ! grid
    integer(kind=ik), intent(in) :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(inout) :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in) :: x0, dx

    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

!---------------------------------------------------------------------------------------------
! variables initialization
    select case(mask_geometry)
    case('moving-shock','sod_shock_tube')
     call draw_sod_shock_tube(mask, x0, dx, Bs, g,'boundary')
    case('vortex_street','cylinder')
      call draw_cylinder(mask, x0, dx, Bs, g )
    case('rhombus','triangle')
      call draw_triangle(mask,x0,dx,Bs,g)
    case('funnel')
      call draw_funnel(mask, x0, dx, Bs, g)
    case default
      call abort(120601,"ERROR: geometry for VPM is unknown"//mask_geometry)
    end select

end subroutine get_mask








!==========================================================================
!> \brief This function computes a penalization term for different geometries
!! in Navier Stokes physics
subroutine add_constraints(rhs ,Bs , g, x0, dx, phi)

    !-------------------------------------------------------
    !> grid parameter
    integer(kind=ik), intent(in)    :: g, Bs
    !> rhs
    real(kind=rk), intent(inout)    :: rhs(Bs+2*g, Bs+2*g,4)
    !> state variables
    real(kind=rk), intent(in)       :: phi(Bs+2*g, Bs+2*g,4)
    !> spacing and origin of block
    real(kind=rk), intent(in)       :: x0(2), dx(2)
    !--------------------------------------------------------
     real(kind=rk)                  :: penalization(Bs+2*g, Bs+2*g,4)

    ! 1. compute volume penalization term for the different case studies

    select case(mask_geometry)
    case('sod_shock_tube')
      call add_sod_shock_tube(penalization, x0, dx, Bs, g, phi )
    case('moving-shock')
      call add_simple_shock(penalization,x0,dx,Bs,g,phi)
    case('vortex_street','cylinder')
      call add_cylinder(penalization, x0, dx, Bs, g, phi )
    case('rhombus','triangle')
      call add_triangle(penalization,x0,dx,Bs,g,phi)
    case('funnel')
      call add_funnel(  penalization, x0, dx, Bs, g, phi )
    case default
      call abort(120401,"ERROR: geometry for VPM is unknown"//mask_geometry)
    end select



    ! 2. add penalty to the right hand side

    ! sqrt(rho) component (density)
    rhs(:,:,1)=rhs(:,:,1) - 0.5_rk/phi(:,:,1)*penalization( :,:,1)
    ! sqrt(rho)u component (momentum)
    rhs(:,:,2)=rhs(:,:,2) - 1.0_rk/phi(:,:,1)*penalization( :,:,2)
    ! sqrt(rho)v component (momentum)
    rhs(:,:,3)=rhs(:,:,3) - 1.0_rk/phi(:,:,1)*penalization( :, :,3)
    ! p component (preasure/energy)
    rhs(:,:,4)=rhs(:,:,4) -                   penalization( :, :,4)

end subroutine add_constraints
!==========================================================================





!==========================================================================
  !> \brief This subroutine returns the value f of a smooth step function \n
  !> The sharp step function would be 1 if delta<=0 and 0 if delta>0 \n
  !> h is the semi-size of the smoothing area, so \n
  !> f is 1 if delta<=0-h \n
  !> f is 0 if delta>0+h \n
  !> f is variable (smooth) in between
  !> \details
  !> \image html maskfunction.bmp "plot of chi(delta)"
  !> \image latex maskfunction.eps "plot of chi(delta)"
    function smoothstep(delta,h)

      implicit none
      real(kind=rk), intent(in)  :: delta,h

      real(kind=rk)              :: smoothstep,f
      !-------------------------------------------------
      ! cos shaped smoothing (compact in phys.space)
      !-------------------------------------------------
      if (delta<=-h) then
        f = 1.0_rk
      elseif ( -h<delta .and. delta<+h  ) then
        f = 0.5_rk * (1.0_rk + dcos((delta+h) * pi / (2.0_rk*h)) )
      else
        f = 0.0_rk
      endif

      smoothstep=f
    end function smoothstep
!==========================================================================





function soft_bump(x,x0,width,h)

  real(kind=rk), intent(in)      :: x, x0, h, width
  real(kind=rk)                  :: soft_bump,d

    d=x-x0

    if (d>=h .and. d<=width-h) then
        soft_bump = 1.0
    elseif ((-h<d) .and.(d<+h)) then
        soft_bump = 0.5 * (1.0 - dcos((d+h) * pi / (2.0*h)) )
    elseif (d>width-h .and. d<width+h) then
        soft_bump = 0.5 * (1.0 - dcos((d-h-width) * pi / (-2.0*(h))) )
    else
        soft_bump = 0.0
    endif

end function soft_bump

function soft_bump2(x,x0,width,h)

  real(kind=rk), intent(in)      :: x, x0, h, width
  real(kind=rk)                  :: soft_bump2,max_R,smooth_width,radius

  soft_bump2=soft_bump(x,x0+h,width-2.0_rk*h,h)

end function soft_bump2

function jet_stream(radius,max_R,smooth_width)

  real(kind=rk), intent(in)      :: radius, smooth_width,max_R
  real(kind=rk)                  :: jet_stream

    jet_stream=0.5_rk*(1-tanh((radius-(max_R-0.5_rk*smooth_width))*2*PI/smooth_width ))
end function jet_stream


!> \brief computes a smooth transition between val_left and val_right
function transition(x,x0,trans_width,val_left,val_rigth)
  !> position x
  real(kind=rk), intent(in)     :: x
  !> beginning point of transition
  real(kind=rk), intent(in)     :: x0
  !> length of transition region
  real(kind=rk), intent(in)     :: trans_width
  !> values at the left and right of the transition region
  real(kind=rk), intent(in)     ::val_left,val_rigth

  !--------------------------------------------
  real(kind=rk)                  :: s,units,transition


    ! convert transition range from 2pi to trans_width
    units  = trans_width/(2*PI)
    ! transition function 0<s<1
    s         = 0.5_rk+0.5_rk * tanh((x-x0-0.5_rk*trans_width)/units)

    transition= s * val_rigth + (1-s) * val_left
end function transition





!> \brief creates the mask of a simple sponge
subroutine simple_sponge(sponge, x0, dx, Bs, g)

    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> sponge term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(out)     :: sponge
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset sponge array
    sponge = 0.0_rk
!---------------------------------------------------------------------------------------------
! main body
    if ( use_sponge .eqv. .false. ) then
      return
    end if


        do ix=1, Bs+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            do iy=1, Bs+2*g

            sponge(ix,iy)=soft_bump2(x,(domain_size(1)-L_sponge),L_sponge,1.5*min(dx(1),dx(2)))


    !        if ((domain_size(1)-x) <= L_sponge) then
    !            sponge(ix,iy) = (x-(domain_size(1)-L_sponge))**2
    !        elseif (x <= L_sponge) then
    !            sponge(ix,iy) = (x-L_sponge)**2
    !        else
    !            sponge(ix,iy) = 0.0_rk
    !        end if
       end do
     end do

end subroutine simple_sponge
















end module module_ns_penalization
