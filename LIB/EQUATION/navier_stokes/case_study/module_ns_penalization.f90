
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
  PUBLIC :: init_mask,add_constraints,get_mask
  !**********************************************************************************************

!  real(kind=rk),    allocatable,     save        :: mask(:,:,:)
  character(len=80),                 save        :: mask_geometry!273.15_rk
  logical           ,                save        :: smooth_mask
  real(kind=rk)                    , save        :: C_eta_inv,C_sp_inv
  real(kind=rk),                     save        :: domain_size(3)
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



  type :: type_funnel_plate
    real(kind=rk) :: x0(1:2)
    real(kind=rk) :: width 
    real(kind=rk) :: r_in
    real(kind=rk) :: r_out   
  end type type_funnel_plate 
 

  type :: type_funnel
      real(kind=rk)       ::outer_diameter         ! outer diameter
      real(kind=rk)       ::max_inner_diameter     ! maximal inner diameter
      real(kind=rk)       ::min_inner_diameter     ! minimal inner diameter
      integer(kind=ik)    ::nr_plates              ! Number of plates
      real(kind=rk)       ::plates_distance        ! distance between origin of plates
      real(kind=rk)       ::plates_thickness       ! 
      real(kind=rk)       ::temperatur             ! temperatur of plates
      

      real(kind=rk)       ::length                 ! total length of funnel
      real(kind=rk)       ::slope                  ! slope of funnel
      real(kind=rk)       ::offset(2)              ! offset of funnel in x and y
      
      ! parameters of flow inlet outlet
      real(kind=rk)       ::pump_diameter          
      real(kind=rk)       ::pump_x_center
      real(kind=rk)       ::jet_radius             ! slope of funnel
      real(kind=rk)       ::wall_thickness       ! 
    
      real(kind=rk)       ::inlet_velocity(2)       ! 
      real(kind=rk)       ::inlet_density       ! 
      real(kind=rk)       ::inlet_pressure       ! 
      real(kind=rk)       ::outlet_pressure       ! 
      real(kind=rk)       ::pump_speed       ! 
      
      type(type_funnel_plate), allocatable:: plate(:)
  end type type_funnel


  !------------------------------------------------
  type(type_funnel)   , save :: funnel
  type(type_cylinder) , save :: cyl
  !------------------------------------------------

contains

include "funnel.f90"
include "vortex_street.f90"



!> \brief reads parameters for mask function from file
subroutine init_mask( filename )

  character(len=*), intent(in) :: filename

    ! inifile structure
    type(inifile) :: FILE
    call read_ini_file_mpi(FILE, filename, .true.)
    
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', smooth_mask, .true.)
    call read_param_mpi(FILE, 'VPM', 'geometry', mask_geometry, "cylinder")
    call read_param_mpi(FILE, 'VPM', 'C_eta', C_eta_inv, 0.01_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Lx', domain_size(1), 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', domain_size(2), 1.0_rk )
    call read_param_mpi(FILE, 'Physics', 'C_sp', C_sp_inv, 0.01_rk )
        ! read adiabatic coefficient
    call read_param_mpi(FILE, 'Navier_Stokes', 'gamma_', gamma_, 0.0_rk )
    ! read specific gas constant
    call read_param_mpi(FILE, 'Navier_Stokes', 'Rs', Rs, 0.0_rk )

    C_eta_inv=1.0_rk/C_eta_inv
    C_sp_inv =1.0_rk/C_sp_inv

    select case(mask_geometry)
    
    case ('vortex_street','cylinder')
      call init_vortex_street(FILE)
    case ('funnel') 
      call init_funnel(FILE)
    case default
      call abort(8546501,"[module_mask.f90] ERROR: geometry for VPM is unknown"//mask_geometry)
    
    end select
  
end subroutine init_mask


 
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
    case('vortex_street','cylinder')
      call draw_cylinder(mask, x0, dx, Bs, g )
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
    case('vortex_street','cylinder')
      call add_cylinder(penalization, x0, dx, Bs, g, phi )
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
    
    !rhs(:,:,4)=rhs(:,:,4) - (gamma_-1)       *penalization( :, :,4)

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


function jet_stream(radius,max_R,smooth_width)

  real(kind=rk), intent(in)      :: radius, smooth_width,max_R
  real(kind=rk)                  :: jet_stream 
    
    jet_stream=0.5_rk*(1-tanh((radius-(max_R-0.5_rk*smooth_width))*2*PI/smooth_width ))
end function jet_stream



end module module_ns_penalization
