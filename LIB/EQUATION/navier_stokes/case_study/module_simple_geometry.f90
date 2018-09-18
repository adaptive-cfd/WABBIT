
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
  use module_helpers


  implicit none

  !**********************************************************************************************
  ! only this functions are visible outside this module
  PUBLIC :: read_params_geometry,add_geometry2D,draw_geometry
  !**********************************************************************************************
  ! make everything private if not explicitly marked public
  PRIVATE
  !**********************************************************************************************
  character(len=80),save  :: MASK_GEOMETRY
  logical,save            :: FREE_OUTLET_WALL=.false.
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


!=========================================================================================
! INITIALIZATIONs
!=========================================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \brief reads parameters for mask function from file
subroutine read_params_geometry( params,FILE )
    implicit none
    !--------------------------------------------
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params
    !---------------------------------------------
    ! available initial conditions of this case
    character(len=*), parameter, dimension(5)  :: INI_LIST  = (/Character(len=14) ::  'mask', 'moving-shock', &
    'pressure-blob', 'zeros', 'simple-shock' /)

    ! add sponge parameters for the in and outflow of the domain
    call init_simple_sponge(params,FILE)

    if (params%mpirank==0) then
     write(*,*)
     write(*,*)
     write(*,*) "PARAMS: simple geometry!"
     write(*,'(" ------------------------")')
   endif

   ! geometry of the solid obstacle
    call read_param_mpi(FILE, 'simple_geometry', 'geometry', mask_geometry, '--' )
    ! free outlet sponge
    call read_param_mpi(FILE, 'simple_geometry', 'free_outlet_wall', FREE_OUTLET_WALL, .false. )


    select case(mask_geometry)
    case ('triangle','rhombus')
      call init_triangle(params,FILE)
    case ('vortex_street','cylinder')
      call init_vortex_street(FILE)
    case default
      call abort(8546501," ERROR: geometry for VPM is unknown"//mask_geometry)
    end select

    ! check if the initial condition is in the list of availabel initial conditions for this case!
    if(  list_contains_name(INI_LIST,params%inicond)==0 ) then
      call abort(111109,"Error inicond ["//trim(params%inicond)//"] not available for this case!")
    endif
end subroutine read_params_geometry
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \brief reads parameters for mask function from file
subroutine init_vortex_street(FILE)

  implicit none
 ! character(len=*), intent(in) :: filename
  type(inifile) , intent(inout) :: FILE
  character(len=*),parameter         :: section='simple_geometry'

  call read_param_mpi(FILE, section, 'x_cntr', cyl%x_cntr,(/ 0.25_rk , 0.5_rk/) )
  call read_param_mpi(FILE, section, 'length', cyl%radius,0.05_rk*min(params_ns%domain_size(1),params_ns%domain_size(2)) )
  cyl%radius=cyl%radius*0.5_rk
  cyl%x_cntr(1)=cyl%x_cntr(1)*params_ns%domain_size(1)
  cyl%x_cntr(2)=cyl%x_cntr(2)*params_ns%domain_size(2)
end subroutine init_vortex_street
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \brief reads parameters for mask function from file
subroutine init_triangle(params,FILE)

 implicit none
 !> params structure of navier stokes
  type(type_params_ns),intent(inout)  :: params
! character(len=*), intent(in) :: filename
 type(inifile) , intent(inout) :: FILE
 character(len=*),parameter         :: section='simple_geometry'

 ! read in the angle between the symmetry axis and the triangle side
 call read_param_mpi(FILE, section, 'angle', triangle%angle, 30.0_rk )
 call read_param_mpi(FILE, section, 'x_cntr', triangle%x_cntr, (/0.2_rk, 0.5_rk, 0.0_rk/) )
 call read_param_mpi(FILE, section, 'length', triangle%length, 0.1*params%domain_size(1) )

 ! convert degrees to radians
 if ( triangle%angle>0.0_rk .and. triangle%angle<90.0_rk ) then
   triangle%angle=triangle%angle*PI/180.0_rk
 else
   call abort(45756,"somebody has to go back to preeshool! 0< angle <90")
 end if

 !convert from height to length
 triangle%length=triangle%length/(2*tan(triangle%angle))

 ! a rhombus is only a double triangle
 if ( mask_geometry=="rhombus" ) then
   triangle%rhombus    =.true.
 else
   triangle%rhombus    =.false.
 end if

 if (triangle%x_cntr(1)>1.0_rk .or. triangle%x_cntr(1)<0.0_rk ) then
   triangle%x_cntr(1) =0.2_rk
 end if

 if (triangle%x_cntr(2)>1.0_rk .or. triangle%x_cntr(2)<0.0_rk ) then
     triangle%x_cntr(2) =0.5_rk
 end if
 triangle%x_cntr(1)=triangle%x_cntr(1)*params%domain_size(1)
 triangle%x_cntr(2)=triangle%x_cntr(2)*params%domain_size(2)
 triangle%x_cntr(3)=triangle%x_cntr(3)*params%domain_size(3)

end subroutine init_triangle
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!=========================================================================================



!==========================================================================
!> This function adds a penalization term
!> to navier stokes equations
subroutine add_geometry2D(penalization, x0, dx, Bs, g, phi)
      implicit none
      !-------------------------------------------------------
      !> grid parameter
      integer(kind=ik), intent(in)    :: g, Bs
      !> rhs
      real(kind=rk), intent(inout)    :: penalization(:,:,:)
      !> state variables
      real(kind=rk), intent(in)       :: phi(:,:,:)
      !> spacing and origin of block
      real(kind=rk), intent(in)       :: x0(3), dx(3)
      !--------------------------------------------------------
      real(kind=rk),allocatable ,save, dimension(:,:)  :: mask,p,rho,u,v
      real(kind=rk),save  :: T0,rho0,u0,p0,Rs
      logical, save       :: tmp_fields_allocated=.false.


      if (size(penalization,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

      if (.not. tmp_fields_allocated) then
        allocate(mask(Bs+2*g,Bs+2*g))
        allocate(p(Bs+2*g,Bs+2*g))
        allocate(rho(Bs+2*g,Bs+2*g))
        allocate(u(Bs+2*g,Bs+2*g))
        allocate(v(Bs+2*g,Bs+2*g))
        tmp_fields_allocated=.true.
      endif

      ! reset mask array
      mask    = 0.0_rk

      T0      = params_ns%initial_temp
      rho0    = params_ns%initial_density
      p0      =  params_ns%initial_pressure
      u0      =  params_ns%initial_velocity(1)
      Rs      =  params_ns%Rs

      rho         = phi(:,:,rhoF)**2
      u           = phi(:,:,UxF)/phi(:,:,rhoF)
      v           = phi(:,:,UyF)/phi(:,:,rhoF)
      p           = phi(:,:,pF)

      ! geometry
      !----------
      select case( MASK_GEOMETRY )
      case('vortex_street','cylinder')
        call draw_cylinder(mask, x0, dx, Bs, g )
      case('rhombus','triangle')
        call draw_triangle(mask, x0, dx, Bs, g )
      case default
        call abort(120401,"ERROR: geometry for VPM is unknown"//mask_geometry)
      end select

      penalization (:,:,rhoF) = 0.0_rk
      ! x-velocity
      penalization(:,:,UxF)=C_eta_inv*mask * ( rho*u )
      ! y-velocity
      penalization(:,:,UyF)=C_eta_inv*mask * ( rho*v )
      ! preasure
      penalization(:,:,pF)=C_eta_inv*mask *( p- rho*Rs*T0 )

      ! sponge
      !--------
      call draw_simple_sponge(mask, x0, dx, Bs, g)

      penalization(:,:,rhoF)= penalization(:,:,rhoF) + C_sp_inv*mask * ( rho - rho0 )
      ! x-velocity
      penalization(:,:,UxF)= penalization(:,:,UxF) + C_sp_inv*mask * ( rho*u - rho0*u0  )
      ! y-velocity
      penalization(:,:,UyF)= penalization(:,:,UyF) + C_sp_inv*mask * ( rho*v )
      ! preasure
      penalization(:,:,pF)= penalization(:,:,pF) + C_sp_inv*mask *( p - p0 )

      ! free outelt wall
      ! ----------------
      ! add a free outlet wall as a sponge in north and south if necessary
      if ( FREE_OUTLET_WALL ) then
        call draw_free_outlet_wall(mask, x0, dx, Bs, g )

        penalization(:,:,1)= penalization(:,:,rhoF) + C_sp_inv*mask * ( rho - rho0 )
        ! x-velocity
        !penalization(:,:,UxF)= penalization(:,:,UxF) + C_sp_inv*mask * ( rho*u - rho0*u0  )
        ! y-velocity
        penalization(:,:,UyF)= penalization(:,:,UyF) + C_sp_inv*mask * ( rho*v )
        ! preasure
        penalization(:,:,pF)= penalization(:,:,pF) + C_sp_inv*mask *( p - p0 )
      endif


end subroutine add_geometry2D
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





!==========================================================================
subroutine draw_cylinder(mask, x0, dx, Bs, g )


    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, y, r, h
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2) - cyl%x_cntr(2)
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1) - cyl%x_cntr(1)
           ! distance from center of cylinder
           r = dsqrt(x*x + y*y)
           mask(ix,iy) = smoothstep( r - cyl%radius, h)
       end do
    end do

end subroutine draw_cylinder
!==========================================================================




!==========================================================================
subroutine draw_triangle(mask, x0, dx, Bs, g )

   implicit none
   ! grid
   integer(kind=ik), intent(in)                              :: Bs, g
   !> mask term for every grid point of this block
   real(kind=rk), dimension(:,:), intent(out)              :: mask
   !> spacing and origin of block
   real(kind=rk), dimension(2), intent(in)                   :: x0, dx

   ! auxiliary variables
   real(kind=rk)                                             :: x, y,h,length,tan_theta,height
   ! loop variables
   integer(kind=ik)                                          :: ix, iy
   !
   logical                                                   :: triangle_is_rhombus

!---------------------------------------------------------------------------------------------
! variables initialization
   if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

   ! reset mask array
   mask  = 0.0_rk

   triangle_is_rhombus =triangle%rhombus
   length              =triangle%length
   tan_theta           =tan(triangle%angle)

!---------------------------------------------------------------------------------------------
! main body


   ! parameter for smoothing function (width)
   h = 1.5_rk*max(dx(1), dx(2))

   do ix=1, Bs+2*g
       x = dble(ix-(g+1)) * dx(1) + x0(1)
       x = x-triangle%x_cntr(1)
       !YES a rhombus is made from two isosceles triangles
       if ( triangle_is_rhombus .and. x>length ) then
         ! reflect triangle at the line at x=length parallel to y axis: <|>
         x = 2 * length - x
       end if
       height = tan_theta*x
       do iy=1, Bs+2*g
          y = dble(iy-(g+1)) * dx(2) + x0(2)
          y = abs(y-triangle%x_cntr(2))

          mask(ix,iy) = soft_bump(x,0.0_rk,length+h,h)*smoothstep(y-height,h)

      end do
   end do
end subroutine draw_triangle
!==========================================================================





end module module_simple_geometry
