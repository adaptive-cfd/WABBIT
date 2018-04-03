
!-----------------------------------------------------------------
!> \file
!> \brief
!! Module of mask function \f$\chi(x,t)\f$
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

module module_mask

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
  PUBLIC :: init_mask,get_mask,get_sponge
  !**********************************************************************************************

!  real(kind=rk),    allocatable,     save        :: mask(:,:,:)
  character(len=80),                 save        :: mask_geometry!273.15_rk
  logical           ,                save        :: smooth_mask,penalization
  real(kind=rk),                     save        :: C_eta
  real(kind=rk),                     save        :: domain_size(3)
  
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
      
      real(kind=rk)       ::length                 ! total length of funnel
      real(kind=rk)       ::slope                  ! slope of funnel
      real(kind=rk)       ::offset(2)              ! offset of funnel in x and y
      
      ! parameters of flow inlet outlet
      real(kind=rk)       ::pump_diameter          
      real(kind=rk)       ::pump_x_center
      real(kind=rk)       ::jet_radius             ! slope of funnel
      real(kind=rk)       ::wall_thickness       ! 
      
      type(type_funnel_plate), allocatable:: plate(:)
  end type type_funnel


  !------------------------------------------------
  type(type_funnel) , save  :: funnel
  type(type_cylinder), save :: cyl
  !------------------------------------------------

contains










!> \brief reads parameters for mask function from file
subroutine init_mask(filename)

  character(len=*), intent(in) :: filename


    ! inifile structure
    type(inifile) :: FILE
    call read_ini_file_mpi(FILE, filename, .true.)

    call read_param_mpi(FILE, 'VPM', 'penalization', penalization, .true.)
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', smooth_mask, .true.)
    call read_param_mpi(FILE, 'VPM', 'geometry', mask_geometry, "cylinder")
    call read_param_mpi(FILE, 'DomainSize', 'Lx', domain_size(1), 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', domain_size(2), 1.0_rk )
    select case(mask_geometry)
    
    case ('cylinder')
    
      call read_param_mpi(FILE, 'VPM', 'x_cntr', cyl%x_cntr,(/0.5_rk, 0.5_rk /) )
      call read_param_mpi(FILE, 'VPM', 'radius', cyl%radius, 0.5_rk )
    
    case ('funnel')
    
      call read_param_mpi(FILE, 'funnel', 'outer_diameter'        , funnel%outer_diameter, domain_size(2)/2.0_rk )
      call read_param_mpi(FILE, 'funnel', 'maximal_inner_diameter', funnel%max_inner_diameter, domain_size(2)/3.0_rk )
      call read_param_mpi(FILE, 'funnel', 'minimal_inner_diameter', funnel%min_inner_diameter, domain_size(2)/4.0_rk )
      call read_param_mpi(FILE, 'funnel', 'Number_of_plates'      , funnel%nr_plates, 30 )
      call read_param_mpi(FILE, 'funnel', 'diameter_per_plate'    , funnel%slope, domain_size(2)/3.0_rk)     
      call read_param_mpi(FILE, 'funnel', 'plates_thickness'      , funnel%plates_thickness, domain_size(1)/100.0_rk)     
      call read_param_mpi(FILE, 'funnel', 'pump_diameter'         , funnel%pump_diameter, domain_size(1)/5.0_rk)     
      call read_param_mpi(FILE, 'funnel', 'jet_diameter'          , funnel%jet_radius, domain_size(2)/20.0_rk)     
      call read_param_mpi(FILE, 'funnel', 'pump_x_center'         , funnel%pump_x_center, domain_size(1)*0.5_rk)     
      
      funnel%length               = domain_size(1)*0.85_rk
      funnel%plates_distance      = (funnel%length-funnel%plates_thickness)/(funnel%nr_plates-1)
      funnel%slope                = funnel%slope/(funnel%plates_distance)
      funnel%wall_thickness       = 0.02*domain_size(1)
      funnel%jet_radius           = funnel%jet_radius/2.0_rk
      write(*,*) "s=", funnel%plates_distance
      write(*,*) "ds=", funnel%slope

      call init_plates(funnel)
    
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
    ! usually, the routine should not be called with no penalization, but if it still
    ! happens, do nothing.
    if ( penalization .eqv. .false.) return

    select case(mask_geometry)
    case('cylinder')
      call draw_cylinder( mask, x0, dx, Bs, g )
    case('funnel')
      call draw_funnel( mask, x0, dx, Bs, g )
    case default
      call abort(120001,"ERROR: geometry for VPM is unknown"//mask_geometry)
    end select

end subroutine get_mask

!--------------------------------------------------------------------------!
!               phils sponge stuff
!--------------------------------------------------------------------------!


!==========================================================================
!> \brief This function computes a 2d sponge term
!!
!! \details The sponge term is
!!   \f{eqnarray*}{
!!           s_q(x,y)=\frac{\chi_{\mathrm sp}(x,y)}{C_{\mathrm sp}}(q(x,y)-q_{\mathrm ref})
!!                          \quad \forall x,y \in \Omega_{\mathrm block}
!!     \f}
!! Where we addopt the notation from <a href="https://arxiv.org/abs/1506.06513">Thomas Engels (2015)</a>
!! - the mask function is \f$\chi_{\rm sp}\f$
!! - \f$C_{\rm sp}\f$ is the sponge coefficient (normaly \f$10^{-1}\f$)

subroutine get_sponge(sponge,Bs, g, x0,dx, rho, u , v , p , C_sp)
    !-------------------------------------------------------
    !> grid parameter
    integer(kind=ik), intent(in)    :: g, Bs
    !> spacing and origin of block
    real(kind=rk), intent(in)       :: x0(2), dx(2)
    !> Sponge coefficient
    real(kind=rk), intent(in)       :: C_sp
    !> quantity \f$q\f$ (veolcity \f$u\f$, preasure \f$p\f$,etc.)
    real(kind=rk), intent(in)       :: rho(Bs+2*g, Bs+2*g),u(Bs+2*g, Bs+2*g),v(Bs+2*g, Bs+2*g),p(Bs+2*g, Bs+2*g)
    !> sponge term \f$s(x,y)\f$
    real(kind=rk)                   :: sponge(Bs+2*g, Bs+2*g,4),mask
    !--------------------------------------------------------
    ! loop variables
    integer                         :: i, n,ix,iy
    ! inverse C_sp
    real(kind=rk)                   :: C_sp_inv,x,y,h

    ! outlets and inlets 
    real(kind=rk)                   :: v_pump, rho_capillary,u_capillary,v_capillary,p_capillary,p_2nd_pump_stage

    C_sp_inv=1.0_rk/C_sp
    sponge=0.0_rk

    ! test!
    v_pump          = 50.0_rk
    rho_capillary   = 1.645_rk
    u_capillary     = 50.0_rk
    v_capillary     = 0.0_rk
    p_capillary     = 101330.0_rk
    p_2nd_pump_stage= 50000.0_rk

    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)
            

            !wall in north and south with 2 pumps
            !------------------------------------
            if (abs(x-funnel%pump_x_center)< funnel%pump_diameter/2) then
                   if (y<funnel%wall_thickness+h) then
                    ! wall in south
                    mask=soft_bump(y,0+h,funnel%wall_thickness-h,h)
                    ! velocity in negative direction (v(ix,iy)-(-v_pump))
                    sponge(ix,iy,3)=C_sp_inv*mask*(v(ix,iy)+v_pump)
                   else if (y>domain_size(2)-funnel%wall_thickness-h) then
                    ! wall in north
                    mask=soft_bump(y,domain_size(2)-funnel%wall_thickness,funnel%wall_thickness-h,h)
                    sponge(ix,iy,3)=C_sp_inv*mask*(v(ix,iy)-v_pump) 
                   else

                   endif

            endif         

            ! wall in EAST with capillary
            !----------------------------
            if (abs(y-domain_size(2)/2)< funnel%jet_radius .and. x<funnel%wall_thickness + h) then
                   ! wall in EAST
                  
                  ! density
                  sponge(ix,iy,1)=(rho(ix,iy)-rho_capillary)
                  ! x-velocity
                  sponge(ix,iy,2)=(u(ix,iy)-u_capillary)
                  ! y-velocity
                  sponge(ix,iy,3)=(v(ix,iy)-v_capillary)
                  ! preasure
                 ! sponge(ix,iy,4)=(p(ix,iy)-p_capillary)
                  
                  sponge(ix,iy,:)=C_sp_inv*smoothstep(x-funnel%wall_thickness,h)* sponge(ix,iy,:)
            endif          

            ! wall in West transition to second pumping stage
            !-------------------------------------------------
            if (abs(y-domain_size(2)/2)< funnel%min_inner_diameter/2 .and. x>domain_size(1)-funnel%wall_thickness - h) then
                  ! only pressure is constraint here 
                  sponge(ix,iy,4)=(p(ix,iy)-p_2nd_pump_stage)
                  
                  sponge(ix,iy,:)=C_sp_inv*smoothstep(domain_size(1)-x-funnel%wall_thickness,h)* sponge(ix,iy,:)
                   ! wall in WEST
                   
            endif



       end do
    end do
end subroutine get_sponge
!==========================================================================



!==========================================================================
! !> \brief This function f(x) implements \n
! !> f(x) is 1 if x(1)<=Lsponge \n
! !> f(x) is 0 else  \n
! function inside_sponge(x)
! !> coordinate vector \f$\vec{x}=(x,y,z)\f$ (real 3d or 2d array)
! real(kind=rk), intent(in)       :: x(:)
! !> logical
! logical                         :: inside_sponge
! ! dimension of array x
! integer                         :: dim,i
! ! size of sponge
! real(kind=rk)                   :: length_sponge

! !> \todo read in length_sponge or thing of something intelligent here
!         length_sponge=params_ns%Lx*0.05
!         if (x(1)<=length_sponge .and. x(1)>=0) then
!             inside_sponge=.true.
!         else
!             inside_sponge=.false.
!         endif


! end function inside_sponge
! !==========================================================================



!==========================================================================
!> \brief This function computes a penalization term
!!
!! \details The penalization term is
!!   \f{eqnarray*}{
!!           p_q(x,y)=\frac{\chi(x,y)}{C_{\eta}}(q(x,y)-q_{\mathrm s})
!!                          \quad \forall x,y \in \Omega_{\mathrm block}
!!     \f}
!! Where we addopt the notation from <a href="https://arxiv.org/abs/1506.06513">Thomas Engels (2015)</a>
!! - the mask function is \f$\chi_{\rm \eta}(\vec{x},t)\f$
!! - \f$C_{\rm \eta}\f$ is the sponge coefficient (normaly \f$10^{-1}\f$)

! elemental function penalization( mask, q, qref, C_eta)

!     !-------------------------------------------------------
!     !> grid parameter
!     real(kind=rk), intent(in)       :: C_eta
!     !> reference value of quantity \f$q\f$ (veolcity \f$u\f$, preasure \f$p\f$,etc.)
!     real(kind=rk), intent(in)       :: qref
!     !> quantity \f$q\f$ (veolcity \f$u\f$, preasure \f$p\f$,etc.)
!     real(kind=rk), intent(in)       :: q
!     !> mask \f$p_q(x,y)\f$
!     real(kind=rk), intent(in)       :: mask
!     !--------------------------------------------------------
!     real(kind=rk)                   :: penalization
!     real(kind=rk)                   :: C_eta_inv

!     ! inverse C_eta
!     C_eta_inv=1.0_rk/C_eta

!     penalization=mask*C_eta_inv*(q-qref)
! end function penalization
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

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body


    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2) - cyl%x_cntr(2)
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1) - cyl%x_cntr(1)
           ! distance from center of cylinder
           r = dsqrt(x*x + y*y)
           if (smooth_mask) then
               mask(ix,iy) = smoothstep( r - cyl%radius, h)
           else
               ! if point is inside the cylinder, set mask to 1
               if (r <= cyl%radius) then
                   mask(ix,iy) = 1.0_rk
               else
                   mask(ix,iy) = 0.0_rk
               end if
           end if
       end do
    end do
end subroutine draw_cylinder





    
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










subroutine draw_funnel(mask, x0, dx, Bs, g )


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
    integer(kind=ik)                                          :: ix, iy,n

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body


    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)
           ! distance frommin_inner_diameterter of cylinder
            do n=1,funnel%nr_plates
                if (x>(funnel%plate(n)%x0(1)-3*h) .and. &
                    x<(funnel%plate(n)%x0(1)+funnel%plates_thickness+3*h)) then
                
                    mask(ix,iy)=mask(ix,iy)+draw_plate(x,y,funnel%plate(n),h)
               endif
            enddo

            ! draw the walls arround the funnel
            mask(ix,iy) =mask(ix,iy) + draw_walls(x,y,funnel,h)

       end do
    end do


end subroutine draw_funnel




function draw_plate(x,y,plate,h)
 
  real(kind=rk), intent(in)          :: x, y, h
  type(type_funnel_plate),intent(in) ::plate

  real(kind=rk)                      ::draw_plate,delta_r

  delta_r     = plate%r_out-plate%r_in
  draw_plate  = soft_bump(x,plate%x0(1),plate%width,h)*(soft_bump(y,plate%x0(2)+plate%r_in,delta_r,h) &
              + soft_bump(y,plate%x0(2)-(plate%r_out),delta_r,h))

end function draw_plate


function draw_walls(x,y,funnel,h)
 
  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::  mask, draw_walls

  mask=0

  if (abs(x-funnel%pump_x_center)> funnel%pump_diameter/2) then
         ! wall in south
         mask=mask+smoothstep(y-funnel%wall_thickness,h)
         ! wall in north
         mask=mask+smoothstep(domain_size(2)-y-funnel%wall_thickness,h)
  endif

  if (abs(y-domain_size(2)/2)> funnel%jet_radius) then
         ! wall in EAST
         mask=mask+smoothstep(x-funnel%wall_thickness,h)
  endif

  if (abs(y-domain_size(2)/2)> funnel%min_inner_diameter/2) then
         ! wall in WEST
         mask=mask+smoothstep(domain_size(1)-x-funnel%wall_thickness,h)
  endif

   ! is needed because mask off walls overlap
  if (mask>1) then
         mask=1
  endif

  draw_walls=mask
end function draw_walls




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













!==========================================================================
    !> \brief initialization of all plates in the funnel
    !> \todo insert picture 
    subroutine init_plates(funnel)
      
      implicit none
      !> geometric parameters of ion funnel
      type(type_funnel), intent(inout) :: funnel

      real(kind=rk)                   :: distance,length,length_focus,width
      integer(kind=ik)                :: n
      type(type_funnel_plate)         :: plate

      allocate(funnel%plate(funnel%nr_plates)) 
      
      distance    =funnel%plates_distance
      length      =funnel%length
      width       =funnel%plates_thickness
      ! length of focus area in funnel
      length_focus           = length - (funnel%max_inner_diameter-funnel%min_inner_diameter)/funnel%slope
      ! origin of funnel
      funnel%offset=(/ domain_size(1)-length-funnel%wall_thickness-distance+width, &
                       domain_size(2)/2/)
      if(funnel%offset(1)<funnel%wall_thickness) then
       call abort(13457,'Error [module_mask.f90]: your funnel is to long')
      endif

      ! initialicd all plates
      do n=1,funnel%nr_plates
        plate%x0(1)     = funnel%offset(1)+(n-1)*distance
        plate%x0(2)     = funnel%offset(2)
        plate%width     = width
        if (plate%x0(1)-funnel%offset(1)<=length_focus) then
           plate%r_in    = funnel%max_inner_diameter/2
        else
           plate%r_in = plate%r_in - funnel%slope*distance/2
        endif
        plate%r_out   =funnel%outer_diameter/2
        funnel%plate(n)=plate   
      enddo
      
    end subroutine init_plates
!==========================================================================


end module module_mask
