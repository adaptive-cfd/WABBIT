
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
  PUBLIC :: init_mask,get_mask
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
      real(kind=rk)       ::dout     ! outer diameter
      real(kind=rk)       ::dmax     ! maximal inner diameter
      real(kind=rk)       ::dmin     ! minimal inner diameter
      integer(kind=ik)    ::Nplate  ! Number of plates
      real(kind=rk)       ::s        ! seperation between plates
      real(kind=rk)       ::ds       ! slope of funnel
      real(kind=rk)       ::offset(2)   ! offset of funnel in x and y
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
    
      call read_param_mpi(FILE, 'funnel', 'outer_diameter'        , funnel%dout, domain_size(2)/2 )
      call read_param_mpi(FILE, 'funnel', 'maximal_inner_diameter', funnel%dmax, domain_size(2)/3 )
      call read_param_mpi(FILE, 'funnel', 'minimal_inner_diameter', funnel%dmin, domain_size(2)/4 )
      call read_param_mpi(FILE, 'funnel', 'Number_of_plates'      , funnel%Nplate, 30 )
      call read_param_mpi(FILE, 'funnel', 'diameter_per_plate'    , funnel%ds, domain_size(2)/30 )
     
      funnel%s      = domain_size(1)*0.9_rk/funnel%Nplate
      funnel%ds     = funnel%ds/funnel%s
      funnel%offset = (/0.15_rk*domain_size(1), domain_size(2)/2.0_rk/)
      write(*,*) "s=", funnel%s
      write(*,*) "ds=", funnel%ds
      write(*,*) "offset=", funnel%offset 
      
      call init_plates(funnel);
    
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
               call smoothstep(mask(ix,iy), r, cyl%radius, h)
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
    subroutine smoothstep(f,x,t,h)
      
      implicit none
      real(kind=rk), intent(out) :: f
      real(kind=rk), intent(in)  :: x,t,h

      !-------------------------------------------------
      ! cos shaped smoothing (compact in phys.space)
      !-------------------------------------------------
      if (x<=t-h) then
        f = 1.0_rk 
      elseif (((t-h)<x).and.(x<(t+h))) then
        f = 0.5_rk * (1.0_rk + dcos((x-t+h) * pi / (2.0_rk*h)) )
      else
        f = 0.0_rk
      endif

    end subroutine smoothstep
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
           ! distance from center of cylinder
            do n=1,funnel%Nplate
                if (x>(funnel%plate(n)%x0(1)-funnel%s/2) .and. &
                    x<(funnel%plate(n)%x0(1)+funnel%plate(n)%width+funnel%s/2)) then
                    mask(ix,iy)=mask(ix,iy)+draw_plate(x,y,funnel%plate(n),h)
                
               endif
            enddo
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






function soft_bump(x,x0,width,h)

  real(kind=rk), intent(in)      :: x, x0, h, width
  real(kind=rk)                  :: soft_bump,d 
    
    d=x-x0;

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

      real(kind=rk)                    :: s,LF,L1,width
      integer(kind=ik)                 :: n
      type(type_funnel_plate)          :: plate


      allocate(funnel%plate(funnel%Nplate)) 
      !-------------------------------------------------
      s       =funnel%s
      ! total length of funnel
      LF      =funnel%Nplate*s 
      L1      = LF - (funnel%dmax-funnel%dmin)/funnel%ds
      width   =s/2
      

      do n=1,funnel%Nplate
        plate%x0(1)     = funnel%offset(1)+(n-1)*s
        plate%x0(2)     = funnel%offset(2)
        plate%width     = width
        if (plate%x0(1)-funnel%offset(1)<=L1) then
           plate%r_in    = funnel%dmax/2
        else
           plate%r_in = plate%r_in - funnel%ds*s/2
        endif
        plate%r_out   =funnel%dout/2
        funnel%plate(n)=plate      
      enddo
      
    end subroutine init_plates
!==========================================================================


end module module_mask
