
!-----------------------------------------------------------------
!> \file
!> \brief
!! Module of mask function \f$\chi(x,t)\f$
!> \details
!> This module implements the mask function on each Block
!!          \f[
!!                           \chi(x,t)\quad \forall\;  x\in\mathcal{B}^l
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
  !------------------------------------------------
  ! available geometrys
  type :: cylinder
      real(kind=rk), dimension(2) :: x_cntr 
      real(kind=rk)               :: radius  
  end type cylinder


  type(cylinder), save :: cyl
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
    select case(mask_geometry)
    case ('cylinder')
      call read_param_mpi(FILE, 'VPM', 'x_cntr', cyl%x_cntr,(/0.5_rk, 0.5_rk /) )
      call read_param_mpi(FILE, 'VPM', 'radius', cyl%radius, 0.5_rk )
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

end module module_mask
