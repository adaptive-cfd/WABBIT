!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_gauss_blob.f90
!> \version 0.5
!> \author engels, msr
!
!> \brief initialize gauss pulse for 2D case \n
!> \note field phi is 3D, but third dimension is not used
!
!> \details
!! input:    - params \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 04/11/16
!!          - switch to v0.4
!!
!! 26/01/17
!!          - use process rank from params struct
!!          - use v0.5 hvy data array
!!
!! 04/04/17
!!          - rewrite to work only on blocks, no large datafield required
!! 03/02/18
!!          - included into new physics modules
! ********************************************************************************************


subroutine inicond_gauss_blob(width, Bs, g, L, u, x0, dx )

    implicit none
    !> actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)   :: u(:,:,:)
    !> spacing and origin of block, domainlength, width of gauss_blob   
    real(kind=rk), intent(in)      :: x0(1:3),dx(1:3),L(1:3),width
    ! grid
    integer(kind=ik),intent(in)    :: Bs, g

    ! auxiliary variable for gauss pulse
    real(kind=rk)                           :: mux, muy, muz, x, z ,y, sigma
    ! loop variables
    integer(kind=ik)                        :: ix, iy, iz


!---------------------------------------------------------------------------------------------
! main body

    ! place pulse in the center of the domain
    mux = 0.5_rk * L(1)
    muy = 0.5_rk * L(2)
    muz = 0.5_rk * L(3)

    ! pulse width
    sigma = width * L(1) * L(2)

    if (size(u,3)==1) then
      ! 2D case
      ! create gauss pulse
      do ix = g+1,Bs+g
        do iy = g+1,Bs+g
          ! compute x,y coordinates from spacing and origin
          x = dble(ix-(g+1)) * dx(1) + x0(1)
          y = dble(iy-(g+1)) * dx(2) + x0(2)
          ! shift to new gauss blob center
          call shift_x_y( x, y, L(1),L(2) )
          ! set actual inicond gauss blob
          u(ix,iy,1) = dexp( -( (x-mux)**2 + (y-muy)**2 ) / sigma )
        end do
      end do

    else
      sigma = width*L(1)
      ! 3D case
      ! create gauss pulse
      do ix = g+1,Bs+g
        do iy = g+1,Bs+g
          do iz = g+1,Bs+g
            ! compute x,y coordinates from spacing and origin
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            y = dble(iy-(g+1)) * dx(2) + x0(2)
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            ! shift to new gauss blob center
            ! call shift_x_y( x, y, params%Lx,params%Ly )
            ! set actual inicond gauss blob
            u(ix,iy,iz) = dexp( -( (x-mux)**2 + (y-muy)**2 +(z-muz)**2 ) / sigma )
          end do
        end do
      end do

  endif

end subroutine inicond_gauss_blob

! function to ensure periodicity:
! shift center of gauss blob to center of computational domain
! if then point(x,y) outside domain -> set coordinates to (periodic) interior point
subroutine shift_x_y( x, y, Lx, Ly )

    implicit none

    ! coordinates
    real(kind=rk), intent(inout)   :: x, y
    ! domain size
    real(kind=rk), intent(in)       :: Lx, Ly

    x = x
    y = y - 0.25_rk

    ! check boundary
    if ( y < 0.0_rk ) then
        y = Ly + y
    end if

    if ( x < 0.0_rk ) then
        x = Lx + x
    end if

    if ( y > Ly ) then
        y = Ly - y
    end if

    if ( x > Lx ) then
        x = Lx - x
    end if

end subroutine shift_x_y
