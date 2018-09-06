
!> \brief initialize gauss pulse for 2D case \n
!> \note field phi is 3D, but third dimension is not used
!
!> \details
!!
!! \date 04/11/16
!!          - switch to v0.4
!!
!! \date 26/01/17
!!          - use process rank from params struct
!!          - use v0.5 hvy data array
!!
!! \date 04/04/17
!!          - rewrite to work only on blocks, no large datafield required
!! \date 03/02/18
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
    real(kind=rk)                           :: x_cntr(3), x, z ,y, sigma
    ! loop variables
    integer(kind=ik)                        :: ix, iy, iz


!---------------------------------------------------------------------------------------------
! main body

    x_cntr = 0.5_rk * L


    if (size(u,3)==1) then
      ! pulse width
      sigma = width * L(1) * L(2)
      ! 2D case
      ! create gauss pulse
      do ix = g+1,Bs+g
        do iy = g+1,Bs+g
          ! compute x,y coordinates from spacing and origin
          x = dble(ix-(g+1)) * dx(1) + x0(1)
          y = dble(iy-(g+1)) * dx(2) + x0(2)
          ! shift to new gauss blob center
          call continue_periodic(x,L(1))
          call continue_periodic(y,L(2))
          ! set actual inicond gauss blob
          u(ix,iy,1) = dexp( -( (x-x_cntr(1))**2 + (y-x_cntr(2))**2 ) / sigma )
        end do
      end do

    else
      sigma = width*minval(L)
      ! 3D case
      ! create gauss pulse
      do ix = g+1,Bs+g
        do iy = g+1,Bs+g
          do iz = g+1,Bs+g

            x = dble(ix-(g+1)) * dx(1) + x0(1)
            y = dble(iy-(g+1)) * dx(2) + x0(2)
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            ! shift to new gauss blob center
            ! set actual inicond gauss blob
            u(ix,iy,iz) = dexp( -( (x-x_cntr(1))**2 + (y-x_cntr(2))**2 +(z-x_cntr(3))**2 ) / (2*sigma)**2 )
          end do


        end do
      end do

  endif

end subroutine inicond_gauss_blob
