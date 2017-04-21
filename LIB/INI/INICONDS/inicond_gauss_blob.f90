! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: inicond_gauss_blob.f90
! version: 0.5
! author: engels, msr
!
! initialize gauss pulse for 2D case
! note: field phi is 3D, but third dimension is not used
!
! input:    - params
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! 26/01/17 - use process rank from params struct
!          - use v0.5 hvy data array
! 04/04/17 - rewrite to wokr only on blocks, no large datafield required
!
! ********************************************************************************************

subroutine inicond_gauss_blob( params, u, x0, dx )

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)    :: params
    ! actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout) :: u(:,:,:,:)
    ! spacing and origin of block
    real(kind=rk), intent(in) :: x0(1:3),dx(1:3)

    ! auxiliary variable for gauss pulse
    real(kind=rk)                           :: mux, muy, x ,y, sigma
    ! loop variables
    integer(kind=ik)                        :: ix, iy
    ! grid
    integer(kind=ik)                        :: Bs, g

!---------------------------------------------------------------------------------------------
! variables initialization
    Bs   = params%number_block_nodes
    g    = params%number_ghost_nodes


!---------------------------------------------------------------------------------------------
! main body

    ! place pulse in the center of the domain
    mux = 0.5_rk * params%Lx;
    !muy = 0.5_rk * params%Ly;
    muy = 0.75_rk * params%Ly;

    ! pulse width
    !sigma     = 0.1e-2_rk * params%Lx * params%Ly
    sigma     = 0.01

    ! create gauss pulse
    do ix = g+1,Bs+g
      do iy = g+1,Bs+g
        ! compute x,y coordinates from spacing and origin
        x = dble(ix-(g+1)) * dx(1) + x0(1)
        y = dble(iy-(g+1)) * dx(2) + x0(2)
        ! set actual inicond gauss blob
        ! FIXME: df=2 ...
        u(ix,iy,1,2) = dexp( -( (x-mux)**2 + (y-muy)**2 ) / sigma )
      end do
    end do

    ! ! it sometimes causes bizarre effects not to delete extremely small numbers:
    ! ! so we do that now.
    ! where ( u<1.0e-13_rk )
    !     u = 0.0_rk
    ! end where

end subroutine inicond_gauss_blob
