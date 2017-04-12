! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: inicond_sphere.f90
! version: 0.5
! author: msr
!
! initialize sphere
!
! input:    - params
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 02/02/17 - create
! 12/04/17 - rewrite to work on one block and initialize with spacing and origin of block
!
! ********************************************************************************************

subroutine inicond_sphere( params, u, x0, dx)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)       :: params

    ! actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)            :: u(:,:,:,:)

    ! spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    ! auxiliary variable for gauss pulse
    real(kind=rk)                           :: mux, muy, muz, sigma, x ,y, z, r, w
    ! loop variables
    integer(kind=ik)                        :: ix, iy, iz

    ! grid
    integer(kind=ik)                        :: Bs, g

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs   = params%number_block_nodes
    g    = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

    ! place sphere in the center of the domain
    mux = 0.5_rk * params%Lx;
    muy = 0.5_rk * params%Ly;
    muz = 0.5_rk * params%Lz;

    ! sphere boundary (layer) width
    sigma     = 20.0_rk / params%Lx

    ! sphere width
    w = 0.0004_rk * params%Lx

    ! create gauss pulse
    do ix = g+1, Bs+g
        do iy = g+1, Bs+g
            do iz = g+1, Bs+g

                ! compute x,y,z coordinates from spacing and origin
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                z = dble(iz-(g+1)) * dx(3) + x0(3)

                ! r coordinate (sphere)
                r = sqrt( (x-mux)**2 + (y-muy)**2 + (z-muz)**2 )

                ! set actual inicond sphere
                ! FIXME: df=2 ...
                u(ix,iy,iz,2) = abs( 0.5_rk * (tanh( sigma * (r-w) ) - 1.0_rk) )

            end do
        end do
    end do

end subroutine inicond_sphere
