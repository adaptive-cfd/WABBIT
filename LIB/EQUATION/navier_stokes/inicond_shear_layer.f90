!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_shear_layer.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize shear layer setup, setup for all datafields here, \todo: get setup values from ini file
!> setup [1]:
!
!>
!! input:    - params \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 16/02/17 - create
!! 08/08/17 - rework to use new heavy data (coordinates) structure
!
! ********************************************************************************************

subroutine inicond_shear_layer(  u, x0, dx ,Bs, g)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)            :: u(:,:,:,:)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3),dx(1:3)

    ! grid
    integer(kind=ik),intent(in)             :: Bs, g

    ! variable for shear layer position
    real(kind=rk)                           :: mux1, mux2, muy, x, y, sigma, w

    ! loop variables
    integer(kind=ik)                        :: ix, iy


    ! p0 value \todo get from ini file
    real(kind=rk)                           :: p0, rho0, u0

!---------------------------------------------------------------------------------------------
! variables initialization


    rho0 = params_ns%initial_density
    p0   = params_ns%initial_pressure
    u0   = params_ns%initial_velocity(1)

!---------------------------------------------------------------------------------------------
! main body

    ! place layer
    mux1 = params_ns%Lx/2.0_rk - 0.25_rk
    mux2 = params_ns%Lx/2.0_rk + 0.25_rk

    muy  = 0.5_rk * params_ns%Ly

    ! boundary layer width
    sigma = params_ns%inicond_width * params_ns%Lx

    ! shear layer width
    w = params_ns%inicond_width

    if (size(u,3)>1) then
        ! nothing to do, set heavy data to zero
        u(:, :, :, :) = 0.0_rk

    else
        ! 2D case
        ! create shear layer, Uy field
        do ix = 1,Bs+2*g
            do iy = 1,Bs+2*g

                ! compute x,y coordinates from spacing and origin
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                ! shear layer, setup [1]
                ! Uy
                if ( x <= 0.5_rk*params_ns%Lx ) then
                    u(ix, iy, 1, UyF) = dtanh( w/params_ns%Lx * ( x - mux1 ) ) + u0
                    u(ix, iy, 1, rhoF) = ( rho0 + dtanh( w * ( x - mux1 ) ) ) / 2.0_rk + rho0
                else
                    u(ix, iy, 1, UyF) = dtanh( w/params_ns%Lx * ( mux2 - x ) ) + u0
                    u(ix, iy, 1, rhoF) = ( rho0 + dtanh( w * ( mux2 - x ) ) ) / 2.0_rk + rho0
                end if

                ! Ux
                !u(ix, iy, 1, UxF) = 0.05_rk * dsin( 8.0_rk * pi * ( y - muy  ) / params%Ly )
                u(ix, iy, 1, UxF) = 0.1_rk * dsin( 2.0_rk * pi * ( y - muy  ) )
            end do
        end do

        ! set p field
        u(:, :, 1, pF)   = p0

        !
        !u(:, :, 1, rhoF) = u(:, :, 1, rhoF)
        u(:, :, 1, rhoF) = dsqrt(u(:, :, 1, rhoF))
        u(:, :, 1, UxF)  = u(:, :, 1, UxF) * (u(:, :, 1, rhoF))
        u(:, :, 1, UyF)  = u(:, :, 1, UyF) * (u(:, :, 1, rhoF))

    endif

end subroutine inicond_shear_layer
