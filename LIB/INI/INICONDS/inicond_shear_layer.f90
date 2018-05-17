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

subroutine inicond_shear_layer( params, u, x0, dx )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)       :: params

    !> actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)            :: u(:,:,:,:)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3),dx(1:3)

    ! variable for shear layer position
    real(kind=rk)                           :: mux1, mux2, muy, x, y, sigma, w

    ! loop variables
    integer(kind=ik)                        :: ix, iy, dF, pF, rhoF, UxF, UyF, UzF

    ! grid
    integer(kind=ik)                        :: Bs, g

    ! p0 value \todo get from ini file
    real(kind=rk)                           :: p0, rho0, u0

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs   = params%number_block_nodes
    g    = params%number_ghost_nodes

    pF   = 0
    rhoF = 0
    UxF  = 0
    UyF  = 0
    UzF  = 0

    rho0 = 1.0_rk
    p0   = 2.5_rk
    u0   = 0.0_rk

    ! find fields
    do dF = 1, params%number_data_fields
        if ( params%physics_ns%names(dF) == "p" ) pF = dF
        if ( params%physics_ns%names(dF) == "rho" ) rhoF = dF
        if ( params%physics_ns%names(dF) == "Ux" ) UxF = dF
        if ( params%physics_ns%names(dF) == "Uy" ) UyF = dF
        if ( params%physics_ns%names(dF) == "Uz" ) UzF = dF
    end do

!---------------------------------------------------------------------------------------------
! main body

    ! place layer
    mux1 = params%Lx/2.0_rk - 0.25_rk
    mux2 = params%Lx/2.0_rk + 0.25_rk

    muy  = 0.5_rk * params%Ly

    ! boundary layer width
    sigma = params%inicond_width * params%Lx

    ! shear layer width
    w = params%inicond_width

    if (params%threeD_case) then
        ! nothing to do, set heavy data to zero
        u(:, :, :, :) = 0.0_rk

    else
        ! 2D case
        ! create shear layer, Uy field
        do ix = g+1,Bs+g
            do iy = g+1,Bs+g

                ! compute x,y coordinates from spacing and origin
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                ! shear layer, setup [1]
                ! Uy
                if ( x <= 0.5_rk*params%Lx ) then
                    u(ix, iy, 1, UyF) = dtanh( w/params%Lx * ( x - mux1 ) ) + u0
                    u(ix, iy, 1, rhoF) = ( rho0 + dtanh( w * ( x - mux1 ) ) ) / 2.0_rk + rho0
                else
                    u(ix, iy, 1, UyF) = dtanh( w/params%Lx * ( mux2 - x ) ) + u0
                    u(ix, iy, 1, rhoF) = ( rho0 + dtanh( w * ( mux2 - x ) ) ) / 2.0_rk + rho0
                end if

                ! Ux
                u(ix, iy, 1, UxF) = 0.01_rk * dsin( 8.0_rk * pi * ( y - muy  ) / params%Ly )

            end do
        end do

        ! set p field
        u(:, :, 1, pF)   = p0

        u(:, :, 1, rhoF) = dsqrt(u(:, :, 1, rhoF))
        u(:, :, 1, UxF)  = u(:, :, 1, UxF) * (u(:, :, 1, rhoF))
        u(:, :, 1, UyF)  = u(:, :, 1, UyF) * (u(:, :, 1, rhoF))

    endif

end subroutine inicond_shear_layer

