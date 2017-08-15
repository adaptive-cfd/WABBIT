!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_shear_layer.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize shear layer setup
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
    real(kind=rk)                           :: mux, muy, x, y, sigma, w

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

    ! place layer in the center of the domain
    mux = 0.5_rk * params%Lx
    muy = 0.5_rk * params%Ly

    ! boundary layer width
    sigma = params%inicond_width * params%Lx

    ! shear layer width
    ! \todo set shear layer width in ini file
    w = 0.15_rk * params%Lx

    if (params%threeD_case) then
        ! nothing to do, set heavy data to zero
        u(:, :, :, 1) = 0.0_rk

    else
        ! 2D case
        ! create shear layer
        do ix = g+1,Bs+g
            do iy = g+1,Bs+g

                ! compute x,y coordinates from spacing and origin
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                ! set actual inicond gauss blob
                u(ix, iy, 1, 1) = abs( 0.5_rk * (tanh( sigma * ( abs(x-mux) - w ) ) - 1.0_rk) )

            end do
        end do

        ! start disturbance
        do ix = g+1,Bs+g
            do iy = g+1,Bs+g
              ! compute x,y coordinates from spacing and origin
              x = dble(ix-(g+1)) * dx(1) + x0(1)
              y = dble(iy-(g+1)) * dx(2) + x0(2)
              ! set actual inicond gauss blob
              u(ix, iy, 1, 1) = u(ix, iy, 1, 1) + 0.001_rk * dexp( -( (x-mux)**2 + (y-muy)**2 ) / (0.01_rk*params%Lx) )
            end do
        end do

    endif

end subroutine inicond_shear_layer
