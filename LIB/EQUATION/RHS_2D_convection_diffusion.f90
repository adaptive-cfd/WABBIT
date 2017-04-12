! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: RHS_2D_convection_diffusion.f90
! version: 0.4
! author: msr, engels
!
! rhs for 2D convection diffusion equation
!
! input:    - datafield, grid parameter, velocity, diffusion coefficient, derivative order
! output:   - RHS(datafield)
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine RHS_2D_convection_diffusion(phi, dx, dy, g, Bs, u01, u02, nu, order_discretization)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! grid parameter
    integer(kind=ik), intent(in)                                :: g, Bs
    ! rhs parameter
    real(kind=rk), intent(in)                                   :: dx, dy, u01, u02, nu
    ! datafield
    real(kind=rk), dimension(Bs+2*g, Bs+2*g), intent(inout)     :: phi
    ! order discretization
    character(len=80)                                           :: order_discretization

    ! auxiliary fields for rhs calculation
    real(kind=rk), dimension(Bs+2*g, Bs+2*g)                    :: rhs
    ! auxiliary variables
    real(kind=rk)                                               :: phi_dx, phi_dy, phi_dxdx, phi_dydy, dx_inv, dy_inv, dx2_inv, dy2_inv
    real(kind=rk)                                               :: a(-3:+3), b1, b2, b3 ,b4 ,b5
    ! loop variables
    integer                                                     :: ix, iy

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    rhs = 0.0_rk

    ! division is expensive, multiplication is cheap, so here we save some time
    dx_inv = 1.0_rk / (dx)
    dy_inv = 1.0_rk / (dy)
    dx2_inv = 1.0_rk / dx**2.0_rk
    dy2_inv = 1.0_rk / dy**2.0_rk

    ! Tam & Webb, 4th order optimized (for first derivative)
    a=(/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
         0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    ! 4th order coefficients for second derivative
    b1 = -1.0_rk/12.0_rk
    b2 = 4.0_rk/3.0_rk
    b3 = -5.0_rk/2.0_rk
    b4 = 4.0_rk/3.0_rk
    b5 = -1.0_rk/12.0_rk

!---------------------------------------------------------------------------------------------
! main body

    if (order_discretization == "FD_2nd_central" ) then
      !-----------------------------------------------------------------------
      ! 2nd order
      !-----------------------------------------------------------------------
      ! loop over interior points only (EXCLUDE ghost nodes)
      do ix = g+1, Bs+g
        do iy = g+1, Bs+g
          phi_dx = (phi(ix+1,iy)-phi(ix-1,iy))/(2.0_rk*dx)
          phi_dy = (phi(ix,iy+1)-phi(ix,iy-1))/(2.0_rk*dy)
          phi_dxdx = (phi(ix-1,iy)-2.d0*phi(ix,iy)+phi(ix+1,iy))*dy2_inv
          phi_dydy = (phi(ix,iy-1)-2.d0*phi(ix,iy)+phi(ix,iy+1))*dy2_inv
          ! compute (assemble) final right hand side
          rhs(ix,iy) = - u01 * phi_dx - u02 * phi_dy &
                       + nu * ( phi_dxdx + phi_dydy )
        end do
      end do

    elseif (order_discretization == "FD_4th_central_optimized") then
      !-----------------------------------------------------------------------
      ! 4th order
      !-----------------------------------------------------------------------
      ! loop over interior points only (EXCLUDE ghost nodes)
      do ix = g+1, Bs+g
        do iy = g+1, Bs+g
          phi_dx = (a(-3)*phi(ix-3,iy) + a(-2)*phi(ix-2,iy) + a(-1)*phi(ix-1,iy) + a(0)*phi(ix,iy)&
                 +  a(+3)*phi(ix+3,iy) + a(+2)*phi(ix+2,iy) + a(+1)*phi(ix+1,iy))*dx_inv
          phi_dy = (a(-3)*phi(ix,iy-3) + a(-2)*phi(ix,iy-2) + a(-1)*phi(ix,iy-1) + a(0)*phi(ix,iy)&
                 +  a(+3)*phi(ix,iy+3) + a(+2)*phi(ix,iy+2) + a(+1)*phi(ix,iy+1))*dy_inv
          phi_dxdx = (b1*phi(ix-2,iy) + b2*phi(ix-1,iy) + b3*phi(ix,iy)&
                   +  b4*phi(ix+1,iy) + b5*phi(ix+2,iy))*dx2_inv
          phi_dydy = (b1*phi(ix,iy-2) + b2*phi(ix,iy-1) + b3*phi(ix,iy)&
                   +  b4*phi(ix,iy+1) + b5*phi(ix,iy+2))*dy2_inv

          ! x =
          ! y =
          ! u01 =
          ! u02 =
          ! compute (assemble) final right hand side
          rhs(ix,iy) = - u01 * phi_dx - u02 * phi_dy &
                       + nu * ( phi_dxdx + phi_dydy )
        end do
      end do

    else
      write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
      write(*,*) order_discretization
      stop
    end if

    ! return (TODO: DO NOT OVERWRITE?)
    phi = rhs

    do ix = 1, Bs+2*g
      do iy = 1, Bs+2*g
        if (abs(phi(ix,iy)) > 1.0e+6_rk) then
          call error_msg("very large values of phi -> we stop here to let you think.")
        endif
      end do
    end do

end subroutine RHS_2D_convection_diffusion
