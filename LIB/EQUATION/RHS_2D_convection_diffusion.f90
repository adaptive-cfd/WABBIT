! ********************************
! WABBIT
! --------------------------------
!
! RHS, 2D, convection-diffusion equation
!
! name: RHS_2D_convection_diffusion.f90
! date: 27.10.2016
! author: msr, engels
! version: 0.3
!
! ********************************

subroutine RHS_2D_convection_diffusion(phi, dx, dy, g, N)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)					            :: g, N
    real(kind=rk), intent(in)                                   :: dx, dy
    real(kind=rk), dimension(N+2*g, N+2*g), intent(inout)       :: phi

    real(kind=rk), dimension(N+2*g, N+2*g)			            :: grad_phi, laplace_phi, rhs
    real(kind=rk)                                               :: phi_dx, phi_dy, phi_dxdx, phi_dydy, dx_inv, dy_inv, dx2_inv, dy2_inv
    real(kind=rk)                                               :: a(-3:+3),b1,b2,b3,b4,b5

    integer                                                     :: ix,iy

    grad_phi 		= 0.0_rk
    laplace_phi	    = 0.0_rk
    rhs             = 0.0_rk

    ! division is expensive, multiplication is cheap, so here we save some time
    dx_inv = 1.0_rk / (dx)
    dy_inv = 1.0_rk / (dy)
    dx2_inv = 1.0_rk / dx**2
    dy2_inv = 1.0_rk / dy**2

    ! Tam & Webb, 4th order optimized (for first derivative)
    a=(/-0.02651995d0, +0.18941314d0, -0.79926643d0, 0.0d0, &
         0.79926643d0, -0.18941314d0, 0.02651995d0/)

    ! 4th order coefficients for second derivative
    b1 = -1.0_rk/12.0_rk
    b2 = 4.0_rk/3.0_rk
    b3 = -5.0_rk/2.0_rk
    b4 = 4.0_rk/3.0_rk
    b5 = -1.0_rk/12.0_rk

    if (params%order_discretization == "FD_2nd_central" ) then
      !-----------------------------------------------------------------------
      ! 2nd order
      !-----------------------------------------------------------------------
      ! loop over interior points only (EXCLUDE ghost nodes)
      do ix = g+1, N+g
        do iy = g+1, N+g
          phi_dx = (phi(ix+1,iy)-phi(ix-1,iy))/(2.0_rk*dx)
          phi_dy = (phi(ix,iy+1)-phi(ix,iy-1))/(2.0_rk*dy)
          phi_dxdx = (phi(ix-1,iy)-2.d0*phi(ix,iy)+phi(ix+1,iy))*dy2_inv
          phi_dydy = (phi(ix,iy-1)-2.d0*phi(ix,iy)+phi(ix,iy+1))*dy2_inv
          ! compute (assemble) final right hand side
          rhs(ix,iy) = - params%u0(1) * phi_dx - params%u0(2) * phi_dy &
                       + params%nu * ( phi_dxdx + phi_dydy )
        end do
      end do

    elseif (params%order_discretization == "FD_4th_central_optimized") then
      !-----------------------------------------------------------------------
      ! 4th order
      !-----------------------------------------------------------------------
      ! loop over interior points only (EXCLUDE ghost nodes)
      do ix = g+1, N+g
        do iy = g+1, N+g
          phi_dx = (a(-3)*phi(ix-3,iy) + a(-2)*phi(ix-2,iy) + a(-1)*phi(ix-1,iy) + a(0)*phi(ix,iy)&
                 +  a(+3)*phi(ix+3,iy) + a(+2)*phi(ix+2,iy) + a(+1)*phi(ix+1,iy))*dx_inv
          phi_dy = (a(-3)*phi(ix,iy-3) + a(-2)*phi(ix,iy-2) + a(-1)*phi(ix,iy-1) + a(0)*phi(ix,iy)&
                 +  a(+3)*phi(ix,iy+3) + a(+2)*phi(ix,iy+2) + a(+1)*phi(ix,iy+1))*dy_inv
          phi_dxdx = (b1*phi(ix-2,iy) + b2*phi(ix-1,iy) + b3*phi(ix,iy)&
                   +  b4*phi(ix+1,iy) + b5*phi(ix+2,iy))*dx2_inv
          phi_dydy = (b1*phi(ix,iy-2) + b2*phi(ix,iy-1) + b3*phi(ix,iy)&
                   +  b4*phi(ix,iy+1) + b5*phi(ix,iy+2))*dy2_inv

          ! compute (assemble) final right hand side
          rhs(ix,iy) = - params%u0(1) * phi_dx - params%u0(2) * phi_dy &
                       + params%nu * ( phi_dxdx + phi_dydy )
        end do
      end do

    else
      write(*,*) "discretization method in params%order_discretization is unknown"
      write(*,*) params%order_discretization
      stop
    end if

    ! return (TODO: DO NOT OVERWRITE?)
    phi = rhs

end subroutine RHS_2D_convection_diffusion
