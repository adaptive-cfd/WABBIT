!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_convection_diffusion.f90
!> \version 0.5
!> \author msr, engels
!
!> \brief rhs for 3D convection diffusion equation
!
!>
!! input:    - datafield, grid parameter, velocity, diffusion coefficient, derivative order \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 02/02/17 - create
!
! ********************************************************************************************

subroutine RHS_3D_convection_diffusion(phi, dx, dy, dz, g, Bs, u01, u02, u03, nu, order_discretization)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)                                        :: g, Bs
    !> rhs parameter
    real(kind=rk), intent(in)                                           :: dx, dy, dz, u01, u02, u03, nu
    !> datafield
    real(kind=rk), dimension(Bs+2*g, Bs+2*g, Bs+2*g), intent(inout)     :: phi
    ! order discretization
    character(len=80)                                                   :: order_discretization

    ! auxiliary fields for rhs calculation
    real(kind=rk), dimension(Bs+2*g, Bs+2*g, Bs+2*g)                    :: rhs
    ! auxiliary variables
    real(kind=rk)                                                       :: phi_dx, phi_dy, phi_dz, phi_dxdx, phi_dydy, phi_dzdz, dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv
    real(kind=rk)                                                       :: a(-3:+3), b1, b2, b3 ,b4 ,b5
    ! loop variables
    integer                                                             :: ix, iy, iz

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    rhs = 0.0_rk

    ! division is expensive, multiplication is cheap, so here we save some time
    dx_inv = 1.0_rk / (dx)
    dy_inv = 1.0_rk / (dy)
    dz_inv = 1.0_rk / (dz)

    dx2_inv = 1.0_rk / dx**2.0_rk
    dy2_inv = 1.0_rk / dy**2.0_rk
    dz2_inv = 1.0_rk / dz**2.0_rk

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
            do iz = g+1, Bs+g

              phi_dx = (phi(ix+1,iy,iz)-phi(ix-1,iy,iz))/(2.0_rk*dx)
              phi_dy = (phi(ix,iy+1,iz)-phi(ix,iy-1,iz))/(2.0_rk*dy)
              phi_dz = (phi(ix,iy,iz+1)-phi(ix,iy,iz-1))/(2.0_rk*dz)

              phi_dxdx = (phi(ix-1,iy,iz)-2.d0*phi(ix,iy,iz)+phi(ix+1,iy,iz))*dy2_inv
              phi_dydy = (phi(ix,iy-1,iz)-2.d0*phi(ix,iy,iz)+phi(ix,iy+1,iz))*dy2_inv
              phi_dzdz = (phi(ix,iy,iz-1)-2.d0*phi(ix,iy,iz)+phi(ix,iy,iz+1))*dz2_inv

              ! compute (assemble) final right hand side
              rhs(ix,iy,iz) = - u01 * phi_dx - u02 * phi_dy - u03 * phi_dz &
                           + nu * ( phi_dxdx + phi_dydy + phi_dzdz )
          end do
        end do
      end do

    elseif (order_discretization == "FD_4th_central_optimized") then
      !-----------------------------------------------------------------------
      ! 4th order
      !-----------------------------------------------------------------------
      ! loop over interior points only (EXCLUDE ghost nodes)
      do ix = g+1, Bs+g
        do iy = g+1, Bs+g
            do iz = g+1, Bs+g

              phi_dx = (a(-3)*phi(ix-3,iy,iz) + a(-2)*phi(ix-2,iy,iz) + a(-1)*phi(ix-1,iy,iz) + a(0)*phi(ix,iy,iz)&
                     +  a(+3)*phi(ix+3,iy,iz) + a(+2)*phi(ix+2,iy,iz) + a(+1)*phi(ix+1,iy,iz))*dx_inv
              phi_dy = (a(-3)*phi(ix,iy-3,iz) + a(-2)*phi(ix,iy-2,iz) + a(-1)*phi(ix,iy-1,iz) + a(0)*phi(ix,iy,iz)&
                     +  a(+3)*phi(ix,iy+3,iz) + a(+2)*phi(ix,iy+2,iz) + a(+1)*phi(ix,iy+1,iz))*dy_inv
              phi_dz = (a(-3)*phi(ix,iy,iz-3) + a(-2)*phi(ix,iy,iz-2) + a(-1)*phi(ix,iy,iz-1) + a(0)*phi(ix,iy,iz)&
                     +  a(+3)*phi(ix,iy,iz+3) + a(+2)*phi(ix,iy,iz+2) + a(+1)*phi(ix,iy,iz+1))*dz_inv

              phi_dxdx = (b1*phi(ix-2,iy,iz) + b2*phi(ix-1,iy,iz) + b3*phi(ix,iy,iz)&
                       +  b4*phi(ix+1,iy,iz) + b5*phi(ix+2,iy,iz))*dx2_inv
              phi_dydy = (b1*phi(ix,iy-2,iz) + b2*phi(ix,iy-1,iz) + b3*phi(ix,iy,iz)&
                       +  b4*phi(ix,iy+1,iz) + b5*phi(ix,iy+2,iz))*dy2_inv
              phi_dzdz = (b1*phi(ix,iy,iz-2) + b2*phi(ix,iy,iz-1) + b3*phi(ix,iy,iz)&
                       +  b4*phi(ix,iy,iz+1) + b5*phi(ix,iy,iz+2))*dz2_inv

              ! compute (assemble) final right hand side
              rhs(ix,iy,iz) = - u01 * phi_dx - u02 * phi_dy - u03 * phi_dz &
                           + nu * ( phi_dxdx + phi_dydy + phi_dzdz )

            end do
        end do
      end do

    else
      write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
      write(*,*) order_discretization
      stop
    end if

    ! return 
    !> \todo DO NOT OVERWRITE?
    phi = rhs

end subroutine RHS_3D_convection_diffusion
