!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_2D_advection.f90
!> \version 0.5
!> \author msr
!
!> \brief rhs for 2D advection equation, calculate velocity field with f_x_y_t subroutine
!
!> \details
!> input:    - datafield, grid parameter, time, derivative order \n
!! output:   - RHS(datafield)
!! \n
!! = log ======================================================================================
!! \n
!! 20/04/17 - create
!
! ********************************************************************************************

subroutine RHS_2D_advection(phi, xx0, ddx, g, Bs, time, order_discretization)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)                                :: g, Bs
    !> rhs parameter
    real(kind=rk), intent(in)                                   :: xx0(1:2), ddx(1:2), time
    !> datafield
    real(kind=rk), dimension(Bs+2*g, Bs+2*g), intent(inout)     :: phi
    ! order discretization
    character(len=80)                                           :: order_discretization

    ! auxiliary fields for rhs calculation
    real(kind=rk), dimension(Bs+2*g, Bs+2*g)                    :: rhs
    ! auxiliary variables
    real(kind=rk)                                               :: phi_dx, phi_dy, dx_inv, dy_inv
    real(kind=rk)                                               :: a(-3:+3), b1, b2, b3 ,b4 ,b5
    ! loop variables
    integer                                                     :: ix, iy

    ! velocity components
    real(kind=rk)                                               :: u01, u02
    ! coordinates
    real(kind=rk)                                               :: x, y

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    rhs = 0.0_rk

    ! division is expensive, multiplication is cheap, so here we save some time
    dx_inv = 1.0_rk / (ddx(1))
    dy_inv = 1.0_rk / (ddx(2))

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

            ! calculate coordinates form ix, iy
            call coords_ix_iy( x, y, ix-g, iy-g, ddx, xx0 )

            ! calculate velocity, depends on x,y,t
            call f_x_y_t( u01, u02, x, y, time )

            ! gradient
            phi_dx = (phi(ix+1,iy)-phi(ix-1,iy))/(2.0_rk*ddx(1))
            phi_dy = (phi(ix,iy+1)-phi(ix,iy-1))/(2.0_rk*ddx(2))
            ! compute (assemble) final right hand side
            rhs(ix,iy) = - u01 * phi_dx - u02 * phi_dy

        end do
      end do

    elseif (order_discretization == "FD_4th_central_optimized") then
      !-----------------------------------------------------------------------
      ! 4th order
      !-----------------------------------------------------------------------
      ! loop over interior points only (EXCLUDE ghost nodes)
      do ix = g+1, Bs+g
        do iy = g+1, Bs+g

            ! calculate coordinates form ix, iy
            call coords_ix_iy( x, y, ix-g, iy-g, ddx, xx0 )

            ! calculate velocity, depends on x,y,t
            call f_x_y_t( u01, u02, x, y, time )

            ! gradient
            phi_dx = (a(-3)*phi(ix-3,iy) + a(-2)*phi(ix-2,iy) + a(-1)*phi(ix-1,iy) + a(0)*phi(ix,iy)&
                   +  a(+3)*phi(ix+3,iy) + a(+2)*phi(ix+2,iy) + a(+1)*phi(ix+1,iy))*dx_inv
            phi_dy = (a(-3)*phi(ix,iy-3) + a(-2)*phi(ix,iy-2) + a(-1)*phi(ix,iy-1) + a(0)*phi(ix,iy)&
                   +  a(+3)*phi(ix,iy+3) + a(+2)*phi(ix,iy+2) + a(+1)*phi(ix,iy+1))*dy_inv

            ! compute (assemble) final right hand side
            rhs(ix,iy) = - u01 * phi_dx - u02 * phi_dy

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

end subroutine RHS_2D_advection

! --------------------------------------------------------------------------------------------
! function to calculate x,y coordinates from loop variables and block origin/spacing
subroutine coords_ix_iy( x, y, ix, iy, ddx, xx0 )

    use module_params

    implicit none

    ! coordinates
    real(kind=rk), intent(inout)   :: x, y
    ! loop parameter
    integer(kind=ik), intent(in)   :: ix, iy
    ! block origin/spacing
    real(kind=rk), intent(in)      :: xx0(1:2), ddx(1:2)

    ! calculate coordinates
    x = xx0(1) + real(ix-1, kind=rk) * ddx(1)
    y = xx0(2) + real(iy-1, kind=rk) * ddx(2)

end subroutine coords_ix_iy

! --------------------------------------------------------------------------------------------
! function to calculate velocity from x,y coordinates and time t
subroutine f_x_y_t( u01, u02, x, y, time )

    use module_params

    implicit none

    ! coordinates
    real(kind=rk), intent(inout)   :: u01, u02
    ! coordinates and time
    real(kind=rk), intent(in)      :: x, y, time

    ! swirl t end
    real(kind=rk)                  :: t_a

    ! set t end
    t_a = 3.0_rk !10.0_rk

    ! calculate velocity
    !u01 = 1.0_rk
    !u02 = 0.5_rk
    u01 = cos((pi*time)/t_a) * (sin(pi*x))**2 * sin(2*pi*y)
    u02 = cos((pi*time)/t_a) * (sin(pi*y))**2 * (-sin(2*pi*x))

end subroutine f_x_y_t
