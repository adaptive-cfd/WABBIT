!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_advection.f90
!> TEST!!!
!
! ********************************************************************************************

subroutine RHS_3D_advection(rhs, phi, xx0, ddx, g, Bs, time, order_discretization)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)                    :: g, Bs
    !> rhs parameter
    real(kind=rk), intent(in)                       :: xx0(1:3), ddx(1:3), time
    !> datafield
    real(kind=rk), intent(in)                       :: phi(Bs+2*g, Bs+2*g, Bs+2*g)
    !> order discretization
    character(len=80)                               :: order_discretization
    !> rhs field
    real(kind=rk), intent(inout)                    :: rhs(Bs+2*g, Bs+2*g, Bs+2*g)

    ! auxiliary variables
    real(kind=rk)                                   :: dx_inv, dy_inv, dz_inv
    real(kind=rk)                                   :: a(-3:+3), b1, b2, b3 ,b4 ,b5
    ! loop variables
    integer                                         :: ix, iy, iz
    ! coordinates
    real(kind=rk)                                   :: x, y, z
    ! swirl t end
    real(kind=rk)                                   :: t_a

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set t end, note: in order to use multiplication we actually set 1/t_a
    t_a = 1.0_rk/3.0_rk

    if (order_discretization == "FD_2nd_central" ) then
        ! division is expensive, multiplication is cheap, so here we save some time
        dx_inv = 1.0_rk / (2.0_rk*ddx(1))
        dy_inv = 1.0_rk / (2.0_rk*ddx(2))
        dz_inv = 1.0_rk / (2.0_rk*ddx(3))
    else
        ! division is expensive, multiplication is cheap, so here we save some time
        dx_inv = 1.0_rk / (ddx(1))
        dy_inv = 1.0_rk / (ddx(2))
        dz_inv = 1.0_rk / (ddx(3))
    end if

    if (order_discretization == "FD_4th_central_optimized") then
        ! Tam & Webb, 4th order optimized (for first derivative)
        a=(/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
             0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

        ! 4th order coefficients for second derivative
        b1 = -1.0_rk/12.0_rk
        b2 = 4.0_rk/3.0_rk
        b3 = -5.0_rk/2.0_rk
        b4 = 4.0_rk/3.0_rk
        b5 = -1.0_rk/12.0_rk
    end if

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

                ! calculate coordinates form ix, iy
                call coords_ix_iy_iz( x, y, z, ix-g, iy-g, iz-g, ddx, xx0 )

                ! compute (assemble) final right hand side
                rhs(ix,iy,iz) = - (dcos((pi*time)*t_a) * (dsin(pi*(x)))**2 &
                                * dsin(2.0_rk*pi*(y))) &
                                * (phi(ix+1,iy,iz)-phi(ix-1,iy,iz))*dx_inv &
                                - (dcos((pi*time)*t_a) * (dsin(pi*(y)))**2 &
                                * (-dsin(2.0_rk*pi*(x)))) &
                                * (phi(ix,iy+1,iz)-phi(ix,iy-1,iz))*dy_inv &
                                - (-dsin((2.0_rk*pi*time)*t_a)) &
                                * (phi(ix,iy,iz+1)-phi(ix,iy,iz-1))*dz_inv

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

                ! calculate coordinates form ix, iy
                call coords_ix_iy_iz( x, y, z, ix-g, iy-g, iz-g, ddx, xx0 )

                ! compute (assemble) final right hand side
                rhs(ix,iy,iz) = - (dcos((pi*time)*t_a) * (dsin( pi*( x ) ))**2 &
                                * dsin(2.0_rk*pi*( y ))) &
                                * (a(-3)*phi(ix-3,iy,iz) + a(-2)*phi(ix-2,iy,iz) + a(-1)*phi(ix-1,iy,iz) + a(0)*phi(ix,iy,iz) &
                                +  a(+3)*phi(ix+3,iy,iz) + a(+2)*phi(ix+2,iy,iz) + a(+1)*phi(ix+1,iy,iz))*dx_inv &
                                - (dcos((pi*time)*t_a) * (dsin(pi*( y )))**2 &
                                * (-dsin(2.0_rk*pi*( x )))) &
                                * (a(-3)*phi(ix,iy-3,iz) + a(-2)*phi(ix,iy-2,iz) + a(-1)*phi(ix,iy-1,iz) + a(0)*phi(ix,iy,iz) &
                                +  a(+3)*phi(ix,iy+3,iz) + a(+2)*phi(ix,iy+2,iz) + a(+1)*phi(ix,iy+1,iz))*dy_inv &
                                - (-dsin((2.0_rk*pi*time)*t_a)) &
                                * (a(-3)*phi(ix,iy,iz-3) + a(-2)*phi(ix,iy,iz-2) + a(-1)*phi(ix,iy,iz-1) + a(0)*phi(ix,iy,iz) &
                                +  a(+3)*phi(ix,iy,iz+3) + a(+2)*phi(ix,iy,iz+2) + a(+1)*phi(ix,iy,iz+1))*dz_inv

            end do
        end do
      end do

    else
      write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
      write(*,*) order_discretization
      stop
    end if

end subroutine RHS_3D_advection

! --------------------------------------------------------------------------------------------
! function to calculate x,y coordinates from loop variables and block origin/spacing
subroutine coords_ix_iy_iz( x, y, z, ix, iy, iz, ddx, xx0 )

    use module_params

    implicit none

    ! coordinates
    real(kind=rk), intent(inout)   :: x, y, z
    ! loop parameter
    integer(kind=ik), intent(in)   :: ix, iy, iz
    ! block origin/spacing
    real(kind=rk), intent(in)      :: xx0(1:3), ddx(1:3)

    ! calculate coordinates
    x = xx0(1) + real(ix-1, kind=rk) * ddx(1)
    y = xx0(2) + real(iy-1, kind=rk) * ddx(2)
    z = xx0(3) + real(iz-1, kind=rk) * ddx(3)

end subroutine coords_ix_iy_iz

!! --------------------------------------------------------------------------------------------
!! function to calculate velocity from x,y coordinates and time t
!subroutine f_x_y_z_t( u01, u02, u03, x, y, z, time )
!
!    use module_params
!
!    implicit none
!
!    ! coordinates
!    real(kind=rk), intent(inout)   :: u01, u02, u03
!    ! coordinates and time
!    real(kind=rk), intent(in)      :: x, y, z, time
!
!    ! swirl t end
!    real(kind=rk)                  :: t_a
!
!    ! set t end
!    t_a = 3.0_rk
!
!    ! calculate velocity
!    u01 = dcos((pi*time)/t_a) * (dsin(pi*x))**2 * dsin(2.0_rk*pi*y)
!    u02 = dcos((pi*time)/t_a) * (dsin(pi*y))**2 * (-dsin(2.0_rk*pi*x))
!    u03 = -dsin((2.0_rk*pi*time)/t_a) !0.0_rk
!
!end subroutine f_x_y_z_t
