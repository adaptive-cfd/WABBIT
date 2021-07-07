!> \brief computation of the divergence from a given velocity field
!-----------------------------------------------------------------------------------------------------
subroutine divergence(u, v, w, dx, Bs, g, discretization, div)
    implicit none
    real(kind=rk), dimension(3), intent(in)        :: dx                        !> origin and spacing of the block
    real(kind=rk), dimension(:,:,:), intent(in)    :: u, v, w                   !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: div                       !> vorticity
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    real(kind=rk)                                  :: u_dx, v_dy, w_dz          !> derivatives
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables
    ! coefficients for Tam&Webb (4th order 1st derivative)
    real(kind=rk), parameter :: a(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)

    div = 0.0_rk
    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)


    if (size(u,3)>2) then
        !-----------------------------------------------------------------------
        ! 3D case 3D case 3D case 3D case
        !-----------------------------------------------------------------------
        dz_inv = 1.0_rk / dx(3)
        select case(discretization)
        case("FD_2nd_central")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
                        u_dx = (u(ix+1,iy,iz)-u(ix-1,iy,iz))*dx_inv*0.5_rk
                        v_dy = (v(ix,iy+1,iz)-v(ix,iy-1,iz))*dy_inv*0.5_rk
                        w_dz = (w(ix,iy,iz+1)-w(ix,iy,iz-1))*dz_inv*0.5_rk

                        div(ix,iy,iz) = u_dx + v_dy + w_dz
                    end do
                end do
            end do
        case("FD_4th_central")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
                        u_dx = (a_FD4(-2)*u(ix-2,iy,iz) +a_FD4(-1)*u(ix-1,iy,iz) +a_FD4(+1)*u(ix+1,iy,iz) + a_FD4(+2)*u(ix+2,iy,iz))*dx_inv
                        v_dy = (a_FD4(-2)*v(ix,iy-2,iz) +a_FD4(-1)*v(ix,iy-1,iz) +a_FD4(+1)*v(ix,iy+1,iz) + a_FD4(+2)*v(ix,iy+2,iz))*dy_inv
                        w_dz = (a_FD4(-2)*w(ix,iy,iz-2) +a_FD4(-1)*w(ix,iy,iz-1) +a_FD4(+1)*w(ix,iy,iz+1) + a_FD4(+2)*w(ix,iy,iz+2))*dz_inv

                        div(ix,iy,iz) = u_dx + v_dy + w_dz
                    end do
                end do
            end do
        case("FD_4th_central_optimized")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
                        u_dx = (a(-3)*u(ix-3,iy,iz) + a(-2)*u(ix-2,iy,iz) &
                        + a(-1)*u(ix-1,iy,iz) + a(+1)*u(ix+1,iy,iz) + a(+2)*u(ix+2,iy,iz) &
                        + a(+3)*u(ix+3,iy,iz))*dx_inv

                        v_dy = (a(-3)*v(ix,iy-3,iz) + a(-2)*v(ix,iy-2,iz) &
                        + a(-1)*v(ix,iy-1,iz) + a(+1)*v(ix,iy+1,iz) + a(+2)*v(ix,iy+2,iz) &
                        + a(+3)*v(ix,iy+3,iz))*dy_inv

                        w_dz = (a(-3)*w(ix,iy,iz-3) + a(-2)*w(ix,iy,iz-2) &
                        + a(-1)*w(ix,iy,iz-1) + a(+1)*w(ix,iy,iz+1) + a(+2)*w(ix,iy,iz+2) &
                        + a(+3)*w(ix,iy,iz+3))*dz_inv

                        div(ix,iy,iz) = u_dx + v_dy + w_dz
                    end do
                end do
            end do
        case default
            call abort(0112181, "ERROR: discretization method in discretization is unknown")
        end select
    else
        !-----------------------------------------------------------------------
        ! 2D case 2D case 2D case 2D case
        !-----------------------------------------------------------------------
        select case(discretization)
        case("FD_2nd_central")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dx = (u(ix+1,iy,1)-u(ix-1,iy,1))*dx_inv*0.5_rk
                    v_dy = (v(ix,iy+1,1)-v(ix,iy-1,1))*dy_inv*0.5_rk
                    div(ix,iy,1) = u_dx + v_dy
                end do
            end do
        case("FD_4th_central")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dx = (a_FD4(-2)*u(ix-2,iy,1) +a_FD4(-1)*u(ix-1,iy,1) +a_FD4(+1)*u(ix+1,iy,1) +a_FD4(+2)*u(ix+2,iy,1))*dx_inv
                    v_dy = (a_FD4(-2)*v(ix,iy-2,1) +a_FD4(-1)*v(ix,iy-1,1) +a_FD4(+1)*v(ix,iy+1,1) +a_FD4(+2)*v(ix,iy+2,1))*dy_inv
                    div(ix,iy,1) = u_dx + v_dy
                end do
            end do
        case("FD_4th_central_optimized")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dx = (a(-3)*u(ix-3,iy,1) + a(-2)*u(ix-2,iy,1) + &
                    a(-1)*u(ix-1,iy,1) + a(0)*u(ix,iy,1)&
                    +  a(+1)*u(ix+1,iy,1) + a(+2)*u(ix+2,iy,1) + a(+3)*u(ix+3,iy,1))*dx_inv

                    v_dy = (a(-3)*v(ix,iy-3,1) + a(-2)*v(ix,iy-2,1) + &
                    a(-1)*v(ix,iy-1,1) + a(0)*v(ix,iy,1)&
                    +  a(+1)*v(ix,iy+1,1) + a(+2)*v(ix,iy+2,1) + a(+3)*v(ix,iy+3,1))*dy_inv

                    div(ix,iy,1) = u_dx + v_dy
                end do
            end do
        case default
            call abort(0112182, "ERROR: discretization method in discretization is unknown")
        end select
    end if

end subroutine divergence
