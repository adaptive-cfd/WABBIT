subroutine compute_Qcriterion(u, v, w, dx, Bs, g, discretization, Qcrit)

    implicit none
    real(kind=rk), dimension(3), intent(in)        :: dx              !> spacing of the block
    real(kind=rk), dimension(:,:,:), intent(in)    :: u, v, w         !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: Qcrit           !> vorticity
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g               !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: uxdx,uxdy,uxdz,uydx,uydy,uydz,uzdx,uzdy,uzdz
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz      ! loop variables
    real(kind=rk)                                  :: a(-3:3)         ! coefficients for Tam&Webb
    real(kind=rk) :: Amatrix(1:3,1:3)

    Qcrit = 0.0_rk

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
        0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then
        !-----------------------------------------------------------------------
        ! 3D case 3D case 3D case 3D case
        !-----------------------------------------------------------------------
        dz_inv = 1.0_rk / dx(3)
        if (discretization == "FD_2nd_central" ) then
            call abort(0112181, "ERROR: Q-criterion using 2nd order not implemented.")

        else if (discretization == "FD_4th_central_optimized") then
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
                        uxdx = (a(-3)*u(ix-3,iy,iz) + a(-2)*u(ix-2,iy,iz) + a(-1)*u(ix-1,iy,iz) + a(0)*u(ix,iy,iz) &
                             +  a(+1)*u(ix+1,iy,iz) + a(+2)*u(ix+2,iy,iz) + a(+3)*u(ix+3,iy,iz))*dx_inv
                        uxdy = (a(-3)*u(ix,iy-3,iz) + a(-2)*u(ix,iy-2,iz) + a(-1)*u(ix,iy-1,iz) + a(0)*u(ix,iy,iz) &
                             +  a(+1)*u(ix,iy+1,iz) + a(+2)*u(ix,iy+2,iz) + a(+3)*u(ix,iy+3,iz))*dy_inv
                        uxdz = (a(-3)*u(ix,iy,iz-3) + a(-2)*u(ix,iy,iz-2) + a(-1)*u(ix,iy,iz-1) + a(0)*u(ix,iy,iz) &
                             +  a(+1)*u(ix,iy,iz+1) + a(+2)*u(ix,iy,iz+2) + a(+3)*u(ix,iy,iz+3))*dz_inv

                        uydx = (a(-3)*v(ix-3,iy,iz) + a(-2)*v(ix-2,iy,iz) + a(-1)*v(ix-1,iy,iz) + a(0)*v(ix,iy,iz) &
                             +  a(+1)*v(ix+1,iy,iz) + a(+2)*v(ix+2,iy,iz) + a(+3)*v(ix+3,iy,iz))*dx_inv
                        uydy = (a(-3)*v(ix,iy-3,iz) + a(-2)*v(ix,iy-2,iz) + a(-1)*v(ix,iy-1,iz) + a(0)*v(ix,iy,iz) &
                             +  a(+1)*v(ix,iy+1,iz) + a(+2)*v(ix,iy+2,iz) + a(+3)*v(ix,iy+3,iz))*dy_inv
                        uydz = (a(-3)*v(ix,iy,iz-3) + a(-2)*v(ix,iy,iz-2) + a(-1)*v(ix,iy,iz-1) + a(0)*v(ix,iy,iz) &
                             +  a(+1)*v(ix,iy,iz+1) + a(+2)*v(ix,iy,iz+2) + a(+3)*v(ix,iy,iz+3))*dz_inv

                        uzdx = (a(-3)*w(ix-3,iy,iz) + a(-2)*w(ix-2,iy,iz) + a(-1)*w(ix-1,iy,iz) + a(0)*w(ix,iy,iz) &
                             +  a(+1)*w(ix+1,iy,iz) + a(+2)*w(ix+2,iy,iz) + a(+3)*w(ix+3,iy,iz))*dx_inv
                        uzdy = (a(-3)*w(ix,iy-3,iz) + a(-2)*w(ix,iy-2,iz) + a(-1)*w(ix,iy-1,iz) + a(0)*w(ix,iy,iz) &
                             +  a(+1)*w(ix,iy+1,iz) + a(+2)*w(ix,iy+2,iz) + a(+3)*w(ix,iy+3,iz))*dy_inv
                        uzdz = (a(-3)*w(ix,iy,iz-3) + a(-2)*w(ix,iy,iz-2) + a(-1)*w(ix,iy,iz-1) + a(0)*w(ix,iy,iz) &
                             +  a(+1)*w(ix,iy,iz+1) + a(+2)*w(ix,iy,iz+2) + a(+3)*w(ix,iy,iz+3))*dz_inv

                        Amatrix(1,:) =  (/ uxdx, uxdy, uxdz/)
                        Amatrix(2,:) =  (/ uydx, uydy, uydz/)
                        Amatrix(3,:) =  (/ uzdx, uzdy, uzdz/)

                        Qcrit(ix,iy,iz) = -0.5d0*( sum( Amatrix*transpose(Amatrix) ) )

                    end do
                end do
            end do
        else
            call abort(0112181, "ERROR: discretization method in discretization is unknown")
        end if
    else
        call abort(0112181, "ERROR: Q-criterion for 2D not implemented.")
    end if
end subroutine compute_Qcriterion
