subroutine compute_Qcriterion(u, dx, Bs, g, discretization, Qcrit)

    implicit none
    real(kind=rk), dimension(3), intent(in)        :: dx              !> spacing of the block
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: u               !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: Qcrit           !> Q-criterion
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g               !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: uxdx,uxdy,uxdz,uydx,uydy,uydz,uzdx,uzdy,uzdz
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz      ! loop variables
    real(kind=rk) :: Amatrix(1:3,1:3)

    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)

    Qcrit = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then
        !-----------------------------------------------------------------------
        ! 3D case 3D case 3D case 3D case
        !-----------------------------------------------------------------------
        dz_inv = 1.0_rk / dx(3)
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                do iz = g+1, Bs(3)+g
                    uxdx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 1)) * dx_inv
                    uxdy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 1)) * dy_inv
                    uxdz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 1)) * dz_inv

                    uydx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 2)) * dx_inv
                    uydy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 2)) * dy_inv
                    uydz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 2)) * dz_inv

                    uzdx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 3)) * dx_inv
                    uzdy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 3)) * dy_inv
                    uzdz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 3)) * dz_inv

                    Amatrix(1,:) =  (/ uxdx, uxdy, uxdz/)
                    Amatrix(2,:) =  (/ uydx, uydy, uydz/)
                    Amatrix(3,:) =  (/ uzdx, uzdy, uzdz/)

                    Qcrit(ix,iy,iz) = -0.5d0*( sum( Amatrix*transpose(Amatrix) ) )
                end do
            end do
        end do
    else
        call abort(0112181, "ERROR: Q-criterion for 2D not implemented.")
    end if
end subroutine compute_Qcriterion
