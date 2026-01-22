!> \brief computation of the laplacian from a given scalar field
!-----------------------------------------------------------------------------------------------------
subroutine compute_laplacian(u, dx, Bs, g, discretization, laplacian)
    implicit none
    real(kind=rk), dimension(3), intent(in)        :: dx                        !> origin and spacing of the block
    real(kind=rk), dimension(:,:,:), intent(in)    :: u                         !> local scalar field
    real(kind=rk), dimension(:,:,:), intent(out)   :: laplacian                 !> laplacian
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    real(kind=rk)                                  :: u_dxx, u_dyy, u_dzz       !> second derivatives
    real(kind=rk)                                  :: dx_inv2, dy_inv2, dz_inv2 !> inverse squared of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables
    !> parameters for FD2 operator
    real(kind=rk), allocatable, dimension(:) :: FD2
    integer(kind=ik) :: FD2_s, FD2_e

    ! Setup finite difference stencils
    call setup_FD2_stencil(discretization, FD2, FD2_s, FD2_e)

    laplacian = 0.0_rk
    dx_inv2 = 1.0_rk / (dx(1) * dx(1))
    dy_inv2 = 1.0_rk / (dx(2) * dx(2))

    if (size(u,3)>2) then
        !-----------------------------------------------------------------------
        ! 3D case 3D case 3D case 3D case
        !-----------------------------------------------------------------------
        dz_inv2 = 1.0_rk / (dx(3) * dx(3))
        
        ! 3D: laplacian u = u_dxx + u_dyy + u_dzz
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    u_dxx = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e, iy, iz)) * dx_inv2
                    u_dyy = sum(FD2(FD2_s:FD2_e) * u(ix, iy+FD2_s:iy+FD2_e, iz)) * dy_inv2
                    u_dzz = sum(FD2(FD2_s:FD2_e) * u(ix, iy, iz+FD2_s:iz+FD2_e)) * dz_inv2

                    laplacian(ix,iy,iz) = u_dxx + u_dyy + u_dzz
                end do
            end do
        end do
    else
        !-----------------------------------------------------------------------
        ! 2D case 2D case 2D case 2D case
        !-----------------------------------------------------------------------
        ! 2D: laplacian u = u_dxx + u_dyy
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                u_dxx = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e, iy, 1)) * dx_inv2
                u_dyy = sum(FD2(FD2_s:FD2_e) * u(ix, iy+FD2_s:iy+FD2_e, 1)) * dy_inv2
                laplacian(ix,iy,1) = u_dxx + u_dyy
            end do
        end do
    end if

end subroutine compute_laplacian
