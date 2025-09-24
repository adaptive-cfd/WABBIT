!> \brief computation of the gradient of a field
!-----------------------------------------------------------------------------------------------------
subroutine compute_gradient(fld, dx, Bs, g, discretization, grad_fld)

    implicit none
    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:), intent(in)    :: fld                       !> local datafields
    real(kind=rk), dimension(:,:,:,:), intent(out) :: grad_fld                  !> gradient
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    real(kind=rk)                                  :: fld_dx, fld_dy, fld_dz    !> derivatives
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables

    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(fld,3)>2) then
        !-----------------------------------------------------------------------
        ! 3D case 3D case 3D case 3D case
        !-----------------------------------------------------------------------
        dz_inv = 1.0_rk / dx(3)
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                do iz = g+1, Bs(3)+g
                    fld_dx = sum(FD1_l(FD1_ls:FD1_le) * fld(ix+FD1_ls:ix+FD1_le, iy, iz)) * dx_inv
                    fld_dy = sum(FD1_l(FD1_ls:FD1_le) * fld(ix, iy+FD1_ls:iy+FD1_le, iz)) * dy_inv
                    fld_dz = sum(FD1_l(FD1_ls:FD1_le) * fld(ix, iy, iz+FD1_ls:iz+FD1_le)) * dz_inv

                    grad_fld(ix,iy,iz,1) = fld_dx
                    grad_fld(ix,iy,iz,2) = fld_dy
                    grad_fld(ix,iy,iz,3) = fld_dz
                end do
            end do
        end do
    else
        !-----------------------------------------------------------------------
        ! 2D case 2D case 2D case 2D case
        !-----------------------------------------------------------------------
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                fld_dx = sum(FD1_l(FD1_ls:FD1_le) * fld(ix+FD1_ls:ix+FD1_le, iy, 1)) * dx_inv
                fld_dy = sum(FD1_l(FD1_ls:FD1_le) * fld(ix, iy+FD1_ls:iy+FD1_le, 1)) * dy_inv
                grad_fld(ix,iy,1,1) = fld_dx
                grad_fld(ix,iy,1,2) = fld_dy
            end do
        end do
    end if

end subroutine compute_gradient
