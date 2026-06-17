!> \brief Subtract the gradient of a scalar field from a vector field
!> \details
!! This routine computes the gradient of a scalar field (e.g., pressure, potential)
!! and subtracts it from a vector field (e.g., velocity). This operation is common in:
!!   - Helmholtz projection: u := u - grad(phi)
!!   - Pressure correction: u := u - grad(p)
!!
!! The operation is performed on the interior domain only (ghost nodes are not updated).
!-----------------------------------------------------------------------------------------------------
subroutine subtract_scalar_gradient(u, dx, Bs, g, discretization, p)

    implicit none
    
    real(kind=rk), dimension(:,:,:,:), intent(inout) :: u                       !> vector field to be corrected (e.g., velocity)
    real(kind=rk), dimension(3), intent(in)          :: dx                      !> spacing of the block
    integer(kind=ik), dimension(3), intent(in)       :: Bs
    integer(kind=ik), intent(in)                     :: g                       !> grid parameters
    character(len=*), intent(in)                     :: discretization
    real(kind=rk), dimension(:,:,:), intent(in)      :: p                       !> scalar field (e.g., pressure, potential)
    real(kind=rk)                                    :: p_dx, p_dy, p_dz        !> derivatives
    real(kind=rk)                                    :: dx_inv, dy_inv, dz_inv  !> inverse of dx, dy, dz
    integer(kind=ik)                                 :: ix, iy, iz              ! loop variables

    !> parameters for FD1_r operator
    real(kind=rk), allocatable, dimension(:) :: FD1_r
    integer(kind=ik) :: FD1_rs, FD1_re
    
    ! Setup finite difference stencils
    call setup_FD1_right_stencil(discretization, FD1_r, FD1_rs, FD1_re)
    
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
                    p_dx = sum(FD1_r(FD1_rs:FD1_re) * p(ix+FD1_rs:ix+FD1_re,iy,iz)) * dx_inv
                    p_dy = sum(FD1_r(FD1_rs:FD1_re) * p(ix,iy+FD1_rs:iy+FD1_re,iz)) * dy_inv
                    p_dz = sum(FD1_r(FD1_rs:FD1_re) * p(ix,iy,iz+FD1_rs:iz+FD1_re)) * dz_inv
                    
                    u(ix,iy,iz,1) = u(ix,iy,iz,1) - p_dx
                    u(ix,iy,iz,2) = u(ix,iy,iz,2) - p_dy
                    u(ix,iy,iz,3) = u(ix,iy,iz,3) - p_dz
                end do
            end do
        end do
    else
        !-----------------------------------------------------------------------
        ! 2D case 2D case 2D case 2D case
        !-----------------------------------------------------------------------
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                p_dx = sum(FD1_r(FD1_rs:FD1_re) * p(ix+FD1_rs:ix+FD1_re,iy,1)) * dx_inv
                p_dy = sum(FD1_r(FD1_rs:FD1_re) * p(ix,iy+FD1_rs:iy+FD1_re,1)) * dy_inv
                
                u(ix,iy,1,1) = u(ix,iy,1,1) - p_dx
                u(ix,iy,1,2) = u(ix,iy,1,2) - p_dy
            end do
        end do
    end if

end subroutine subtract_scalar_gradient
