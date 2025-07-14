!> \brief computation of the divergence from a given velocity field
!-----------------------------------------------------------------------------------------------------
subroutine compute_divergence(u, dx, Bs, g, discretization, div)
    implicit none
    real(kind=rk), dimension(3), intent(in)        :: dx                        !> origin and spacing of the block
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: u                         !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: div                       !> vorticity
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    real(kind=rk)                                  :: u_dx, v_dy, w_dz          !> derivatives
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables
    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)

    div = 0.0_rk
    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)


    if (size(u,3)>2) then
        !-----------------------------------------------------------------------
        ! 3D case 3D case 3D case 3D case
        !-----------------------------------------------------------------------
        dz_inv = 1.0_rk / dx(3)
        
        ! 3D: div u = u_dx + v_dy + w_dz
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    u_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz,1)) * dx_inv
                    v_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz,2)) * dy_inv
                    w_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le,3)) * dz_inv

                    div(ix,iy,iz) = u_dx + v_dy + w_dz
                end do
            end do
        end do
    else
        !-----------------------------------------------------------------------
        ! 2D case 2D case 2D case 2D case
        !-----------------------------------------------------------------------
        ! 2D: div u = u_dx + v_dy
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                u_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, 1,1)) * dx_inv
                v_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, 1,2)) * dy_inv
                div(ix,iy,1) = u_dx + v_dy
            end do
        end do
    end if

end subroutine compute_divergence