!> \brief compute derivative for time step t (for saving it on disk)
! ************************************************************************************
subroutine compute_derivative(u, dx, Bs, g, derivative_dim, derivative_order, discretization, u_out)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:), intent(in)    :: u                         !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: u_out
    integer(kind=ik), intent(in)                   :: derivative_dim            !> Dimension to apply the derivative - 1:x, 2:y, 3:z
    integer(kind=ik), intent(in)                   :: derivative_order          !> Order to apply the derivative - 1:first order (dx), 2:second order (dxdx)
    character(len=*), intent(in)                   :: discretization            !> Convergence order of stencil to be used
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz, dim, gs(3), s
    !> parameters for FD1L, FD2, FD1X operators
    real(kind=rk), allocatable, dimension(:) :: FD1_l, FD2, FD1X
    integer(kind=ik) :: FD1_ls, FD1_le, FD2_s, FD2_e, FD1X_s, FD1X_e

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)
    call setup_FD2_stencil(discretization, FD2, FD2_s, FD2_e)
    call setup_FD1X_stencil(discretization, FD1X, FD1X_s, FD1X_e)

    dim = 2
    if (size(u,3)> 1) dim = 3
    ! define g as 3D vector, so that we can merge 2D and 3D loops
    gs = 0
    gs(1:dim) = g

    u_out(:,:,:) = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    if (dim == 3) then
        dz_inv = 1.0_rk / dx(3)
        dz2_inv = 1.0_rk / (dx(3)**2)
    endif

    ! Unified derivative computation using module variables
    select case(derivative_order)
    case(1) ! first order derivative
        select case(derivative_dim)
        case(1) ! derivative in x-direction
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        u_out(ix,iy,iz) = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz)) * dx_inv
                    end do
                end do
            end do
        case(2) ! derivative in y-direction
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        u_out(ix,iy,iz) = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz)) * dy_inv
                    end do
                end do
            end do
        case(3) ! derivative in z-direction
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        u_out(ix,iy,iz) = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le)) * dz_inv
                    end do
                end do
            end do
        end select
    case(2) ! second order derivative
        select case(derivative_dim)
        case(1) ! derivative in x-direction
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        u_out(ix,iy,iz) = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e, iy, iz)) * dx2_inv
                    end do
                end do
            end do
        case(2) ! derivative in y-direction
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        u_out(ix,iy,iz) = sum(FD2(FD2_s:FD2_e) * u(ix, iy+FD2_s:iy+FD2_e, iz)) * dy2_inv
                    end do
                end do
            end do
        case(3) ! derivative in z-direction
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        u_out(ix,iy,iz) = sum(FD2(FD2_s:FD2_e) * u(ix, iy, iz+FD2_s:iz+FD2_e)) * dz2_inv
                    end do
                end do
            end do
        end select
    case(11) ! cross derivatives (using FD1X stencils in cross-shaped pattern)
        select case(derivative_dim)
        case(12) ! derivative in xy-direction (mixed partial d²/dxdy)
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        ! Apply cross-shaped stencil for mixed partial: d²u/dxdy
                        ! FD1X stencil is applied twice: once for each diagonal
                        ! First diagonal (i,j): from (-s,-s) to (s,s) with positive sign
                        ! Second diagonal (i,-j): from (-s,s) to (s,-s) with negative sign
                        u_out(ix,iy,iz) = 0.0_rk
                        do s = FD1X_s, FD1X_e
                            ! First diagonal: positive contribution
                            u_out(ix,iy,iz) = u_out(ix,iy,iz) + FD1X(s) * u(ix+s, iy+s, iz)
                            ! Second diagonal: negative contribution
                            u_out(ix,iy,iz) = u_out(ix,iy,iz) - FD1X(s) * u(ix+s, iy-s, iz)
                        end do
                        u_out(ix,iy,iz) = u_out(ix,iy,iz) * dx_inv * dy_inv
                    end do
                end do
            end do
        case(13) ! derivative in xz-direction (mixed partial d²/dxdz)
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        ! Apply cross-shaped stencil for mixed partial: d²u/dxdz
                        ! FD1X stencil is applied twice: once for each diagonal
                        ! First diagonal (i,k): from (-s,-s) to (s,s) with positive sign
                        ! Second diagonal (i,-k): from (-s,s) to (s,-s) with negative sign
                        u_out(ix,iy,iz) = 0.0_rk
                        do s = FD1X_s, FD1X_e
                            ! First diagonal: positive contribution
                            u_out(ix,iy,iz) = u_out(ix,iy,iz) + FD1X(s) * u(ix+s, iy, iz+s)
                            ! Second diagonal: negative contribution
                            u_out(ix,iy,iz) = u_out(ix,iy,iz) - FD1X(s) * u(ix+s, iy, iz-s)
                        end do
                        u_out(ix,iy,iz) = u_out(ix,iy,iz) * dx_inv * dz_inv
                    end do
                end do
            end do
        case(23) ! derivative in yz-direction (mixed partial d²/dydz)
            do iz = gs(3)+1, Bs(3)+gs(3)
                do iy = gs(2)+1, Bs(2)+gs(2)
                    do ix = gs(1)+1, Bs(1)+gs(1)
                        ! Apply cross-shaped stencil for mixed partial: d²u/dydz
                        ! FD1X stencil is applied twice: once for each diagonal
                        ! First diagonal (j,k): from (-s,-s) to (s,s) with positive sign
                        ! Second diagonal (j,-k): from (-s,s) to (s,-s) with negative sign
                        u_out(ix,iy,iz) = 0.0_rk
                        do s = FD1X_s, FD1X_e
                            ! First diagonal: positive contribution
                            u_out(ix,iy,iz) = u_out(ix,iy,iz) + FD1X(s) * u(ix, iy+s, iz+s)
                            ! Second diagonal: negative contribution
                            u_out(ix,iy,iz) = u_out(ix,iy,iz) - FD1X(s) * u(ix, iy+s, iz-s)
                        end do
                        u_out(ix,iy,iz) = u_out(ix,iy,iz) * dy_inv * dz_inv
                    end do
                end do
            end do
        end select
    case default
        write(*,*) "ERROR: derivative_order", derivative_order, "is not supported"
        call abort(240321, "ERROR: derivative_order in compute_derivative is unknown")
    end select

end subroutine compute_derivative