!> \brief compute derivative for time step t (for saving it on disk)
! ********************************************************************************************
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
    integer(kind=ik)                               :: ix, iy, iz, dim, gs(3)
    ! coefficients for Tam&Webb (4th order 1st derivative)
    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! 4th order central FD scheme
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    real(kind=rk), parameter :: b_FD4(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! 6th order central FD scheme
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/) ! 1st derivative
    real(kind=rk), parameter :: b_FD6(-3:3) = (/ 1.0_rk/90.0_rk, -3.0_rk/20.0_rk, 3.0_rk/2.0_rk, -49.0_rk/18.0_rk, 3.0_rk/2.0_rk, -3.0_rk/20.0_rk, 1.0_rk/90.0_rk/) ! 2nd derivative

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

    select case(discretization)
    case("FD_2nd_central")
        select case(derivative_order)
        case(1) ! first order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (u(ix+1, iy, iz) - u(ix-1, iy, iz)) * dx_inv * 0.5_rk
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (u(ix, iy+1, iz) - u(ix, iy-1, iz)) * dy_inv * 0.5_rk
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (u(ix, iy, iz+1) - u(ix, iy, iz-1)) * dx_inv * 0.5_rk
                end do; end do; end do
            end select
        case(2) ! second order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (u(ix-1, iy  , iz  ) -2.0_rk*u(ix, iy, iz) + u(ix+1, iy  , iz  ))*dx2_inv
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (u(ix  , iy+1, iz  ) -2.0_rk*u(ix, iy, iz) + u(ix  , iy+1, iz  ))*dy2_inv
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (u(ix  , iy  , iz-1) -2.0_rk*u(ix, iy, iz) + u(ix  , iy  , iz+1))*dz2_inv
                end do; end do; end do
            end select
        case(11) ! cross derivatives
            select case(derivative_dim)
            case(12) ! derivative in xy-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = 1.0_rk/4.0_rk*(u(ix+1, iy+1, iz  ) - u(ix-1, iy+1, iz  ) - u(ix+1, iy-1, iz  ) + u(ix-1, iy-1, iz  ))*dx_inv*dy_inv
                end do; end do; end do
            case(13) ! derivative in xz-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = 1.0_rk/4.0_rk*(u(ix+1, iy  , iz+1) - u(ix-1, iy  , iz+1) - u(ix+1, iy  , iz-1) + u(ix-1, iy  , iz-1))*dx_inv*dz_inv
                end do; end do; end do
            case(23) ! derivative in yz-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = 1.0_rk/4.0_rk*(u(ix  , iy+1, iz+1) - u(ix  , iy-1, iz+1) - u(ix  , iy+1, iz-1) + u(ix  , iy-1, iz-1))*dy_inv*dz_inv
                end do; end do; end do
            end select
        end select
    case("FD_4th_central")
        select case(derivative_order)
        case(1) ! first order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_FD4(-2)*u(ix-2,iy,iz) +a_FD4(-1)*u(ix-1,iy,iz) +a_FD4(+1)*u(ix+1,iy,iz) +a_FD4(+2)*u(ix+2,iy,iz))*dx_inv
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_FD4(-2)*u(ix,iy-2,iz) +a_FD4(-1)*u(ix,iy-1,iz) +a_FD4(+1)*u(ix,iy+1,iz) +a_FD4(+2)*u(ix,iy+2,iz))*dy_inv
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_FD4(-2)*u(ix,iy,iz-2) +a_FD4(-1)*u(ix,iy,iz-1) +a_FD4(+1)*u(ix,iy,iz+1) +a_FD4(+2)*u(ix,iy,iz+2))*dz_inv
                end do; end do; end do
            end select
        case(2) ! second order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz)   = (b_FD4(-2)*u(ix-2,iy,iz) &
                            + b_FD4(-1)*u(ix-1,iy,iz) &
                            + b_FD4(0)*u(ix,iy,iz) &
                            + b_FD4(+1)*u(ix+1,iy,iz) &
                            + b_FD4(+2)*u(ix+2,iy,iz))*dx2_inv
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz)   = (b_FD4(-2)*u(ix,iy-2,iz) &
                            + b_FD4(-1)*u(ix,iy-1,iz) &
                            + b_FD4(0)*u(ix,iy,iz) &
                            + b_FD4(+1)*u(ix,iy+1,iz) &
                            + b_FD4(+2)*u(ix,iy+2,iz))*dy2_inv
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz)   = (b_FD4(-2)*u(ix,iy,iz-2) &
                            + b_FD4(-1)*u(ix,iy,iz-1) &
                            + b_FD4(0)*u(ix,iy,iz) &
                            + b_FD4(+1)*u(ix,iy,iz+1) &
                            + b_FD4(+2)*u(ix,iy,iz+2))*dz2_inv
                end do; end do; end do
            end select
        case(11) ! cross derivatives
            select case(derivative_dim)
            case(12) ! derivative in xy-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (1.0_rk/3.0_rk*(u(ix+1, iy+1, iz  ) - u(ix-1, iy+1, iz  ) - u(ix+1, iy-1, iz  ) + u(ix-1, iy-1, iz  )) + \
                                     -1.0_rk/48.0_rk*(u(ix+2, iy+2, iz  ) - u(ix-2, iy+2, iz  ) - u(ix+2, iy-2, iz  ) + u(ix-2, iy-2, iz  )))*dx_inv*dy_inv
                end do; end do; end do
            case(13) ! derivative in xz-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (1.0_rk/3.0_rk*(u(ix+1, iy  , iz+1) - u(ix-1, iy  , iz+1) - u(ix+1, iy  , iz-1) + u(ix-1, iy  , iz-1)) + \
                                     -1.0_rk/48.0_rk*(u(ix+2, iy  , iz+2) - u(ix-2, iy  , iz+2) - u(ix+2, iy  , iz-2) + u(ix-2, iy  , iz-2)))*dx_inv*dz_inv
                end do; end do; end do
            case(23) ! derivative in yz-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (1.0_rk/3.0_rk*(u(ix  , iy+1, iz+1) - u(ix  , iy-1, iz+1) - u(ix  , iy+1, iz-1) + u(ix  , iy-1, iz-1)) + \
                                     -1.0_rk/48.0_rk*(u(ix  , iy+2, iz+2) - u(ix  , iy-2, iz+2) - u(ix  , iy+2, iz-2) + u(ix  , iy-2, iz-2)))*dy_inv*dz_inv
                end do; end do; end do
            end select
        end select
    case("FD_6th_central")
        select case(derivative_order)
        case(1) ! first order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_FD6(-3)*u(ix-3,iy,iz) &
                           +a_FD6(-2)*u(ix-2,iy,iz) &
                           +a_FD6(-1)*u(ix-1,iy,iz) &
                           +a_FD6(+1)*u(ix+1,iy,iz) &
                           +a_FD6(+2)*u(ix+2,iy,iz) &
                           +a_FD6(+3)*u(ix+3,iy,iz))*dx_inv
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_FD6(-3)*u(ix,iy-3,iz) &
                           +a_FD6(-2)*u(ix,iy-2,iz) &
                           +a_FD6(-1)*u(ix,iy-1,iz) &
                           +a_FD6(+1)*u(ix,iy+1,iz) &
                           +a_FD6(+2)*u(ix,iy+2,iz) &
                           +a_FD6(+3)*u(ix,iy+3,iz))*dy_inv
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_FD6(-3)*u(ix,iy,iz-3) &
                           +a_FD6(-2)*u(ix,iy,iz-2) &
                           +a_FD6(-1)*u(ix,iy,iz-1) &
                           +a_FD6(+1)*u(ix,iy,iz+1) &
                           +a_FD6(+2)*u(ix,iy,iz+2) &
                           +a_FD6(+3)*u(ix,iy,iz+3))*dz_inv
                end do; end do; end do
            end select
        case(2) ! second order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz)  = (b_FD6(-3)*u(ix-3,iy,iz) &
                            + b_FD6(-2)*u(ix-2,iy,iz) &
                            + b_FD6(-1)*u(ix-1,iy,iz) &
                            + b_FD6(0)*u(ix,iy,iz) &
                            + b_FD6(+1)*u(ix+1,iy,iz) &
                            + b_FD6(+2)*u(ix+2,iy,iz) &
                            + b_FD6(+3)*u(ix+3,iy,iz))*dx2_inv
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz)  = (b_FD6(-3)*u(ix,iy-3,iz) &
                            + b_FD6(-2)*u(ix,iy-2,iz) &
                            + b_FD6(-1)*u(ix,iy-1,iz) &
                            + b_FD6(0)*u(ix,iy,iz) &
                            + b_FD6(+1)*u(ix,iy+1,iz) &
                            + b_FD6(+2)*u(ix,iy+2,iz) &
                            + b_FD6(+3)*u(ix,iy+3,iz))*dy2_inv
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz)  = (b_FD6(-3)*u(ix,iy,iz-3) &
                            + b_FD6(-2)*u(ix,iy,iz-2) &
                            + b_FD6(-1)*u(ix,iy,iz-1) &
                            + b_FD6(0)*u(ix,iy,iz) &
                            + b_FD6(+1)*u(ix,iy,iz+1) &
                            + b_FD6(+2)*u(ix,iy,iz+2) &
                            + b_FD6(+3)*u(ix,iy,iz+3))*dz2_inv
                end do; end do; end do
            end select
        end select
    case("FD_4th_central_optimized")
        select case(derivative_order)
        case(1) ! first order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_TW4(-3)*u(ix-3,iy,iz) &
                           +a_TW4(-2)*u(ix-2,iy,iz) &
                           +a_TW4(-1)*u(ix-1,iy,iz) &
                           +a_TW4(+1)*u(ix+1,iy,iz) &
                           +a_TW4(+2)*u(ix+2,iy,iz) &
                           +a_TW4(+3)*u(ix+3,iy,iz))*dx_inv
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_TW4(-3)*u(ix,iy-3,iz) &
                           +a_TW4(-2)*u(ix,iy-2,iz) &
                           +a_TW4(-1)*u(ix,iy-1,iz) &
                           +a_TW4(+1)*u(ix,iy+1,iz) &
                           +a_TW4(+2)*u(ix,iy+2,iz) &
                           +a_TW4(+3)*u(ix,iy+3,iz))*dy_inv
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (a_TW4(-3)*u(ix,iy,iz-3) &
                           +a_TW4(-2)*u(ix,iy,iz-2) &
                           +a_TW4(-1)*u(ix,iy,iz-1) &
                           +a_TW4(+1)*u(ix,iy,iz+1) &
                           +a_TW4(+2)*u(ix,iy,iz+2) &
                           +a_TW4(+3)*u(ix,iy,iz+3))*dz_inv
                end do; end do; end do
            end select
        case(2) ! second order derivative
            select case(derivative_dim)
            case(1) ! derivative in x-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (b_FD4(-2)*u(ix-2,iy,iz) &
                           + b_FD4(-1)*u(ix-1,iy,iz) &
                           + b_FD4(0)*u(ix,iy,iz) &
                           + b_FD4(+1)*u(ix+1,iy,iz) &
                           + b_FD4(+2)*u(ix+2,iy,iz))*dx2_inv
                end do; end do; end do
            case(2) ! derivative in y-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (b_FD4(-2)*u(ix,iy-2,iz) &
                           + b_FD4(-1)*u(ix,iy-1,iz) &
                           + b_FD4(0)*u(ix,iy,iz) &
                           + b_FD4(+1)*u(ix,iy+1,iz) &
                           + b_FD4(+2)*u(ix,iy+2,iz))*dy2_inv
                end do; end do; end do
            case(3) ! derivative in z-direction
                do iz = gs(3)+1, Bs(3)+gs(3); do iy = gs(2)+1, Bs(2)+gs(2); do ix = gs(1)+1, Bs(1)+gs(1)
                    u_out(ix,iy,iz) = (b_FD4(-2)*u(ix,iy,iz-2) &
                           + b_FD4(-1)*u(ix,iy,iz-1) &
                           + b_FD4(0)*u(ix,iy,iz) &
                           + b_FD4(+1)*u(ix,iy,iz+1) &
                           + b_FD4(+2)*u(ix,iy,iz+2))*dz2_inv
                end do; end do; end do
            end select
        end select
    case default
        write(*,*) discretization
        call abort(240321, "ERROR: discretization method in discretization is unknown")
    end select

end subroutine compute_derivative