
!> on a block compute the dissipation rate (local point-wise),
!> i.e. \varepsilon = u_j * \nabla^2 u_j (Einstein summation convention)
!> This computes the velocity-weighted Laplacian for each component
subroutine compute_dissipation(u, dx, Bs, g, discretization, dissipation_rate)
    implicit none

    real(kind=rk), intent(in)                  :: u(:,:,:,:)               !< datafields
    real(kind=rk), dimension(3), intent(in)    :: dx                       !< block spacing
    integer(kind=ik), dimension(3), intent(in) :: Bs                       !< block size
    integer(kind=ik), intent(in)               :: g                        !< ghost points
    character(len=*), intent(in)               :: discretization           !< discretization order for second order FD stencil
    real(kind=rk), intent(inout)               :: dissipation_rate(:,:,:)  !< local dissipation rate

    !> parameters for FD2 operator
    real(kind=rk), allocatable, dimension(:) :: FD2
    integer(kind=ik) :: FD2_s, FD2_e

    !> inverse of dx, dy, dz
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, u_dxdx, u_dydy, u_dzdz, &
    v_dxdx, v_dydy, v_dzdz, w_dxdx, w_dydy, w_dzdz
    ! loop variables
    integer(kind=ik) :: ix, iy, iz

    ! Setup stencils using the unified interface from module_operators
    call setup_FD2_stencil(discretization, FD2, FD2_s, FD2_e)

    ! Dissipation rate is u_j*laplace(u_j) (einstein summation convention)

    if ( size(u,3)==1) then
        dx2_inv = 1.0_rk / dx(1)**2
        dy2_inv = 1.0_rk / dx(2)**2

        iz = 1

        ! Use the unified stencil infrastructure for all discretization orders
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! Second derivatives using module stencils
                u_dxdx = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e,iy,iz,1)) * dx2_inv
                v_dxdx = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e,iy,iz,2)) * dx2_inv
                
                u_dydy = sum(FD2(FD2_s:FD2_e) * u(ix,iy+FD2_s:iy+FD2_e,iz,1)) * dy2_inv
                v_dydy = sum(FD2(FD2_s:FD2_e) * u(ix,iy+FD2_s:iy+FD2_e,iz,2)) * dy2_inv

                dissipation_rate(ix,iy,iz) = u(ix,iy,iz,1)*(u_dxdx+u_dydy) + u(ix,iy,iz,2)*(v_dxdx+v_dydy)
            end do
        end do

    else
        dx2_inv = 1.0_rk / dx(1)**2
        dy2_inv = 1.0_rk / dx(2)**2
        dz2_inv = 1.0_rk / dx(3)**2

        ! Use the unified stencil infrastructure for all discretization orders
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! Second derivatives using module stencils
                    u_dxdx = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e,iy,iz,1)) * dx2_inv
                    u_dydy = sum(FD2(FD2_s:FD2_e) * u(ix,iy+FD2_s:iy+FD2_e,iz,1)) * dy2_inv
                    u_dzdz = sum(FD2(FD2_s:FD2_e) * u(ix,iy,iz+FD2_s:iz+FD2_e,1)) * dz2_inv

                    v_dxdx = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e,iy,iz,2)) * dx2_inv
                    v_dydy = sum(FD2(FD2_s:FD2_e) * u(ix,iy+FD2_s:iy+FD2_e,iz,2)) * dy2_inv
                    v_dzdz = sum(FD2(FD2_s:FD2_e) * u(ix,iy,iz+FD2_s:iz+FD2_e,2)) * dz2_inv

                    w_dxdx = sum(FD2(FD2_s:FD2_e) * u(ix+FD2_s:ix+FD2_e,iy,iz,3)) * dx2_inv
                    w_dydy = sum(FD2(FD2_s:FD2_e) * u(ix,iy+FD2_s:iy+FD2_e,iz,3)) * dy2_inv
                    w_dzdz = sum(FD2(FD2_s:FD2_e) * u(ix,iy,iz+FD2_s:iz+FD2_e,3)) * dz2_inv

                    dissipation_rate(ix,iy,iz) = u(ix,iy,iz,1)*(u_dxdx+u_dydy+u_dzdz) + u(ix,iy,iz,2)*(v_dxdx+v_dydy+v_dzdz) + u(ix,iy,iz,3)*(w_dxdx+w_dydy+w_dzdz) 
                end do
            end do
        end do
    endif

end subroutine compute_dissipation
