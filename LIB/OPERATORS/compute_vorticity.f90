!> \brief compute vorticity for time step t (for saving it on disk)
! ********************************************************************************************
subroutine compute_vorticity(u, dx, Bs, g, discretization, vorticity)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: u                         !> local datafields
    real(kind=rk), dimension(:,:,:,:), intent(out) :: vorticity
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: u_dy, u_dz, v_dx, v_dz, w_dx, w_dy
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables

    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)


    vorticity = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then ! 3D case
        dz_inv = 1.0_rk / dx(3)

        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    u_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 1)) * dy_inv
                    u_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 1)) * dz_inv
                    v_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 2)) * dx_inv
                    v_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 2)) * dz_inv
                    w_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 3)) * dx_inv
                    w_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 3)) * dy_inv

                    vorticity(ix,iy,iz,1) = w_dy - v_dz
                    vorticity(ix,iy,iz,2) = u_dz - w_dx
                    vorticity(ix,iy,iz,3) = v_dx - u_dy
                end do
            end do
        end do

    else
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                u_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, 1, 1)) * dy_inv
                v_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, 1, 2)) * dx_inv

                vorticity(ix,iy,1,1) = v_dx - u_dy
            end do
        end do
    end if

    ! deallocate the stencils if they were allocated
    if (allocated(FD1_l)) deallocate(FD1_l)

end subroutine compute_vorticity


subroutine compute_vorticity_abs(u, dx, Bs, g, discretization, vor_abs)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    implicit none

    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: u                         !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: vor_abs                   !> vorticity
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: u_dy, u_dz, v_dx, v_dz, w_dx, w_dy
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables

    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)


    vor_abs = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then ! 3D case
        dz_inv = 1.0_rk / dx(3)
        
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                do iz = g+1, Bs(3)+g
                    u_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 1)) * dy_inv
                    u_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 1)) * dz_inv
                    v_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 2)) * dx_inv
                    v_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 2)) * dz_inv
                    w_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 3)) * dx_inv
                    w_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 3)) * dy_inv

                    vor_abs(ix,iy,iz) = sqrt( (w_dy - v_dz)**2 + (u_dz - w_dx)**2 + (v_dx - u_dy)**2)
                end do
            end do
        end do
    else
        call abort(23321, "ERROR: vor-abs makes not much sense for 2D data.")
    end if

    ! deallocate the stencils if they were allocated
    if (allocated(FD1_l)) deallocate(FD1_l)

end subroutine compute_vorticity_abs


subroutine compute_helicity(u, dx, Bs, g, discretization, helicity)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    implicit none

    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: u                         !> local datafields
    real(kind=rk), dimension(:,:,:,:), intent(out) :: helicity                 !> helicity vector field
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: u_dy, u_dz, v_dx, v_dz, w_dx, w_dy
    real(kind=rk)                                  :: vort_x, vort_y, vort_z    !> vorticity components
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables

    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)

    helicity = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then ! 3D case
        dz_inv = 1.0_rk / dx(3)
        
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                do iz = g+1, Bs(3)+g
                    ! Compute vorticity components
                    u_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 1)) * dy_inv
                    u_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 1)) * dz_inv
                    v_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 2)) * dx_inv
                    v_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 2)) * dz_inv
                    w_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 3)) * dx_inv
                    w_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 3)) * dy_inv

                    vort_x = w_dy - v_dz
                    vort_y = u_dz - w_dx
                    vort_z = v_dx - u_dy

                    ! Compute helicity vector: H = u × ω
                    helicity(ix,iy,iz,1) = u(ix,iy,iz,1) * vort_x
                    helicity(ix,iy,iz,2) = u(ix,iy,iz,2) * vort_y
                    helicity(ix,iy,iz,3) = u(ix,iy,iz,3) * vort_z
                end do
            end do
        end do
    else
        call abort(23322, "ERROR: helicity makes not much sense for 2D data.")
    end if

    ! deallocate the stencils if they were allocated
    if (allocated(FD1_l)) deallocate(FD1_l)

end subroutine compute_helicity


subroutine compute_helicity_abs(u, dx, Bs, g, discretization, hel_abs)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    implicit none

    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: u                         !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: hel_abs                  !> helicity magnitude
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: u_dy, u_dz, v_dx, v_dz, w_dx, w_dy
    real(kind=rk)                                  :: vort_x, vort_y, vort_z    !> vorticity components
    real(kind=rk)                                  :: helicity_val              !> helicity value
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables

    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)

    hel_abs = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then ! 3D case
        dz_inv = 1.0_rk / dx(3)
        
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                do iz = g+1, Bs(3)+g
                    ! Compute vorticity components
                    u_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 1)) * dy_inv
                    u_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 1)) * dz_inv
                    v_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 2)) * dx_inv
                    v_dz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 2)) * dz_inv
                    w_dx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 3)) * dx_inv
                    w_dy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 3)) * dy_inv

                    vort_x = w_dy - v_dz
                    vort_y = u_dz - w_dx
                    vort_z = v_dx - u_dy

                    ! Compute helicity: H = u · ω
                    helicity_val = u(ix,iy,iz,1) * vort_x + u(ix,iy,iz,2) * vort_y + u(ix,iy,iz,3) * vort_z
                    
                    ! Compute helicity magnitude
                    hel_abs(ix,iy,iz) = abs(helicity_val)
                end do
            end do
        end do
    else
        call abort(23323, "ERROR: helicity-abs makes not much sense for 2D data.")
    end if

    ! deallocate the stencils if they were allocated
    if (allocated(FD1_l)) deallocate(FD1_l)

end subroutine compute_helicity_abs


!> \brief Compute vorticity stretching: alpha = omega_hat_i * e_ij * omega_hat_j
!! where omega_hat = omega / |omega| is normalized vorticity
!! and e_ij = 0.5*(u_i,j + u_j,i) is the symmetric strain rate tensor
!! This quantity measures the alignment of vorticity with the principal stretching direction
subroutine compute_vorticity_stretching(u, dx, Bs, g, discretization, vor_stretch)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    implicit none

    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: u                         !> local datafields
    real(kind=rk), dimension(:,:,:), intent(out)   :: vor_stretch               !> vorticity stretching (scalar)
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    
    !> derivatives
    real(kind=rk)                                  :: uxdx, uxdy, uxdz, uydx, uydy, uydz, uzdx, uzdy, uzdz
    real(kind=rk)                                  :: u_dy, u_dz, v_dx, v_dz, w_dx, w_dy
    real(kind=rk)                                  :: vort_x, vort_y, vort_z    !> vorticity components
    real(kind=rk)                                  :: omega_mag, omega_mag_inv  !> vorticity magnitude and its inverse
    real(kind=rk)                                  :: omega_hat(3)              !> normalized vorticity
    real(kind=rk)                                  :: e_tensor(3,3)             !> strain rate tensor
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv    !> inverse of dx, dy, dz
    real(kind=rk), parameter                       :: eps = 1.0e-10_rk          !> regularization parameter for |omega| -> 0
    integer(kind=ik)                               :: ix, iy, iz, i, j          !> loop variables

    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l
    integer(kind=ik) :: FD1_ls, FD1_le

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(discretization, FD1_l, FD1_ls, FD1_le)

    vor_stretch = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then ! 3D case only
        dz_inv = 1.0_rk / dx(3)
        
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                do iz = g+1, Bs(3)+g
                    ! ===================================================================
                    ! Compute all 9 velocity gradient components
                    ! ===================================================================
                    uxdx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 1)) * dx_inv
                    uxdy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 1)) * dy_inv
                    uxdz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 1)) * dz_inv

                    uydx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 2)) * dx_inv
                    uydy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 2)) * dy_inv
                    uydz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 2)) * dz_inv

                    uzdx = sum(FD1_l(FD1_ls:FD1_le) * u(ix+FD1_ls:ix+FD1_le, iy, iz, 3)) * dx_inv
                    uzdy = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy+FD1_ls:iy+FD1_le, iz, 3)) * dy_inv
                    uzdz = sum(FD1_l(FD1_ls:FD1_le) * u(ix, iy, iz+FD1_ls:iz+FD1_le, 3)) * dz_inv

                    ! ===================================================================
                    ! Compute symmetric strain rate tensor: e_ij = 0.5 * (u_i,j + u_j,i)
                    ! ===================================================================
                    e_tensor(1,1) = uxdx                      ! e_xx = du/dx
                    e_tensor(2,2) = uydy                      ! e_yy = dv/dy
                    e_tensor(3,3) = uzdz                      ! e_zz = dw/dz
                    e_tensor(1,2) = 0.5_rk * (uxdy + uydx)    ! e_xy = 0.5*(du/dy + dv/dx)
                    e_tensor(1,3) = 0.5_rk * (uxdz + uzdx)    ! e_xz = 0.5*(du/dz + dw/dx)
                    e_tensor(2,3) = 0.5_rk * (uydz + uzdy)    ! e_yz = 0.5*(dv/dz + dw/dy)
                    e_tensor(2,1) = e_tensor(1,2)             ! symmetry
                    e_tensor(3,1) = e_tensor(1,3)             ! symmetry
                    e_tensor(3,2) = e_tensor(2,3)             ! symmetry

                    ! ===================================================================
                    ! Compute vorticity components
                    ! ===================================================================
                    vort_x = uzdy - uydz   ! omega_x = dw/dy - dv/dz
                    vort_y = uxdz - uzdx   ! omega_y = du/dz - dw/dx
                    vort_z = uydx - uxdy   ! omega_z = dv/dx - du/dy

                    ! ===================================================================
                    ! Compute vorticity magnitude and normalize
                    ! ===================================================================
                    omega_mag = sqrt(vort_x**2 + vort_y**2 + vort_z**2)
                    
                    ! Singularity handling: avoid division by zero when |omega| -> 0
                    ! Use regularization: omega_hat = omega / (|omega| + eps)
                    if (omega_mag > eps) then
                        omega_mag_inv = 1.0_rk / omega_mag
                        omega_hat(1) = vort_x * omega_mag_inv
                        omega_hat(2) = vort_y * omega_mag_inv
                        omega_hat(3) = vort_z * omega_mag_inv
                    else
                        ! In regions of negligible vorticity, set stretching to zero
                        omega_hat = 0.0_rk
                    end if

                    ! ===================================================================
                    ! Compute vorticity stretching: alpha = omega_hat_i * e_ij * omega_hat_j
                    ! ===================================================================
                    vor_stretch(ix,iy,iz) = 0.0_rk
                    do i = 1, 3
                        do j = 1, 3
                            vor_stretch(ix,iy,iz) = vor_stretch(ix,iy,iz) + omega_hat(i) * e_tensor(i,j) * omega_hat(j)
                        end do
                    end do
                end do
            end do
        end do
    else
        ! 2D case: vorticity stretching is not meaningful in 2D
        call abort(23324, "ERROR: vorticity stretching is only defined for 3D flows.")
    end if

    ! deallocate the stencils if they were allocated
    if (allocated(FD1_l)) deallocate(FD1_l)

end subroutine compute_vorticity_stretching
