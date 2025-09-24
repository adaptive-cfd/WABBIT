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

end subroutine compute_helicity_abs
