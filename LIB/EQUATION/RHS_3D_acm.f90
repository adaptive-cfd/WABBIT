!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_acm.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief RHS for 3D artificial compressibility method
!
!>
!! input:    - datafield, grid parameter, derivative order \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 27/07/17 - create
! ********************************************************************************************

subroutine RHS_3D_acm(params, g, Bs, dx, x0, N_dF, phi, order_discretization, volume_int, time)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

    use module_operators

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> physics parameter structure
    type (type_params), intent(in)                         :: params

    !> grid parameter
    integer(kind=ik), intent(in)                           :: g, Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in)                :: x0, dx

    !> number of datafields
    integer(kind=ik), intent(in)                           :: N_dF
    !> datafields
    real(kind=rk), intent(inout)                           :: phi(Bs+2*g, Bs+2*g, Bs+2*g, N_dF)
    !> discretization order
    character(len=80), intent(in)                          :: order_discretization
    !> global volume integral of the last time step
    real(kind=rk), dimension(3), intent(in)                :: volume_int
    !> time
    real(kind=rk), intent(in)                              :: time

    !> RHS
    real(kind=rk), dimension(Bs+2*g, Bs+2*g, Bs+2*g, N_dF) :: rhs
    
    !> mask term for every grid point in this block
    real(kind=rk), dimension(Bs+2*g, Bs+2*g, Bs+2*g)       :: mask
    !> velocity of the solid
    real(kind=rk), dimension(Bs+2*g, Bs+2*g ,Bs+2*g, 3)    :: us
    !> forcing term
    real(kind=rk), dimension(3)                            :: forcing
    !> local datafields
    real(kind=rk), dimension(Bs+2*g, Bs+2*g, Bs+2*g )      :: u, v, w, p

    !> inverse dx, physics/acm parameters
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, c_0,&
                                                      nu, eps, eps_inv, gamma
    !> derivatives
    real(kind=rk)                                  :: div_U, u_dx, u_dy, u_dz, u_dxdx, u_dydy, u_dzdz, &
                                                     v_dx, v_dy, v_dz, v_dxdx, v_dydy, v_dzdz, &
                                                     w_dx, w_dy, w_dz, w_dxdx, w_dydy, w_dzdz, & 
                                                     p_dx, p_dy, p_dz, penalx, penaly, penalz
    !> loop variables
    integer(kind=rk)                               :: ix, iy, iz
    !> coefficients for Tam&Webb
    real(kind=rk)                                  :: a(-3:3)
    real(kind=rk)                                  :: b(-2:2)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set parameters for readability
    c_0         = params%physics_acm%c_0
    nu          = params%physics_acm%nu
    eps         = params%eps_penal
    gamma       = params%physics_acm%gamma_p

    
    u = phi(:,:,:,1)
    v = phi(:,:,:,2)
    p = phi(:,:,:,3)

    rhs  = 0.0_rk
    mask = 0.0_rk
    us   = 0.0_rk

    dx_inv = 1.0_rk / (2.0_rk*dx(1))
    dy_inv = 1.0_rk / (2.0_rk*dx(2))
    dz_inv = 1.0_rk / (2.0_rk*dx(3))

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    dz2_inv = 1.0_rk / (dx(3)**2)

    eps_inv = 1.0_rk / eps

    ! Tam & Webb, 4th order optimized (for first derivative)
    a=(/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
         0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    ! 4th order coefficients for second derivative
    b = (/ -1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, &
            4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)

!---------------------------------------------------------------------------------------------
! main body

    if (params%penalization) then
        ! create mask term for every grid point in this block
        call create_mask_3D(params, mask, x0, dx, Bs, g)
        mask = mask*eps_inv
    end if

    call  compute_forcing(forcing, volume_int, params%Lx, params%Ly, params%Lz, time)

   if (order_discretization == "FD_2nd_central" ) then
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        do ix = g+1, Bs+g
            do iy = g+1, Bs+g
                do iz = g+1, Bs+g

                    ! first and second derivatives of u,v,w
                    u_dx = (u(ix+1,iy,iz)-u(ix-1,iy,iz))*dx_inv
                    u_dy = (u(ix,iy+1,iz)-u(ix,iy-1,iz))*dy_inv
                    u_dz = (u(ix,iy,iz+1)-u(ix,iy,iz-1))*dz_inv

                    u_dxdx = (u(ix-1,iy,iz)-2.0_rk*u(ix,iy,iz)+u(ix+1,iy,iz))*dx2_inv
                    u_dydy = (u(ix,iy-1,iz)-2.0_rk*u(ix,iy,iz)+u(ix,iy+1,iz))*dy2_inv
                    u_dzdz = (u(ix,iy,iz-1)-2.0_rk*u(ix,iy,iz)+u(ix,iy,iz+1))*dz2_inv

                    v_dx = (v(ix+1,iy,iz)-v(ix-1,iy,iz))*dx_inv
                    v_dy = (v(ix,iy+1,iz)-v(ix,iy-1,iz))*dy_inv
                    v_dz = (v(ix,iy,iz+1)-v(ix,iy,iz-1))*dz_inv

                    v_dxdx = (v(ix-1,iy,iz)-2.0_rk*v(ix,iy,iz)+v(ix+1,iy,iz))*dx2_inv
                    v_dydy = (v(ix,iy-1,iz)-2.0_rk*v(ix,iy,iz)+v(ix,iy+1,iz))*dy2_inv
                    v_dzdz = (v(ix,iy,iz-1)-2.0_rk*v(ix,iy,iz)+v(ix,iy,iz+1))*dz2_inv

                    w_dx = (w(ix+1,iy,iz)-w(ix-1,iy,iz))*dx_inv
                    w_dy = (w(ix,iy+1,iz)-w(ix,iy-1,iz))*dy_inv
                    w_dz = (w(ix,iy,iz+1)-w(ix,iy,iz-1))*dz_inv

                    w_dxdx = (w(ix-1,iy,iz)-2.0_rk*w(ix,iy,iz)+w(ix+1,iy,iz))*dx2_inv
                    w_dydy = (w(ix,iy-1,iz)-2.0_rk*w(ix,iy,iz)+w(ix,iy+1,iz))*dy2_inv
                    w_dzdz = (w(ix,iy,iz-1)-2.0_rk*w(ix,iy,iz)+w(ix,iy,iz+1))*dz2_inv
                    
                    ! first derivative of p
                    p_dx = (p(ix+1,iy,iz)-p(ix-1,iy,iz))*dx_inv
                    p_dy = (p(ix,iy+1,iz)-p(ix,iy-1,iz))*dy_inv
                    p_dy = (p(ix,iy,iz+1)-p(ix,iy,iz-1))*dz_inv

                    div_U = u_dx + v_dy + w_dz

                    penalx = -mask(ix,iy,iz)*(u(ix,iy,iz)-us(ix,iy,iz,1))
                    penaly = -mask(ix,iy,iz)*(v(ix,iy,iz)-us(ix,iy,iz,2))
                    penalz = -mask(ix,iy,iz)*(w(ix,iy,iz)-us(ix,iy,iz,3))

                    rhs(ix,iy,iz,1) = -u(ix,iy,iz)*u_dx - v(ix,iy,iz)*u_dy - w(ix,iy,iz)*u_dz - p_dx &
                                   + nu*(u_dxdx + u_dydy + u_dzdz) + penalx + forcing(1)
                    rhs(ix,iy,iz,2) = -u(ix,iy,iz)*v_dx - v(ix,iy,iz)*v_dy - w(ix,iy,iz)*v_dz - p_dy &
                                    + nu*(v_dxdx + v_dydy + v_dzdz) + penaly + forcing(2)
                    rhs(ix,iy,iz,3) = -u(ix,iy,iz)*w_dx - v(ix,iy,iz)*w_dy - w(ix,iy,iz)*w_dz - p_dz &
                                    + nu*(w_dxdx + w_dydy + w_dzdz) + penalz + forcing(3)
                    rhs(ix,iy,iz,4) = -(c_0**2)*div_U - gamma*p(ix,iy,iz)
                end do
            end do
        end do

    else if (order_discretization == "FD_4th_central_optimized") then
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do ix = g+1, Bs+g
            do iy = g+1, Bs+g
                do iz = g+1, Bs+g
                    ! first derivatives of u, v, p
                    u_dx = (a(-3)*u(ix-3,iy,iz) + a(-2)*u(ix-2,iy,iz) + a(-1)*u(ix-1,iy,iz) + a(0)*u(ix,iy,iz)&
                     +  a(+1)*u(ix+1,iy,iz) + a(+2)*u(ix+2,iy,iz) + a(+3)*u(ix+3,iy,iz))*dx_inv
                    u_dy = (a(-3)*u(ix,iy-3,iz) + a(-2)*u(ix,iy-2,iz) + a(-1)*u(ix,iy-1,iz) + a(0)*u(ix,iy,iz)&
                     +  a(+1)*u(ix,iy+1,iz) + a(+2)*u(ix,iy+2,iz) + a(+3)*u(ix,iy+3,iz))*dy_inv
                    u_dz = (a(-3)*u(ix,iy,iz-3) + a(-2)*u(ix,iy,iz-2) + a(-1)*u(ix,iy,iz-1) + a(0)*u(ix,iy,iz)&
                     +  a(+1)*u(ix,iy,iz+1) + a(+2)*u(ix,iy,iz+2) + a(+3)*u(ix,iy,iz+3))*dz_inv

                    v_dx = (a(-3)*v(ix-3,iy,iz) + a(-2)*v(ix-2,iy,iz) + a(-1)*v(ix-1,iy,iz) + a(0)*v(ix,iy,iz)&
                     +  a(+1)*v(ix+1,iy,iz) + a(+2)*v(ix+2,iy,iz) + a(+3)*v(ix+3,iy,iz))*dx_inv
                    v_dy = (a(-3)*v(ix,iy-3,iz) + a(-2)*v(ix,iy-2,iz) + a(-1)*v(ix,iy-1,iz) + a(0)*v(ix,iy,iz)&
                     +  a(+1)*v(ix,iy+1,iz) + a(+2)*v(ix,iy+2,iz) + a(+3)*v(ix,iy+3,iz))*dy_inv
                    v_dz = (a(-3)*v(ix,iy,iz-3) + a(-2)*v(ix,iy,iz-2) + a(-1)*v(ix,iy,iz-1) + a(0)*v(ix,iy,iz)&
                     +  a(+1)*v(ix,iy,iz+1) + a(+2)*v(ix,iy,iz+2) + a(+3)*v(ix,iy,iz+3))*dz_inv

                    w_dx = (a(-3)*w(ix-3,iy,iz) + a(-2)*w(ix-2,iy,iz) + a(-1)*w(ix-1,iy,iz) + a(0)*w(ix,iy,iz)&
                     +  a(+1)*w(ix+1,iy,iz) + a(+2)*w(ix+2,iy,iz) + a(+3)*w(ix+3,iy,iz))*dx_inv
                    w_dy = (a(-3)*w(ix,iy-3,iz) + a(-2)*w(ix,iy-2,iz) + a(-1)*w(ix,iy-1,iz) + a(0)*w(ix,iy,iz)&
                     +  a(+1)*w(ix,iy+1,iz) + a(+2)*w(ix,iy+2,iz) + a(+3)*w(ix,iy+3,iz))*dy_inv
                    w_dz = (a(-3)*w(ix,iy,iz-3) + a(-2)*w(ix,iy,iz-2) + a(-1)*w(ix,iy,iz-1) + a(0)*w(ix,iy,iz)&
                     +  a(+1)*w(ix,iy,iz+1) + a(+2)*w(ix,iy,iz+2) + a(+3)*w(ix,iy,iz+3))*dz_inv

                    p_dx = (a(-3)*p(ix-3,iy,iz) + a(-2)*p(ix-2,iy,iz) + a(-1)*p(ix-1,iy,iz) + a(0)*p(ix,iy,iz)&
                     +  a(+1)*p(ix+1,iy,iz) + a(+2)*p(ix+2,iy,iz) + a(+3)*p(ix+3,iy,iz))*dx_inv
                    p_dy = (a(-3)*p(ix,iy-3,iz) + a(-2)*p(ix,iy-2,iz) + a(-1)*p(ix,iy-1,iz) + a(0)*p(ix,iy,iz)&
                     +  a(+1)*p(ix,iy+1,iz) + a(+2)*p(ix,iy+2,iz) + a(+3)*u(ix,iy+3,iz))*dy_inv
                    p_dz = (a(-3)*p(ix,iy,iz-3) + a(-2)*p(ix,iy,iz-2) + a(-1)*p(ix,iy,iz-1) + a(0)*p(ix,iy,iz)&
                     +  a(+1)*p(ix,iy,iz+1) + a(+2)*p(ix,iy,iz+2) + a(+3)*p(ix,iy,iz+3))*dx_inv

                   ! second derivatives of u, v and w
                    u_dxdx = (b(-2)*u(ix-2,iy,iz) + b(-1)*u(ix-1,iy,iz) + b(0)*u(ix,iy,iz)&
                      +  b(1)*u(ix+1,iy,iz) + b(2)*u(ix+2,iy,iz))*dx2_inv
                    u_dydy = (b(-2)*u(ix,iy-2,iz) + b(-1)*u(ix,iy-1,iz) + b(0)*u(ix,iy,iz)&
                      +  b(1)*u(ix,iy+1,iz) + b(2)*u(ix,iy+2,iz))*dy2_inv
                    u_dzdz = (b(-2)*u(ix,iy,iz-2) + b(-1)*u(ix,iy,iz-1) + b(0)*u(ix,iy,iz)&
                      +  b(1)*u(ix,iy,iz+1) + b(2)*u(ix,iy,iz+2))*dz2_inv
                    v_dxdx = (b(-2)*v(ix-2,iy,iz) + b(-1)*v(ix-1,iy,iz) + b(0)*v(ix,iy,iz)&
                      +  b(1)*v(ix+1,iy,iz) + b(2)*v(ix+2,iy,iz))*dx2_inv
                    v_dydy = (b(-2)*v(ix,iy-2,iz) + b(-1)*v(ix,iy-1,iz) + b(0)*v(ix,iy,iz)&
                      +  b(1)*v(ix,iy+1,iz) + b(2)*v(ix,iy+2,iz))*dy2_inv
                    v_dzdz = (b(-2)*v(ix,iy,iz-2) + b(-1)*v(ix,iy,iz-1) + b(0)*v(ix,iy,iz)&
                      +  b(1)*v(ix,iy,iz+1) + b(2)*v(ix,iy,iz+2))*dz2_inv
                    w_dxdx = (b(-2)*w(ix-2,iy,iz) + b(-1)*w(ix-1,iy,iz) + b(0)*w(ix,iy,iz)&
                      +  b(1)*w(ix+1,iy,iz) + b(2)*w(ix+2,iy,iz))*dx2_inv
                    w_dydy = (b(-2)*w(ix,iy-2,iz) + b(-1)*w(ix,iy-1,iz) + b(0)*w(ix,iy,iz)&
                      +  b(1)*w(ix,iy+1,iz) + b(2)*w(ix,iy+2,iz))*dy2_inv
                    w_dzdz = (b(-2)*w(ix,iy,iz-2) + b(-1)*w(ix,iy,iz-1) + b(0)*w(ix,iy,iz)&
                      +  b(1)*w(ix,iy,iz+1) + b(2)*w(ix,iy,iz+2))*dz2_inv

                    div_U = u_dx + v_dy + w_dz

                    penalx = -mask(ix,iy,iz)*(u(ix,iy,iz)-us(ix,iy,iz,1))
                    penaly = -mask(ix,iy,iz)*(v(ix,iy,iz)-us(ix,iy,iz,2))
                    penalz = -mask(ix,iy,iz)*(w(ix,iy,iz)-us(ix,iy,iz,3))

                    rhs(ix,iy,iz,1) = -u(ix,iy,iz)*u_dx - v(ix,iy,iz)*u_dy - w(ix,iy,iz)*u_dz - p_dx &
                                   + nu*(u_dxdx + u_dydy + u_dzdz) + penalx + forcing(1)
                    rhs(ix,iy,iz,2) = -u(ix,iy,iz)*v_dx - v(ix,iy,iz)*v_dy - w(ix,iy,iz)*v_dz - p_dy &
                                    + nu*(v_dxdx + v_dydy + v_dzdz) + penaly + forcing(2)
                    rhs(ix,iy,iz,3) = -u(ix,iy,iz)*w_dx - v(ix,iy,iz)*w_dy - w(ix,iy,iz)*w_dz - p_dz &
                                    + nu*(w_dxdx + w_dydy + w_dzdz) + penalz + forcing(3)
                    rhs(ix,iy,iz,4) = -(c_0**2)*div_U - gamma*p(ix,iy,iz)
                end do
            end do
        end do

    else
      write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
      write(*,*) order_discretization
      stop
    end if

    !> \todo DO NOT OVERWRITE?
    phi = rhs

end subroutine RHS_3D_acm
