!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_2D_acm.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief RHS for 2D artificial compressibility method
!
!>
!! input:    - datafield, grid parameter, derivative order \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 27/06/17 - create
! ********************************************************************************************

subroutine RHS_2D_acm(params, g, Bs, dx, x0, N_dF, phi, order_discretization, volume_int, time)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params
    use module_operators

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> physics parameter structure
    type (type_params), intent(in)                 :: params

    !> grid parameter
    integer(kind=ik), intent(in)                   :: g, Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx

    !> number of datafields
    integer(kind=ik), intent(in)                   :: N_dF
    !> datafields
    real(kind=rk), intent(inout)                   :: phi(Bs+2*g, Bs+2*g, N_dF)
    !> discretization order
    character(len=80), intent(in)                  :: order_discretization
    !> global volume integral
    real(kind=rk), dimension(3), intent(in)        :: volume_int
    !> time
    real(kind=rk), intent(in)                      :: time

    !> RHS
    real(kind=rk), dimension(Bs+2*g, Bs+2*g, N_dF) :: rhs

    !> mask term for every grid point in this block
    real(kind=rk), dimension(Bs+2*g, Bs+2*g)       :: mask, sponge
    !> velocity of the solid
    real(kind=rk), dimension(Bs+2*g, Bs+2*g, 2)    :: us
    !> forcing term
    real(kind=rk), dimension(3)                    :: forcing


    !> local datafields
    real(kind=rk), dimension(Bs+2*g, Bs+2*g)       :: u, v, p
    !>
    real(kind=rk)                                  :: dx_inv, dy_inv, dx2_inv, dy2_inv, c_0, nu, eps, eps_inv, gamma
    real(kind=rk)                                  :: div_U, u_dx, u_dy, u_dxdx, u_dydy, v_dx, v_dy, v_dxdx, &
                                                      v_dydy, p_dx, p_dy, penalx, penaly, alpha
    ! loop variables
    integer(kind=rk)                               :: ix, iy
    ! coefficients for Tam&Webb
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

    u = phi(:,:,1)
    v = phi(:,:,2)
    p = phi(:,:,3)

    rhs    = 0.0_rk
    mask   = 0.0_rk
    us     = 0.0_rk
    sponge = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    eps_inv = 1.0_rk / eps

    alpha = 100.0_rk

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! 4th order coefficients for second derivative
    b = (/ -1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
!---------------------------------------------------------------------------------------------
! main body

    if (params%penalization) then
        ! create mask term for every grid point in this block
        call create_mask_2D(params, mask, x0, dx, Bs, g)
        mask = mask*eps_inv
    end if

    if (params%physics_acm%forcing) then
        call compute_forcing(forcing, volume_int, params%Lx, params%Ly, 1.0_rk, time)
    else
        forcing = 0.0_rk
    end if

    call sponge_2D(params, sponge, x0, dx, Bs, g)
    sponge=alpha*sponge

    if (order_discretization == "FD_2nd_central" ) then
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        do ix = g+1, Bs+g
            do iy = g+1, Bs+g

                u_dx = (u(ix+1,iy)-u(ix-1,iy))*dx_inv*0.5_rk
                u_dy = (u(ix,iy+1)-u(ix,iy-1))*dy_inv*0.5_rk
                u_dxdx = (u(ix-1,iy)-2.0_rk*u(ix,iy)+u(ix+1,iy))*dx2_inv
                u_dydy = (u(ix,iy-1)-2.0_rk*u(ix,iy)+u(ix,iy+1))*dy2_inv

                v_dx = (v(ix+1,iy)-v(ix-1,iy))*dx_inv*0.5_rk
                v_dy = (v(ix,iy+1)-v(ix,iy-1))*dy_inv*0.5_rk
                v_dxdx = (v(ix-1,iy)-2.0_rk*v(ix,iy)+v(ix+1,iy))*dx2_inv
                v_dydy = (v(ix,iy-1)-2.0_rk*v(ix,iy)+v(ix,iy+1))*dy2_inv

                p_dx = (p(ix+1,iy)-p(ix-1,iy))*dx_inv*0.5_rk
                p_dy = (p(ix,iy+1)-p(ix,iy-1))*dy_inv*0.5_rk

                div_U = u_dx + v_dy

                penalx = -mask(ix,iy)*(u(ix,iy)-us(ix,iy,1))-alpha*sponge(ix,iy)*(u(ix,iy)-1.0_rk)
                penaly = -mask(ix,iy)*(v(ix,iy)-us(ix,iy,2))

                rhs(ix,iy,1) = -u(ix,iy)*u_dx - v(ix,iy)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx + forcing(1)
                rhs(ix,iy,2) = -u(ix,iy)*v_dx - v(ix,iy)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly + forcing(2)
                rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*p(ix,iy)

            end do
        end do

    else if (order_discretization == "FD_4th_central_optimized") then
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do ix = g+1, Bs+g
          do iy = g+1, Bs+g

            ! first derivatives of u, v, p
            u_dx = (a(-3)*u(ix-3,iy) + a(-2)*u(ix-2,iy) + a(-1)*u(ix-1,iy) + a(0)*u(ix,iy) + a(+1)*u(ix+1,iy) + a(+2)*u(ix+2,iy) + a(+3)*u(ix+3,iy))*dx_inv
            u_dy = (a(-3)*u(ix,iy-3) + a(-2)*u(ix,iy-2) + a(-1)*u(ix,iy-1) + a(0)*u(ix,iy) + a(+1)*u(ix,iy+1) + a(+2)*u(ix,iy+2) + a(+3)*u(ix,iy+3))*dy_inv
            v_dx = (a(-3)*v(ix-3,iy) + a(-2)*v(ix-2,iy) + a(-1)*v(ix-1,iy) + a(0)*v(ix,iy) + a(+1)*v(ix+1,iy) + a(+2)*v(ix+2,iy) + a(+3)*v(ix+3,iy))*dx_inv
            v_dy = (a(-3)*v(ix,iy-3) + a(-2)*v(ix,iy-2) + a(-1)*v(ix,iy-1) + a(0)*v(ix,iy) + a(+1)*v(ix,iy+1) + a(+2)*v(ix,iy+2) + a(+3)*v(ix,iy+3))*dy_inv
            p_dx = (a(-3)*p(ix-3,iy) + a(-2)*p(ix-2,iy) + a(-1)*p(ix-1,iy) + a(0)*p(ix,iy) + a(+1)*p(ix+1,iy) + a(+2)*p(ix+2,iy) + a(+3)*p(ix+3,iy))*dx_inv
            p_dy = (a(-3)*p(ix,iy-3) + a(-2)*p(ix,iy-2) + a(-1)*p(ix,iy-1) + a(0)*p(ix,iy) + a(+1)*p(ix,iy+1) + a(+2)*p(ix,iy+2) + a(+3)*p(ix,iy+3))*dy_inv

            ! second derivatives of u and v
            u_dxdx = (b(-2)*u(ix-2,iy) + b(-1)*u(ix-1,iy) + b(0)*u(ix,iy) + b(+1)*u(ix+1,iy) + b(+2)*u(ix+2,iy))*dx2_inv
            u_dydy = (b(-2)*u(ix,iy-2) + b(-1)*u(ix,iy-1) + b(0)*u(ix,iy) + b(+1)*u(ix,iy+1) + b(+2)*u(ix,iy+2))*dy2_inv
            v_dxdx = (b(-2)*v(ix-2,iy) + b(-1)*v(ix-1,iy) + b(0)*v(ix,iy) + b(+1)*v(ix+1,iy) + b(+2)*v(ix+2,iy))*dx2_inv
            v_dydy = (b(-2)*v(ix,iy-2) + b(-1)*v(ix,iy-1) + b(0)*v(ix,iy) + b(+1)*v(ix,iy+1) + b(+2)*v(ix,iy+2))*dy2_inv

            div_U = u_dx + v_dy

            penalx = -mask(ix,iy)*(u(ix,iy)-us(ix,iy,1))
            penaly = -mask(ix,iy)*(v(ix,iy)-us(ix,iy,2))

            rhs(ix,iy,1) = -u(ix,iy)*u_dx - v(ix,iy)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx + forcing(1)
            rhs(ix,iy,2) = -u(ix,iy)*v_dx - v(ix,iy)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly + forcing(2)
            rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*p(ix,iy)
          end do
        end do

    else
      write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
      write(*,*) order_discretization
      stop
    end if

    !> \todo DO NOT OVERWRITE?
    phi = rhs

end subroutine RHS_2D_acm
