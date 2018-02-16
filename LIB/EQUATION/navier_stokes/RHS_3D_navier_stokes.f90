!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_navier_stokes.f90
!> \version 0.5
!> \author msr
!
!> \brief RHS for 3D navier stokes equation
!
!>
!! input:    - datafield, grid parameter, velocity, diffusion coefficient, derivative order \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 14/02/17 - create
!
! ********************************************************************************************

subroutine RHS_3D_navier_stokes(g, Bs, dx, dy, dz, phi,rhs)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> grid parameter
    integer(kind=ik), intent(in)                            :: g, Bs
    !> rhs parameter
    real(kind=rk), intent(in)                               :: dx, dy, dz
    !> datafields
    real(kind=rk), intent(in)                              :: phi(:, :, :, :)
    ! rhs array
    real(kind=rk),intent(out)                               :: rhs(:, :, :,:)

     ! adiabatic coefficien t
    real(kind=rk)                                           :: gamma_
    ! specific gas constant
    real(kind=rk)                                           :: Rs
    ! isochoric heat capacity
    real(kind=rk)                                           :: Cv
    ! isobaric heat capacity
    real(kind=rk)                                           :: Cp
    ! prandtl number
    real(kind=rk)                                           :: Pr
    ! dynamic viscosity
    real(kind=rk)                                           :: mu0
    ! dissipation switch
    logical                                                 :: dissipation

    ! variables
    real(kind=rk)                                           :: rho(Bs+2*g, Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g, Bs+2*g), w(Bs+2*g, Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               T(Bs+2*g, Bs+2*g, Bs+2*g), mu(Bs+2*g, Bs+2*g, Bs+2*g), mu_d(Bs+2*g, Bs+2*g, Bs+2*g), lambda(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               fric_p(Bs+2*g, Bs+2*g, Bs+2*g), fric_u(Bs+2*g, Bs+2*g, Bs+2*g), fric_v(Bs+2*g, Bs+2*g, Bs+2*g), fric_w(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               fric_T1(Bs+2*g, Bs+2*g, Bs+2*g), fric_T2(Bs+2*g, Bs+2*g, Bs+2*g), fric_T3(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               tau11(Bs+2*g, Bs+2*g, Bs+2*g), tau22(Bs+2*g, Bs+2*g, Bs+2*g), tau33(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               tau12(Bs+2*g, Bs+2*g, Bs+2*g), tau13(Bs+2*g, Bs+2*g, Bs+2*g), tau23(Bs+2*g, Bs+2*g, Bs+2*g)
    ! derivatives
    real(kind=rk)                                           :: u_x(Bs+2*g, Bs+2*g, Bs+2*g), u_y(Bs+2*g, Bs+2*g, Bs+2*g), u_z(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               v_x(Bs+2*g, Bs+2*g, Bs+2*g), v_y(Bs+2*g, Bs+2*g, Bs+2*g), v_z(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               w_x(Bs+2*g, Bs+2*g, Bs+2*g), w_y(Bs+2*g, Bs+2*g, Bs+2*g), w_z(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               p_x(Bs+2*g, Bs+2*g, Bs+2*g), p_y(Bs+2*g, Bs+2*g, Bs+2*g), p_z(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               T_x(Bs+2*g, Bs+2*g, Bs+2*g), T_y(Bs+2*g, Bs+2*g, Bs+2*g), T_z(Bs+2*g, Bs+2*g, Bs+2*g), &
                                                               div_U(Bs+2*g, Bs+2*g, Bs+2*g)

    ! dummy field
    real(kind=rk)                                           :: dummy(Bs+2*g, Bs+2*g, Bs+2*g)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set physics parameters for readability
    gamma_      = params_ns%gamma_
    Rs          = params_ns%Rs
    Cv          = params_ns%Cv
    Cp          = params_ns%Cp
    Pr          = params_ns%Pr
    mu0         = params_ns%mu0
    dissipation = params_ns%dissipation

    ! variables
    rho         = phi(:,:,:,1)**2
    u           = phi(:,:,:,2)/phi(:,:,:,1)
    v           = phi(:,:,:,3)/phi(:,:,:,1)
    w           = phi(:,:,:,4)/phi(:,:,:,1)
    p           = phi(:,:,:,5)

    ! rhs
    rhs         = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! derivatives
    call grad_zentral_3D( Bs, g, dx, dy, dz, u, u_x, u_y, u_z)
    call grad_zentral_3D( Bs, g, dx, dy, dz, v, v_x, v_y, v_z)
    call grad_zentral_3D( Bs, g, dx, dy, dz, w, w_x, w_y, w_z)

    call diff1x_zentral_3D( Bs, g, dx, p, p_x)
    call diff1y_zentral_3D( Bs, g, dy, p, p_y)
    call diff1z_zentral_3D( Bs, g, dz, p, p_z)

    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    call diff1x_zentral_3D( Bs, g, dx, rho*u, dummy)
    rhs(:,:,:,1) = -dummy
    call diff1y_zentral_3D( Bs, g, dy, rho*v, dummy)
    rhs(:,:,:,1) = rhs(:,:,:,1) - dummy
    call diff1z_zentral_3D( Bs, g, dz, rho*w, dummy)
    rhs(:,:,:,1) = rhs(:,:,:,1) - dummy

    rhs(:,:,:,1) = rhs(:,:,:,1) * 0.5_rk/phi(:,:,:,1)

    ! friction
    if (dissipation) then

        ! Compute mu
        T    = p / (rho*Rs)
        mu   = mu0
        mu_d = 0.0_rk

        ! thermal conductivity
        lambda  = Cp * mu/Pr

        ! tau11
        tau11 = mu * 2.0_rk * u_x

        call diff1x_zentral_3D( Bs, g, dx, u, dummy)
        div_U = dummy
        call diff1y_zentral_3D( Bs, g, dy, v, dummy)
        div_U = div_U + dummy
        call diff1z_zentral_3D( Bs, g, dz, w, dummy)
        div_U = div_U + dummy

        tau11 = tau11 + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U

        ! tau22
        tau22 = mu * 2.0_rk * v_y
        tau22 = tau22 + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U

        ! tau33
        tau33 = mu * 2.0_rk * w_z
        tau33 = tau33 + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U

        ! tau12
        tau12 = mu * ( v_x + u_y )

        ! tau13
        tau13 = mu * ( w_x + u_z )

        ! tau23
        tau23 = mu * ( w_y + v_z )

        ! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
        call diff1x_zentral_3D( Bs, g, dx, tau11, dummy)
        fric_u = dummy
        call diff1y_zentral_3D( Bs, g, dy, tau12, dummy)
        fric_u = fric_u + dummy
        call diff1z_zentral_3D( Bs, g, dz, tau13, dummy)
        fric_u = fric_u + dummy

        fric_u = fric_u / phi(:,:,:,1)

        call diff1x_zentral_3D( Bs, g, dx, tau12, dummy)
        fric_v = dummy
        call diff1y_zentral_3D( Bs, g, dy, tau22, dummy)
        fric_v = fric_v + dummy
        call diff1z_zentral_3D( Bs, g, dz, tau23, dummy)
        fric_v = fric_v + dummy

        fric_v = fric_v / phi(:,:,:,1)

        call diff1x_zentral_3D( Bs, g, dx, tau13, dummy)
        fric_w = dummy
        call diff1y_zentral_3D( Bs, g, dy, tau23, dummy)
        fric_w = fric_w + dummy
        call diff1z_zentral_3D( Bs, g, dz, tau33, dummy)
        fric_w = fric_w + dummy

        fric_w = fric_w / phi(:,:,:,1)

        ! Friction terms for the energy equation
        ! Heat Flux
        call grad_zentral_3D( Bs, g, dx, dy, dz, T, T_x, T_y, T_z)

        fric_T1 = lambda * T_x
        fric_T2 = lambda * T_y
        fric_T3 = lambda * T_z

        ! All simple divergence terms for u_i*tau_ik and phi_k
        call diff1x_zentral_3D( Bs, g, dx, ( u*tau11 + v*tau12 + w*tau13 + fric_T1 ), dummy)
        fric_p = dummy
        call diff1y_zentral_3D( Bs, g, dy, ( u*tau12 + v*tau22 + w*tau23 + fric_T2 ), dummy)
        fric_p = fric_p + dummy
        call diff1z_zentral_3D( Bs, g, dz, ( u*tau12 + v*tau22 + w*tau33 + fric_T3 ), dummy)
        fric_p = fric_p + dummy

        ! u_i*dx_k (tau_ik) terms
        call diff1x_zentral_3D( Bs, g, dx, tau11, dummy)
        fric_p = fric_p - u*dummy
        call diff1y_zentral_3D( Bs, g, dy, tau12, dummy)
        fric_p = fric_p - u*dummy
        call diff1z_zentral_3D( Bs, g, dz, tau13, dummy)
        fric_p = fric_p - u*dummy

        call diff1x_zentral_3D( Bs, g, dx, tau12, dummy)
        fric_p = fric_p - v*dummy
        call diff1y_zentral_3D( Bs, g, dy, tau22, dummy)
        fric_p = fric_p - v*dummy
        call diff1z_zentral_3D( Bs, g, dz, tau23, dummy)
        fric_p = fric_p - v*dummy

        call diff1x_zentral_3D( Bs, g, dx, tau13, dummy)
        fric_p = fric_p - w*dummy
        call diff1y_zentral_3D( Bs, g, dy, tau23, dummy)
        fric_p = fric_p - w*dummy
        call diff1z_zentral_3D( Bs, g, dz, tau33, dummy)
        fric_p = fric_p - w*dummy

        fric_p = ( gamma_ - 1 ) * fric_p

    else

        fric_p = 0.0_rk
        fric_u = 0.0_rk
        fric_v = 0.0_rk
        fric_w = 0.0_rk

    end if

    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    call diff1x_zentral_3D( Bs, g, dx, (u * p), dummy)
    rhs(:,:,:,5) = - dummy
    call diff1y_zentral_3D( Bs, g, dy, (v * p), dummy)
    rhs(:,:,:,5) = rhs(:,:,:,5) - dummy
    call diff1z_zentral_3D( Bs, g, dz, (w * p), dummy)
    rhs(:,:,:,5) = rhs(:,:,:,5) - dummy

    rhs(:,:,:,5) = rhs(:,:,:,5) * gamma_

    rhs(:,:,:,5) = rhs(:,:,:,5) + (gamma_ - 1.0_rk) * (u*p_x + v*p_y + w*p_z)

    rhs(:,:,:,5) = rhs(:,:,:,5) + fric_p

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    call diff1x_zentral_3D( Bs, g, dx, (u * rho * u), dummy)
    rhs(:,:,:,2) = - 0.5_rk * dummy
    call diff1y_zentral_3D( Bs, g, dy, (v * rho * u), dummy)
    rhs(:,:,:,2) = rhs(:,:,:,2) - 0.5_rk * dummy
    call diff1z_zentral_3D( Bs, g, dz, (w * rho * u), dummy)
    rhs(:,:,:,2) = rhs(:,:,:,2) - 0.5_rk * dummy

    rhs(:,:,:,2) = rhs(:,:,:,2) - 0.5_rk * rho * u * u_x
    rhs(:,:,:,2) = rhs(:,:,:,2) - 0.5_rk * rho * v * u_y
    rhs(:,:,:,2) = rhs(:,:,:,2) - 0.5_rk * rho * w * u_z

    rhs(:,:,:,2) = rhs(:,:,:,2) - p_x

    rhs(:,:,:,2) = rhs(:,:,:,2) / phi(:,:,:,1)

    rhs(:,:,:,2) = rhs(:,:,:,2) + fric_u

    ! RHS of  momentum equation for v
    call diff1x_zentral_3D( Bs, g, dx, (u * rho * v), dummy)
    rhs(:,:,:,3) = - 0.5_rk * dummy
    call diff1y_zentral_3D( Bs, g, dy, (v * rho * v), dummy)
    rhs(:,:,:,3) = rhs(:,:,:,3) - 0.5_rk * dummy
    call diff1z_zentral_3D( Bs, g, dz, (w * rho * v), dummy)
    rhs(:,:,:,3) = rhs(:,:,:,3) - 0.5_rk * dummy

    rhs(:,:,:,3) = rhs(:,:,:,3) - 0.5_rk * rho * u * v_x
    rhs(:,:,:,3) = rhs(:,:,:,3) - 0.5_rk * rho * v * v_y
    rhs(:,:,:,3) = rhs(:,:,:,3) - 0.5_rk * rho * w * v_z

    rhs(:,:,:,3) = rhs(:,:,:,3) - p_y

    rhs(:,:,:,3) = rhs(:,:,:,3) / phi(:,:,:,1)

    rhs(:,:,:,3) = rhs(:,:,:,3) + fric_v

    ! RHS of  momentum equation for v
    call diff1x_zentral_3D( Bs, g, dx, (u * rho * w), dummy)
    rhs(:,:,:,4) = - 0.5_rk * dummy
    call diff1y_zentral_3D( Bs, g, dy, (v * rho * w), dummy)
    rhs(:,:,:,4) = rhs(:,:,:,4) - 0.5_rk * dummy
    call diff1z_zentral_3D( Bs, g, dz, (w * rho * w), dummy)
    rhs(:,:,:,4) = rhs(:,:,:,4) - 0.5_rk * dummy

    rhs(:,:,:,4) = rhs(:,:,:,4) - 0.5_rk * rho * u * v_x
    rhs(:,:,:,4) = rhs(:,:,:,4) - 0.5_rk * rho * v * v_y
    rhs(:,:,:,4) = rhs(:,:,:,4) - 0.5_rk * rho * w * v_z

    rhs(:,:,:,4) = rhs(:,:,:,4) - p_z

    rhs(:,:,:,4) = rhs(:,:,:,4) / phi(:,:,:,1)

    rhs(:,:,:,4) = rhs(:,:,:,4) + fric_w

end subroutine RHS_3D_navier_stokes

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

subroutine grad_zentral_3D(Bs, g, dx, dy, dz, q, qx, qy, qz)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dx, dy, dz
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qx(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qy(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qz(Bs+2*g, Bs+2*g, Bs+2*g)

    !! XXX !!
    call diffx_c_3D( Bs, g, dx, q, qx)

    !! YYY !!
    call diffy_c_3D( Bs, g, dy, q, qy)

    !! ZZZ !!
    call diffz_c_3D( Bs, g, dz, q, qz)

end subroutine grad_zentral_3D

!---------------------------------------------------------------------------------------------

subroutine diff1x_zentral_3D(Bs, g, dx, q, qx)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qx(Bs+2*g, Bs+2*g, Bs+2*g)

    !! XXX !!
    call diffx_c_3D( Bs, g, dx, q, qx)

end subroutine diff1x_zentral_3D

!---------------------------------------------------------------------------------------------

subroutine diff1y_zentral_3D(Bs, g, dy, q, qy)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qy(Bs+2*g, Bs+2*g, Bs+2*g)

    !! YYY !!
    call diffy_c_3D( Bs, g, dy, q, qy)

end subroutine diff1y_zentral_3D

!---------------------------------------------------------------------------------------------

subroutine diff1z_zentral_3D(Bs, g, dz, q, qz)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qz(Bs+2*g, Bs+2*g, Bs+2*g)

    !! ZZZ !!
    call diffz_c_3D( Bs, g, dz, q, qz)

end subroutine diff1z_zentral_3D

!---------------------------------------------------------------------------------------------

subroutine  diffx_c_3D( Bs, g, dx, u, dudx)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs+2*g, Bs+2*g, Bs+2*g)

    integer                         :: i, n

    n = size(u,1)

    dudx(1,:,:) = ( u(n-1,:,:) - 8.0_rk*u(n,:,:) + 8.0_rk*u(2,:,:) - u(3,:,:) ) / (12.0_rk*dx)
    dudx(2,:,:) = ( u(n,:,:)   - 8.0_rk*u(1,:,:) + 8.0_rk*u(3,:,:) - u(4,:,:) ) / (12.0_rk*dx)

    forall ( i = 3:n-2 )
       dudx(i,:,:) = ( u(i-2,:,:) - 8.0_rk*u(i-1,:,:) + 8.0_rk*u(i+1,:,:) - u(i+2,:,:) ) / (12.0_rk*dx)
    end forall

    dudx(n-1,:,:) = ( u(n-3,:,:) - 8.0_rk*u(n-2,:,:) + 8.0_rk*u(n,:,:) - u(1,:,:) ) / (12.0_rk*dx)
    dudx(n,:,:)   = ( u(n-2,:,:) - 8.0_rk*u(n-1,:,:) + 8.0_rk*u(1,:,:) - u(2,:,:) ) / (12.0_rk*dx)

end subroutine diffx_c_3D

!---------------------------------------------------------------------------------------------

subroutine  diffy_c_3D( Bs, g, dy, u, dudy)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs+2*g, Bs+2*g, Bs+2*g)

    integer                         :: i, n

    n = size(u,1)

    dudy(:,1,:) = ( u(:,n-1,:) - 8.0_rk*u(:,n,:) + 8.0_rk*u(:,2,:) - u(:,3,:) ) / (12.0_rk*dy)
    dudy(:,2,:) = ( u(:,n,:)   - 8.0_rk*u(:,1,:) + 8.0_rk*u(:,3,:) - u(:,4,:) ) / (12.0_rk*dy)

    forall ( i = 3:n-2 )
       dudy(:,i,:) = ( u(:,i-2,:) - 8.0_rk*u(:,i-1,:) + 8.0_rk*u(:,i+1,:) - u(:,i+2,:) ) / (12.0_rk*dy)
    end forall

    dudy(:,n-1,:) = ( u(:,n-3,:) - 8.0_rk*u(:,n-2,:) + 8.0_rk*u(:,n,:) - u(:,1,:) ) / (12.0_rk*dy)
    dudy(:,n,:)   = ( u(:,n-2,:) - 8.0_rk*u(:,n-1,:) + 8.0_rk*u(:,1,:) - u(:,2,:) ) / (12.0_rk*dy)

end subroutine diffy_c_3D

!---------------------------------------------------------------------------------------------

subroutine  diffz_c_3D( Bs, g, dz, u, dudz)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs+2*g, Bs+2*g, Bs+2*g)

    integer                         :: i, n

    n = size(u,1)

    dudz(:,:,1) = ( u(:,:,n-1) - 8.0_rk*u(:,:,n) + 8.0_rk*u(:,:,2) - u(:,:,3) ) / (12.0_rk*dz)
    dudz(:,:,2) = ( u(:,:,n)   - 8.0_rk*u(:,:,1) + 8.0_rk*u(:,:,3) - u(:,:,4) ) / (12.0_rk*dz)

    forall ( i = 3:n-2 )
       dudz(:,:,i) = ( u(:,:,i-2) - 8.0_rk*u(:,:,i-1) + 8.0_rk*u(:,:,i+1) - u(:,:,i+2) ) / (12.0_rk*dz)
    end forall

    dudz(:,:,n-1) = ( u(:,:,n-3) - 8.0_rk*u(:,:,n-2) + 8.0_rk*u(:,:,n) - u(:,:,1) ) / (12.0_rk*dz)
    dudz(:,:,n)   = ( u(:,:,n-2) - 8.0_rk*u(:,:,n-1) + 8.0_rk*u(:,:,1) - u(:,:,2) ) / (12.0_rk*dz)

end subroutine diffz_c_3D
