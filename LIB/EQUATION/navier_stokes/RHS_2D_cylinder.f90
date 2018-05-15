
!--------------------------------------------------------------------------------------------------------------------------------------------------------
!> \file
!> \brief Right hand side for 2D navier stokes equation in cylindrical coordinates
!>        -------------------------------------------------------------------------
!> The right hand side of navier stokes in the skew symmetric form is implemented as follows:
!>\f{eqnarray*}{
!!     \partial_t \sqrt{\rho} &=& -\frac{1}{2J\sqrt{\rho}} \nabla \cdot (\rho \vec{u})-\frac{1}{\sqrt{\rho}}\frac{1}{C_{\rm SP} } (\rho-\rho^{\rm SP}) \\
!!    \partial_t (\sqrt{\rho} u_\alpha) &=& -\frac{1}{2J \sqrt{\rho}}
!!                                          \left[
!!                                                       (u_\alpha \partial_\beta (\rho u_\beta)+
!!                                                        u_\beta \rho \partial_\beta u_\alpha)
!!                                            \right]
!!                                           -\frac{1}{J \sqrt{\rho}} \partial_\beta \tau_{\alpha\beta}
!!                                           -\frac{1}{\sqrt{\rho}} \partial_\alpha p
!!                                            -\frac{1}{\sqrt{\rho}} \frac{1}{C_{\rm SP} }(\rho u_\alpha-\rho^{\rm SP} u_\alpha^{\rm SP})
!!                                           -\frac{\chi}{2\sqrt{\rho}C_\eta} (\rho u_\alpha)\\
!!    \partial_t p &=& -\frac{\gamma}{J} \partial_\beta( u_\beta p) + (\gamma-1)(u_\alpha \partial_\alpha p)
!!                                      +\frac{\gamma -1}{J}
!!                                           \left[
!!                                                       \partial_\alpha(u_\beta \tau_{\alpha\beta}+\phi_\alpha)
!!                                                       - u_\alpha\partial_\beta \tau_{\alpha\beta}
!!                                            \right]   -\frac{\gamma-1}{C_{\rm SP} } (p-p^{\rm SP})
!!                                             -\frac{\chi}{C_\eta} (p -\rho R_s T)
!!\f}
!> \version 0.5
!! \date 08/12/16 - create \n
!!  \date 13/2/18 - include mask and sponge terms (commit 1cf9d2d53ea76e3fa52f887d593fad5826afec88)
!> \author msr
!--------------------------------------------------------------------------------------------------------------------------------------------------------










!>\brief main function of RHS_2D_cylinders
subroutine RHS_2D_cylinder( g, Bs, x0, dx, phi, rhs)
!---------------------------------------------------------------------------------------------
!
    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)                            :: g, Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(2), intent(in)                  :: x0, dx
    !> datafields
    real(kind=rk), intent(in)                            :: phi(:, :, :)
    !> rhs array
    real(kind=rk), intent(inout)                            :: rhs(:, :, :)

    ! adiabatic coefficient
    real(kind=rk)                                           :: gamma_
    ! specific gas constant
    real(kind=rk)                                           :: Rs
    ! isochoric heat capacity
    real(kind=rk)                                           :: Cv
    ! isobaric heat capacity
    real(kind=rk)                                           :: Cp
    ! prandtl number
    real(kind=rk)                                           :: Pr
    ! dx(2)namic viscosity
    real(kind=rk)                                           :: mu0
    ! dissipation switch
    logical                                                 :: dissipation

    ! variables
    real(kind=rk)                                           :: rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g), &
                                                               T(Bs+2*g, Bs+2*g), mu(Bs+2*g, Bs+2*g), mu_d(Bs+2*g, Bs+2*g), lambda(Bs+2*g, Bs+2*g), &
                                                               fric_p(Bs+2*g, Bs+2*g), fric_u(Bs+2*g, Bs+2*g), fric_v(Bs+2*g, Bs+2*g), &
                                                               fric_T1(Bs+2*g, Bs+2*g), fric_T2(Bs+2*g, Bs+2*g), &
                                                               tau11(Bs+2*g, Bs+2*g), tau22(Bs+2*g, Bs+2*g), tau33(Bs+2*g, Bs+2*g), &
                                                               tau12(Bs+2*g, Bs+2*g), tau13(Bs+2*g, Bs+2*g), tau23(Bs+2*g, Bs+2*g)
                                                                 ! derivatives
    real(kind=rk)                                           :: u_x(Bs+2*g, Bs+2*g), u_y(Bs+2*g, Bs+2*g), v_x(Bs+2*g, Bs+2*g), v_y(Bs+2*g, Bs+2*g), &
                                                               p_x(Bs+2*g, Bs+2*g), p_y(Bs+2*g, Bs+2*g), T_x(Bs+2*g, Bs+2*g), T_y(Bs+2*g, Bs+2*g),&
                                                               div_U(Bs+2*g, Bs+2*g)
    real(kind=rk)                                           :: r(Bs+2*g),y

    ! dummy field
    real(kind=rk)                                           :: dummy(Bs+2*g, Bs+2*g)
    integer(kind=ik)                                        :: iy
    integer(kind=ik)   ,save                                     :: k=0

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
    rho         = phi(:,:,1)**2
    u           = phi(:,:,2)/phi(:,:,1)
    v           = phi(:,:,3)/phi(:,:,1)
    p           = phi(:,:,4)


    ! rhs
    rhs         = 0.0_rk
    r=0
    do iy=g+1, Bs+g
      k=k+1
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r(iy) = abs(y-params_ns%R_max)
       write(*,*) k, y,r(iy)
    end do


    call abort(73563,"Error [RHS_2D_cylinder.f90]: implement me!")
!---------------------------------------------------------------------------------------------
! main bodx(2)

    ! derivatives
    call grad_zentral( Bs, g, dx(1), dx(2), u, u_x, u_y)
    call grad_zentral( Bs, g, dx(1), dx(2), v, v_x, v_y)

    call diff1x_zentral( Bs, g, dx(1), p, p_x)
    call diff1y_zentral( Bs, g, dx(2), p, p_y)

    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    call diff1x_zentral( Bs, g, dx(1), rho*u, dummy)
    rhs(:,:,1) = -dummy
    call diff1y_zentral( Bs, g, dx(2), rho*v, dummy)
    rhs(:,:,1) = rhs(:,:,1) - dummy
    rhs(:,:,1) = rhs(:,:,1) * 0.5_rk/phi(:,:,1)

    ! friction
    if (dissipation) then

        ! Compute mu
        T    = p / (rho*Rs) ! ideal gas
        mu   = mu0
        mu_d = 0.0_rk

        ! thermal conductivity
        lambda  = Cp * mu/Pr

        ! tau11
        tau11 = mu * 2.0_rk * u_x

        !> \todo why not just simply div_U= u_x + v_y ?
        call diff1x_zentral( Bs, g, dx(1), u, dummy)
        div_U = dummy
        call diff1y_zentral( Bs, g, dx(2), v, dummy)
        div_U = div_U + dummy

        tau11 = tau11 + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U

        ! tau22
        tau22 = mu * 2.0_rk * v_y
        tau22 = tau22 + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U

        ! tau33
        tau33 = 0.0_rk
        tau33 = tau33 + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U

        ! tau12
        tau12 = mu * ( v_x + u_y )

        ! tau13
        tau13 = 0.0_rk

        ! tau23
        tau23 = 0.0_rk

        ! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
        call diff1x_zentral( Bs, g, dx(1), tau11, dummy)
        fric_u = dummy
        call diff1y_zentral( Bs, g, dx(2), tau12, dummy)
        fric_u = fric_u + dummy

        fric_u = fric_u / phi(:,:,1)

        call diff1x_zentral( Bs, g, dx(1), tau12, dummy)
        fric_v = dummy
        call diff1y_zentral( Bs, g, dx(2), tau22, dummy)
        fric_v = fric_v + dummy

        fric_v = fric_v / phi(:,:,1)

        ! Friction terms for the energy equation
        ! Heat Flux
        call grad_zentral( Bs, g, dx(1), dx(2), T, T_x, T_y)

        fric_T1 = lambda * T_x
        fric_T2 = lambda * T_y

        ! All simple divergence terms for u_i*tau_ik and phi_k
        call diff1x_zentral( Bs, g, dx(1), ( u*tau11 + v*tau12 + fric_T1 ), dummy)
        fric_p = dummy
        call diff1y_zentral( Bs, g, dx(2), ( u*tau12 + v*tau22 + fric_T2 ), dummy)
        fric_p = fric_p + dummy

        ! u_i*dx(1)_k (tau_ik) terms
        call diff1x_zentral( Bs, g, dx(1), tau11, dummy)
        fric_p = fric_p - u*dummy
        call diff1y_zentral( Bs, g, dx(2), tau12, dummy)
        fric_p = fric_p - u*dummy

        call diff1x_zentral( Bs, g, dx(1), tau12, dummy)
        fric_p = fric_p - v*dummy
        call diff1y_zentral( Bs, g, dx(2), tau22, dummy)
        fric_p = fric_p - v*dummy

        fric_p = ( gamma_ - 1 ) * fric_p

    else

        fric_p = 0.0_rk
        fric_u = 0.0_rk
        fric_v = 0.0_rk

    end if

    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    call diff1x_zentral( Bs, g, dx(1), (u * p), dummy)
    rhs(:,:,4) = - dummy
    call diff1y_zentral( Bs, g, dx(2), (v * p), dummy)
    rhs(:,:,4) = rhs(:,:,4) - dummy

    rhs(:,:,4) = rhs(:,:,4) * gamma_

    rhs(:,:,4) = rhs(:,:,4) + (gamma_ - 1.0_rk) * (u*p_x + v*p_y)

    rhs(:,:,4) = rhs(:,:,4) + fric_p

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    call diff1x_zentral( Bs, g, dx(1), (u * rho * u), dummy)
    rhs(:,:,2) = - 0.5_rk * dummy
    call diff1y_zentral( Bs, g, dx(2), (v * rho * u), dummy)
    rhs(:,:,2) = rhs(:,:,2) - 0.5_rk * dummy

    rhs(:,:,2) = rhs(:,:,2) - 0.5_rk * rho * u * u_x
    rhs(:,:,2) = rhs(:,:,2) - 0.5_rk * rho * v * u_y

    rhs(:,:,2) = rhs(:,:,2) - p_x

    rhs(:,:,2) = rhs(:,:,2) / phi(:,:,1)

    rhs(:,:,2) = rhs(:,:,2) + fric_u

    ! RHS of  momentum equation for v
    call diff1x_zentral( Bs, g, dx(1), (u * rho * v), dummy)
    rhs(:,:,3) = - 0.5_rk * dummy
    call diff1y_zentral( Bs, g, dx(2), (v * rho * v), dummy)
    rhs(:,:,3) = rhs(:,:,3) - 0.5_rk * dummy

    rhs(:,:,3) = rhs(:,:,3) - 0.5_rk * rho * u * v_x
    rhs(:,:,3) = rhs(:,:,3) - 0.5_rk * rho * v * v_y

    rhs(:,:,3) = rhs(:,:,3) - p_y

    rhs(:,:,3) = rhs(:,:,3) / phi(:,:,1)

    rhs(:,:,3) = rhs(:,:,3) + fric_v


end subroutine RHS_2D_cylinder
