
!--------------------------------------------------------------------------------------------------------------------------------------------------------
!> \file
!> \brief Right hand side for 2D navier stokes equation
!>        ---------------------------------------------
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
!! \details
!! With the head flux \f$\phi=\lambda T(x,y,t)\f$ and \f$\lambda \in \mathrm{R}\f$
!! The friction term is further expanded:
!>\f{eqnarray*}{
!!
!!    \partial_\alpha(u_\beta \tau_{\alpha\beta}+\phi_\alpha) - u_\alpha\partial_\beta \tau_{\alpha\beta}
!!    =
!!\f}
!> \version 0.5
!> \date 05/03/18 - create (commit beb3fa44e) \n
!> \author Pkrah
!--------------------------------------------------------------------------------------------------------------------------------------------------------










!>\brief main function of RHS_2D_navier_stokes
subroutine rhs_ns_2D( g, Bs, x0, delta_x, phi, rhs)
!---------------------------------------------------------------------------------------------
!
    implicit none
    !> grid parameter
    integer(kind=ik), intent(in)                            :: g, Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(2), intent(in)                  :: x0, delta_x
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
    ! dynamic viscosity
    real(kind=rk)                                           :: dx,dy,dx2_12_inv,dy2_12_inv,sqrt_rho_inv
    ! dissipation switch
    logical                                                 :: dissipation

    ! variables
    real(kind=rk)                                           :: rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g),vRho(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g), &
                                                               uRho(Bs+2*g, Bs+2*g),uRhou(Bs+2*g, Bs+2*g),uRhov(Bs+2*g, Bs+2*g),vRhov(Bs+2*g, Bs+2*g), &
                                                               T(Bs+2*g, Bs+2*g),tau(1:3, 1:3),&
                                                               mask(Bs+2*g, Bs+2*g), sponge(Bs+2*g, Bs+2*g),up(Bs+2*g, Bs+2*g), vp(Bs+2*g, Bs+2*g)
    ! derivatives
    real(kind=rk)                                           :: u_x, u_xx, u_xy, u_y, u_yy, v_x, v_xx, v_xy, v_y, v_yy, p_x, p_y, T_xx, T_yy, &
                                                               div_U, uRho_x, uRhou_x, vRhov_y, vRhou_y, vRho_y,uRhov_x, vp_y, up_x
    !friction parameters
    real(kind=rk)                                           ::  mu_d, mu, mu0, lambda

    ! friction_terms
    real(kind=rk)                                           ::  fric_p, fric_u, fric_v

    !loop variables
    integer(kind=ik)                                        :: ix,iy
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
    rho         = phi(:,:,1)**2
    u           = phi(:,:,2)/phi(:,:,1)
    v           = phi(:,:,3)/phi(:,:,1)
    p           = phi(:,:,4)
    uRho        = u*rho
    uRhou       = u*uRho
    uRhov       = v*uRho
    vRho        = v*rho
    vRhov       = v*vRho
    uP          = u*p
    vP          = v*p
    ! discretization constant
    dx=delta_x(1)
    dy=delta_x(2)
    dx2_12_inv=1.0_rk/(12.0_rk*dx**2)
    dy2_12_inv=1.0_rk/(12.0_rk*dy**2)

    ! rhs
    rhs         = 0.0_rk

    if (dissipation) then
        mu   = mu0
        mu_d = 0.0_rk

        ! thermal conductivity
        lambda  = Cp * mu/Pr
        T    = p / (rho*Rs) ! ideal gas
    end if

    if (params_ns%sponge_layer) then
        call get_sponge(sponge, x0, delta_x, Bs, g)
        ! sponge coefficient
        sponge = 1.0_rk/params_ns%C_sp * sponge
    end if

    if (params_ns%penalization) then
        call get_mask(mask, x0, delta_x, Bs, g )
        ! mask coefficient
        mask = 1.0_rk/params_ns%C_eta*mask
    endif
!---------------------------------------------------------------------------------------------
! main body

       !-----------------------------------------------------------------------
        do ix = g+1, Bs+g
            do iy = g+1, Bs+g

                ! Derivatives
                !------------

                ! velocity (u,v)    u_x=du/dx u_xx=d^2u/dx^2
                u_x    = df_dx_central(u,dx,ix,iy)
                u_xx  = ( -u(ix+2,iy) + 16.0_rk*u(ix+1,iy) - 30.0_rk*u(ix,iy) + 16.0_rk*u(ix-1,iy) - u(ix-2,iy) ) * (dx2_12_inv)
                u_y    = df_dy_central(u,dy,ix,iy)
                u_yy  = ( -u(ix,iy+2) + 16.0_rk*u(ix,iy+1) - 30.0_rk*u(ix,iy) + 16.0_rk*u(ix,iy-1) - u(ix,iy-2) ) * (dy2_12_inv)
                u_xy  = df_dxdy_central(u,dx,dy,ix,iy)

                v_x    = df_dx_central(v,dx,ix,iy)
                v_xx  = ( -v(ix+2,iy) + 16.0_rk*v(ix+1,iy) - 30.0_rk*v(ix,iy) + 16.0_rk*v(ix-1,iy) - v(ix-2,iy) ) * (dx2_12_inv)
                v_y    = df_dy_central(v,dy,ix,iy)
                v_yy  = ( -v(ix,iy+2) + 16.0_rk*v(ix,iy+1) - 30.0_rk*v(ix,iy) + 16.0_rk*v(ix,iy-1) - v(ix,iy-2) ) * (dy2_12_inv)
                v_xy  = df_dxdy_central(v,dx,dy,ix,iy)

                !preasure
                p_x    = df_dx_central(p,dx,ix,iy)
                p_y    = df_dy_central(p,dy,ix,iy)

                ! massflux
                uRho_x = df_dx_central(uRho,dx,ix,iy)
                vRho_y = df_dy_central(vRho,dy,ix,iy)

                uRhou_x = df_dx_central(uRhou,dx,ix,iy)
                vRhou_y = df_dy_central(uRhov,dy,ix,iy)

                uRhov_x = df_dx_central(uRhov,dx,ix,iy)
                vRhov_y = df_dy_central(vRhov,dy,ix,iy)

                uP_x    = df_dx_central(uP,dx,ix,iy)
                vP_y    = df_dy_central(vP,dy,ix,iy)

                div_U = u_x + v_y
                !------------------------------------------------------------------------------------------------------------------


                sqrt_rho_inv = 1/phi(ix,iy,1)


                ! Friction
                !---------
                ! this if should go outside the 4 loop
                if (dissipation) then

                    ! friction tensor
                    tau(1,1)   = mu * 2.0_rk * u_x
                    tau(1,1)   = tau(1,1) + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U
                    tau(2,2)   = mu * 2.0_rk * v_y
                    tau(2,2)   = tau(2,2) + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U
                    tau(3,3)   = 0.0_rk
                    tau(3,3)   = tau(3,3) + ( mu_d - 2.0_rk/3.0_rk * mu ) * div_U
                    tau(1,2)   = mu * ( v_x + u_y )
                    tau(1,3)   = 0.0_rk
                    tau(2,3)   = 0.0_rk

                    ! laplace T
                    T_xx  = ( -T(ix+2,iy) + 16.0_rk*T(ix+1,iy) - 30.0_rk*T(ix,iy) + 16.0_rk*T(ix-1,iy) - T(ix-2,iy) ) * (dx2_12_inv)
                    T_yy  = ( -T(ix,iy+2) + 16.0_rk*T(ix,iy+1) - 30.0_rk*T(ix,iy) + 16.0_rk*T(ix,iy-1) - T(ix,iy-2) ) * (dy2_12_inv)

                    ! friction term of velocity
                    fric_u=mu*u_yy + (1.0_rk/3.0_rk*mu+mu_d)*v_xy+(4.0_rk/3.0_rk*mu+mu_d)*u_xx
                    fric_u=fric_u/sqrt_rho_inv
                    fric_v=mu*v_xx + (1.0_rk/3.0_rk*mu+mu_d)*u_xy+(4.0_rk/3.0_rk*mu+mu_d)*v_yy
                    fric_u=fric_v/sqrt_rho_inv

                    ! friction term of preasure
                    fric_p= (gamma_-1)*(u_x * tau(1,1) + v_x * tau(1,2) + u_y * tau(1,2) + v_y * tau(2,2) + lambda* (T_xx + T_yy))
                else

                    fric_p = 0.0_rk
                    fric_u = 0.0_rk
                    fric_v = 0.0_rk

                end if


                ! Governing Equations
                !--------------------

                ! conservation of mass
                rhs(ix,iy,1) = -(uRho_x + vRho_y)*0.5_rk*sqrt_rho_inv
                ! conservation of momentum
                rhs(ix,iy,2) = -sqrt_rho_inv*(0.5_rk*                                  &
                                ( rho(ix,iy)*u(ix,iy)*u_x + rho(ix,iy)*v(ix,iy)*u_y +  &
                                  uRhou_x                 + vRhou_y                  ) &
                                + p_x) + fric_u

                rhs(ix,iy,3) = -sqrt_rho_inv*(0.5_rk*                                  &
                                ( rho(ix,iy)*u(ix,iy)*v_x + rho(ix,iy)*v(ix,iy)*v_y +  &
                                  uRhov_x                 + vRhov_y                  ) &
                                + p_y) + fric_v
                ! conservation of energy
                rhs(ix,iy,4) = -gamma_ * ( uP_x + vP_y ) + (gamma_-1)*( u(ix,iy)*p_x + v(ix,iy)*p_y ) + fric_p



                ! Add Sponge and penalization terms
                !----------------------------------

               if (params_ns%sponge_layer) then
                   ! add sponge
                   ! rho component (density)
                   rhs(ix,iy,1)=rhs(ix,iy,1) - 0.5_rk* sqrt_rho_inv*sponge(ix,iy) * ( rho(ix,iy)  -   rho0_)
                   ! sqrt(rho)u component (momentum)
                   rhs(ix,iy,2)=rhs(ix,iy,2) -         sqrt_rho_inv*sponge(ix,iy) * ( uRho(ix,iy) - rho0_*u0_)
                   ! sqrt(rho)v component (momentum)
                   rhs(ix,iy,3)=rhs(ix,iy,3) -         sqrt_rho_inv*sponge(ix,iy) * ( vRho(ix,iy) - 0.0_rk   )
                   ! p component (preasure/energy)
                   rhs(ix,iy,4)=rhs(ix,iy,4) - (gamma_-1)          *sponge(ix,iy) * ( p(ix,iy)    - p0_ )
               endif
               if (params_ns%penalization) then
                   ! add mask
                   ! sqrt(rho)u component (momentum)
                   rhs(ix,iy,2)=rhs(ix,iy,2) - sqrt_rho_inv *mask(ix,iy) * ( uRho(ix,iy) - 0.0_rk )
                   ! sqrt(rho)v component (momentum)
                   rhs(ix,iy,3)=rhs(ix,iy,3) - sqrt_rho_inv *mask(ix,iy) * ( vRho(ix,iy) - 0.0_rk )
                   ! p component (preasure/energy)
                   rhs(ix,iy,4)=rhs(ix,iy,4) -               mask(ix,iy) * ( p(ix,iy)    - rho(ix,iy)*Rs*T0_ )
                   
               endif


            end do
        end do




end subroutine rhs_ns_2d


!---------------------------------------------------------------------------------------------
pure function  df_dx_central( f,dx,ix,iy)

    integer(kind=ik), intent(in)    :: ix,iy
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: f(:,:)
    real(kind=rk)                   :: df_dx_central

    df_dx_central = ( f(ix-2,iy) - 8.0_rk*f(ix-1,iy) + 8.0_rk*f(ix+1,iy) - f(ix+2,iy) ) / (12.0_rk*dx)

end function df_dx_central


pure function  df_dy_central(f,dy,ix,iy)

    integer(kind=ik), intent(in)       :: ix,iy
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: f(:,:)
    real(kind=rk)                   :: df_dy_central

    df_dy_central = ( f(ix,iy-2) - 8.0_rk*f(ix,iy-1) + 8.0_rk*f(ix,iy+1) - f(ix,iy+2) ) / (12.0_rk*dy)

end function df_dy_central

pure function  df_dxdy_central(f,dx,dy,ix,iy)

    integer(kind=ik), intent(in)    :: ix,iy
    real(kind=rk), intent(in)       :: dx,dy
    real(kind=rk), intent(in)       :: f(:,:)
    real(kind=rk)                   :: df_dxdy_central

    df_dxdy_central      = f(ix-2,iy-2) - 8.0_rk    * f(ix-1,iy-2) + 8.0_rk   * f(ix+1,iy-2) - f(ix+2,iy-2) &
                         - f(ix-2,iy-1) + 64.0_rk   * f(ix-1,iy-1) - 64.0_rk  * f(ix+1,iy-1) + f(ix+2,iy-1) &
                         + f(ix-2,iy+1) - 64.0_rk   * f(ix-1,iy+1) + 64.0_rk  * f(ix+1,iy+1) - f(ix+2,iy+1) &
                         - f(ix-2,iy+2) + 8.0_rk    * f(ix-1,iy+2) - 8.0_rk   * f(ix+1,iy+2) + f(ix+2,iy+2)
    df_dxdy_central      = df_dxdy_central/(144.0_rk*dx*dy)

end function df_dxdy_central




!--------------------------------------------------------------------------!
!               phils sponge stuff
!--------------------------------------------------------------------------!


!==========================================================================
!> \brief This function computes a 2d sponge term
!!
!! \details The sponge term is
!!   \f{eqnarray*}{
!!           s_q(x,y)=\frac{\chi_{\mathrm sp}(x,y)}{C_{\mathrm sp}}(q(x,y)-q_{\mathrm ref})
!!                          \quad \forall x,y \in \Omega_{\mathrm block}
!!     \f}
!! Where we addopt the notation from <a href="https://arxiv.org/abs/1506.06513">Thomas Engels (2015)</a>
!! - the mask function is \f$\chi_{\rm sp}\f$
!! - \f$C_{\rm sp}\f$ is the sponge coefficient (normaly \f$10^{-1}\f$)

subroutine get_sponge(sponge, x0, dx, Bs, g)

    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> sponge term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(out)     :: sponge
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, ddx
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset sponge array
    sponge = 0.0_rk
!---------------------------------------------------------------------------------------------
! main body
    ddx = 0.1_rk*params_ns%Lx

    do iy=1, Bs+2*g
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)
           if ((params_ns%Lx-x) <= ddx) then
               sponge(ix,iy) = (x-(params_ns%Lx-ddx))**2
           elseif (x <= ddx) then
               sponge(ix,iy) = (x-ddx)**2
           else
               sponge(ix,iy) = 0.0_rk
           end if
       end do
    end do

end subroutine get_sponge
