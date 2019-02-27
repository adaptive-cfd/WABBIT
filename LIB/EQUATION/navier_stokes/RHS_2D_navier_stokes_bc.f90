
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
!> \version 0.5
!! \date 08/12/16 - create \n
!!  \date 13/2/18 - include mask and sponge terms (commit 1cf9d2d53ea76e3fa52f887d593fad5826afec88)
!> \author msr
!--------------------------------------------------------------------------------------------------------------------------------------------------------

!>\brief main function of RHS_2D_navier_stokes
subroutine RHS_2D_navier_stokes_BC( g, Bs, x0, delta_x, phi, rhs, boundary_flag)
!---------------------------------------------------------------------------------------------
!

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)                            :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(2), intent(in)                 :: x0, delta_x
    !> datafields
    real(kind=rk), intent(in)                               :: phi(:, :, :)
    !> rhs array
    real(kind=rk), intent(inout)                            :: rhs(:, :, :)
    ! when implementing boundary conditions, it is necessary to now if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
    ! boundary_flag(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    integer(kind=2), intent(in)                             :: boundary_flag(3)

    ! adiabatic coefficient
    real(kind=rk)                                           :: gamma_
    ! specific gas constant
    real(kind=rk)                                           :: Rs_inv
    ! isochoric heat capacity
    real(kind=rk)                                           :: Cv
    ! isobaric heat capacity
    real(kind=rk)                                           :: Cp
    ! prandtl number
    real(kind=rk)                                           :: Pr
    ! dynamic viscosity
    real(kind=rk)                                           :: mu0, mu_d, mu, lambda
    ! dissipation switch
    logical                                                 :: dissipation
    ! spacing
    real(kind=rk)                                           :: dx, dy

    ! variables
    real(kind=rk) :: rho(Bs(1)+2*g, Bs(2)+2*g), u(Bs(1)+2*g, Bs(2)+2*g), v(Bs(1)+2*g, Bs(2)+2*g), p(Bs(1)+2*g, Bs(2)+2*g), T(Bs(1)+2*g, Bs(2)+2*g), &
                      tau11(Bs(1)+2*g, Bs(2)+2*g), tau22(Bs(1)+2*g, Bs(2)+2*g), tau33(Bs(1)+2*g, Bs(2)+2*g), tau12(Bs(1)+2*g, Bs(2)+2*g)
    ! dummy field
    real(kind=rk) :: dummy(Bs(1)+2*g, Bs(2)+2*g), dummy2(Bs(1)+2*g, Bs(2)+2*g), dummy3(Bs(1)+2*g, Bs(2)+2*g), dummy4(Bs(1)+2*g, Bs(2)+2*g)


    ! inverse sqrt(rho) field
    real(kind=rk)                                           :: phi1_inv(Bs(1)+2*g, Bs(2)+2*g)

    ! loop variables
    integer(kind=ik)                                        :: i, j,n_eqn
    ! inverse sqrt(rho) field
    real(kind=rk),save,allocatable   :: phi_prime(:,:,:),phi_ref(:,:,:),mask(:,:,:)
    logical ,save :: allocated_penal_fields=.false.

    ! optimization
    ! - do not use all ghost nodes (note: use only two ghost nodes to get correct second derivatives)
    ! - write loops explicitly,
    ! - use multiplication instead of division
    ! - access array in column-major order
    ! - reduce number of additionaly variables -> lead to direct calculation of rhs terms after derivation

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set physics parameters for readability
    gamma_      = params_ns%gamma_
    Rs_inv      = 1.0_rk/params_ns%Rs
    Cv          = params_ns%Cv
    Cp          = params_ns%Cp
    Pr          = params_ns%Pr
    mu0         = params_ns%mu0
    dissipation = params_ns%dissipation

    dummy=0
    dummy2=0
    dummy3=0
    dummy4=0
    ! primitive variables
    do j = 1, Bs(2)+2*g
        do i = 1, Bs(1)+2*g
            rho(i,j)       = phi(i,j,1) * phi(i,j,1)
            phi1_inv(i,j)  = 1.0_rk / phi(i,j,1)
            u(i,j)         = phi(i,j,2) * phi1_inv(i,j)
            v(i,j)         = phi(i,j,3) * phi1_inv(i,j)
            p(i,j)         = phi(i,j,4)
        end do
    end do

    ! Compute mu and T
    if (dissipation) then
        do j = 1, Bs(2)+2*g
            do i = 1, Bs(1)+2*g
                T(i,j) = p(i,j) * phi1_inv(i,j) * phi1_inv(i,j) * Rs_inv
            end do
        end do
        mu   = mu0
        mu_d = 0.0_rk
        ! thermal conductivity
        lambda= Cp * mu/Pr
    end if

    ! discretization constant
    dx = delta_x(1)
    dy = delta_x(2)

!---------------------------------------------------------------------------------------------
! main body

    ! derivatives
    ! u_x, u_y
    !---------------------------------------------------------------------------------------------
    call diffxy( Bs, g, dx, dy, u, dummy, dummy2,boundary_flag)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,2) = - 0.5_rk * rho(i,j) * ( u(i,j) * dummy(i,j) + v(i,j) * dummy2(i,j))
        end do
    end do

    if (dissipation) then
        ! u_x
        tau11 = ( mu * 2.0_rk +  mu_d - 2.0_rk/3.0_rk * mu ) * dummy
        tau22 = ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy
        tau33 = ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy
        ! u_y
        tau12 = mu * dummy2
    end if

    ! v_x, v_y
    !---------------------------------------------------------------------------------------------
    call diffxy( Bs, g, dx, dy, v, dummy, dummy2,boundary_flag)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,3) = - 0.5_rk * rho(i,j) * ( u(i,j) * dummy(i,j) + v(i,j) * dummy2(i,j))
        end do
    end do

    if (dissipation) then
        ! v_x
        tau12 = tau12 + mu * dummy
        ! v_y
        tau11 = tau11 + ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy2
        tau22 = tau22 + ( mu * 2.0_rk + mu_d - 2.0_rk/3.0_rk * mu ) * dummy2
        tau33 = tau33 + ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy2
    end if

    ! p_x, p_y
    !---------------------------------------------------------------------------------------------
    call diffxy( Bs, g, dx, dy, p, dummy, dummy2,boundary_flag)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,2) = rhs(i,j,2) - dummy(i,j)
            rhs(i,j,3) = rhs(i,j,3) - dummy2(i,j)
            rhs(i,j,4) = (gamma_ - 1.0_rk) * ( u(i,j) * dummy(i,j) + v(i,j) * dummy2(i,j) )
        end do
    end do

    ! friction
    if (dissipation) then

        ! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
        ! tau11_x
        !---------------------------------------------------------------------------------------------
        call diffx( Bs, g, dx, tau11, dummy,boundary_flag(1))

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,2) = rhs(i,j,2) + dummy(i,j)
                rhs(i,j,4) = rhs(i,j,4) - ( gamma_ - 1.0_rk ) * u(i,j) * dummy(i,j)
            end do
        end do

        ! tau12_y
        !---------------------------------------------------------------------------------------------
        call diffy( Bs, g, dy, tau12, dummy, boundary_flag(2))

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,2) = rhs(i,j,2) + dummy(i,j)
                rhs(i,j,4) = rhs(i,j,4) - ( gamma_ - 1.0_rk ) * u(i,j) * dummy(i,j)
            end do
        end do

        ! tau12_x
        !---------------------------------------------------------------------------------------------
        call diffx( Bs, g, dx, tau12, dummy, boundary_flag(1))

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,3) = rhs(i,j,3) + dummy(i,j)
                rhs(i,j,4) = rhs(i,j,4) - ( gamma_ - 1.0_rk ) * v(i,j) * dummy(i,j)
            end do
        end do

        ! tau22_y
        !---------------------------------------------------------------------------------------------
        call diffy( Bs, g, dy, tau22, dummy, boundary_flag(2))

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,3) = rhs(i,j,3) + dummy(i,j)
                rhs(i,j,4) = rhs(i,j,4) - ( gamma_ - 1.0_rk ) * v(i,j) * dummy(i,j)
            end do
        end do

        ! Friction terms for the energy equation
        ! Heat Flux
        call diffxy( Bs, g, dx, dy, T, dummy, dummy2, boundary_flag)

        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy3(i,j)  = u(i,j)*tau11(i,j) + v(i,j)*tau12(i,j) + lambda * dummy(i,j)
                dummy4(i,j) = u(i,j)*tau12(i,j) + v(i,j)*tau22(i,j) + lambda * dummy2(i,j)
            end do
        end do
        call diffx( Bs, g, dx, dummy3, dummy, boundary_flag(1))
        call diffy( Bs, g, dy, dummy4, dummy2, boundary_flag(2))

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,4) = rhs(i,j,4) + ( gamma_ - 1.0_rk ) * ( dummy(i,j) + dummy2(i,j) )
            end do
        end do

    end if

    ! EQUATIONS
    ! --------------------------------------------------------------------------------------------------------------
    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = rho(i,j)*u(i,j)
            dummy2(i,j) = rho(i,j)*v(i,j)
        end do
    end do
    call diffx( Bs, g, dx, dummy,  dummy3, boundary_flag(1))
    call diffy( Bs, g, dy, dummy2, dummy4, boundary_flag(2))

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1) = (-dummy3(i,j) - dummy4(i,j)) * 0.5_rk * phi1_inv(i,j)
        end do
    end do

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = u(i,j)*rho(i,j)*u(i,j)
            dummy2(i,j) = v(i,j)*rho(i,j)*u(i,j)
        end do
    end do
    call diffx( Bs, g, dx, dummy,  dummy3, boundary_flag(1))
    call diffy( Bs, g, dy, dummy2, dummy4, boundary_flag(2))

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,2) = ( rhs(i,j,2) - 0.5_rk * ( dummy3(i,j) + dummy4(i,j) ) ) * phi1_inv(i,j)
        end do
    end do

    ! RHS of  momentum equation for v
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = u(i,j)*rho(i,j)*v(i,j)
            dummy2(i,j) = v(i,j)*rho(i,j)*v(i,j)
        end do
    end do
    call diffx( Bs, g, dx, dummy,  dummy3, boundary_flag(1))
    call diffy( Bs, g, dy, dummy2, dummy4, boundary_flag(2))

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,3) = ( rhs(i,j,3) - 0.5_rk * ( dummy3(i,j) + dummy4(i,j) ) ) * phi1_inv(i,j)
        end do
    end do

    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = u(i,j)*p(i,j)
            dummy2(i,j) = v(i,j)*p(i,j)
        end do
    end do
    call diffx( Bs, g, dx, dummy,  dummy3, boundary_flag(1))
    call diffy( Bs, g, dy, dummy2, dummy4, boundary_flag(2))

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,4) = rhs(i,j,4) - gamma_ * ( dummy3(i,j) + dummy4(i,j) )
        end do
    end do

    if (params_ns%penalization) then
        ! add volume penalization
        if (.not. allocated_penal_fields) then
          allocated_penal_fields=.true.
          n_eqn=params_ns%n_eqn
          allocate( mask(Bs(1)+2*g,Bs(2)+2*g,n_eqn), &
                    phi_prime(Bs(1)+2*g,Bs(2)+2*g,n_eqn),&
                    phi_ref(Bs(1)+2*g,Bs(2)+2*g,n_eqn))
        endif
        phi_prime(:,:,rhoF)= rho
        phi_prime(:,:,UxF )= u
        phi_prime(:,:,UyF )= v
        phi_prime(:,:,pF  )= p

        call compute_mask_and_ref2D(params_ns, Bs, g, x0, delta_x, phi_prime, mask, phi_ref)
        do j = g+1, Bs(2)+g
          do i = g+1, Bs(1)+g
            ! density
            rhs(i,j,rhoF)=rhs(i,j,rhoF) -0.5_rk*phi1_inv(i,j)*mask(i,j,rhoF)*(rho(i,j)-  Phi_ref(i,j,rhoF) )
            ! x-velocity
            rhs(i,j,UxF)=rhs(i,j,UxF) -1.0_rk*phi1_inv(i,j)*mask(i,j,UxF)*(rho(i,j)*u(i,j)-  Phi_ref(i,j,UxF) )
            ! y-velocity
            rhs(i,j,UyF)=rhs(i,j,UyF) -1.0_rk*phi1_inv(i,j)*mask(i,j,UyF)*(rho(i,j)*v(i,j)-  Phi_ref(i,j,UyF) )
            ! preasure
            rhs(i,j,pF)=rhs(i,j,pF)                        -mask(i,j,pF)*(p(i,j)- Phi_ref(i,j,pF) )
          end do
        end do
    endif

end subroutine RHS_2D_navier_stokes_BC

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

subroutine  diffxy(Bs, g, dx, dy, u, dudx, dudy,boundary_flag)
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    integer(kind=2), intent(in)     :: boundary_flag(3)
    real(kind=rk), intent(in)       :: dx, dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g), dudy(Bs(1)+2*g, Bs(2)+2*g)

    integer                         :: i, j
    real(kind=rk)                   :: dx_inv, dy_inv

    ! - do not use all ghost nodes (note: use only two ghost nodes to get correct second derivatives)
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    call DUS_Dx(dudx,dx,u,boundary_flag(1))
    call DUS_Dy(dudy,dy,u,boundary_flag(2))

    ! dudx=0
    ! dudy=0
    ! dx_inv = 1.0_rk/(12.0_rk*dx)
    ! dy_inv = 1.0_rk/(12.0_rk*dy)
    !
    ! do j = g-1, Bs+g+2
    !     do i = g-1, Bs+g+2
    !         dudx(i,j) = ( u(i-2,j) - 8.0_rk*u(i-1,j) + 8.0_rk*u(i+1,j) - u(i+2,j) ) * dx_inv
    !         dudy(i,j) = ( u(i,j-2) - 8.0_rk*u(i,j-1) + 8.0_rk*u(i,j+1) - u(i,j+2) ) * dy_inv
    !     end do
    ! end do

end subroutine diffxy

subroutine  diffx( Bs, g, dx, u, dudx,boundary_flag)
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    integer(kind=2), intent(in)     :: boundary_flag
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g)

    integer                         :: i, j
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order
    call DUS_Dx(dudx,dx,u,boundary_flag)
    ! dx_inv = 1.0_rk/(12.0_rk*dx)
    !
    ! do j = g+1, Bs+g
    !     do i = g+1, Bs+g
    !         dudx(i,j) = ( u(i-2,j) - 8.0_rk*u(i-1,j) + 8.0_rk*u(i+1,j) - u(i+2,j) ) * dx_inv
    !     end do
    ! end do
    !

end subroutine diffx

subroutine  diffy( Bs, g, dy, u, dudy, boundary_flag)

    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    integer(kind=2), intent(in)     :: boundary_flag
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g)

    integer                         :: i, j
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order
    call DUS_Dy(dudy,dy,u,boundary_flag)
    ! dudy=0
    ! dy_inv = 1.0_rk/(12.0_rk*dy)
    !
    ! do j = g+1, Bs+g
    !     do i = g+1, Bs+g
    !         dudy(i,j) = ( u(i,j-2) - 8.0_rk*u(i,j-1) + 8.0_rk*u(i,j+1) - u(i,j+2) ) * dy_inv
    !     end do
    ! end do

end subroutine diffy
