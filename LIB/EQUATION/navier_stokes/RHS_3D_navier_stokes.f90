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

subroutine RHS_3D_navier_stokes(g, Bs, x0, delta_x, phi, rhs)!, boundary_flag)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> grid parameter
    integer(kind=ik), intent(in)                            :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> rhs parameter
    real(kind=rk), dimension(3), intent(in)                 :: x0,delta_x
    !> datafields
    real(kind=rk), intent(in)                              :: phi(:, :, :, :)
    ! rhs array
    real(kind=rk),intent(inout)                               :: rhs(:, :, :,:)
    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
    ! boundary_flag(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    !integer(kind=2), intent(in)                             :: boundary_flag(3)

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
    real(kind=rk)                                           :: mu0, mu_d, mu, lambda

    ! spacing
    real(kind=rk)                                           :: dx, dy, dz

    ! dissipation switch
    logical                                                 :: dissipation

    ! variables
    real(kind=rk)                                           :: rho(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), v(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), w(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               p(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), T(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               tau11(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau22(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau33(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               tau12(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau13(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau23(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    ! dummy field
    real(kind=rk)                                           :: dummy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy2(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy3(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               dummy4(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy5(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy6(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    ! inverse sqrt(rho) field
    real(kind=rk)                                           :: phi1_inv(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    ! loop variables
    integer(kind=ik)                                        :: i, j, k,n_eqn
    ! penalization fields
    real(kind=rk), allocatable,save   :: phi_prime(:,:,:,:), phi_ref(:,:,:,:), mask(:,:,:,:)
    logical ,save :: allocated_penal_fields=.false.
!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set physics parameters for readability
    gamma_      = params_ns%gamma_
    ! only use 1/Rs
    Rs          = 1.0_rk/params_ns%Rs
    Cv          = params_ns%Cv
    Cp          = params_ns%Cp
    Pr          = params_ns%Pr
    mu0         = params_ns%mu0
    dissipation = params_ns%dissipation

    ! variables
    do k = 1, Bs(3)+2*g
       do j = 1, Bs(2)+2*g
            do i = 1, Bs(1)+2*g
                rho(i,j,k)       = phi(i,j,k,1) * phi(i,j,k,1)
                phi1_inv(i,j,k)  = 1.0_rk / phi(i,j,k,1)
                u(i,j,k)         = phi(i,j,k,2) * phi1_inv(i,j,k)
                v(i,j,k)         = phi(i,j,k,3) * phi1_inv(i,j,k)
                w(i,j,k)         = phi(i,j,k,4) * phi1_inv(i,j,k)
                p(i,j,k)         = phi(i,j,k,5)
            end do
        end do
    end do

    ! Compute mu and T
    if (dissipation) then
        do k = 1, Bs(3)+2*g
            do j = 1, Bs(2)+2*g
                do i = 1, Bs(1)+2*g
                    T(i,j,k) = p(i,j,k) * phi1_inv(i,j,k) * phi1_inv(i,j,k) * Rs
                end do
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
    dz = delta_x(3)

!---------------------------------------------------------------------------------------------
! main body

    ! derivatives
    ! u_x, u_y, u_z
    !---------------------------------------------------------------------------------------------
    call diffxyz_c_3D_opt( Bs, g, dx, dy, dz, u, dummy, dummy2, dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,2) = - 0.5_rk * rho(i,j,k) * ( u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k))
            end do
        end do
    end do

    if (dissipation) then
        ! u_x
        tau11 = ( mu * 2.0_rk +  mu_d - 2.0_rk/3.0_rk * mu ) * dummy
        tau22 = ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy
        tau33 = ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy
        ! u_y
        tau12 = mu * dummy2
        ! u_z
        tau13 = mu * dummy3
    end if

    ! v_x, v_y, v_z
    !---------------------------------------------------------------------------------------------
    call diffxyz_c_3D_opt( Bs, g, dx, dy, dz, v, dummy, dummy2, dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,3) = - 0.5_rk * rho(i,j,k) * ( u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k))
            end do
        end do
    end do

    if (dissipation) then
        ! v_x
        tau12 = tau12 + mu * dummy
        ! v_y
        tau11 = tau11 + ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy2
        tau22 = tau22 + ( mu * 2.0_rk + mu_d - 2.0_rk/3.0_rk * mu ) * dummy2
        tau33 = tau33 + ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy2
        ! v_z
        tau23 = mu * dummy3
    end if

    ! w_x, w_y, w_z
    !---------------------------------------------------------------------------------------------
    call diffxyz_c_3D_opt( Bs, g, dx, dy, dz, w, dummy, dummy2, dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,4) = - 0.5_rk * rho(i,j,k) * ( u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k))
            end do
        end do
    end do

    if (dissipation) then
        ! w_x
        tau13 = tau13 + mu * dummy
        ! w_z
        tau11 = tau11 + ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy3
        tau22 = tau22 + ( mu_d - 2.0_rk/3.0_rk * mu ) * dummy3
        tau33 = tau33 + ( mu * 2.0_rk + mu_d - 2.0_rk/3.0_rk * mu ) * dummy3
        ! w_y
        tau23 = tau23 + mu * dummy2
    end if

    ! p_x, p_y, p_z
    !---------------------------------------------------------------------------------------------
    call diffxyz_c_3D_opt( Bs, g, dx, dy, dz, p, dummy, dummy2, dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,2) = rhs(i,j,k,2) - dummy(i,j,k)
                rhs(i,j,k,3) = rhs(i,j,k,3) - dummy2(i,j,k)
                rhs(i,j,k,4) = rhs(i,j,k,4) - dummy3(i,j,k)
                rhs(i,j,k,5) = (gamma_ - 1.0_rk) * ( u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k) )
            end do
        end do
    end do

    ! friction
    if (dissipation) then

        ! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
        ! tau11_x
        !---------------------------------------------------------------------------------------------
        call diffx_c_3D_opt( Bs, g, dx, tau11, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,2) = rhs(i,j,k,2) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * u(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau12_x
        !---------------------------------------------------------------------------------------------
        call diffx_c_3D_opt( Bs, g, dx, tau12, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,3) = rhs(i,j,k,3) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * v(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau12_y
        !---------------------------------------------------------------------------------------------
        call diffy_c_3D_opt( Bs, g, dy, tau12, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,2) = rhs(i,j,k,2) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * u(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau13_x
        !---------------------------------------------------------------------------------------------
        call diffx_c_3D_opt( Bs, g, dx, tau13, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,4) = rhs(i,j,k,4) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * w(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau13_z
        !---------------------------------------------------------------------------------------------
        call diffz_c_3D_opt( Bs, g, dz, tau13, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,2) = rhs(i,j,k,2) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * u(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau22_y
        !---------------------------------------------------------------------------------------------
        call diffy_c_3D_opt( Bs, g, dy, tau22, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,3) = rhs(i,j,k,3) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * v(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau23_y
        !---------------------------------------------------------------------------------------------
        call diffy_c_3D_opt( Bs, g, dy, tau23, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,4) = rhs(i,j,k,4) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * w(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau23_z
        !---------------------------------------------------------------------------------------------
        call diffz_c_3D_opt( Bs, g, dz, tau23, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,3) = rhs(i,j,k,3) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * v(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau33_z
        !---------------------------------------------------------------------------------------------
        call diffz_c_3D_opt( Bs, g, dz, tau33, dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,4) = rhs(i,j,k,4) + dummy(i,j,k)
                    rhs(i,j,k,5) = rhs(i,j,k,5) - ( gamma_ - 1.0_rk ) * w(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! Friction terms for the energy equation
        ! Heat Flux
        call diffxyz_c_3D_opt( Bs, g, dx, dy, dz, T, dummy, dummy2, dummy3)

        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2
                    dummy4(i,j,k) = u(i,j,k)*tau11(i,j,k) + v(i,j,k)*tau12(i,j,k) + w(i,j,k)*tau13(i,j,k) + lambda * dummy(i,j,k)
                    dummy5(i,j,k) = u(i,j,k)*tau12(i,j,k) + v(i,j,k)*tau22(i,j,k) + w(i,j,k)*tau23(i,j,k) + lambda * dummy2(i,j,k)
                    dummy6(i,j,k) = u(i,j,k)*tau13(i,j,k) + v(i,j,k)*tau23(i,j,k) + w(i,j,k)*tau33(i,j,k) + lambda * dummy3(i,j,k)
                end do
            end do
        end do
        call diffx_c_3D_opt( Bs, g, dx, dummy4, dummy)
        call diffy_c_3D_opt( Bs, g, dy, dummy5, dummy2)
        call diffz_c_3D_opt( Bs, g, dz, dummy6, dummy3)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,5) = rhs(i,j,k,5) + ( gamma_ - 1.0_rk ) * ( dummy(i,j,k) + dummy2(i,j,k) + dummy3(i,j,k) )
                end do
            end do
        end do

    end if

    ! EQUATIONS
    ! --------------------------------------------------------------------------------------------------------------
    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = rho(i,j,k)*u(i,j,k)
                dummy2(i,j,k) = rho(i,j,k)*v(i,j,k)
                dummy3(i,j,k) = rho(i,j,k)*w(i,j,k)
            end do
        end do
    end do
    call diffx_c_3D_opt( Bs, g, dx, dummy,  dummy4)
    call diffy_c_3D_opt( Bs, g, dy, dummy2, dummy5)
    call diffz_c_3D_opt( Bs, g, dz, dummy3, dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,1) = (-dummy4(i,j,k) - dummy5(i,j,k) - dummy6(i,j,k)) * 0.5_rk * phi1_inv(i,j,k)
            end do
        end do
    end do

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*rho(i,j,k)*u(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*rho(i,j,k)*u(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*rho(i,j,k)*u(i,j,k)
            end do
        end do
    end do
    call diffx_c_3D_opt( Bs, g, dx, dummy,  dummy4)
    call diffy_c_3D_opt( Bs, g, dy, dummy2, dummy5)
    call diffz_c_3D_opt( Bs, g, dz, dummy3, dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,2) = ( rhs(i,j,k,2) - 0.5_rk * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) ) ) * phi1_inv(i,j,k)
            end do
        end do
    end do

    ! RHS of  momentum equation for v
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*rho(i,j,k)*v(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*rho(i,j,k)*v(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*rho(i,j,k)*v(i,j,k)
            end do
        end do
    end do
    call diffx_c_3D_opt( Bs, g, dx, dummy,  dummy4)
    call diffy_c_3D_opt( Bs, g, dy, dummy2, dummy5)
    call diffz_c_3D_opt( Bs, g, dz, dummy3, dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,3) = ( rhs(i,j,k,3) - 0.5_rk * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k)) ) * phi1_inv(i,j,k)
            end do
        end do
    end do

    ! RHS of  momentum equation for w
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*rho(i,j,k)*w(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*rho(i,j,k)*w(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*rho(i,j,k)*w(i,j,k)
            end do
        end do
    end do
    call diffx_c_3D_opt( Bs, g, dx, dummy,  dummy4)
    call diffy_c_3D_opt( Bs, g, dy, dummy2, dummy5)
    call diffz_c_3D_opt( Bs, g, dz, dummy3, dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,4) = ( rhs(i,j,k,4) - 0.5_rk * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k)) ) * phi1_inv(i,j,k)
            end do
        end do
    end do

    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*p(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*p(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*p(i,j,k)
            end do
        end do
    end do
    call diffx_c_3D_opt( Bs, g, dx, dummy,  dummy4)
    call diffy_c_3D_opt( Bs, g, dy, dummy2, dummy5)
    call diffz_c_3D_opt( Bs, g, dz, dummy3, dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,5) = rhs(i,j,k,5) - gamma_ * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
            end do
        end do
    end do

    if (params_ns%penalization) then
        ! add volume penalization
        if (.not. allocated_penal_fields) then
          allocated_penal_fields=.true.
          n_eqn=params_ns%n_eqn
          allocate( mask(Bs(1)+2*g,Bs(2)+2*g,Bs(3)+2*g,n_eqn), &
                    phi_prime(Bs(1)+2*g,Bs(2)+2*g,Bs(3)+2*g,n_eqn),&
                    phi_ref(Bs(1)+2*g,Bs(2)+2*g,Bs(3)+2*g,n_eqn))
        endif
        phi_prime(:,:,:,rhoF)= rho
        phi_prime(:,:,:,UxF )= u
        phi_prime(:,:,:,UyF )= v
        phi_prime(:,:,:,UzF )= w
        phi_prime(:,:,:,pF  )= p

        call compute_mask_and_ref3D(params_ns, Bs, g, x0, delta_x, phi_prime, mask, phi_ref)
        do k = 1, Bs(3)+2*g
          do j = 1, Bs(2)+2*g
            do i = 1, Bs(1)+2*g
              ! density
              rhs(i,j,k,rhoF)=rhs(i,j,k,rhoF) -0.5_rk*phi1_inv(i,j,k)*mask(i,j,k,rhoF)*(phi_prime(i,j,k,rhoF)-  Phi_ref(i,j,k,rhoF) )
              ! x-velocity
              rhs(i,j,k,UxF)=rhs(i,j,k,UxF) -1.0_rk*phi1_inv(i,j,k)*mask(i,j,k,UxF)*(phi_prime(i,j,k,rhoF)*phi_prime(i,j,k,UxF)-  Phi_ref(i,j,k,UxF) )
              ! y-velocity
              rhs(i,j,k,UyF)=rhs(i,j,k,UyF) -1.0_rk*phi1_inv(i,j,k)*mask(i,j,k,UyF)*(phi_prime(i,j,k,rhoF)*phi_prime(i,j,k,UyF)-  Phi_ref(i,j,k,UyF) )
              ! z-velocity
              rhs(i,j,k,UzF)=rhs(i,j,k,UzF) -1.0_rk*phi1_inv(i,j,k)*mask(i,j,k,UzF)*(phi_prime(i,j,k,rhoF)*phi_prime(i,j,k,UzF)-  Phi_ref(i,j,k,UzF) )
              ! preasure
              rhs(i,j,k,pF)=rhs(i,j,k,pF)                        -mask(i,j,k,pF)*(phi_prime(i,j,k,pF)- Phi_ref(i,j,k,pF) )
            end do
          end do
        end do
    endif

end subroutine RHS_3D_navier_stokes

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

subroutine  diffxyz_c_3D_opt( Bs, g, dx, dy, dz, u, dudx, dudy, dudz)

    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dx, dy, dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv, dy_inv, dz_inv

    ! - do not use all ghost nodes (note: use only two ghost nodes to get correct second derivatives)
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)
    dy_inv = 1.0_rk/(12.0_rk*dy)
    dz_inv = 1.0_rk/(12.0_rk*dz)

    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do

end subroutine diffxyz_c_3D_opt

subroutine  diffx_c_3D_opt( Bs, g, dx, u, dudx)

    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
            end do
        end do
    end do

end subroutine diffx_c_3D_opt

subroutine  diffy_c_3D_opt( Bs, g, dy, u, dudy)

    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dy_inv = 1.0_rk/(12.0_rk*dy)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
            end do
        end do
    end do

end subroutine diffy_c_3D_opt

subroutine  diffz_c_3D_opt( Bs, g, dz, u, dudz)

    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dz_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dz_inv = 1.0_rk/(12.0_rk*dz)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do

end subroutine diffz_c_3D_opt

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!--- OLD -------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

subroutine grad_zentral_3D(Bs, g, dx, dy, dz, q, qx, qy, qz)
    use module_params
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dx, dy, dz
    real(kind=rk), intent(in)       :: q(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: qx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: qy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: qz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

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
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: q(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: qx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    !! XXX !!
    call diffx_c_3D( Bs(1), g, dx, q, qx)

end subroutine diff1x_zentral_3D

!---------------------------------------------------------------------------------------------

subroutine diff1y_zentral_3D(Bs, g, dy, q, qy)
    use module_params
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: q(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: qy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    !! YYY !!
    call diffy_c_3D( Bs, g, dy, q, qy)

end subroutine diff1y_zentral_3D

!---------------------------------------------------------------------------------------------

subroutine diff1z_zentral_3D(Bs, g, dz, q, qz)
    use module_params
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: q(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: qz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    !! ZZZ !!
    call diffz_c_3D( Bs, g, dz, q, qz)

end subroutine diff1z_zentral_3D

!---------------------------------------------------------------------------------------------

subroutine  diffx_c_3D( Bs, g, dx, u, dudx)
    use module_params
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, n

    n = size(u,2)

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
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, n

    n = size(u,3)

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
    integer(kind=ik), intent(in)    :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

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
