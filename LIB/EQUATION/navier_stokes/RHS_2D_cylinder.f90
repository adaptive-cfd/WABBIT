!--------------------------------------------------------------------------------------------------------------------------------------------------------
!> \file
!> \brief Right hand side for 2D navier stokes equation in cylindrical coordinates
!>        -------------------------------------------------------------------------
!>\f{eqnarray*}{
!! \frac{\partial \sqrt{\rho}}{\partial t} &=- \left(\frac{1}{r} \frac{\partial (\rho r u_r)}{\partial r} + \frac{\partial (\rho u_z)}{\partial z}\right)  \frac{1}{2\sqrt{\rho}}	\\
!! \frac{\partial (\sqrt{\rho} u_r)}{\partial t}   &= \left[-\frac{1}{2}\left(   \frac{1}{r}\frac{\partial (r \rho u_r u_r)}{\partial r} + \frac{\partial (\rho u_r u_z)}{\partial z} +  \rho u_r \frac{\partial ( u_r)}{\partial r} + \rho u_z \frac{\partial (u_r)}{\partial z}  \right)  -\frac{\partial p}{\partial r} + F_r\right]\frac{1}{\sqrt{\rho}}\\
!! \frac{\partial (\sqrt{\rho} u_z)}{\partial t} &=  \left[-\frac{1}{2}\left(\frac{1}{r}\frac{\partial (r \rho u_r u_z)}			{\partial r} + \frac{\partial (\rho u_z u_z)}{\partial z} + \rho u_r \frac{\partial ( u_z)}{\partial r} + \rho u_z \frac{\partial (u_z)}{\partial z}  \right)  -\frac{\partial p}{\partial z} + F_z\right]\frac{1}{\sqrt{\rho}}\\
!! \frac{\partial p}{\partial t}  &= - \gamma \left(  \frac{1}{r} \frac{\partial(r u_r P)}{\partial r}   +  \ \frac{\partial( u_z P)}{\partial z}     \right) \\
!! &+ (\gamma-1)\left(
!! \frac{1}{r}\frac{\partial}{\partial r} \left[ \lambda r \frac{\partial T}{\partial r}  \right] + \frac{\partial}{\partial z} \left[ \lambda  \frac{\partial T}{\partial z}  \right] + \phi +  u_r\frac{\partial p}{\partial r} + u_z\frac{\partial p}{\partial z} - u_r F_r - u_z F_z \right)
!!\f}
!! with:
!>\f{eqnarray*}{
!!   F_r &= \frac{1}{r}\frac{\partial (r\tau_{rr})}{\partial r}  +
!!   \frac{\partial \tau_{rz}}{\partial z} - \frac{\tau_{\theta \theta}}{r}\\
!!   F_z &= \frac{1}{r}\frac{\partial (r\tau_{rz})}{\partial r}  +
!!   \frac{\partial \tau_{zz}}{\partial z} \\
!!   \phi &= \frac{1}{r}\frac{\partial}{\partial r} \left( r\,(u_r\tau_{rr} + u_z\tau_{rz}\,) \right) 	+
!!   \frac{\partial}{\partial z} \left( u_r\tau_{rz} + u_z\tau_{zz} \right)\\
!!   \tau_{rr} &= 2\mu\left(\frac{\partial u_r}{\partial r} - \frac{1}{3}(\nabla \cdot \vec{u})\right)\\
!!   \tau_{\theta \theta} &=2\mu \left( \frac{u_r}{r} -  \frac{1}{3}(\nabla \cdot \vec{u})\right)    \\
!!   \tau_{zz} &= 2 \mu \left( \frac{\partial u_z}{\partial z } -  \frac{1}{3}(\nabla \cdot \vec{u})\right) \\
!!   \tau_{rz} &=\mu \left[ \frac{\partial u_z}{\partial r} + \frac{\partial u_r}{\partial z}  \right]\\
!!   \nabla \cdot \vec{u} &= \frac{1}{r}\frac{\partial (ru_r)}{\partial r} + \frac{\partial u_z}{\partial z}
!!\f}
!! \date 08.12.18 - creation
!> \author Pkrah
!--------------------------------------------------------------------------------------------------------------------------------------------------------
!>\brief main function of RHS_2D_cylinders
subroutine RHS_2D_cylinder( g, Bs, x0, dx, phi, rhs, boundary_flag)

    implicit none
    !--------------------------------------------------------
    integer(kind=ik), intent(in)            :: g, Bs         !< # ghost and bulk points
    real(kind=rk), dimension(2), intent(in) :: x0, dx        !< grid coordinates
    real(kind=rk), intent(in)               :: phi(:, :, :)  !< statevector
    real(kind=rk), intent(inout)            :: rhs(:, :, :)  !< rhs array
    ! when implementing boundary conditions, it is necessary to now if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
    ! boundary_flag(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! curently only acessible in the local stage
    integer(kind=2), intent(in)            :: boundary_flag(3)
    !--------------------------------------------------------
    real(kind=rk)     :: gamma_, Rs, Cv, Cp, Pr, mu0 ! physics constants
    real(kind=rk)     :: dr, dz ! lattice spacing cylinder coordinates
    logical           :: dissipation
    ! variables
    real(kind=rk)  :: rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g), &
                   T(Bs+2*g, Bs+2*g), mu(Bs+2*g, Bs+2*g), mu_d(Bs+2*g, Bs+2*g), lambda(Bs+2*g, Bs+2*g), &
                   fric_p(Bs+2*g, Bs+2*g), fric_u(Bs+2*g, Bs+2*g), fric_v(Bs+2*g, Bs+2*g), &
                   lambdaT_r(Bs+2*g, Bs+2*g), lambdaT_z(Bs+2*g, Bs+2*g), &
                   tau_rr(Bs+2*g, Bs+2*g), tau_zz(Bs+2*g, Bs+2*g), &
                   tau_tt(Bs+2*g, Bs+2*g), tau_rz(Bs+2*g, Bs+2*g)
    ! derivatives
    real(kind=rk)  :: rho_v(Bs+2*g, Bs+2*g), rho_u(Bs+2*g, Bs+2*g), &
                   u_z(Bs+2*g, Bs+2*g), u_r(Bs+2*g, Bs+2*g), &
                   v_z(Bs+2*g, Bs+2*g), v_r(Bs+2*g, Bs+2*g), &
                   p_r(Bs+2*g, Bs+2*g), p_z(Bs+2*g, Bs+2*g), sqrt_rho_inv(Bs+2*g, Bs+2*g)

    real(kind=rk)  :: r(Bs+2*g,Bs+2*g), r_inv(Bs+2*g, Bs+2*g), r0
    ! tmp1 field
    real(kind=rk)       :: tmp1(Bs+2*g, Bs+2*g)
    integer(kind=ik)    :: i,j
    !----------------------------------------------------------

    ! pysical constants
    gamma_      = params_ns%gamma_
    Rs          = params_ns%Rs
    Cv          = params_ns%Cv
    Cp          = params_ns%Cp
    Pr          = params_ns%Pr
    mu0         = params_ns%mu0
    dissipation = params_ns%dissipation
    dz          = dx(1) ! lattice spacing in axial direction
    dr          = dx(2) ! lattice spacing of the radial component
    ! The total grid is shifted by R_min, which accounts for the
    ! infinitesimal cylinder centered arround the symmetrie axis.
    r0 = x0(2) + params_ns%R_min

    ! preperation of physical fields and arrays
    do j = 1, Bs+2*g  ! index of radial component
        do i = 1, Bs+2*g  ! index of axial component
            r(i,j)     = dble(j-(g+1)) * dr + r0
            r_inv(i,j) = 1/r(i,j)
            rho(i,j)       = phi(i,j,rhoF) * phi(i,j,rhoF)
            sqrt_rho_inv(i,j)  = 1.0_rk / phi(i,j,rhoF)
            u(i,j)         = phi(i,j,UyF) * sqrt_rho_inv(i,j)
            v(i,j)         = phi(i,j,UxF) * sqrt_rho_inv(i,j)
            p(i,j)         = phi(i,j,pF)
            rho_v(i,j)     = rho(i,j)*v(i,j)
            rho_u(i,j)     = rho(i,j)*u(i,j)
        end do
    end do

    ! definitely needed multiple times
    call D_r(u, u_r)
    call D_z(u, u_z)
    call D_r(v, v_r)
    call D_z(v, v_z)

    ! ------------ continuums equation ---------------!
    call D_z( rho_v, tmp1)
    rhs(:,:,rhoF) = tmp1
    call D_r( rho_u, tmp1)
    rhs(:,:,rhoF) =  rhs(:,:,rhoF) + tmp1 + rho_u*r_inv
    rhs(:,:,rhoF) = -rhs(:,:,rhoF) * 0.5_rk * sqrt_rho_inv


    ! ------------- momentum equation ---------------- !
    if (dissipation) then
         ! Compute mu
         T    = p / (rho*Rs) ! ideal gas
         mu   = mu0
         ! thermal conductivity
         lambda  = Cp * mu/Pr
         ! stress tensor
         tau_rr = 2.0_rk/3.0_rk * mu * (-v_z + 2.0_rk*u_r - u*r_inv)
         tau_tt = 2.0_rk/3.0_rk * mu * (-v_z -        u_r + 2.0_rk*u*r_inv)
         tau_zz = 2.0_rk/3.0_rk * mu * (2.0_rk*v_z -  u_r - u*r_inv)
         tau_rz =                 mu * ( v_r + u_z)
         ! Friction terms for the momentum equation
         call D_r(tau_rr,tmp1)
         fric_u = tmp1 + r_inv * (tau_rr - tau_tt)
         call D_z(tau_rz,tmp1)
         fric_u = fric_u + tmp1

         call D_r(tau_rz,tmp1)
         fric_v = tmp1 + tau_rz*r_inv
         call D_z(tau_zz, tmp1)
         fric_v = fric_v + tmp1
         ! Friction terms for the energy equation
         ! Heat Flux
         call D_r( T, tmp1)
         lambdaT_r = r*lambda*tmp1
         call D_z( T, tmp1)
         lambdaT_z = lambda * tmp1

         call D_r(lambdaT_r, tmp1)
         fric_p = r_inv * tmp1
         call D_z(lambdaT_z, tmp1)
         fric_p = fric_p + tmp1

         ! Friction terms because of stress (at work)
         call D_r(u*tau_rr + v*tau_rz, tmp1)
         fric_p = fric_p + tmp1
         call D_z(u*tau_rz + v*tau_zz, tmp1)
         fric_p = fric_p + tmp1

         fric_p = fric_p + r_inv * (u*tau_rr + v*tau_rz)

    else
         fric_p = 0.0_rk
         fric_u = 0.0_rk
         fric_v = 0.0_rk
    end if

    ! momentum in radial direction
    call D_r( rho_u*u, tmp1)
    rhs(:, :, UyF) = tmp1 + rho_u*u*r_inv
    call D_z(rho_v*u, tmp1)
    rhs(:, :, UyF) = rhs(:, :, UyF) + tmp1 + u_z * rho_v + u_r * rho_u
    call D_r(p, p_r)
    rhs(:, :, UyF) = -(rhs(:, :, UyF)*0.5_rk + p_r - fric_u)*sqrt_rho_inv

    ! momentum in axial direction
    call D_r( rho_u*v, tmp1)
    rhs(:, :, UxF) = tmp1 + rho_u * v * r_inv
    call D_z(rho_v*v, tmp1)
    rhs(:, :, UxF) = rhs(:, :, UxF) + tmp1 + v_z * rho_v + v_r * rho_u
    call D_z(p, p_z)
    rhs(:, :, UxF) = -(rhs(:, :, UxF)*0.5_rk + p_z - fric_v)*sqrt_rho_inv

    ! ---------- energy equation -------------- !
    call D_r(u*p, tmp1)
    rhs(:, :, pF) = (gamma_ -1) * (-u*p_r - v*p_z + u*fric_u + v*fric_v - fric_p ) &
                  + gamma_*tmp1
    call D_z(v*p, tmp1)
    rhs(:, :, pF) = -(rhs(:, :, pF) + gamma_*(u*p*r_inv + tmp1))

    contains

        ! ----------- inline functions ------------!
        ! for better readabillity we redefine the
        ! derivatives:
        !   # D_r - radial derivative
        !   # D_z - axial derivative

        !> Derivative in the radial direction
        subroutine  D_r( u, dudr)
            real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g)
            real(kind=rk), intent(out)      :: dudr(Bs+2*g, Bs+2*g)
            !> \details Note Bs, g, dz, boundary_flag are defined in the supfunction!
            call diffy_c_opt( Bs, g, dr, u, dudr, boundary_flag(2))

        end subroutine D_r

        !> Derivative in the axial direction
        subroutine  D_z( u, dudz)
            real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g)
            real(kind=rk), intent(out)      :: dudz(Bs+2*g, Bs+2*g)
            !> \details Note Bs, g, dz, boundary_flag are defined in the supfunction!
            call diffx_c_opt( Bs, g, dz, u, dudz, boundary_flag(1))

        end subroutine D_z

end subroutine RHS_2D_cylinder
