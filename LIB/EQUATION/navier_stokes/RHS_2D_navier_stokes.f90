!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_2D_navier_stokes.f90
!> \version 0.4
!> \author msr
!
!> \brief RHS for 2D navier stokes equation
!
!>
!! input:    - datafield, grid parameter, velocity, diffusion coefficient, derivative order \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 08/12/16 - create \n
!! 08/1/2018- include mask and sponge terms
! ********************************************************************************************
!>\details
!> We implement the right hand side of navier stokes in the skew symmetric form:
!>\f{eqnarray*}{
!!     \partial_t \sqrt{\rho} &=& -\frac{1}{2J\sqrt{\rho}} \nabla \cdot (\rho \vec{u})-\frac{1}{\sqrt{\rho}}\frac{1}{C_{\rm SP} } (\rho-\rho^{\rm SP}) \\
!!    \partial_t (\sqrt{\rho} u_\alpha) &=& -\frac{1}{2J \sqrt{\rho}} 
!!                                          \left[  
!!                                                       (u_\alpha \partial_\beta (\rho u_\beta)+
!!                                                        u_\beta \rho \partial_\beta u_\alpha)
!!                                            \right] 
!!                                           -\frac{1}{J \sqrt{\rho}} \partial_\beta \tau_{\alpha\beta}
!!                                           -\frac{1}{\sqrt{\rho}} \partial_\alpha p
!!                                            -\frac{1}{\sqrt{\rho}} \frac{1}{C_{\rm SP} }(\rho u_\alpha-\rho^{\rm SP} u_\alpha^{\rm SP}) \\
!!    \partial_t p &=& -\frac{\gamma}{J} \partial_\beta( u_\beta p) + (\gamma-1)(u_\alpha \partial_\alpha p)
!!                                      +\frac{\gamma -1}{J}
!!                                           \left[ 
!!                                                       \partial_\alpha(u_\beta \tau_{\alpha\beta}+\phi_\alpha) 
!!                                                       - u_\alpha\partial_\beta \tau_{\alpha\beta}
!!                                            \right]    
!!\f}

!! where the friction terms \f$ r_1,r_2,r_3 \f$ are
subroutine RHS_2D_navier_stokes( g, Bs, x0, delta_x, phi, rhs)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !< grid parameter
    integer(kind=ik), intent(in)                            :: g, Bs
    !< origin and spacing of the block
    real(kind=rk), dimension(2), intent(in)                  :: x0, delta_x
    !> datafields
    real(kind=rk), intent(in)                            :: phi(:, :, :)
    ! rhs array
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
    real(kind=rk)                                           :: mu0,dx,dy
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

    ! dummy field
    real(kind=rk)                                           :: dummy(Bs+2*g, Bs+2*g)

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
    v           = phi(:,:,3)/phi(:,:,2)
    p           = phi(:,:,4)
    ! discretization constant
    dx=delta_x(1)
    dy=delta_x(2)

    ! rhs
    rhs         = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! derivatives
    call grad_zentral( Bs, g, dx, dy, u, u_x, u_y)
    call grad_zentral( Bs, g, dx, dy, v, v_x, v_y)

    call diff1x_zentral( Bs, g, dx, p, p_x)
    call diff1y_zentral( Bs, g, dy, p, p_y)

    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    call diff1x_zentral( Bs, g, dx, rho*u, dummy)
    rhs(:,:,1) = -dummy
    call diff1y_zentral( Bs, g, dy, rho*v, dummy)
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
        call diff1x_zentral( Bs, g, dx, u, dummy)
        div_U = dummy
        call diff1y_zentral( Bs, g, dy, v, dummy)
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
        call diff1x_zentral( Bs, g, dx, tau11, dummy)
        fric_u = dummy
        call diff1y_zentral( Bs, g, dy, tau12, dummy)
        fric_u = fric_u + dummy

        fric_u = fric_u / phi(:,:,1)

        call diff1x_zentral( Bs, g, dx, tau12, dummy)
        fric_v = dummy
        call diff1y_zentral( Bs, g, dy, tau22, dummy)
        fric_v = fric_v + dummy

        fric_v = fric_v / phi(:,:,1)

        ! Friction terms for the energy equation
        ! Heat Flux
        call grad_zentral( Bs, g, dx, dy, T, T_x, T_y)

        fric_T1 = lambda * T_x
        fric_T2 = lambda * T_y

        ! All simple divergence terms for u_i*tau_ik and phi_k
        call diff1x_zentral( Bs, g, dx, ( u*tau11 + v*tau12 + fric_T1 ), dummy)
        fric_p = dummy
        call diff1y_zentral( Bs, g, dy, ( u*tau12 + v*tau22 + fric_T2 ), dummy)
        fric_p = fric_p + dummy

        ! u_i*dx_k (tau_ik) terms
        call diff1x_zentral( Bs, g, dx, tau11, dummy)
        fric_p = fric_p - u*dummy
        call diff1y_zentral( Bs, g, dy, tau12, dummy)
        fric_p = fric_p - u*dummy

        call diff1x_zentral( Bs, g, dx, tau12, dummy)
        fric_p = fric_p - v*dummy
        call diff1y_zentral( Bs, g, dy, tau22, dummy)
        fric_p = fric_p - v*dummy

        fric_p = ( gamma_ - 1 ) * fric_p

    else

        fric_p = 0.0_rk
        fric_u = 0.0_rk
        fric_v = 0.0_rk

    end if

    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    call diff1x_zentral( Bs, g, dx, (u * p), dummy)
    rhs(:,:,4) = - dummy
    call diff1y_zentral( Bs, g, dy, (v * p), dummy)
    rhs(:,:,4) = rhs(:,:,4) - dummy

    rhs(:,:,4) = rhs(:,:,4) * gamma_

    rhs(:,:,4) = rhs(:,:,4) + (gamma_ - 1.0_rk) * (u*p_x + v*p_y)

    rhs(:,:,4) = rhs(:,:,4) + fric_p

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    call diff1x_zentral( Bs, g, dx, (u * rho * u), dummy)
    rhs(:,:,2) = - 0.5_rk * dummy
    call diff1y_zentral( Bs, g, dy, (v * rho * u), dummy)
    rhs(:,:,2) = rhs(:,:,2) - 0.5_rk * dummy

    rhs(:,:,2) = rhs(:,:,2) - 0.5_rk * rho * u * u_x
    rhs(:,:,2) = rhs(:,:,2) - 0.5_rk * rho * v * u_y

    rhs(:,:,2) = rhs(:,:,2) - p_x

    rhs(:,:,2) = rhs(:,:,2) / phi(:,:,1)

    rhs(:,:,2) = rhs(:,:,2) + fric_u

    ! RHS of  momentum equation for v
    call diff1x_zentral( Bs, g, dx, (u * rho * v), dummy)
    rhs(:,:,3) = - 0.5_rk * dummy
    call diff1y_zentral( Bs, g, dy, (v * rho * v), dummy)
    rhs(:,:,3) = rhs(:,:,3) - 0.5_rk * dummy

    rhs(:,:,3) = rhs(:,:,3) - 0.5_rk * rho * u * v_x
    rhs(:,:,3) = rhs(:,:,3) - 0.5_rk * rho * v * v_y

    rhs(:,:,3) = rhs(:,:,3) - p_y

    rhs(:,:,3) = rhs(:,:,3) / phi(:,:,1)

    rhs(:,:,3) = rhs(:,:,3) + fric_v


    ! SPONGE (bob)
    if (params_ns%sponge_layer) then
        ! add spnge 
        ! rho component (density)  
        rhs(:,:,1)=rhs(:,:,1) - 0.5_rk/phi(:,:,1)*sponge( Bs, g, x0,delta_x, rho,   rho0_       , params_ns%C_sp)
        ! sqrt(rho)u component (momentum)
        rhs(:,:,2)=rhs(:,:,2) - 1.0_rk/phi(:,:,1)*sponge( Bs, g, x0,delta_x, rho*u, rho0_*u0_   , params_ns%C_sp)
        ! sqrt(rho)v component (momentum)
        rhs(:,:,3)=rhs(:,:,3) - 1.0_rk/phi(:,:,1)*sponge( Bs, g, x0,delta_x, rho*v, rho0_*v0_   , params_ns%C_sp)
        ! p component (preasure/energy)
        rhs(:,:,4)=rhs(:,:,4) - (gamma_-1)       *sponge( Bs, g, x0,delta_x, p    , p0_         , params_ns%C_sp)
        
    endif



end subroutine RHS_2D_navier_stokes

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

subroutine grad_zentral(Bs, g, dx, dy, q, qx, qy)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dx, dy
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qx(Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qy(Bs+2*g, Bs+2*g)

    !! XXX !!
    call diffx_c( Bs, g, dx, q, qx)

    !! YYY !!
    call diffy_c( Bs, g, dy, q, qy)

end subroutine grad_zentral

!---------------------------------------------------------------------------------------------

subroutine diff1x_zentral(Bs, g, dx, q, qx)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qx(Bs+2*g, Bs+2*g)

    !! XXX !!
    call diffx_c( Bs, g, dx, q, qx)

end subroutine diff1x_zentral

!---------------------------------------------------------------------------------------------

subroutine diff1y_zentral(Bs, g, dy, q, qy)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: qy(Bs+2*g, Bs+2*g)

    !! XXX !!
    call diffy_c( Bs, g, dy, q, qy)

end subroutine diff1y_zentral

!---------------------------------------------------------------------------------------------

subroutine  diffx_c( Bs, g, dx, u, dudx)
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs+2*g, Bs+2*g)

    integer                         :: i, n

    n = size(u,1)

    dudx(1,:) = ( u(n-1,:) - 8.0_rk*u(n,:) + 8.0_rk*u(2,:) - u(3,:) ) / (12.0_rk*dx)
    dudx(2,:) = ( u(n,:)   - 8.0_rk*u(1,:) + 8.0_rk*u(3,:) - u(4,:) ) / (12.0_rk*dx)

    forall ( i = 3:n-2 )
       dudx(i,:) = ( u(i-2,:) - 8.0_rk*u(i-1,:) + 8.0_rk*u(i+1,:) - u(i+2,:) ) / (12.0_rk*dx)
    end forall

    dudx(n-1,:) = ( u(n-3,:) - 8.0_rk*u(n-2,:) + 8.0_rk*u(n,:) - u(1,:) ) / (12.0_rk*dx)
    dudx(n,:)   = ( u(n-2,:) - 8.0_rk*u(n-1,:) + 8.0_rk*u(1,:) - u(2,:) ) / (12.0_rk*dx)

end subroutine diffx_c

!---------------------------------------------------------------------------------------------

subroutine  diffy_c( Bs, g, dy, u, dudy)
    use module_params
    integer(kind=ik), intent(in)    :: g, Bs
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs+2*g, Bs+2*g)

    integer                         :: i, n

    n = size(u,1)

    dudy(:,1) = ( u(:,n-1) - 8.0_rk*u(:,n) + 8.0_rk*u(:,2) - u(:,3) ) / (12.0_rk*dy)
    dudy(:,2) = ( u(:,n)   - 8.0_rk*u(:,1) + 8.0_rk*u(:,3) - u(:,4) ) / (12.0_rk*dy)

    forall ( i = 3:n-2 )
       dudy(:,i) = ( u(:,i-2) - 8.0_rk*u(:,i-1) + 8.0_rk*u(:,i+1) - u(:,i+2) ) / (12.0_rk*dy)
    end forall

    dudy(:,n-1) = ( u(:,n-3) - 8.0_rk*u(:,n-2) + 8.0_rk*u(:,n) - u(:,1) ) / (12.0_rk*dy)
    dudy(:,n)   = ( u(:,n-2) - 8.0_rk*u(:,n-1) + 8.0_rk*u(:,1) - u(:,2) ) / (12.0_rk*dy)

end subroutine diffy_c







!--------------------------------------------------------------------------!
!               phils sponge stuff
!--------------------------------------------------------------------------!


!==========================================================================
!> \brief This function computes a 2d sponge term 
!!   \f{equation}{
!!           s(x,y)=\frac{\chi_{\mathrm sp}(x,y)}{C_{\mathrm sp}}(q(x,y)-q_{\mathrm ref})
!!                          \quad \forall x,y \in \Omega_{\mathrm block} 
!!     \f}
!! Where we addopt the notation from <a href="https://arxiv.org/abs/1506.06513">Thomas Engels (2015)</a> 
!! - the mask function is \f$\chi_{\rm sp}\f$ 
!! - \f$C_{\rm sp}\f$ is the sponge coefficient (normaly \f$10^{-1}\f$)

function sponge( Bs, g, x0,dx, q, qref, C_sp)
    !-------------------------------------------------------
    !> grid parameter
    integer(kind=ik), intent(in)    :: g, Bs
    !> spacing and origin of block
    real(kind=rk), intent(in)       :: x0(2), dx(2)                    
    !> Sponge coefficient
    real(kind=rk), intent(in)       :: C_sp
    !> reference value of quantity \f$q\f$ (veolcity \f$u\f$, preasure \f$p\f$,etc.)
    real(kind=rk), intent(in)       :: qref
    !> quantity \f$q\f$ (veolcity \f$u\f$, preasure \f$p\f$,etc.)   
    real(kind=rk), intent(in)       :: q(Bs+2*g, Bs+2*g)
    !> sponge term \f$s(x,y)\f$
    real(kind=rk)                   :: sponge(Bs+2*g, Bs+2*g)
    !--------------------------------------------------------
    ! loop variables
    integer                         :: i, n,ix,iy
    ! inverse C_sp
    real(kind=rk)                   :: C_sp_inv,x,y
    
    C_sp_inv=1.0_rk/C_sp

     do ix=1, Bs+2*g
        x = dble(ix-(g+1)) * dx(1) + x0(1)
        
        do iy=1, Bs+2*g
            y = dble(iy-(g+1)) * dx(2) + x0(2)

            if (inside_sponge((/x,y/))) then
                   sponge(ix,iy) = C_sp_inv*(q(ix,iy)-qref)
            else
                   sponge(ix,iy) = 0.0_rk
            end if

       end do
    end do

end function sponge
!==========================================================================



!==========================================================================
!> \brief This function f(x) implements \n
!> f(x) is 1 if x(1)<=Lsponge \n
!> f(x) is 0 else  \n
function inside_sponge(x)
!> coordinate vector \f$\vec{x}=(x,y,z)\f$ (real 3d or 2d array)
real(kind=rk), intent(in)       :: x(:)    
!> logical 
logical                         :: inside_sponge
! dimension of array x
integer                         :: dim,i
! size of sponge
real(kind=rk)                   :: length_sponge

!> \todo read in length_sponge or thing of something intelligent here
        length_sponge=params_ns%Lx*0.05 
        if (x(1)<=length_sponge .and. x(1)>=0) then
            inside_sponge=.true.
        else
            inside_sponge=.false.
        endif


end function inside_sponge
!==========================================================================




subroutine create_mask_2D_NEW(mask, x0, dx, Bs, g,geometry )

    ! use module_params
    ! use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in)                    :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(inout)    :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)         :: x0, dx
    !< geometry to of body \f$\Omega_s\f$
     character(len=80), intent(in)                  :: geometry


    real(kind=rk) :: cx ,cy,R_cyl,x,y,r,h
    integer :: iy,ix

    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")


    select case(geometry)
    !case('cylinder')
      !call draw_cylinder( mask, x0, dx, Bs, g )
    !case('two-cylinders')
      !call draw_two_cylinders( mask, x0, dx, Bs, g )
   ! case('rectangle')

    case default
      call abort(120001,"ERROR: geometry for VPM is unknown"//geometry)
    end select

end subroutine create_mask_2D_NEW


! subroutine draw_cylinder(mask, x0, dx, Bs, g )

!     use module_params
!     use module_precision

!     implicit none

!     ! grid
!     integer(kind=ik), intent(in)                              :: Bs, g
!     !> mask term for every grid point of this block
!     real(kind=rk), dimension(:,:), intent(out)     :: mask
!     !> spacing and origin of block
!     real(kind=rk), dimension(2), intent(in)                   :: x0, dx

!     ! auxiliary variables
!     real(kind=rk)                                             :: x, y, r, h
!     ! loop variables
!     integer(kind=ik)                                          :: ix, iy

! !---------------------------------------------------------------------------------------------
! ! variables initialization
!     if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

!     ! reset mask array
!     mask = 0.0_rk

! !---------------------------------------------------------------------------------------------
! ! main body


!     ! parameter for smoothing function (width)
!     h = 1.5_rk*max(dx(1), dx(2))

!     do ix=1, Bs+2*g
!        x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%x_cntr(1)
!        do iy=1, Bs+2*g
!            y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%x_cntr(2)
!            ! distance from center of cylinder
!            r = dsqrt(x*x + y*y)
!            if (params_acm%smooth_mask) then
!                call smoothstep(mask(ix,iy), r, params_acm%R_cyl, h)
!            else
!                ! if point is inside the cylinder, set mask to 1
!                if (r <= params_acm%R_cyl) then
!                    mask(ix,iy) = 1.0_rk
!                else
!                    mask(ix,iy) = 0.0_rk
!                end if
!            end if
!        end do
!     end do

! end subroutine draw_cylinder



! subroutine draw_two_cylinders( mask, x0, dx, Bs, g)

!   use module_params
!   use module_precision

!   implicit none

!   ! grid
!   integer(kind=ik), intent(in)                              :: Bs, g
!   !> mask term for every grid point of this block
!   real(kind=rk), dimension(:,:), intent(out)     :: mask
!   !> spacing and origin of block
!   real(kind=rk), dimension(2), intent(in)                   :: x0, dx

!   ! auxiliary variables
!   real(kind=rk)                                             :: x1, x2, y1, y2, R, cx1, cx2, cy1,&
!   cy2, r_1, r_2, h, mask1, mask2
!   ! loop variables
!   integer(kind=ik)                                          :: ix, iy

!   !---------------------------------------------------------------------------------------------
!   ! variables initialization
!   if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

!   ! reset mask array
!   mask = 0.0_rk
!   mask1 = 0.0_rk
!   mask2 = 0.0_rk

!   !---------------------------------------------------------------------------------------------
!   ! main body

!   ! center of the first cylinder
!   cx1 = 0.5884_rk*params_acm%Lx
!   cy1 = 0.4116_rk*params_acm%Ly

!   ! center of the second cylinder
!   cx2 = 0.4116_rk*params_acm%Lx
!   cy2 = 0.5884_rk*params_acm%Ly

!   ! radius of the cylinders
!   R = params_acm%R_cyl
!   ! parameter for smoothing function (width)
!   h = 1.5_rk*max(dx(1), dx(2))

!   do ix=1, Bs+2*g
!     x1 = dble(ix-(g+1)) * dx(1) + x0(1) - cx1
!     x2 = dble(ix-(g+1)) * dx(1) + x0(1) - cx2
!     do iy=1, Bs+2*g
!       y1 = dble(iy-(g+1)) * dx(2) + x0(2) - cy1
!       y2 = dble(iy-(g+1)) * dx(2) + x0(2) - cy2
!       ! distance from center of cylinder 1
!       r_1 = dsqrt(x1*x1 + y1*y1)
!       ! distance from center of cylinder 2
!       r_2 = dsqrt(x2*x2 + y2*y2)
!       if (params_acm%smooth_mask) then
!         call smoothstep(mask1, r_1, R, h)
!         call smoothstep(mask2, r_2, R, h)
!         mask(ix,iy) = mask1 + mask2
!       else
!         ! if point is inside one of the cylinders, set mask to 1
!         if (r_1 <= R) then
!           mask(ix,iy) = 1.0_rk
!         elseif ( r_2 <= R) then
!           mask(ix,iy) = 1.0_rk
!         else
!           mask(ix,iy) = 0.0_rk
!         end if
!       end if
!     end do
!   end do


! end subroutine draw_two_cylinders