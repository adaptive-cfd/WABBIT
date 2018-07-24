!-----------------------------------------------------------------------------
! main level wrapper to set the right hand side on a block. Note this is completely
! independent of the grid any an MPI formalism, neighboring relations and the like.
! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
! from that. Ghost nodes are assumed to be sync'ed.
!-----------------------------------------------------------------------------
subroutine RHS_ACM( time, u, g, x0, dx, rhs, stage )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! output. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in) :: stage

    ! local variables
    integer(kind=ik) :: Bs, mpierr
    real(kind=rk) :: tmp(1:3), tmp2

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.
        ! performs initializations in the RHS module, such as resetting integrals

        params_acm%mean_flow = 0.0_rk
        params_acm%mean_p = 0.0_rk

        if (params_acm%geometry == "Insect") call Update_Insect(time, Insect)

    case ("integral_stage")
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! For some RHS, the eqn depend not only on local, block based qtys, such as
        ! the state vector, but also on the entire grid, for example to compute a
        ! global forcing term (e.g. in FSI the forces on bodies). As the physics
        ! modules cannot see the grid, (they only see blocks), in order to encapsulate
        ! them nicer, two RHS stages have to be defined: integral / local stage.
        !
        ! called for each block.
        if (maxval(abs(u))>1.0e5) then
            call abort(6661,"ACM fail: very very large values in state vector.")
        endif

        if (params_acm%dim == 2) then
            params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1))*dx(1)*dx(2)
            params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 2))*dx(1)*dx(2)
            params_acm%mean_p = params_acm%mean_p + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 3))*dx(1)*dx(2)
        else
            params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 1))*dx(1)*dx(2)*dx(3)
            params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 2))*dx(1)*dx(2)*dx(3)
            params_acm%mean_flow(3) = params_acm%mean_flow(3) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 3))*dx(1)*dx(2)*dx(3)
            params_acm%mean_p = params_acm%mean_p + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 4))*dx(1)*dx(2)*dx(3)
        endif ! NOTE: MPI_SUM is perfomed in the post_stage.

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.

        tmp = params_acm%mean_flow
        call MPI_ALLREDUCE(tmp, params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        tmp2 = params_acm%mean_p
        call MPI_ALLREDUCE(tmp2, params_acm%mean_p, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        if (params_acm%dim == 2) then
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly)
            params_acm%mean_p = params_acm%mean_p / (params_acm%Lx*params_acm%Ly)
        else
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
            params_acm%mean_p = params_acm%mean_p / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
        endif

    case ("local_stage")
        !-------------------------------------------------------------------------
        ! 4th stage: local evaluation of RHS on all blocks
        !-------------------------------------------------------------------------
        ! the second stage then is what you would usually do: evaluate local differential
        ! operators etc.
        !
        ! called for each block.

        if (params_acm%dim == 2) then
            ! this is a 2d case (ux,uy,p)
            call RHS_2D_acm(g, Bs, dx(1:2), x0(1:2), u(:,:,1,:), params_acm%discretization, time, rhs(:,:,1,:))

        else
            ! this is a 3d case (ux,uy,uz,p)
            call RHS_3D_acm(g, Bs, dx, x0, u, params_acm%discretization, time, rhs)

        endif

    case default
        call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")

    end select


end subroutine RHS_ACM









!> \file
!> \brief Right hand side for 2D artificial compressibility equations
!>        ---------------------------------------------
!> The right hand side for the artificial compressibility equations reads:
!>\f{eqnarray*}{
!! \partial_t u &=& -u \nabla \cdot u - \nu \nabla^2 u +  \frac{1}{\rho} \nabla p - \chi \frac{1}{C_\eta} (u-u_s)  \\
!! \partial_t p &=& -c_0^2 \nabla \cdot u - \gamma p
!!\f}
!> \version 0.5
!> \version  27/06/17 - create
!!
!> \author sm, engels
!--------------------------------------------------------------------------------------------

subroutine RHS_2D_acm(g, Bs, dx, x0, phi, order_discretization, time, rhs)

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)            :: g, Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(2), intent(in) :: x0, dx
    !> datafields
    real(kind=rk), intent(inout)            :: phi(:,:,:)
    real(kind=rk), intent(inout)            :: rhs(:,:,:)
    !> discretization order
    character(len=80), intent(in)           :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    !> mask term for every grid point in this block
    real(kind=rk), allocatable, save :: mask(:, :), sponge(:, :)
    !> velocity of the solid
    real(kind=rk), allocatable, save :: us(:, :, :)
    !> forcing term
    real(kind=rk), dimension(3) :: forcing
    !>
    real(kind=rk) :: dx_inv, dy_inv, dx2_inv, dy2_inv, c_0, nu, eps, eps_inv, gamma
    real(kind=rk) :: div_U, u_dx, u_dy, u_dxdx, u_dydy, v_dx, v_dy, v_dxdx, &
                     v_dydy, p_dx, p_dy, penalx, penaly, x, y, term_2
    ! loop variables
    integer(kind=rk) :: ix, iy, idir
    ! coefficients for Tam&Webb
    real(kind=rk) :: a(-3:3)
    real(kind=rk) :: b(-2:2)

    !---------------------------------------------------------------------------------------------
    ! variables initialization
    if (.not. allocated(sponge)) allocate(sponge(1:Bs+2*g, 1:Bs+2*g))
    if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g))
    if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1:2))

    ! set parameters for readability
    c_0         = params_acm%c_0
    nu          = params_acm%nu
    eps         = params_acm%C_eta
    gamma       = params_acm%gamma_p

    mask   = 0.0_rk
    us     = 0.0_rk
    sponge = 0.0_rk
    rhs = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    eps_inv = 1.0_rk / eps

    if (size(phi,1)/=Bs+2*g .or. size(phi,2)/=Bs+2*g .or. size(phi,3)/=3) then
        call abort(66233,"wrong size, I go for a walk instead.")
    endif

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    ! 4th order coefficients for second derivative
    b = (/ -1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)

    if (maxval(abs(params_acm%mean_flow)) > 1.0e3_rk ) then
        call abort(887, "ACM: meanflow out of bounds")
    endif

!---------------------------------------------------------------------------------------------
! main body

    if (params_acm%penalization) then
        ! create mask term for every grid point in this block
        call create_mask_2D(time, x0, dx, Bs, g, mask, us)
    end if


    if (order_discretization == "FD_2nd_central" ) then
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        do iy = g+1, Bs+g
            do ix = g+1, Bs+g

                u_dx   = (phi(ix+1,iy,1) - phi(ix-1,iy,1))*dx_inv*0.5_rk
                u_dy   = (phi(ix,iy+1,1) - phi(ix,iy-1,1))*dy_inv*0.5_rk
                u_dxdx = (phi(ix-1,iy,1) -2.0_rk*phi(ix,iy,1) +phi(ix+1,iy,1))*dx2_inv
                u_dydy = (phi(ix,iy-1,1) -2.0_rk*phi(ix,iy,1) +phi(ix,iy+1,1))*dy2_inv

                v_dx   = (phi(ix+1,iy,2) -phi(ix-1,iy,2))*dx_inv*0.5_rk
                v_dy   = (phi(ix,iy+1,2) -phi(ix,iy-1,2))*dy_inv*0.5_rk
                v_dxdx = (phi(ix-1,iy,2) -2.0_rk*phi(ix,iy,2) +phi(ix+1,iy,2))*dx2_inv
                v_dydy = (phi(ix,iy-1,2) -2.0_rk*phi(ix,iy,2) +phi(ix,iy+1,2))*dy2_inv

                p_dx = (phi(ix+1,iy,3) -phi(ix-1,iy,3))*dx_inv*0.5_rk
                p_dy = (phi(ix,iy+1,3) -phi(ix,iy-1,3))*dy_inv*0.5_rk

                div_U = u_dx + v_dy

                penalx = -mask(ix,iy) * eps_inv * (phi(ix,iy,1) -us(ix,iy,1))
                penaly = -mask(ix,iy) * eps_inv * (phi(ix,iy,2) -us(ix,iy,2))

                ! actual RHS. note mean flow forcing is just a constant and added at the end of the routine
                rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)

            end do
        end do

    else if (order_discretization == "FD_4th_central_optimized") then
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do iy = g+1, Bs+g
            do ix = g+1, Bs+g
                ! first derivatives of u, v, p
                u_dx = (a(-3)*phi(ix-3,iy,1) + a(-2)*phi(ix-2,iy,1) + a(-1)*phi(ix-1,iy,1) &
                     +  a(0 )*phi(ix,iy,1) + a(+1)*phi(ix+1,iy,1) + a(+2)*phi(ix+2,iy,1) + a(+3)*phi(ix+3,iy,1))*dx_inv
                u_dy = (a(-3)*phi(ix,iy-3,1) + a(-2)*phi(ix,iy-2,1) + a(-1)*phi(ix,iy-1,1) &
                     +  a(0 )*phi(ix,iy,1) + a(+1)*phi(ix,iy+1,1) + a(+2)*phi(ix,iy+2,1) + a(+3)*phi(ix,iy+3,1))*dy_inv

                v_dx = (a(-3)*phi(ix-3,iy,2) + a(-2)*phi(ix-2,iy,2) + a(-1)*phi(ix-1,iy,2) &
                     +  a(0 )*phi(ix,iy,2) + a(+1)*phi(ix+1,iy,2) + a(+2)*phi(ix+2,iy,2) + a(+3)*phi(ix+3,iy,2))*dx_inv
                v_dy = (a(-3)*phi(ix,iy-3,2) + a(-2)*phi(ix,iy-2,2) + a(-1)*phi(ix,iy-1,2) &
                     +  a(0 )*phi(ix,iy,2) + a(+1)*phi(ix,iy+1,2) + a(+2)*phi(ix,iy+2,2) + a(+3)*phi(ix,iy+3,2))*dy_inv

                p_dx = (a(-3)*phi(ix-3,iy,3) + a(-2)*phi(ix-2,iy,3) + a(-1)*phi(ix-1,iy,3) &
                     +  a(0 )*phi(ix,iy,3) + a(+1)*phi(ix+1,iy,3) + a(+2)*phi(ix+2,iy,3) + a(+3)*phi(ix+3,iy,3))*dx_inv
                p_dy = (a(-3)*phi(ix,iy-3,3) + a(-2)*phi(ix,iy-2,3) + a(-1)*phi(ix,iy-1,3) &
                     +  a(0 )*phi(ix,iy,3) + a(+1)*phi(ix,iy+1,3) + a(+2)*phi(ix,iy+2,3) + a(+3)*phi(ix,iy+3,3))*dy_inv

                ! second derivatives of u and v
                u_dxdx = (b(-2)*phi(ix-2,iy,1) + b(-1)*phi(ix-1,iy,1) + b(0)*phi(ix,iy,1) &
                       +  b(+1)*phi(ix+1,iy,1) + b(+2)*phi(ix+2,iy,1))*dx2_inv
                u_dydy = (b(-2)*phi(ix,iy-2,1) + b(-1)*phi(ix,iy-1,1) + b(0)*phi(ix,iy,1) &
                       +  b(+1)*phi(ix,iy+1,1) + b(+2)*phi(ix,iy+2,1))*dy2_inv

                v_dxdx = (b(-2)*phi(ix-2,iy,2) + b(-1)*phi(ix-1,iy,2) + b(0)*phi(ix,iy,2) &
                       +  b(+1)*phi(ix+1,iy,2) + b(+2)*phi(ix+2,iy,2))*dx2_inv
                v_dydy = (b(-2)*phi(ix,iy-2,2) + b(-1)*phi(ix,iy-1,2) + b(0)*phi(ix,iy,2) &
                       +  b(+1)*phi(ix,iy+1,2) + b(+2)*phi(ix,iy+2,2))*dy2_inv

                div_U = u_dx + v_dy

                penalx = -mask(ix,iy)*eps_inv*(phi(ix,iy,1)-us(ix,iy,1))
                penaly = -mask(ix,iy)*eps_inv*(phi(ix,iy,2)-us(ix,iy,2))

                rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
            end do
        end do

    else
        call abort(441166, "Discretization unkown "//order_discretization//", I ll walk into the light now." )
    end if

    ! --------------------------------------------------------------------------
    ! forcing term.
    ! --------------------------------------------------------------------------
    ! is forcing used at all?
    if (params_acm%forcing) then
        do idir = 1, 2
            select case (params_acm%forcing_type(idir))
            case('accelerate')
                forcing(idir) = max(0.0_rk, params_acm%u_mean_set(idir)-params_acm%mean_flow(idir)) &
                * startup_conditioner(time, 0.0_rk, 0.5_rk)
                ! add forcing to right hand side
                rhs(:,:,idir) = rhs(:,:,idir) + forcing(idir)

            case('fixed')
                ! note fixed forcing directly modifies the state vetor and is not a
                ! forcing term in the conventional sense.
                forcing(idir) = 0.0_rk
                phi(:,:,idir) = phi(:,:,idir) - params_acm%mean_flow(idir) + params_acm%u_mean_set(idir)
                ! add forcing to right hand side
                rhs(:,:,idir) = rhs(:,:,idir) + forcing(idir)

            case('none')
                ! do nothing in this direction.

            case('taylor_green')
                if (idir==1) then
                    do iy = g+1, Bs+g
                        do ix = g+1, Bs+g
                            x = x0(1) + dble(ix-g-1) * dx(1)
                            y = x0(2) + dble(iy-g-1) * dx(2)
                            call continue_periodic(x,params_acm%Lx)
                            call continue_periodic(y,params_acm%Ly)
                            term_2 = 2.0_rk*nu*dcos(time) - dsin(time)
                            forcing(1) = dsin(x - params_acm%u_mean_set(1)*time) * dcos(y - params_acm%u_mean_set(2)*time) *term_2
                            forcing(2) = -dcos(x - params_acm%u_mean_set(1)*time) * dsin(y - params_acm%u_mean_set(2) * time) *term_2
                            rhs(ix,iy,1:2) = rhs(ix,iy,1:2) + forcing(1:2)
                        end do
                    end do
                end if

            case default
                call abort(7710, "ACM::rhs.f90::meanflow forcing type unkown")

            end select
        end do
    end if

    ! remove mean pressure. NOTE: there is some oddities here, as the code modifies the
    ! state vector and not the RHS.
    if (params_acm%p_mean_zero) then
        phi(:,:,3) = phi(:,:,3) - params_acm%mean_p
    end if

    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    if (params_acm%use_sponge) then
        call sponge_2D(sponge, x0, dx, Bs, g)
        eps_inv = 1.0_rk / params_acm%C_sponge

        ! NOTE: the sponge term acts, if active, on ALL components, ux,uy,p
        ! which is different from the penalization term, which acts only on ux,uy and not p
        rhs(:,:,1) = rhs(:,:,1) - (phi(:,:,1)-params_acm%u_mean_set(1)) * sponge * eps_inv
        rhs(:,:,2) = rhs(:,:,2) - (phi(:,:,2)-params_acm%u_mean_set(2)) * sponge * eps_inv
        rhs(:,:,3) = rhs(:,:,3) - phi(:,:,3)*sponge*eps_inv
    end if


end subroutine RHS_2D_acm




subroutine RHS_3D_acm(g, Bs, dx, x0, phi, order_discretization, time, rhs)
    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)            :: g, Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in) :: x0, dx
    !> datafields
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)
    real(kind=rk), intent(inout)            :: rhs(:,:,:,:)
    !> discretization order
    character(len=80), intent(in)           :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    !> temporary, persistent arrays
    !> mask term for every grid point in this block
    real(kind=rk), allocatable, save :: mask(:, :, :), sponge(:, :, :)
    !> velocity of the solid
    real(kind=rk), allocatable, save :: us(:, :, :, :)

    !> forcing term
    real(kind=rk), dimension(3) :: forcing

    !> inverse dx, physics/acm parameters
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, c_0, &
                     nu, eps, eps_inv, gamma
    !> derivatives
    real(kind=rk) :: div_U, u_dx, u_dy, u_dz, u_dxdx, u_dydy, u_dzdz, &
                     v_dx, v_dy, v_dz, v_dxdx, v_dydy, v_dzdz, &
                     w_dx, w_dy, w_dz, w_dxdx, w_dydy, w_dzdz, &
                     p_dx, p_dy, p_dz, penalx, penaly, penalz
    !> loop variables
    integer(kind=rk) :: ix, iy, iz
    !> coefficients for Tam&Webb
    real(kind=rk) :: a(-3:3)
    real(kind=rk) :: b(-2:2)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    if (.not. allocated(sponge)) allocate(sponge(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
    if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
    if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g, 1:3))

    ! set parameters for readability
    c_0         = params_acm%c_0
    nu          = params_acm%nu
    eps         = params_acm%C_eta
    gamma       = params_acm%gamma_p

    mask = 0.0_rk
    us   = 0.0_rk
    rhs = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dz_inv = 1.0_rk / dx(3)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    dz2_inv = 1.0_rk / (dx(3)**2)

    eps_inv = 1.0_rk / eps

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    ! 4th order coefficients for second derivative
    b = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)

!---------------------------------------------------------------------------------------------
! main body

    if (params_acm%penalization) then
        ! create mask term for every grid point in this block
        call create_mask_3D(time, x0, dx, Bs, g, mask, us)
        mask = mask * eps_inv
    end if

    if (order_discretization == "FD_2nd_central" ) then
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        do iz = g+1, Bs+g
            do iy = g+1, Bs+g
                do ix = g+1, Bs+g

                    ! first and second derivatives of u,v,w
                    u_dx = (phi(ix+1,iy,iz,1)-phi(ix-1,iy,iz,1))*dx_inv*0.5_rk
                    u_dy = (phi(ix,iy+1,iz,1)-phi(ix,iy-1,iz,1))*dy_inv*0.5_rk
                    u_dz = (phi(ix,iy,iz+1,1)-phi(ix,iy,iz-1,1))*dz_inv*0.5_rk

                    u_dxdx = (phi(ix-1,iy,iz,1)-2.0_rk*phi(ix,iy,iz,1)+phi(ix+1,iy,iz,1))*dx2_inv
                    u_dydy = (phi(ix,iy-1,iz,1)-2.0_rk*phi(ix,iy,iz,1)+phi(ix,iy+1,iz,1))*dy2_inv
                    u_dzdz = (phi(ix,iy,iz-1,1)-2.0_rk*phi(ix,iy,iz,1)+phi(ix,iy,iz+1,1))*dz2_inv

                    v_dx = (phi(ix+1,iy,iz,2)-phi(ix-1,iy,iz,2))*dx_inv*0.5_rk
                    v_dy = (phi(ix,iy+1,iz,2)-phi(ix,iy-1,iz,2))*dy_inv*0.5_rk
                    v_dz = (phi(ix,iy,iz+1,2)-phi(ix,iy,iz-1,2))*dz_inv*0.5_rk

                    v_dxdx = (phi(ix-1,iy,iz,2)-2.0_rk*phi(ix,iy,iz,2)+phi(ix+1,iy,iz,2))*dx2_inv
                    v_dydy = (phi(ix,iy-1,iz,2)-2.0_rk*phi(ix,iy,iz,2)+phi(ix,iy+1,iz,2))*dy2_inv
                    v_dzdz = (phi(ix,iy,iz-1,2)-2.0_rk*phi(ix,iy,iz,2)+phi(ix,iy,iz+1,2))*dz2_inv

                    w_dx = (phi(ix+1,iy,iz,3)-phi(ix-1,iy,iz,3))*dx_inv*0.5_rk
                    w_dy = (phi(ix,iy+1,iz,3)-phi(ix,iy-1,iz,3))*dy_inv*0.5_rk
                    w_dz = (phi(ix,iy,iz+1,3)-phi(ix,iy,iz-1,3))*dz_inv*0.5_rk

                    w_dxdx = (phi(ix-1,iy,iz,3)-2.0_rk*phi(ix,iy,iz,3)+phi(ix+1,iy,iz,3))*dx2_inv
                    w_dydy = (phi(ix,iy-1,iz,3)-2.0_rk*phi(ix,iy,iz,3)+phi(ix,iy+1,iz,3))*dy2_inv
                    w_dzdz = (phi(ix,iy,iz-1,3)-2.0_rk*phi(ix,iy,iz,3)+phi(ix,iy,iz+1,3))*dz2_inv

                    ! first derivative of p
                    p_dx = (phi(ix+1,iy,iz,4)-phi(ix-1,iy,iz,4))*dx_inv*0.5_rk
                    p_dy = (phi(ix,iy+1,iz,4)-phi(ix,iy-1,iz,4))*dy_inv*0.5_rk
                    p_dz = (phi(ix,iy,iz+1,4)-phi(ix,iy,iz-1,4))*dz_inv*0.5_rk

                    div_U = u_dx + v_dy + w_dz

                    penalx = -mask(ix,iy,iz)*(phi(ix,iy,iz,1)-us(ix,iy,iz,1))
                    penaly = -mask(ix,iy,iz)*(phi(ix,iy,iz,2)-us(ix,iy,iz,2))
                    penalz = -mask(ix,iy,iz)*(phi(ix,iy,iz,3)-us(ix,iy,iz,3))

                    rhs(ix,iy,iz,1) = -phi(ix,iy,iz,1)*u_dx - phi(ix,iy,iz,2)*u_dy - phi(ix,iy,iz,3)*u_dz - p_dx &
                                    + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                    rhs(ix,iy,iz,2) = -phi(ix,iy,iz,1)*v_dx - phi(ix,iy,iz,2)*v_dy - phi(ix,iy,iz,3)*v_dz - p_dy &
                                    + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                    rhs(ix,iy,iz,3) = -phi(ix,iy,iz,1)*w_dx - phi(ix,iy,iz,2)*w_dy - phi(ix,iy,iz,3)*w_dz - p_dz &
                                    + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                    rhs(ix,iy,iz,4) = -(c_0**2)*div_U - gamma*phi(ix,iy,iz,4)
                end do
            end do
        end do

    else if (order_discretization == "FD_4th_central_optimized") then
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do iz = g+1, Bs+g
            do iy = g+1, Bs+g
                do ix = g+1, Bs+g
                    ! first derivatives of u, v, p
                    u_dx = (a(-3)*phi(ix-3,iy,iz,1) + a(-2)*phi(ix-2,iy,iz,1) + a(-1)*phi(ix-1,iy,iz,1) + a(0)*phi(ix,iy,iz,1)&
                         +  a(+1)*phi(ix+1,iy,iz,1) + a(+2)*phi(ix+2,iy,iz,1) + a(+3)*phi(ix+3,iy,iz,1))*dx_inv
                    u_dy = (a(-3)*phi(ix,iy-3,iz,1) + a(-2)*phi(ix,iy-2,iz,1) + a(-1)*phi(ix,iy-1,iz,1) + a(0)*phi(ix,iy,iz,1)&
                         +  a(+1)*phi(ix,iy+1,iz,1) + a(+2)*phi(ix,iy+2,iz,1) + a(+3)*phi(ix,iy+3,iz,1))*dy_inv
                    u_dz = (a(-3)*phi(ix,iy,iz-3,1) + a(-2)*phi(ix,iy,iz-2,1) + a(-1)*phi(ix,iy,iz-1,1) + a(0)*phi(ix,iy,iz,1)&
                         +  a(+1)*phi(ix,iy,iz+1,1) + a(+2)*phi(ix,iy,iz+2,1) + a(+3)*phi(ix,iy,iz+3,1))*dz_inv

                    v_dx = (a(-3)*phi(ix-3,iy,iz,2) + a(-2)*phi(ix-2,iy,iz,2) + a(-1)*phi(ix-1,iy,iz,2) + a(0)*phi(ix,iy,iz,2)&
                         +  a(+1)*phi(ix+1,iy,iz,2) + a(+2)*phi(ix+2,iy,iz,2) + a(+3)*phi(ix+3,iy,iz,2))*dx_inv
                    v_dy = (a(-3)*phi(ix,iy-3,iz,2) + a(-2)*phi(ix,iy-2,iz,2) + a(-1)*phi(ix,iy-1,iz,2) + a(0)*phi(ix,iy,iz,2)&
                         +  a(+1)*phi(ix,iy+1,iz,2) + a(+2)*phi(ix,iy+2,iz,2) + a(+3)*phi(ix,iy+3,iz,2))*dy_inv
                    v_dz = (a(-3)*phi(ix,iy,iz-3,2) + a(-2)*phi(ix,iy,iz-2,2) + a(-1)*phi(ix,iy,iz-1,2) + a(0)*phi(ix,iy,iz,2)&
                         +  a(+1)*phi(ix,iy,iz+1,2) + a(+2)*phi(ix,iy,iz+2,2) + a(+3)*phi(ix,iy,iz+3,2))*dz_inv

                    w_dx = (a(-3)*phi(ix-3,iy,iz,3) + a(-2)*phi(ix-2,iy,iz,3) + a(-1)*phi(ix-1,iy,iz,3) + a(0)*phi(ix,iy,iz,3)&
                         +  a(+1)*phi(ix+1,iy,iz,3) + a(+2)*phi(ix+2,iy,iz,3) + a(+3)*phi(ix+3,iy,iz,3))*dx_inv
                    w_dy = (a(-3)*phi(ix,iy-3,iz,3) + a(-2)*phi(ix,iy-2,iz,3) + a(-1)*phi(ix,iy-1,iz,3) + a(0)*phi(ix,iy,iz,3)&
                         +  a(+1)*phi(ix,iy+1,iz,3) + a(+2)*phi(ix,iy+2,iz,3) + a(+3)*phi(ix,iy+3,iz,3))*dy_inv
                    w_dz = (a(-3)*phi(ix,iy,iz-3,3) + a(-2)*phi(ix,iy,iz-2,3) + a(-1)*phi(ix,iy,iz-1,3) + a(0)*phi(ix,iy,iz,3)&
                         +  a(+1)*phi(ix,iy,iz+1,3) + a(+2)*phi(ix,iy,iz+2,3) + a(+3)*phi(ix,iy,iz+3,3))*dz_inv

                    p_dx = (a(-3)*phi(ix-3,iy,iz,4) + a(-2)*phi(ix-2,iy,iz,4) + a(-1)*phi(ix-1,iy,iz,4) + a(0)*phi(ix,iy,iz,4)&
                         +  a(+1)*phi(ix+1,iy,iz,4) + a(+2)*phi(ix+2,iy,iz,4) + a(+3)*phi(ix+3,iy,iz,4))*dx_inv
                    p_dy = (a(-3)*phi(ix,iy-3,iz,4) + a(-2)*phi(ix,iy-2,iz,4) + a(-1)*phi(ix,iy-1,iz,4) + a(0)*phi(ix,iy,iz,4)&
                         +  a(+1)*phi(ix,iy+1,iz,4) + a(+2)*phi(ix,iy+2,iz,4) + a(+3)*phi(ix,iy+3,iz,4))*dy_inv
                    p_dz = (a(-3)*phi(ix,iy,iz-3,4) + a(-2)*phi(ix,iy,iz-2,4) + a(-1)*phi(ix,iy,iz-1,4) + a(0)*phi(ix,iy,iz,4)&
                         +  a(+1)*phi(ix,iy,iz+1,4) + a(+2)*phi(ix,iy,iz+2,4) + a(+3)*phi(ix,iy,iz+3,4))*dz_inv

                    ! second derivatives of u, v and w
                    u_dxdx = (b(-2)*phi(ix-2,iy,iz,1) + b(-1)*phi(ix-1,iy,iz,1) + b(0)*phi(ix,iy,iz,1)&
                           +  b(+1)*phi(ix+1,iy,iz,1) + b(+2)*phi(ix+2,iy,iz,1))*dx2_inv
                    u_dydy = (b(-2)*phi(ix,iy-2,iz,1) + b(-1)*phi(ix,iy-1,iz,1) + b(0)*phi(ix,iy,iz,1)&
                           +  b(+1)*phi(ix,iy+1,iz,1) + b(+2)*phi(ix,iy+2,iz,1))*dy2_inv
                    u_dzdz = (b(-2)*phi(ix,iy,iz-2,1) + b(-1)*phi(ix,iy,iz-1,1) + b(0)*phi(ix,iy,iz,1)&
                           +  b(+1)*phi(ix,iy,iz+1,1) + b(+2)*phi(ix,iy,iz+2,1))*dz2_inv

                    v_dxdx = (b(-2)*phi(ix-2,iy,iz,2) + b(-1)*phi(ix-1,iy,iz,2) + b(0)*phi(ix,iy,iz,2)&
                           +  b(+1)*phi(ix+1,iy,iz,2) + b(+2)*phi(ix+2,iy,iz,2))*dx2_inv
                    v_dydy = (b(-2)*phi(ix,iy-2,iz,2) + b(-1)*phi(ix,iy-1,iz,2) + b(0)*phi(ix,iy,iz,2)&
                           +  b(+1)*phi(ix,iy+1,iz,2) + b(+2)*phi(ix,iy+2,iz,2))*dy2_inv
                    v_dzdz = (b(-2)*phi(ix,iy,iz-2,2) + b(-1)*phi(ix,iy,iz-1,2) + b(0)*phi(ix,iy,iz,2)&
                           +  b(+1)*phi(ix,iy,iz+1,2) + b(+2)*phi(ix,iy,iz+2,2))*dz2_inv

                    w_dxdx = (b(-2)*phi(ix-2,iy,iz,3) + b(-1)*phi(ix-1,iy,iz,3) + b(0)*phi(ix,iy,iz,3)&
                           +  b(+1)*phi(ix+1,iy,iz,3) + b(+2)*phi(ix+2,iy,iz,3))*dx2_inv
                    w_dydy = (b(-2)*phi(ix,iy-2,iz,3) + b(-1)*phi(ix,iy-1,iz,3) + b(0)*phi(ix,iy,iz,3)&
                           +  b(+1)*phi(ix,iy+1,iz,3) + b(+2)*phi(ix,iy+2,iz,3))*dy2_inv
                    w_dzdz = (b(-2)*phi(ix,iy,iz-2,3) + b(-1)*phi(ix,iy,iz-1,3) + b(0)*phi(ix,iy,iz,3)&
                           +  b(+1)*phi(ix,iy,iz+1,3) + b(+2)*phi(ix,iy,iz+2,3))*dz2_inv

                    div_U = u_dx + v_dy + w_dz

                    penalx = -mask(ix,iy,iz)*(phi(ix,iy,iz,1)-us(ix,iy,iz,1))
                    penaly = -mask(ix,iy,iz)*(phi(ix,iy,iz,2)-us(ix,iy,iz,2))
                    penalz = -mask(ix,iy,iz)*(phi(ix,iy,iz,3)-us(ix,iy,iz,3))

                    rhs(ix,iy,iz,1) = -phi(ix,iy,iz,1)*u_dx - phi(ix,iy,iz,2)*u_dy - phi(ix,iy,iz,3)*u_dz - p_dx &
                                    + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                    rhs(ix,iy,iz,2) = -phi(ix,iy,iz,1)*v_dx - phi(ix,iy,iz,2)*v_dy - phi(ix,iy,iz,3)*v_dz - p_dy &
                                    + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                    rhs(ix,iy,iz,3) = -phi(ix,iy,iz,1)*w_dx - phi(ix,iy,iz,2)*w_dy - phi(ix,iy,iz,3)*w_dz - p_dz &
                                    + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                    rhs(ix,iy,iz,4) = -(c_0**2)*div_U - gamma*phi(ix,iy,iz,4)
                end do
            end do
        end do

    else
        call abort(441167, "3d Discretization unkown "//order_discretization//", I ll walk into the light now." )
    end if

    ! remove mean pressure. NOTE: there is some oddities here, as the code modifies the
    ! state vector and not the RHS.
    if (params_acm%p_mean_zero) then
        phi(:,:,:,4) = phi(:,:,:,4) - params_acm%mean_p
    end if

    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    if (params_acm%use_sponge) then
        call sponge_3D(sponge, x0, dx, Bs, g)
        eps_inv = 1.0_rk / params_acm%C_sponge
        sponge = sponge * eps_inv

        ! NOTE: the sponge term acts, if active, on ALL components, ux,uy,p
        ! which is different from the penalization term, which acts only on ux,uy and not p
        rhs(:,:,:,1) = rhs(:,:,:,1) - (phi(:,:,:,1)-params_acm%u_mean_set(1)) * sponge
        rhs(:,:,:,2) = rhs(:,:,:,2) - (phi(:,:,:,2)-params_acm%u_mean_set(2)) * sponge
        rhs(:,:,:,3) = rhs(:,:,:,3) - (phi(:,:,:,3)-params_acm%u_mean_set(3)) * sponge
        rhs(:,:,:,4) = rhs(:,:,:,4) - phi(:,:,:,4)*sponge
    end if

end subroutine RHS_3D_acm
