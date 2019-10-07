!-----------------------------------------------------------------------------
! main level wrapper to set the right hand side on a block. Note this is completely
! independent of the grid and any MPI formalism, neighboring relations and the like.
! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
! from that. Ghost nodes are assumed to be sync'ed.
!-----------------------------------------------------------------------------
subroutine RHS_ACM( time, u, g, x0, dx, rhs, mask, stage )
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

    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(in) :: mask(1:,1:,1:,1:)


    ! local variables
    integer(kind=ik) :: mpierr, i, dim
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: tmp(1:3), tmp2

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    dim = params_acm%dim

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
        do i = 1, size(u,4)
            if (maxval(abs(u(:,:,:,i))) > 1.0e4_rk) then
                write(*,'("maxval in u(:,:,:,",i2,") = ", es15.8)') i, maxval(abs(u(:,:,:,i)))
                call abort(0409201933,"ACM fail: very very large values in state vector.")
            endif
        enddo

        if (params_acm%dim == 2) then
            if (params_acm%u_mean_zero .or. params_acm%forcing) then
                params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, 1))*dx(1)*dx(2)
                params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, 2))*dx(1)*dx(2)
            endif
            if (params_acm%p_mean_zero) then
                params_acm%mean_p = params_acm%mean_p + sum(u(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, 3))*dx(1)*dx(2)
            endif

        else
            if (params_acm%u_mean_zero .or. params_acm%forcing) then
                params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 1))*dx(1)*dx(2)*dx(3)
                params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 2))*dx(1)*dx(2)*dx(3)
                params_acm%mean_flow(3) = params_acm%mean_flow(3) + sum(u(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 3))*dx(1)*dx(2)*dx(3)
            endif
            if (params_acm%p_mean_zero) then
                params_acm%mean_p = params_acm%mean_p + sum(u(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 4))*dx(1)*dx(2)*dx(3)
            endif
        endif ! NOTE: MPI_SUM is perfomed in the post_stage.

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.

        if (params_acm%u_mean_zero .or. params_acm%forcing) then
            tmp = params_acm%mean_flow
            call MPI_ALLREDUCE(tmp, params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            params_acm%mean_flow = params_acm%mean_flow / product(params_acm%domain_size(1:dim))
        endif
        if (params_acm%p_mean_zero) then
            tmp2 = params_acm%mean_p
            call MPI_ALLREDUCE(tmp2, params_acm%mean_p, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            params_acm%mean_p = params_acm%mean_p / product(params_acm%domain_size(1:dim))
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
            ! --------------------------------------------------------------------------
            ! flow
            ! --------------------------------------------------------------------------
            if (params_acm%compute_flow) then
                ! this is a 2d case (ux,uy,p)
                call RHS_2D_acm(g, Bs, dx(1:2), x0(1:2), u(:,:,1,:), params_acm%discretization, &
                time, rhs(:,:,1,:), mask(:,:,1,:))
            endif
            ! --------------------------------------------------------------------------
            ! passive scalars
            ! --------------------------------------------------------------------------
            if (params_acm%use_passive_scalar) then
                call RHS_2D_scalar(g, Bs, dx, x0, u, params_acm%discretization, time, rhs, mask)
            endif

        else
            ! --------------------------------------------------------------------------
            ! flow
            ! --------------------------------------------------------------------------
            if (params_acm%compute_flow) then
                ! this is a 3d case (ux,uy,uz,p)
                call RHS_3D_acm(g, Bs, dx, x0, u, params_acm%discretization, time, rhs, mask)
            endif

            ! --------------------------------------------------------------------------
            ! passive scalars
            ! --------------------------------------------------------------------------
            if (params_acm%use_passive_scalar) then
                call RHS_3D_scalar(g, Bs, dx, x0, u, params_acm%discretization, time, rhs, mask)
            endif

        endif

    case default
        call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")

    end select


end subroutine RHS_ACM






subroutine RHS_2D_acm(g, Bs, dx, x0, phi, order_discretization, time, rhs, mask)

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)            :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(2), intent(in) :: x0, dx
    !> datafields
    real(kind=rk), intent(inout)            :: phi(:,:,:)
    real(kind=rk), intent(inout)            :: rhs(:,:,:)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(in)               :: mask(:,:,:)
    !> discretization order
    character(len=80), intent(in)           :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    !> forcing term
    real(kind=rk), dimension(3) :: forcing
    !>
    real(kind=rk) :: dx_inv, dy_inv, dx2_inv, dy2_inv, c_0, nu, eps, eps_inv, gamma
    real(kind=rk) :: div_U, u_dx, u_dy, u_dxdx, u_dydy, v_dx, v_dy, v_dxdx, &
                     v_dydy, p_dx, p_dy, penalx, penaly, x, y, term_2, spo
    ! loop variables
    integer(kind=rk) :: ix, iy, idir
    ! coefficients for Tam&Webb
    real(kind=rk) :: a(-3:3)
    real(kind=rk) :: b(-2:2)


    ! set parameters for readability
    c_0         = params_acm%c_0
    nu          = params_acm%nu
    eps         = params_acm%C_eta
    gamma       = params_acm%gamma_p

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    eps_inv = 1.0_rk / eps

    if (size(phi,1)/=Bs(1)+2*g .or. size(phi,2)/=Bs(2)+2*g .or. size(phi,3)/=params_acm%dim+1+params_acm%N_scalars) then
        call abort(66233,"wrong size, I go for a walk instead.")
    endif

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    ! 4th order coefficients for second derivative
    b = (/ -1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)

    if (maxval(abs(params_acm%mean_flow)) > 1.0e3_rk ) then
        call abort(887, "ACM: meanflow out of bounds")
    endif




    if (order_discretization == "FD_2nd_central" ) then
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g

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

                penalx = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                penaly = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,2) -mask(ix,iy,3))

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
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
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

                penalx = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                penaly = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,2) -mask(ix,iy,3))

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
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            x = x0(1) + dble(ix-g-1) * dx(1)
                            y = x0(2) + dble(iy-g-1) * dx(2)
                            call continue_periodic(x,params_acm%domain_size(1))
                            call continue_periodic(y,params_acm%domain_size(2))
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

    ! remove mean flow
    if (params_acm%u_mean_zero) then
        phi(:,:,1) = phi(:,:,1) - params_acm%mean_flow(1)
        phi(:,:,2) = phi(:,:,2) - params_acm%mean_flow(2)
    end if

    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    if (params_acm%use_sponge) then
       !  ! create the mask for the sponge on this block
       ! call sponge_2D(sponge, x0, dx, Bs, g)

        ! avoid division by multiplying with inverse
        eps_inv = 1.0_rk / params_acm%C_sponge

        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! NOTE: the sponge term acts, if active, on ALL components, ux,uy,p
                ! which is different from the penalization term, which acts only on ux,uy and not p
                spo = mask(ix,iy,6) * eps_inv

                rhs(ix,iy,1) = rhs(ix,iy,1) - (phi(ix,iy,1)-params_acm%u_mean_set(1)) * spo
                rhs(ix,iy,2) = rhs(ix,iy,2) - (phi(ix,iy,2)-params_acm%u_mean_set(2)) * spo
                rhs(ix,iy,3) = rhs(ix,iy,3) - phi(ix,iy,3)*spo
            enddo
        enddo
    end if


end subroutine RHS_2D_acm




subroutine RHS_3D_acm(g, Bs, dx, x0, phi, order_discretization, time, rhs, mask)
    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)            :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in) :: x0, dx
    !> datafields
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)
    real(kind=rk), intent(inout)            :: rhs(:,:,:,:)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(in)               :: mask(:,:,:,:)
    !> discretization order
    character(len=80), intent(in)           :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time


    !> forcing term
    real(kind=rk), dimension(3) :: forcing

    !> inverse dx, physics/acm parameters
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, c_0, &
                     nu, eps, eps_inv, gamma, spo
    !> derivatives
    real(kind=rk) :: div_U, u_dx, u_dy, u_dz, u_dxdx, u_dydy, u_dzdz, &
                     v_dx, v_dy, v_dz, v_dxdx, v_dydy, v_dzdz, &
                     w_dx, w_dy, w_dz, w_dxdx, w_dydy, w_dzdz, &
                     p_dx, p_dy, p_dz, penalx, penaly, penalz, u, v, w, p, chi
    !> loop variables
    integer(kind=rk) :: ix, iy, iz
    !> coefficients for Tam&Webb
    real(kind=rk), parameter :: a(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)


    ! set parameters for readability
    c_0         = params_acm%c_0
    nu          = params_acm%nu
    eps         = params_acm%C_eta
    gamma       = params_acm%gamma_p

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dz_inv = 1.0_rk / dx(3)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    dz2_inv = 1.0_rk / (dx(3)**2)

    eps_inv = 1.0_rk / eps



    if (order_discretization == "FD_2nd_central" ) then
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first and second derivatives of u,v,w
                    u_dx = (phi(ix+1, iy,   iz  , 1) - phi(ix-1, iy,   iz  , 1))*dx_inv*0.5_rk
                    u_dy = (phi(ix  , iy+1, iz  , 1) - phi(ix,   iy-1, iz  , 1))*dy_inv*0.5_rk
                    u_dz = (phi(ix  , iy,   iz+1, 1) - phi(ix,   iy,   iz-1, 1))*dz_inv*0.5_rk

                    v_dx = (phi(ix+1, iy  , iz  , 2) - phi(ix-1, iy  , iz  , 2))*dx_inv*0.5_rk
                    v_dy = (phi(ix  , iy+1, iz  , 2) - phi(ix  , iy-1, iz  , 2))*dy_inv*0.5_rk
                    v_dz = (phi(ix  , iy  , iz+1, 2) - phi(ix  , iy  , iz-1, 2))*dz_inv*0.5_rk

                    w_dx = (phi(ix+1, iy  , iz  , 3) - phi(ix-1, iy  , iz  , 3))*dx_inv*0.5_rk
                    w_dy = (phi(ix  , iy+1, iz  , 3) - phi(ix  , iy-1, iz  , 3))*dy_inv*0.5_rk
                    w_dz = (phi(ix  , iy  , iz+1, 3) - phi(ix  , iy  , iz-1, 3))*dz_inv*0.5_rk

                    u_dxdx = (phi(ix-1, iy  , iz  , 1) -2.0_rk*phi(ix, iy, iz, 1) + phi(ix+1, iy  , iz  , 1))*dx2_inv
                    u_dydy = (phi(ix  , iy-1, iz  , 1) -2.0_rk*phi(ix, iy, iz, 1) + phi(ix  , iy+1, iz  , 1))*dy2_inv
                    u_dzdz = (phi(ix  , iy  , iz-1, 1) -2.0_rk*phi(ix, iy, iz, 1) + phi(ix  , iy  , iz+1, 1))*dz2_inv

                    v_dxdx = (phi(ix-1, iy  , iz  , 2) -2.0_rk*phi(ix, iy, iz, 2) + phi(ix+1, iy  , iz  , 2))*dx2_inv
                    v_dydy = (phi(ix  , iy-1, iz  , 2) -2.0_rk*phi(ix, iy, iz, 2) + phi(ix  , iy+1, iz  , 2))*dy2_inv
                    v_dzdz = (phi(ix  , iy  , iz-1, 2) -2.0_rk*phi(ix, iy, iz, 2) + phi(ix  , iy  , iz+1, 2))*dz2_inv

                    w_dxdx = (phi(ix-1, iy  , iz  , 3) -2.0_rk*phi(ix, iy, iz, 3) + phi(ix+1, iy  , iz  , 3))*dx2_inv
                    w_dydy = (phi(ix  , iy-1, iz  , 3) -2.0_rk*phi(ix, iy, iz, 3) + phi(ix  , iy+1, iz  , 3))*dy2_inv
                    w_dzdz = (phi(ix  , iy  , iz-1, 3) -2.0_rk*phi(ix, iy, iz, 3) + phi(ix  , iy  , iz+1, 3))*dz2_inv

                    ! first derivative of p
                    p_dx = (phi(ix+1, iy  , iz  , 4) -phi(ix-1, iy  , iz  , 4))*dx_inv*0.5_rk
                    p_dy = (phi(ix  , iy+1, iz  , 4) -phi(ix  , iy-1, iz  , 4))*dy_inv*0.5_rk
                    p_dz = (phi(ix  , iy  , iz+1, 4) -phi(ix  , iy  , iz-1, 4))*dz_inv*0.5_rk

                    div_U = u_dx + v_dy + w_dz

                    u = phi(ix, iy, iz, 1)
                    v = phi(ix, iy, iz, 2)
                    w = phi(ix, iy, iz, 3)
                    p = phi(ix, iy, iz, 4)

                    chi = mask(ix,iy,iz,1) * eps_inv
                    penalx = -chi * (u - mask(ix,iy,iz,2))
                    penaly = -chi * (v - mask(ix,iy,iz,3))
                    penalz = -chi * (w - mask(ix,iy,iz,4))

                    rhs(ix,iy,iz,1) = (-u*u_dx - v*u_dy - w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                    rhs(ix,iy,iz,2) = (-u*v_dx - v*v_dy - w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                    rhs(ix,iy,iz,3) = (-u*w_dx - v*w_dy - w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                    rhs(ix,iy,iz,4) = -(c_0**2)*div_U - gamma*p
                end do
            end do
        end do

    else if (order_discretization == "FD_4th_central_optimized") then
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    u_dx = (a(-3)*phi(ix-3,iy,iz,1) + a(-2)*phi(ix-2,iy,iz,1) + a(-1)*phi(ix-1,iy,iz,1) + a(0)*phi(ix,iy,iz,1) &
                         +  a(+1)*phi(ix+1,iy,iz,1) + a(+2)*phi(ix+2,iy,iz,1) + a(+3)*phi(ix+3,iy,iz,1))*dx_inv
                    u_dy = (a(-3)*phi(ix,iy-3,iz,1) + a(-2)*phi(ix,iy-2,iz,1) + a(-1)*phi(ix,iy-1,iz,1) + a(0)*phi(ix,iy,iz,1) &
                         +  a(+1)*phi(ix,iy+1,iz,1) + a(+2)*phi(ix,iy+2,iz,1) + a(+3)*phi(ix,iy+3,iz,1))*dy_inv
                    u_dz = (a(-3)*phi(ix,iy,iz-3,1) + a(-2)*phi(ix,iy,iz-2,1) + a(-1)*phi(ix,iy,iz-1,1) + a(0)*phi(ix,iy,iz,1) &
                         +  a(+1)*phi(ix,iy,iz+1,1) + a(+2)*phi(ix,iy,iz+2,1) + a(+3)*phi(ix,iy,iz+3,1))*dz_inv

                    v_dx = (a(-3)*phi(ix-3,iy,iz,2) + a(-2)*phi(ix-2,iy,iz,2) + a(-1)*phi(ix-1,iy,iz,2) + a(0)*phi(ix,iy,iz,2) &
                         +  a(+1)*phi(ix+1,iy,iz,2) + a(+2)*phi(ix+2,iy,iz,2) + a(+3)*phi(ix+3,iy,iz,2))*dx_inv
                    v_dy = (a(-3)*phi(ix,iy-3,iz,2) + a(-2)*phi(ix,iy-2,iz,2) + a(-1)*phi(ix,iy-1,iz,2) + a(0)*phi(ix,iy,iz,2) &
                         +  a(+1)*phi(ix,iy+1,iz,2) + a(+2)*phi(ix,iy+2,iz,2) + a(+3)*phi(ix,iy+3,iz,2))*dy_inv
                    v_dz = (a(-3)*phi(ix,iy,iz-3,2) + a(-2)*phi(ix,iy,iz-2,2) + a(-1)*phi(ix,iy,iz-1,2) + a(0)*phi(ix,iy,iz,2) &
                         +  a(+1)*phi(ix,iy,iz+1,2) + a(+2)*phi(ix,iy,iz+2,2) + a(+3)*phi(ix,iy,iz+3,2))*dz_inv

                    w_dx = (a(-3)*phi(ix-3,iy,iz,3) + a(-2)*phi(ix-2,iy,iz,3) + a(-1)*phi(ix-1,iy,iz,3) + a(0)*phi(ix,iy,iz,3) &
                         +  a(+1)*phi(ix+1,iy,iz,3) + a(+2)*phi(ix+2,iy,iz,3) + a(+3)*phi(ix+3,iy,iz,3))*dx_inv
                    w_dy = (a(-3)*phi(ix,iy-3,iz,3) + a(-2)*phi(ix,iy-2,iz,3) + a(-1)*phi(ix,iy-1,iz,3) + a(0)*phi(ix,iy,iz,3) &
                         +  a(+1)*phi(ix,iy+1,iz,3) + a(+2)*phi(ix,iy+2,iz,3) + a(+3)*phi(ix,iy+3,iz,3))*dy_inv
                    w_dz = (a(-3)*phi(ix,iy,iz-3,3) + a(-2)*phi(ix,iy,iz-2,3) + a(-1)*phi(ix,iy,iz-1,3) + a(0)*phi(ix,iy,iz,3) &
                         +  a(+1)*phi(ix,iy,iz+1,3) + a(+2)*phi(ix,iy,iz+2,3) + a(+3)*phi(ix,iy,iz+3,3))*dz_inv

                    p_dx = (a(-3)*phi(ix-3,iy,iz,4) + a(-2)*phi(ix-2,iy,iz,4) + a(-1)*phi(ix-1,iy,iz,4) + a(0)*phi(ix,iy,iz,4) &
                         +  a(+1)*phi(ix+1,iy,iz,4) + a(+2)*phi(ix+2,iy,iz,4) + a(+3)*phi(ix+3,iy,iz,4))*dx_inv
                    p_dy = (a(-3)*phi(ix,iy-3,iz,4) + a(-2)*phi(ix,iy-2,iz,4) + a(-1)*phi(ix,iy-1,iz,4) + a(0)*phi(ix,iy,iz,4) &
                         +  a(+1)*phi(ix,iy+1,iz,4) + a(+2)*phi(ix,iy+2,iz,4) + a(+3)*phi(ix,iy+3,iz,4))*dy_inv
                    p_dz = (a(-3)*phi(ix,iy,iz-3,4) + a(-2)*phi(ix,iy,iz-2,4) + a(-1)*phi(ix,iy,iz-1,4) + a(0)*phi(ix,iy,iz,4) &
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

                    u = phi(ix, iy, iz, 1)
                    v = phi(ix, iy, iz, 2)
                    w = phi(ix, iy, iz, 3)
                    p = phi(ix, iy, iz, 4)

                    chi = mask(ix,iy,iz,1) * eps_inv
                    penalx = -chi * (u - mask(ix,iy,iz,2))
                    penaly = -chi * (v - mask(ix,iy,iz,3))
                    penalz = -chi * (w - mask(ix,iy,iz,4))

                    rhs(ix,iy,iz,1) = (-u*u_dx - v*u_dy - w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                    rhs(ix,iy,iz,2) = (-u*v_dx - v*v_dy - w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                    rhs(ix,iy,iz,3) = (-u*w_dx - v*w_dy - w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                    rhs(ix,iy,iz,4) = -(c_0**2)*div_U - gamma*p
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

    ! remove mean flow
    if (params_acm%u_mean_zero) then
        phi(:,:,:,1) = phi(:,:,:,1) - params_acm%mean_flow(1)
        phi(:,:,:,2) = phi(:,:,:,2) - params_acm%mean_flow(2)
        phi(:,:,:,3) = phi(:,:,:,3) - params_acm%mean_flow(3)
    end if

    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    if (params_acm%use_sponge) then
        ! create the mask for the sponge on this block
       ! call sponge_3D(sponge, x0, dx, Bs, g)

        ! avoid division by multiplying with inverse
        eps_inv = 1.0_rk / params_acm%C_sponge

        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! NOTE: the sponge term acts, if active, on ALL components, ux,uy,p
                    ! which is different from the penalization term, which acts only on ux,uy and not p
                    ! NOTE: sponge mask set in grid_qty
                    spo = mask(ix,iy,iz,6) * eps_inv

                    rhs(ix,iy,iz,1) = rhs(ix,iy,iz,1) - (phi(ix,iy,iz,1)-params_acm%u_mean_set(1)) * spo
                    rhs(ix,iy,iz,2) = rhs(ix,iy,iz,2) - (phi(ix,iy,iz,2)-params_acm%u_mean_set(2)) * spo
                    rhs(ix,iy,iz,3) = rhs(ix,iy,iz,3) - (phi(ix,iy,iz,3)-params_acm%u_mean_set(3)) * spo
                    rhs(ix,iy,iz,4) = rhs(ix,iy,iz,4) - (phi(ix,iy,iz,4))*spo
                end do
            end do
        end do
    end if



end subroutine RHS_3D_acm




subroutine RHS_3D_scalar(g, Bs, dx, x0, phi, order_discretization, time, rhs, mask)
    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)            :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in) :: x0, dx
    !> datafields
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)
    real(kind=rk), intent(inout)            :: rhs(:,:,:,:)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(in)               :: mask(:,:,:,:)
    !> discretization order
    character(len=80), intent(in)           :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    integer(kind=rk) :: ix, iy, iz, iscalar, j

    !> coefficients for Tam&Webb
    real(kind=rk), parameter :: a(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    real(kind=rk) :: kappa, x, y, z, masksource, nu, R, R0sq
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv
    real(kind=rk) :: ux,uy,uz,&
    usx,usy,usz,wx,wy,wz,gx,gy,gz,D,chi,chidx,chidz,chidy,D_dx,D_dy,D_dz,gxx,gyy,gzz
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're ngeligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: source(:,:,:)

    if (.not. allocated(source)) allocate(source(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g))

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dz_inv = 1.0_rk / dx(3)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    dz2_inv = 1.0_rk / (dx(3)**2)

    nu = params_acm%nu

    if (order_discretization == "FD_2nd_central" ) then
        call abort(2208191, "passive scalar implemented only with 4th order....sorry, I am lazy")
    else if (order_discretization == "FD_4th_central_optimized") then
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do iscalar = 1, params_acm%N_scalars
            ! actual index of this scalar in the array
            j = iscalar + (params_acm%dim + 1)

            ! compute diffusivity from schmidt number (and fluid viscosity)
            kappa = params_acm%schmidt_numbers(iscalar) * nu

            ! reset source term for each scalar.
            source = 0.0_rk

            ! 1st: compute source terms (note the strcmp needs to be outside the loop)
            select case (params_acm%scalar_source_type(iscalar))
            case ("gaussian")
                do iz = g+1, Bs(3)+g
                    z = (x0(3) + dble(iz-g-1)*dx(3) - params_acm%z0source(iscalar))**2
                    do iy = g+1, Bs(2)+g
                        y = (x0(2) + dble(iy-g-1)*dx(2) - params_acm%y0source(iscalar))**2
                        do ix = g+1, Bs(1)+g
                            x = (x0(1) + dble(ix-g-1)*dx(1) - params_acm%x0source(iscalar))**2

                            R = x + y + z ! note this is (x-x0)**2

                            masksource = dexp( -R / (params_acm%widthsource(iscalar)**2)  )

                            if (masksource > 1.0d-6) then
                                ! for the source term, we use the usual dirichlet C_eta
                                ! to force scalar to 1
                                ! source(ix,iy,iz) = -1.0d0*(phi(ix,iy,iz,j)-masksource-) / params_acm%C_eta
                                source(ix,iy,iz) = (masksource - phi(ix,iy,iz,j)) / params_acm%C_eta
                                ! source(ix,iy,iz) = -masksource*(phi(ix,iy,iz,j)-1.d0) / params_acm%C_eta
                            endif
                        end do
                    end do
                end do
            case ("circular")
                R0sq = params_acm%widthsource(iscalar)**2
                do iz = g+1, Bs(3)+g
                    z = (x0(3) + dble(iz-g-1)*dx(3) - params_acm%z0source(iscalar))**2
                    do iy = g+1, Bs(2)+g
                        y = (x0(2) + dble(iy-g-1)*dx(2) - params_acm%y0source(iscalar))**2
                        do ix = g+1, Bs(1)+g
                            x = (x0(1) + dble(ix-g-1)*dx(1) - params_acm%x0source(iscalar))**2

                            R = x + y + z ! note this is (x-x0)**2

                            if ( R <= R0sq ) then
                                ! for the source term, we use the usual dirichlet C_eta
                                ! to force scalar to 1
                                source(ix,iy,iz) = -(phi(ix,iy,iz,j)-1.d0) / params_acm%C_eta
                            endif
                        end do
                    end do
                end do

            case ("mask_color_emission")
                call abort(26081919,"lazy tommy not done yet")

            case ("none", "empty")
                ! do nothing

            case default
                call abort(2608191,"scalar source is unkown.")

            end select


            ! sponge layer
            if (params_acm%absorbing_sponge) then
                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            ! for the source term, we use the usual dirichlet C_eta
                            ! to force scalar to 1
                            source(ix,iy,iz) = source(ix,iy,iz) - mask(ix,iy,iz,6)*phi(ix,iy,iz,j) / params_acm%C_eta
                        end do
                    end do
                enddo
            endif


            ! 2nd: compute rhs for this scalar.
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ux = phi(ix,iy,iz,1)
                        uy = phi(ix,iy,iz,2)
                        uz = phi(ix,iy,iz,3)

                        usx = mask(ix,iy,iz,2)
                        usy = mask(ix,iy,iz,3)
                        usz = mask(ix,iy,iz,4)

                        ! ATTENTION you need to sync the mask
                        chi = mask(ix,iy,iz,1)

                        ! penalized diffusion coefficient at this point
                        D = kappa*(1.0_rk - chi) + params_acm%scalar_Ceta(iscalar)*chi

                        ! this is the vector in front of the gradient
                        wx = -((1.d0-chi)*ux + chi*usx)
                        wy = -((1.d0-chi)*uy + chi*usy)
                        wz = -((1.d0-chi)*uz + chi*usz)

                        ! gradient of passive scalar
                        gx = (a(-3)*phi(ix-3,iy,iz,j)&
                             +a(-2)*phi(ix-2,iy,iz,j)&
                             +a(-1)*phi(ix-1,iy,iz,j)&
                             +a( 0)*phi(ix  ,iy,iz,j)&
                             +a(+2)*phi(ix+2,iy,iz,j)&
                             +a(+3)*phi(ix+3,iy,iz,j)&
                             +a(+1)*phi(ix+1,iy,iz,j))*dx_inv

                        gy = (a(-3)*phi(ix,iy-3,iz,j)&
                             +a(-2)*phi(ix,iy-2,iz,j)&
                             +a(-1)*phi(ix,iy-1,iz,j)&
                             +a( 0)*phi(ix,iy  ,iz,j)&
                             +a(+2)*phi(ix,iy+2,iz,j)&
                             +a(+3)*phi(ix,iy+3,iz,j)&
                             +a(+1)*phi(ix,iy+1,iz,j))*dy_inv

                        gz = (a(-3)*phi(ix,iy,iz-3,j)&
                             +a(-2)*phi(ix,iy,iz-2,j)&
                             +a(-1)*phi(ix,iy,iz-1,j)&
                             +a( 0)*phi(ix,iy,iz  ,j)&
                             +a(+2)*phi(ix,iy,iz+2,j)&
                             +a(+3)*phi(ix,iy,iz+3,j)&
                             +a(+1)*phi(ix,iy,iz+1,j))*dz_inv

                        ! gradient of mask function ( we need that for the diffusive term)
                        ! since this guy reads div( (kappa(1-mask) + eps*mask) * grad(phi) )
                        ! so this boils down to d/dx (D*gx) = D_dx*gx + D*gxx
                        ! so we need D_dx and this is kappa*(1-mask_dx)+ eps*mask_dx
                        chidx = (a(-3)*mask(ix-3,iy,iz, 1)&
                                +a(-2)*mask(ix-2,iy,iz, 1)&
                                +a(-1)*mask(ix-1,iy,iz, 1)&
                                +a( 0)*mask(ix  ,iy,iz, 1)&
                                +a(+3)*mask(ix+3,iy,iz, 1)&
                                +a(+2)*mask(ix+2,iy,iz, 1)&
                                +a(+1)*mask(ix+1,iy,iz, 1))*dx_inv

                        chidy = (a(-3)*mask(ix,iy-3,iz, 1)&
                                +a(-2)*mask(ix,iy-2,iz, 1)&
                                +a(-1)*mask(ix,iy-1,iz, 1)&
                                +a( 0)*mask(ix,iy  ,iz, 1)&
                                +a(+3)*mask(ix,iy+3,iz, 1)&
                                +a(+2)*mask(ix,iy+2,iz, 1)&
                                +a(+1)*mask(ix,iy+1,iz, 1))*dy_inv

                        chidz = (a(-3)*mask(ix,iy,iz-3, 1)&
                                +a(-2)*mask(ix,iy,iz-2, 1)&
                                +a(-1)*mask(ix,iy,iz-1, 1)&
                                +a( 0)*mask(ix,iy,iz  , 1)&
                                +a(+3)*mask(ix,iy,iz+3, 1)&
                                +a(+2)*mask(ix,iy,iz+2, 1)&
                                +a(+1)*mask(ix,iy,iz+1, 1))*dz_inv

                        D_dx = kappa*(1.0_rk-chidx) + params_acm%scalar_Ceta(iscalar) * chidx
                        D_dy = kappa*(1.0_rk-chidy) + params_acm%scalar_Ceta(iscalar) * chidy
                        D_dz = kappa*(1.0_rk-chidz) + params_acm%scalar_Ceta(iscalar) * chidz

                        ! second derivatives of passive scalar
                        gxx = (b(-2)*phi(ix-2,iy,iz ,j)&
                              +b(-1)*phi(ix-1,iy,iz ,j)&
                              +b( 0)*phi(ix  ,iy,iz ,j)&
                              +b(+1)*phi(ix+1,iy,iz ,j)&
                              +b(+2)*phi(ix+2,iy,iz ,j))*dx2_inv
                        gyy = (b(-2)*phi(ix,iy-2,iz ,j)&
                              +b(-1)*phi(ix,iy-1,iz ,j)&
                              +b( 0)*phi(ix,iy  ,iz ,j)&
                              +b(+1)*phi(ix,iy+1,iz ,j)&
                              +b(+2)*phi(ix,iy+2,iz ,j))*dy2_inv
                        gzz = (b(-2)*phi(ix,iy,iz-2 ,j)&
                              +b(-1)*phi(ix,iy,iz-1 ,j)&
                              +b( 0)*phi(ix,iy,iz   ,j)&
                              +b(+1)*phi(ix,iy,iz+1 ,j)&
                              +b(+2)*phi(ix,iy,iz+2 ,j))*dz2_inv

                        ! assemble everything
                        rhs(ix,iy,iz,j) = wx*gx + wy*gy + wz*gz & ! penalized convection term
                        + D_dx*gx + D*gxx & ! penalized laplacian
                        + D_dy*gy + D*gyy &
                        + D_dz*gz + D*gzz &
                        + source(ix,iy,iz)
                    end do
                end do
            end do
        end do ! loop over scalars


    else
        call abort(441167, "3d Discretization unkown "//order_discretization//", I ll walk into the light now." )
    end if
end subroutine RHS_3D_scalar


subroutine RHS_2D_scalar(g, Bs, dx, x0, phi, order_discretization, time, rhs, mask)
    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)            :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in) :: x0, dx
    !> datafields
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)
    real(kind=rk), intent(inout)            :: rhs(:,:,:,:)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(in)               :: mask(:,:,:,:)
    !> discretization order
    character(len=80), intent(in)           :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    integer(kind=rk) :: ix, iy, iscalar, j

    !> coefficients for Tam&Webb
    real(kind=rk), parameter :: a(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    real(kind=rk) :: kappa, x, y, masksource, nu, R
    real(kind=rk) :: dx_inv, dy_inv, dx2_inv, dy2_inv
    real(kind=rk) :: ux, uy, usx, usy, wx, wy, gx, gy, D, chi, chidx, chidy, D_dx, D_dy, gxx, gyy
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're ngeligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: source(:,:,:)

    if (.not. allocated(source)) allocate(source(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
    source = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    nu = params_acm%nu

    if (order_discretization == "FD_2nd_central" ) then
        call abort(2208199, "passive scalar implemented only with 4th order....sorry, I am lazy")

    else if (order_discretization == "FD_4th_central_optimized") then
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do iscalar = 1, params_acm%N_scalars
            ! actual index of this scalar in the array
            j = iscalar + (params_acm%dim + 1)

            ! compute diffusivity from schmidt number (and of course fluid viscosity)
            kappa = nu / params_acm%schmidt_numbers(iscalar)

            source = 0.0_rk

            ! 1st: compute source terms (note the strcmp needs to be outside the loop)
            select case (params_acm%scalar_source_type(iscalar))
            case ("gaussian")
                do iy = g+1, Bs(2)+g
                    y = (x0(2) + dble(iy-g-1)*dx(2) - params_acm%y0source(iscalar))**2
                    do ix = g+1, Bs(1)+g
                        x = (x0(1) + dble(ix-g-1)*dx(1) - params_acm%x0source(iscalar))**2

                        masksource = dexp( -(x + y) / (params_acm%widthsource(iscalar))**2  )

                        if (masksource > 1.0d-6) then
                            ! for the source term, we use the usual dirichlet C_eta
                            ! to force scalar to 1
                            source(ix,iy,1) = (masksource - phi(ix,iy,1,j)) / params_acm%C_eta
                            ! source(ix,iy,1) = -masksource*(phi(ix,iy,1,j)-1.d0) / params_acm%C_eta
                        endif
                    end do
                end do

            case ("mask_color_emission")
                where ( abs(mask(:,:,:,5) - params_acm%widthsource(iscalar)) <=1.0e-8 )
                    source = -mask(:,:,:,5)*(phi(:,:,:,j)-1.d0) / params_acm%C_eta
                end where

            case ("none", "empty")
                ! do nothing.

            case default
                call abort(2608191,"scalar source is unkown.")

            end select


            ! sponge layer
            if (params_acm%absorbing_sponge) then
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! for the source term, we use the usual dirichlet C_eta
                        ! to force scalar to 0
                        source(ix,iy,1) = source(ix,iy,1) - mask(ix,iy,1,6)*phi(ix,iy,1,j) / params_acm%C_eta
                    end do
                end do
            endif


            ! 2nd: compute rhs for this scalar.
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ux = phi(ix,iy,1,1)
                    uy = phi(ix,iy,1,2)

                    ! ATTENTION you need to sync the mask
                    chi = mask(ix,iy,1,1)
                    usx = mask(ix,iy,1,2)
                    usy = mask(ix,iy,1,3)

                    ! penalized diffusion coefficient at this point
                    D = kappa*(1.0_rk - chi) + params_acm%scalar_Ceta(iscalar)*chi

                    ! this is the vector in front of the gradient
                    wx = -((1.d0-chi)*ux + chi*usx)
                    wy = -((1.d0-chi)*uy + chi*usy)

                    ! gradient of passive scalar
                    gx = (a(-3)*phi(ix-3,iy,1,j)&
                         +a(-2)*phi(ix-2,iy,1,j)&
                         +a(-1)*phi(ix-1,iy,1,j)&
                         +a( 0)*phi(ix  ,iy,1,j)&
                         +a(+2)*phi(ix+2,iy,1,j)&
                         +a(+3)*phi(ix+3,iy,1,j)&
                         +a(+1)*phi(ix+1,iy,1,j))*dx_inv

                    gy = (a(-3)*phi(ix,iy-3,1,j)&
                         +a(-2)*phi(ix,iy-2,1,j)&
                         +a(-1)*phi(ix,iy-1,1,j)&
                         +a( 0)*phi(ix,iy  ,1,j)&
                         +a(+2)*phi(ix,iy+2,1,j)&
                         +a(+3)*phi(ix,iy+3,1,j)&
                         +a(+1)*phi(ix,iy+1,1,j))*dy_inv

                    ! gradient of mask function ( we need that for the diffusive term)
                    ! since this guy reads div( (kappa(1-mask) + eps*mask) * grad(phi) )
                    ! so this boils down to d/dx (D*gx) = D_dx*gx + D*gxx
                    ! so we need D_dx and this is kappa*(1-mask_dx)+ eps*mask_dx
                    chidx = (a(-3)*mask(ix-3,iy,1, 1)&
                            +a(-2)*mask(ix-2,iy,1, 1)&
                            +a(-1)*mask(ix-1,iy,1, 1)&
                            +a( 0)*mask(ix  ,iy,1, 1)&
                            +a(+3)*mask(ix+3,iy,1, 1)&
                            +a(+2)*mask(ix+2,iy,1, 1)&
                            +a(+1)*mask(ix+1,iy,1, 1))*dx_inv

                    chidy = (a(-3)*mask(ix,iy-3,1, 1)&
                            +a(-2)*mask(ix,iy-2,1, 1)&
                            +a(-1)*mask(ix,iy-1,1, 1)&
                            +a( 0)*mask(ix,iy  ,1, 1)&
                            +a(+3)*mask(ix,iy+3,1, 1)&
                            +a(+2)*mask(ix,iy+2,1, 1)&
                            +a(+1)*mask(ix,iy+1,1, 1))*dy_inv

                    D_dx = kappa*(-chidx) + params_acm%scalar_Ceta(iscalar) * chidx
                    D_dy = kappa*(-chidy) + params_acm%scalar_Ceta(iscalar) * chidy

                    ! second derivatives of passive scalar
                    gxx = (b(-2)*phi(ix-2,iy,1 ,j)&
                          +b(-1)*phi(ix-1,iy,1 ,j)&
                          +b( 0)*phi(ix  ,iy,1 ,j)&
                          +b(+1)*phi(ix+1,iy,1 ,j)&
                          +b(+2)*phi(ix+2,iy,1 ,j))*dx2_inv
                    gyy = (b(-2)*phi(ix,iy-2,1 ,j)&
                          +b(-1)*phi(ix,iy-1,1 ,j)&
                          +b( 0)*phi(ix,iy  ,1 ,j)&
                          +b(+1)*phi(ix,iy+1,1 ,j)&
                          +b(+2)*phi(ix,iy+2,1 ,j))*dy2_inv

                    ! assemble everything
                    rhs(ix,iy,1,j) = wx*gx + wy*gy & ! penalized convection term
                                   + D_dx*gx + D*gxx + D_dy*gy + D*gyy & ! penalized laplacian
                                   + source(ix, iy, 1)
                end do
            end do
        end do ! loop over scalars


    else
        call abort(441167, "3d Discretization unkown "//order_discretization//", I ll walk into the light now." )
    end if
end subroutine RHS_2D_scalar
