!-----------------------------------------------------------------------------
! main level wrapper to set the right hand side on a block. Note this is completely
! independent of the grid and any MPI formalism, neighboring relations and the like.
! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
! from that. Ghost nodes are assumed to be sync'ed.
!-----------------------------------------------------------------------------
subroutine RHS_NSPP( time, u, g, x0, dx, rhs, mask, stage, n_domain, discretization_overwrite )
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

    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The n_domain tells you if the local field is adjacent to a domain boundary:
    ! n_domain(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    ! NOTE: NSPP only supports symmetry BC for the moment (which is handled by wabbit and not NSPP)
    integer(kind=2), intent(in) :: n_domain(3)

    ! overwrite the discretization for this call, important for working with different discretization for solving the pressure-poisson-equation
    character(len=cshort), intent(in), optional :: discretization_overwrite


    real(kind=rk), allocatable, save :: vor(:,:,:,:)
    integer(kind=ik) :: mpierr, i, dim, ix, iy, iz
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: tmp(1:3), tmp2, dV, dV2, penal(1:3), C_eta_apply(1:ncolors), x, y, z, f_block(1:3)
    integer(kind=2) :: color
    integer(kind=ik) :: i_insect
    character(len=clong) :: discretization


    if (.not. params_nspp%initialized) write(*,*) "WARNING: RHS_NSPP called but NSPP not initialized"

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    dim = params_nspp%dim
    discretization = params_nspp%discretization
    if (present(discretization_overwrite)) discretization = discretization_overwrite

    ! EXPERIMENTAL. Soft penalization: gently turn on penalization at the beginning (t=0) of a simulation. This makes the initial wave 
    ! travelling through the domain smoother. Its thickness is usually C_eta*C_0 (dimension: L). Slowly decreasing C_eta from C_eta_start=1.0
    ! to its final value (small, say C_eta=1e-4) thus drastically increases the shock width. The shock intensity decays as R^-1 (2D) or 
    ! R^-2 (3D), but even in 3D it can still require resolution of the shock wave that can be expensive.
    if (params_nspp%soft_penalization_startup) then
        if (time < params_nspp%penalization_startup_time) then
            params_nspp%C_eta_temp = params_nspp%C_eta_start
        elseif (time < params_nspp%penalization_startup_time + params_nspp%penalization_startup_tau) then
            ! the formulation with EXP seems to be more efficient in damping the initial shock wave; 2D tests showed that nicely.
            ! The new defaults are C_eta_start = 1.0  and  penalization_startup_tau=0.20. Note relatively long before this time the penalization
            ! is (almost) completely activated because of the exponential decay
            ! Min is used to cap the exponential to 1.0, so that nothing is happening before the startup time
            params_nspp%C_eta_temp = params_nspp%C_eta + (params_nspp%C_eta_start - params_nspp%C_eta)*min(exp(-20.0*(time-params_nspp%penalization_startup_time)/params_nspp%penalization_startup_tau), 1.0_rk)
        else
            ! constant value of C_eta after the startup phase
            params_nspp%C_eta_temp = params_nspp%C_eta
        endif
        C_eta_apply = params_nspp%C_eta
        C_eta_apply(params_nspp%penalization_startup_colors:) = params_nspp%C_eta_temp
    else
        ! constant value of C_eta
        C_eta_apply = params_nspp%C_eta
    endif

    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.
        ! performs initializations in the RHS module, such as resetting integrals

        ! Linear Forcing for HIT (Lundgren) requires us to know kinetic energy and dissipation
        ! rate at all times, so compute that, if we use the forcing.
        if (params_nspp%HIT_linear_forcing) then
            params_nspp%e_kin = 0.0_rk
            params_nspp%enstrophy = 0.0_rk
            params_nspp%mean_flow = 0.0_rk
        endif

        ! let's update all insects. If there is none, then this is just an empty loop, so no problemo
        call Update_All_Insects(time)

        if (params_nspp%use_free_flight_solver) then
            ! reset forces, we compute their value now
            do i_insect = 1, n_insects
                Insects(i_insect)%force_g = 0.0_rk
                Insects(i_insect)%moment_g = 0.0_rk
            enddo
        endif

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
            if (maxval(abs(u(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,i))) > LIM_DIVERGED) then
                write(*,'("RHS: maxval in u(",i2,") = ", es10.3, ", Block with origin", 3(1x,es9.2), " and dx", 3(1x,es9.2))') i, maxval(abs(u(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,i))), x0, dx

                ! call dump_block_fancy(u(:,:,:,i:i), "block_NSPP_diverged_RHS.txt", Bs, g)

                ! done by all ranks but well I hope the cluster can take one for the team.
                ! This (empty) file is for scripting purposes on the supercomputers.
                open (77, file='NSPP_diverged', status='replace')
                close(77)
                
                call abort(0409201933,"NSPP fail: very very large values in state vector.")
            endif
        enddo

        dV = product(dx(1:params_nspp%dim))

        ! Linear Forcing for HIT (Lundgren) requires us to know kinetic energy and dissipation
        ! rate at all times, so compute that, if we use the forcing.
        if (params_nspp%HIT_linear_forcing) then
            ! vorticity work array
            if (.not. allocated(vor) ) allocate(vor(1:size(u,1), 1:size(u,2), 1:size(u,3), 1:3 ))
            ! to compute the current dissipation rate
            call compute_vorticity(u(:,:,:,1:3), dx, Bs, g, discretization, vor(:,:,:,:))



            params_nspp%mean_flow(1) = params_nspp%mean_flow(1) + dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1))
            params_nspp%mean_flow(2) = params_nspp%mean_flow(2) + dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 2))
            if (params_nspp%dim==2) then
                params_nspp%e_kin = params_nspp%e_kin + 0.5_rk*dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g,    :       , 1:params_nspp%dim)**2)
                params_nspp%enstrophy = params_nspp%enstrophy + 0.5_rk*dv*sum(vor(g+1:Bs(1)+g, g+1:Bs(2)+g,    :       , 1)**2)
            else
                params_nspp%e_kin = params_nspp%e_kin + 0.5_rk*dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:params_nspp%dim)**2)
                params_nspp%enstrophy = params_nspp%enstrophy + 0.5_rk*dv*sum(vor(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:3)**2)
                params_nspp%mean_flow(3) = params_nspp%mean_flow(3) + dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 3))
            endif

            ! NOTE: MPI_SUM is perfomed in the post_stage.
        endif

        ! if (params_nspp%geometry == "Insect".and. params_nspp%use_free_flight_solver) then
        if (params_nspp%use_free_flight_solver) then
            dV = product(dx(1:params_nspp%dim))
            f_block = 0.0_rk
            do i_insect = 1, n_insects
                do iz = merge(1, g+1, params_nspp%dim == 2), merge(1, Bs(3)+g, params_nspp%dim == 2)
                    if (params_nspp%dim == 2) then
                        z = 0.0_rk
                    else
                        z = x0(3) + dble(iz-(g+1)) * dx(3) - Insects(i_insect)%xc_body_g(3)
                    endif
                    do iy = g+1, Bs(2)+g
                        y = x0(2) + dble(iy-(g+1)) * dx(2) - Insects(i_insect)%xc_body_g(2)
                        do ix = g+1, Bs(1)+g
                            x = x0(1) + dble(ix-(g+1)) * dx(1) - Insects(i_insect)%xc_body_g(1)

                            ! get this points color
                            color = int( mask(ix, iy, iz, 5), kind=2 )

                            ! only include parts of the insect
                            if (any(color == (/insects(i_insect)%color_body, insects(i_insect)%color_l, insects(i_insect)%color_r, insects(i_insect)%color_l2, insects(i_insect)%color_r2/))) then
                                ! penalization term
                                penal = -mask(ix,iy,iz,1) * (u(ix,iy,iz,1:3) - mask(ix,iy,iz,2:4)) / C_eta_apply(color)

                                f_block(1:params_nspp%dim) = f_block(1:params_nspp%dim) - penal(1:params_nspp%dim)

                                ! x_lev = periodize_coordinate(x_lev, (/xl,yl,zl/))

                                ! moments. For insects, we compute the total moment wrt to the body center
                                insects(i_insect)%moment_g = insects(i_insect)%moment_g - cross((/x, y, z/), penal)*dV
                            endif
                        enddo
                    enddo
                enddo
                insects(i_insect)%force_g = insects(i_insect)%force_g + f_block*dV
            enddo
            ! NOTE: MPI_SUM is perfomed in the post_stage.

        endif

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.

        ! Linear Forcing for HIT (Lundgren) requires us to know kinetic energy and dissipation
        ! rate at all times, so compute that, if we use the forcing.
        if (params_nspp%HIT_linear_forcing) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            params_nspp%mean_flow = params_nspp%mean_flow / product(params_nspp%domain_size(1:params_nspp%dim))
            params_nspp%dissipation = params_nspp%enstrophy * params_nspp%nu
        endif

        if (params_nspp%use_free_flight_solver) then
            do i_insect = 1, n_insects
                call MPI_ALLREDUCE(MPI_IN_PLACE, insects(i_insect)%force_g, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
                call MPI_ALLREDUCE(MPI_IN_PLACE, insects(i_insect)%moment_g, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            enddo
        endif

    case ("local_stage")
        !-------------------------------------------------------------------------
        ! 4th stage: local evaluation of RHS on all blocks
        !-------------------------------------------------------------------------
        ! Evaluate local differential operators for the Navier-Stokes equations
        ! using the generic finite difference operators from module_operators.
        ! This computes: dU/dt = -1/2 (div(U\otimes U) + U\cdot grad(U)) + \nu \Delta U - \chi/C_eta (U - U_s)
        ! where the skew-symmetric formulation ensures energy conservation.
        !
        ! Called for each block individually.

        if (params_nspp%compute_flow) then
            ! Compute velocity RHS using generic finite difference operators
            call RHS_NSPP_Velocity(g, Bs, dx, x0, u, discretization, &
                                   time, rhs, mask, n_domain)
        endif

        ! --------------------------------------------------------------------------
        ! passive scalars
        ! --------------------------------------------------------------------------
        if (params_nspp%use_passive_scalar) then
            if (params_nspp%dim == 2) then
                call RHS_2D_scalar(g, Bs, dx, x0, u, discretization, time, rhs, mask, n_domain)
            else
                call RHS_3D_scalar(g, Bs, dx, x0, u, discretization, time, rhs, mask, n_domain)
            endif
        endif

    case ("divergence_stage")
        !-------------------------------------------------------------------------
        ! 5th stage: divergence computation
        !-------------------------------------------------------------------------
        ! Compute the divergence of the velocity RHS for the pressure Poisson equation.
        ! The divergence field is stored in rhs(:,:,:,1) to serve as the RHS for
        ! the pressure Poisson equation: Laplacian(p) = div(dU/dt)
        ! Uses left-sided finite difference stencil (FD1_l), which is complementary to
        ! the right-sided stencil (FD1_r) used in the pressure gradient computation below.
        call compute_divergence(rhs(:,:,:,1:params_nspp%dim), dx, Bs, g, &
                                discretization, u(:,:,:,1))

    case ("pressure_gradient_stage")
        !-------------------------------------------------------------------------
        ! 6th stage: pressure gradient computation
        !-------------------------------------------------------------------------
        ! Apply the pressure gradient correction to the velocity RHS:
        !   RHS_corrected = RHS_uncorrected - grad(p)
        ! where p (pressure) is stored in u(:,:,:,params_nspp%dim+1).
        ! This corrects the velocity RHS to enforce the divergence-free constraint.
        ! Uses right-sided finite difference stencil (FD1_r), which is complementary to
        ! the left-sided stencil (FD1_l) used in the divergence computation above.
        call subtract_scalar_gradient(rhs(:,:,:,1:params_nspp%dim), dx, Bs, g, &
                                      discretization, u(:,:,:,params_nspp%dim+1))

    case default
        call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")

    end select


end subroutine RHS_NSPP



!-----------------------------------------------------------------------------
! Compute velocity RHS for Navier-Stokes equations with generic operators
!-----------------------------------------------------------------------------
! This subroutine computes the right-hand side of the velocity equations
! using generic finite difference operators. It supports both 2D and 3D flows
! and implements the skew-symmetric formulation for energy conservation.
!
! Mathematical formulation:
!   dU/dt = -1/2 (div(U\otimes U) + U\cdot grad(U)) + \nu \Delta U - \chi/C_eta (U - U_s)
!
! where:
!   - First term: skew-symmetric convection (energy-conserving)
!   - Second term: viscous diffusion
!   - Third term: volume penalization for immersed boundaries
!   - Fourth term: sponge layer (optional)
!   - Fifth term: HIT linear forcing (optional)
!-----------------------------------------------------------------------------
subroutine RHS_NSPP_Velocity(g, Bs, dx, x0, phi, order_discretization, time, rhs, mask, n_domain)
    use module_operators
    implicit none
    
    !> grid parameter
    integer(kind=ik), intent(in) :: g
    !> block dimensions (without ghost nodes)
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> origin and spacing of the block
    real(kind=rk), dimension(1:3), intent(in) :: dx, x0
    !> datafields (velocity field, components 1:dim contain velocity)
    real(kind=rk), intent(inout) :: phi(1:,1:,1:,1:)
    !> discretization order
    character(len=cshort), intent(in) :: order_discretization
    !> time
    real(kind=rk), intent(in) :: time
    !> output: velocity RHS
    real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)
    !> mask array: component 1 = chi, components 2:dim+1 = target velocity, component 6 = sponge mask
    real(kind=rk), intent(in) :: mask(1:,1:,1:,1:)
    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The n_domain tells you if the local field is adjacent to a domain boundary:
    ! n_domain(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    integer(kind=2), intent(in) :: n_domain(3)
    
    ! Local variables for derivatives and velocities
    integer(kind=ik) :: ix, iy, iz, dim, color
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, C_eta_apply(1:ncolors)
    real(kind=rk) :: u_dx, u_dy, u_dz, u_dxdx, u_dydy, u_dzdz, &
                     v_dx, v_dy, v_dz, v_dxdx, v_dydy, v_dzdz, &
                     w_dx, w_dy, w_dz, w_dxdx, w_dydy, w_dzdz, &
                     penalx, penaly, penalz, ux, uy, uz, &
                     uu_dx, uv_dy, uw_dz, vu_dx, vv_dy, vw_dz, wu_dx, wv_dy, ww_dz
    
    ! Variables for sponge term
    real(kind=rk) :: spo
    
    ! Variables for HIT linear forcing
    real(kind=rk) :: A_forcing, G_gain, e_kin_set, t_l_inf

    !> finite difference stencils
    real(kind=rk), allocatable :: FD1_l(:), FD1_r(:), FD2(:)
    integer :: FD1_ls, FD1_le, FD1_rs, FD1_re, FD2_s, FD2_e
    
    ! Extract commonly used parameters from params_nspp for convenience
    dim = params_nspp%dim
    
    ! Size check for safety
    if (dim == 2) then
        if (size(phi,1)/=Bs(1)+2*g .or. size(phi,2)/=Bs(2)+2*g .or. size(phi,3)/=1 .or. size(phi,4)<params_nspp%dim+params_nspp%N_scalars+params_nspp%N_time_statistics) then
            write(*, '("RHS_NSPP_Velocity: size(phi) = (", 4I10, "), expected = (", 4I10, ")")') size(phi,1), size(phi,2), size(phi,3), size(phi,4), Bs(1)+2*g, Bs(2)+2*g, 1, params_nspp%dim+params_nspp%N_scalars+params_nspp%N_time_statistics
            call abort(66233,"RHS_NSPP_Velocity: wrong array size (2D), I go for a walk instead.")
        endif
    else
        if (size(phi,1)/=Bs(1)+2*g .or. size(phi,2)/=Bs(2)+2*g .or. size(phi,3)/=Bs(3)+2*g .or. size(phi,4)<params_nspp%dim+params_nspp%N_scalars+params_nspp%N_time_statistics) then
            write(*, '("RHS_NSPP_Velocity: size(phi) = (", 4I10, "), expected = (", 4I10, ")")') size(phi,1), size(phi,2), size(phi,3), size(phi,4), Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, params_nspp%dim+params_nspp%N_scalars+params_nspp%N_time_statistics
            call abort(66234,"RHS_NSPP_Velocity: wrong array size (3D), I go for a walk instead.")
        endif
    endif
    
    ! =========================================================================
    ! Setup finite difference stencils
    ! =========================================================================
    call setup_FD1_left_stencil(order_discretization, FD1_l, FD1_ls, FD1_le)
    call setup_FD1_right_stencil(order_discretization, FD1_r, FD1_rs, FD1_re)
    
    ! Performance optimization: skip FD2 setup if viscosity is zero
    if (abs(params_nspp%nu) < 1.0e-14_rk) then
        allocate(FD2(0:0))
        FD2(0) = 0.0_rk
        FD2_s = 0
        FD2_e = 0
    else
        call setup_FD2_stencil(order_discretization, FD2, FD2_s, FD2_e)
    endif
    
    ! Precompute inverse grid spacings and penalization parameter
    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    ! for now - c_eta can only vary if the geometry is faded in or not. Later, this can be changed for full flexibility for each color
    C_eta_apply = params_nspp%C_eta
    C_eta_apply(params_nspp%penalization_startup_colors:) = params_nspp%C_eta_temp
    
    if (dim == 3) then
        dz_inv = 1.0_rk / dx(3)
        dz2_inv = 1.0_rk / (dx(3)**2)
    endif
    
    ! =========================================================================
    ! Compute RHS: dU/dt = -1/2 (div(U\otimes U) + U\cdot grad(U)) + \nu \Delta U - \chi/C_eta(U-U_s)
    ! =========================================================================
    
    if (.not. params_nspp%skew_symmetry) then
        ! Convective (non-skew-symmetric) formulation is not yet implemented
        call abort(2501042, "RHS_NSPP_Velocity: Non-skew-symmetric formulation not implemented. Set skew_symmetry=.true.")
    endif
    
    if (dim == 2) then
        ! =====================================================================
        ! 2D Case - Skew-symmetric formulation
        ! =====================================================================
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! Velocity components
                ux = phi(ix, iy, 1, 1)
                uy = phi(ix, iy, 1, 2)

                ! First derivatives for advection form: U\cdot grad(U)
                u_dx = sum(FD1_r(FD1_rs:FD1_re) * phi(ix+FD1_rs:ix+FD1_re,iy,1,1)) * dx_inv
                u_dy = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy+FD1_rs:iy+FD1_re,1,1)) * dy_inv
                
                v_dx = sum(FD1_r(FD1_rs:FD1_re) * phi(ix+FD1_rs:ix+FD1_re,iy,1,2)) * dx_inv
                v_dy = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy+FD1_rs:iy+FD1_re,1,2)) * dy_inv

                ! Divergence form: div(U\otimes U)
                uu_dx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,1,1) * phi(ix+FD1_ls:ix+FD1_le,iy,1,1)) * dx_inv
                uv_dy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,1,1) * phi(ix,iy+FD1_ls:iy+FD1_le,1,2)) * dy_inv

                vu_dx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,1,2) * phi(ix+FD1_ls:ix+FD1_le,iy,1,1)) * dx_inv
                vv_dy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,1,2) * phi(ix,iy+FD1_ls:iy+FD1_le,1,2)) * dy_inv

                ! Second derivatives for diffusion
                u_dxdx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,1,1)) * dx2_inv
                u_dydy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,1,1)) * dy2_inv

                v_dxdx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,1,2)) * dx2_inv
                v_dydy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,1,2)) * dy2_inv

                ! Penalization: -\chi/C_eta * (U - U_target)
                color = int(mask(ix,iy,1,5), kind=2)
                penalx = -mask(ix,iy,1,1) / C_eta_apply( color ) * (ux - mask(ix,iy,1,2))
                penaly = -mask(ix,iy,1,1) / C_eta_apply( color ) * (uy - mask(ix,iy,1,3))

                ! Assemble RHS: -1/2 (div(U\otimes U) + U\cdot grad(U)) + \nu \Delta U + penalization
                rhs(ix,iy,1,1) = -0.5_rk*(uu_dx + uv_dy + ux*u_dx + uy*u_dy) + params_nspp%nu*(u_dxdx + u_dydy) + penalx
                rhs(ix,iy,1,2) = -0.5_rk*(vu_dx + vv_dy + ux*v_dx + uy*v_dy) + params_nspp%nu*(v_dxdx + v_dydy) + penaly
            enddo
        enddo
        
    else
        ! =====================================================================
        ! 3D Case - Skew-symmetric formulation
        ! =====================================================================
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! Velocity components
                    ux = phi(ix, iy, iz, 1)
                    uy = phi(ix, iy, iz, 2)
                    uz = phi(ix, iy, iz, 3)

                    ! First derivatives for advection form: U\cdot grad(U)
                    u_dx = sum(FD1_r(FD1_rs:FD1_re) * phi(ix+FD1_rs:ix+FD1_re,iy,iz,1)) * dx_inv
                    u_dy = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy+FD1_rs:iy+FD1_re,iz,1)) * dy_inv
                    u_dz = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy,iz+FD1_rs:iz+FD1_re,1)) * dz_inv

                    v_dx = sum(FD1_r(FD1_rs:FD1_re) * phi(ix+FD1_rs:ix+FD1_re,iy,iz,2)) * dx_inv
                    v_dy = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy+FD1_rs:iy+FD1_re,iz,2)) * dy_inv
                    v_dz = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy,iz+FD1_rs:iz+FD1_re,2)) * dz_inv

                    w_dx = sum(FD1_r(FD1_rs:FD1_re) * phi(ix+FD1_rs:ix+FD1_re,iy,iz,3)) * dx_inv
                    w_dy = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy+FD1_rs:iy+FD1_re,iz,3)) * dy_inv
                    w_dz = sum(FD1_r(FD1_rs:FD1_re) * phi(ix,iy,iz+FD1_rs:iz+FD1_re,3)) * dz_inv

                    ! Divergence form: div(U\otimes U)
                    uu_dx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,1) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,1)) * dx_inv
                    uv_dy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,1) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,2)) * dy_inv
                    uw_dz = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,1) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,3)) * dz_inv

                    vu_dx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,2) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,1)) * dx_inv
                    vv_dy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,2) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,2)) * dy_inv
                    vw_dz = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,2) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,3)) * dz_inv

                    wu_dx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,3) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,1)) * dx_inv
                    wv_dy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,3) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,2)) * dy_inv
                    ww_dz = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,3) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,3)) * dz_inv

                    ! Second derivatives for diffusion
                    u_dxdx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,iz,1)) * dx2_inv
                    u_dydy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,iz,1)) * dy2_inv
                    u_dzdz = sum(FD2(FD2_s:FD2_e) * phi(ix,iy,iz+FD2_s:iz+FD2_e,1)) * dz2_inv

                    v_dxdx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,iz,2)) * dx2_inv
                    v_dydy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,iz,2)) * dy2_inv
                    v_dzdz = sum(FD2(FD2_s:FD2_e) * phi(ix,iy,iz+FD2_s:iz+FD2_e,2)) * dz2_inv

                    w_dxdx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,iz,3)) * dx2_inv
                    w_dydy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,iz,3)) * dy2_inv
                    w_dzdz = sum(FD2(FD2_s:FD2_e) * phi(ix,iy,iz+FD2_s:iz+FD2_e,3)) * dz2_inv

                    ! Penalization: -\chi/C_eta * (U - U_target)
                    color = int(mask(ix,iy,iz,5), kind=2)
                    penalx = -mask(ix,iy,iz,1) / C_eta_apply(color) * (ux - mask(ix,iy,iz,2))
                    penaly = -mask(ix,iy,iz,1) / C_eta_apply(color) * (uy - mask(ix,iy,iz,3))
                    penalz = -mask(ix,iy,iz,1) / C_eta_apply(color) * (uz - mask(ix,iy,iz,4))

                    ! Assemble RHS: -1/2 (div(U\otimes U) + U\cdot grad(U)) + \nu \Delta U + penalization
                    rhs(ix,iy,iz,1) = -0.5_rk * (uu_dx + uv_dy + uw_dz + ux*u_dx + uy*u_dy + uz*u_dz) &
                        + params_nspp%nu*(u_dxdx + u_dydy + u_dzdz) + penalx

                    rhs(ix,iy,iz,2) = -0.5_rk * (vu_dx + vv_dy + vw_dz + ux*v_dx + uy*v_dy + uz*v_dz) &
                        + params_nspp%nu*(v_dxdx + v_dydy + v_dzdz) + penaly

                    rhs(ix,iy,iz,3) = -0.5_rk * (wu_dx + wv_dy + ww_dz + ux*w_dx + uy*w_dy + uz*w_dz) &
                        + params_nspp%nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                enddo
            enddo
        enddo
    endif
    
    ! =========================================================================
    ! Sponge layer term (optional)
    ! =========================================================================
    ! The sponge term damps velocity fluctuations near domain boundaries to
    ! prevent reflections. It acts on ALL velocity components.
    ! Term: -\chi_sponge/C_sponge * (U - U_target)
    ! 
    ! Special case: Lamballais cylinder geometry uses C_eta_ring instead of C_sponge
    ! Reference: Gautier, R., Biau, D., Lamballais, E.: A reference solution of 
    ! the flow over a circular cylinder at Re = 40, Computers & Fluids 75, 103–111, 2013
    ! =========================================================================
    if (params_nspp%use_sponge) then
        ! Check for Lamballais special treatment (2D only)
        if (dim == 2 .and. (any(params_nspp%geometries(:) == "lamballais") .or. any(params_nspp%geometries(:) == "lamballais-local"))) then
            ! Special treatment for Lamballais: use C_eta_ring for the outer ring penalization
            ! Note: C_eta_ring can differ from C_eta in local variation cases
            ! The ring (stored in sponge mask array) penalizes velocity towards target values            
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! Sponge mask (component 6) contains ring mask
                    ! Components 2,3 contain target velocities in the ring
                    rhs(ix,iy,1,1) = rhs(ix,iy,1,1) - mask(ix,iy,1,6)*(phi(ix,iy,1,1) - mask(ix,iy,1,2))/ params_nspp%C_eta_ring
                    rhs(ix,iy,1,2) = rhs(ix,iy,1,2) - mask(ix,iy,1,6)*(phi(ix,iy,1,2) - mask(ix,iy,1,3))/ params_nspp%C_eta_ring
                    ! Note: Pressure component (if present) would be handled separately in pressure equation
                enddo
            enddo
        else
            ! Standard sponge treatment            
            if (dim == 2) then
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! Sponge mask is stored in component 6 of mask array
                        spo = mask(ix,iy,1,6) / params_nspp%C_sponge
                        
                        rhs(ix,iy,1,1) = rhs(ix,iy,1,1) - (phi(ix,iy,1,1) - params_nspp%u_mean_set(1)) * spo
                        rhs(ix,iy,1,2) = rhs(ix,iy,1,2) - (phi(ix,iy,1,2) - params_nspp%u_mean_set(2)) * spo
                    enddo
                enddo
            else
                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            ! Sponge mask is stored in component 6 of mask array
                            spo = mask(ix,iy,iz,6) / params_nspp%C_sponge
                            
                            rhs(ix,iy,iz,1) = rhs(ix,iy,iz,1) - (phi(ix,iy,iz,1) - params_nspp%u_mean_set(1)) * spo
                            rhs(ix,iy,iz,2) = rhs(ix,iy,iz,2) - (phi(ix,iy,iz,2) - params_nspp%u_mean_set(2)) * spo
                            rhs(ix,iy,iz,3) = rhs(ix,iy,iz,3) - (phi(ix,iy,iz,3) - params_nspp%u_mean_set(3)) * spo
                        enddo
                    enddo
                enddo
            endif
        endif
    endif
    
    ! =========================================================================
    ! HIT (Homogeneous Isotropic Turbulence) linear forcing (optional)
    ! =========================================================================
    ! Linear forcing to maintain constant kinetic energy in HIT simulations.
    ! Based on Lundgren forcing with energy feedback control.
    ! Reference: Bassenne et al. (2016)
    ! 
    ! The forcing is applied to velocity fluctuations (removing mean flow)
    ! to maintain target kinetic energy while accounting for dissipation.
    ! =========================================================================
    if (params_nspp%HIT_linear_forcing) then
        G_gain = params_nspp%HIT_gain
        e_kin_set = params_nspp%HIT_energy * product(params_nspp%domain_size(1:params_nspp%dim))
        t_l_inf = 1.0_rk  ! sqrt(nu/epsilon), adjusted via gain
        
        ! Compute forcing amplitude: A = (epsilon - G*(E-E_target)/t_l) / (2*E)
        A_forcing = (params_nspp%dissipation - G_gain * (params_nspp%e_kin - e_kin_set) / t_l_inf) / (2.0_rk * params_nspp%e_kin)
        
        if (dim == 2) then
            ! Apply forcing to velocity fluctuations (subtract mean flow)
            rhs(:,:,1,1) = rhs(:,:,1,1) + A_forcing * (phi(:,:,1,1) - params_nspp%mean_flow(1))
            rhs(:,:,1,2) = rhs(:,:,1,2) + A_forcing * (phi(:,:,1,2) - params_nspp%mean_flow(2))
        else
            ! Apply forcing to velocity fluctuations (subtract mean flow)
            rhs(:,:,:,1) = rhs(:,:,:,1) + A_forcing * (phi(:,:,:,1) - params_nspp%mean_flow(1))
            rhs(:,:,:,2) = rhs(:,:,:,2) + A_forcing * (phi(:,:,:,2) - params_nspp%mean_flow(2))
            rhs(:,:,:,3) = rhs(:,:,:,3) + A_forcing * (phi(:,:,:,3) - params_nspp%mean_flow(3))
        endif
    endif
    
    ! Clean up allocated stencil arrays
    deallocate(FD1_l, FD1_r, FD2)

end subroutine RHS_NSPP_Velocity



subroutine RHS_3D_scalar(g, Bs, dx, x0, phi, order_discretization, time, rhs, mask, n_domain)
    use module_operators
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
    character(len=cshort), intent(in)       :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time
    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The n_domain tells you if the local field is adjacent to a domain boundary:
    ! n_domain(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    integer(kind=2), intent(in) :: n_domain(3)

    integer(kind=ik) :: ix, iy, iz, iscalar, j

    !> parameters for FD1L, FD2 operators (generalized stencil approach)
    real(kind=rk), allocatable, dimension(:) :: FD1_l, FD2
    integer(kind=ik) :: FD1_ls, FD1_le, FD2_s, FD2_e

    real(kind=rk) :: kappa, x, y, z, masksource, nu, R, R0sq, C_eta_apply(1:ncolors)
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv
    real(kind=rk) :: ux, uy, uz,&
    usx,usy,usz,wx,wy,wz,gx,gy,gz,D,chi,chidx,chidz,chidy,D_dx,D_dy,D_dz,gxx,gyy,gzz
    real(kind=rk) :: phi_dx, phi_dy, phi_dz, phi_dxdx, phi_dydy, phi_dzdz
    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're negligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: source(:,:,:)

    if (.not. allocated(source)) allocate(source(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g))

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dz_inv = 1.0_rk / dx(3)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    dz2_inv = 1.0_rk / (dx(3)**2)

    nu = params_nspp%nu

    ! for now - c_eta can only vary if the geometry is faded in or not. Later, this can be changed for full flexibility for each color
    C_eta_apply = params_nspp%C_eta
    C_eta_apply(params_nspp%penalization_startup_colors:) = params_nspp%C_eta_temp

    !-----------------------------------------------------------------------
    ! passive scalar equations: loop over all scalars and compute RHS
    !-----------------------------------------------------------------------
    ! Setup stencils using the unified interface from module_operators
    call setup_FD1_left_stencil(order_discretization, FD1_l, FD1_ls, FD1_le)
    call setup_FD2_stencil(order_discretization, FD2, FD2_s, FD2_e)

    ! Loop over all scalars
    do iscalar = 1, params_nspp%N_scalars
        ! actual index of this scalar in the array
        j = iscalar + (params_nspp%dim + 1)

        ! compute diffusivity from schmidt number (and fluid viscosity)
        kappa = nu / params_nspp%schmidt_numbers(iscalar)

        ! reset source term for each scalar.
        source = 0.0_rk

        ! 1st: compute source terms (note the strcmp needs to be outside the loop)
        select case (params_nspp%scalar_source_type(iscalar))
        case ("gaussian")
            do iz = g+1, Bs(3)+g
                z = (x0(3) + dble(iz-g-1)*dx(3) - params_nspp%z0source(iscalar))**2
                do iy = g+1, Bs(2)+g
                    y = (x0(2) + dble(iy-g-1)*dx(2) - params_nspp%y0source(iscalar))**2
                    do ix = g+1, Bs(1)+g
                        x = (x0(1) + dble(ix-g-1)*dx(1) - params_nspp%x0source(iscalar))**2

                        R = x + y + z ! note this is (x-x0)**2

                        masksource = dexp( -R / (params_nspp%widthsource(iscalar)**2)  )

                        if (masksource > 1.0d-6) then
                            ! for the source term, we use the usual dirichlet C_eta
                            ! to force scalar to 1
                            ! source(ix,iy,iz) = -1.0d0*(phi(ix,iy,iz,j)-masksource-) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) )
                            source(ix,iy,iz) = (masksource - phi(ix,iy,iz,j)) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) )
                            ! source(ix,iy,iz) = -masksource*(phi(ix,iy,iz,j)-1.d0) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) )
                        endif
                    end do
                end do
            end do
        case ("circular")
            R0sq = params_nspp%widthsource(iscalar)**2
            do iz = g+1, Bs(3)+g
                z = (x0(3) + dble(iz-g-1)*dx(3) - params_nspp%z0source(iscalar))**2
                do iy = g+1, Bs(2)+g
                    y = (x0(2) + dble(iy-g-1)*dx(2) - params_nspp%y0source(iscalar))**2
                    do ix = g+1, Bs(1)+g
                        x = (x0(1) + dble(ix-g-1)*dx(1) - params_nspp%x0source(iscalar))**2

                        R = x + y + z ! note this is (x-x0)**2

                        if ( R <= R0sq ) then
                            ! for the source term, we use the usual dirichlet C_eta
                            ! to force scalar to 1
                            source(ix,iy,iz) = -(phi(ix,iy,iz,j)-1.d0) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) )
                        endif
                    end do
                end do
            end do

        case ("inflow-x")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        x = x0(1) + dble(ix-(g+1))*dx(1)
                        if ( x <= params_nspp%widthsource(iscalar) ) then
                            ! INFLOW
                            source(ix,iy,iz) = (1.0_rk - phi(ix,iy,iz,j)) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) ) ! for the source term, we use the usual dirichlet C_eta
                        endif
                    end do
                end do
            end do

        case ("in+outflow-x")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        x = x0(1) + dble(ix-(g+1))*dx(1)
                        if ( x <= params_nspp%widthsource(iscalar) ) then
                            ! INFLOW
                            source(ix,iy,iz) = (1.0_rk - phi(ix,iy,iz,j)) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) ) ! for the source term, we use the usual dirichlet C_eta
                        endif
                        if ( x >= params_nspp%domain_size(1)-params_nspp%widthsource(iscalar) ) then
                            ! OUTFLOW
                            source(ix,iy,iz) = (0.0_rk - phi(ix,iy,iz,j)) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) ) ! for the source term, we use the usual dirichlet C_eta
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
        if (params_nspp%absorbing_sponge) then
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! for the source term, we use the usual dirichlet C_eta
                        ! to force scalar to 1
                        source(ix,iy,iz) = source(ix,iy,iz) - mask(ix,iy,iz,6)*phi(ix,iy,iz,j) / params_nspp%C_sponge
                    end do
                end do
            enddo
        endif

        if (params_nspp%scalar_BC_type == "neumann") then
            ! 2nd: compute rhs for this scalar using generalized approach
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
                        D = kappa*(1.0_rk - chi) + params_nspp%scalar_Ceta(iscalar)*chi

                        ! this is the vector in front of the gradient
                        wx = -((1.d0-chi)*ux + chi*usx)
                        wy = -((1.d0-chi)*uy + chi*usy)
                        wz = -((1.d0-chi)*uz + chi*usz)

                        ! gradient of passive scalar using generalized stencils
                        gx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,j)) * dx_inv
                        gy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,j)) * dy_inv
                        gz = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,j)) * dz_inv

                        ! gradient of mask function (we need that for the diffusive term)
                        ! since this guy reads div( (kappa(1-mask) + eps*mask) * grad(phi) )
                        ! so this boils down to d/dx (D*gx) = D_dx*gx + D*gxx
                        ! so we need D_dx and this is kappa*(1-mask_dx)+ eps*mask_dx
                        chidx = sum(FD1_l(FD1_ls:FD1_le) * mask(ix+FD1_ls:ix+FD1_le,iy,iz,1)) * dx_inv
                        chidy = sum(FD1_l(FD1_ls:FD1_le) * mask(ix,iy+FD1_ls:iy+FD1_le,iz,1)) * dy_inv
                        chidz = sum(FD1_l(FD1_ls:FD1_le) * mask(ix,iy,iz+FD1_ls:iz+FD1_le,1)) * dz_inv

                        D_dx = kappa*(-chidx) + params_nspp%scalar_Ceta(iscalar) * chidx
                        D_dy = kappa*(-chidy) + params_nspp%scalar_Ceta(iscalar) * chidy
                        D_dz = kappa*(-chidz) + params_nspp%scalar_Ceta(iscalar) * chidz

                        ! second derivatives of passive scalar using generalized stencils
                        gxx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,iz,j)) * dx2_inv
                        gyy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,iz,j)) * dy2_inv
                        gzz = sum(FD2(FD2_s:FD2_e) * phi(ix,iy,iz+FD2_s:iz+FD2_e,j)) * dz2_inv

                        ! assemble everything
                        rhs(ix,iy,iz,j) = wx*gx + wy*gy + wz*gz & ! penalized convection term
                        + D_dx*gx + D*gxx & ! penalized laplacian
                        + D_dy*gy + D*gyy &
                        + D_dz*gz + D*gzz &
                        + source(ix,iy,iz)
                    end do
                end do
            end do
        elseif (params_nspp%scalar_BC_type == "dirichlet") then
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ux = phi(ix,iy,iz,1)
                        uy = phi(ix,iy,iz,2)
                        uz = phi(ix,iy,iz,3)

                        chi = mask(ix,iy,iz,1)

                        ! gradient using generalized stencils
                        phi_dx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,iz,j)) * dx_inv
                        phi_dy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,iz,j)) * dy_inv
                        phi_dz = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy,iz+FD1_ls:iz+FD1_le,j)) * dz_inv

                        ! laplace using generalized stencils
                        phi_dxdx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,iz,j)) * dx2_inv
                        phi_dydy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,iz,j)) * dy2_inv
                        phi_dzdz = sum(FD2(FD2_s:FD2_e) * phi(ix,iy,iz+FD2_s:iz+FD2_e,j)) * dz2_inv

                        ! easy RHS for dirichlet BC
                        rhs(ix,iy,iz,j) = -ux*phi_dx -uy*phi_dy -uz*phi_dz &
                        + kappa*(phi_dxdx + phi_dydy + phi_dzdz) &
                        - chi*(phi(ix,iy,iz,j) - 0.0_rk) / C_eta_apply( int(mask(ix,iy,iz,5), kind=2) ) & ! Dirichlet penalization for obstacle mask (instead of Neumann penalization)
                        + source(ix, iy, iz) ! source term is actually a dirichlet penalization term as well
                    enddo
                enddo
            enddo
        else
            call abort(22092234, "scalar_BC_type unkown"//trim(adjustl(params_nspp%scalar_BC_type))//". Time to fake a smile :)" )
        endif
    end do ! loop over scalars

end subroutine RHS_3D_scalar


subroutine RHS_2D_scalar(g, Bs, dx, x0, phi, order_discretization, time, rhs, mask, n_domain)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_globals
    use module_operators
    
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
    character(len=clong), intent(in)       :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time
    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The n_domain tells you if the local field is adjacent to a domain boundary:
    ! n_domain(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    integer(kind=2), intent(in) :: n_domain(3)

    integer(kind=ik) :: ix, iy, iscalar, j

    real(kind=rk) :: kappa, x, y, masksource, nu, R, C_eta_apply(1:ncolors)
    real(kind=rk) :: dx_inv, dy_inv, dx2_inv, dy2_inv
    real(kind=rk) :: ux, uy, usx, usy, wx, wy, gx, gy, D, chi, chidx, chidy, D_dx, D_dy, gxx, gyy
    real(kind=rk) :: phi_dx, phi_dy, phi_dxdx, phi_dydy

    !> parameters for FD1L, FD2 operators (generalized stencil approach)
    real(kind=rk), allocatable, dimension(:) :: FD1_l, FD2
    integer(kind=ik) :: FD1_ls, FD1_le, FD2_s, FD2_e

    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're negligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: source(:,:,:)

    if (.not. allocated(source)) allocate(source(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
    source = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    nu = params_nspp%nu

    ! for now - c_eta can only vary if the geometry is faded in or not. Later, this can be changed for full flexibility for each color
    C_eta_apply = params_nspp%C_eta
    C_eta_apply(params_nspp%penalization_startup_colors:) = params_nspp%C_eta_temp

    !-----------------------------------------------------------------------
    ! passive scalar equations: loop over all scalars and compute RHS
    !-----------------------------------------------------------------------
    ! Setup stencils using the unified interface from module_operators
    call setup_FD1_left_stencil(order_discretization, FD1_l, FD1_ls, FD1_le)
    call setup_FD2_stencil(order_discretization, FD2, FD2_s, FD2_e)

    ! Loop over all scalars
    do iscalar = 1, params_nspp%N_scalars
        ! actual index of this scalar in the array
        j = iscalar + (params_nspp%dim + 1)

        ! compute diffusivity from schmidt number (and of course fluid viscosity)
        kappa = nu / params_nspp%schmidt_numbers(iscalar)

        source = 0.0_rk

        ! 1st: compute source terms (note the strcmp needs to be outside the loop)
        select case (params_nspp%scalar_source_type(iscalar))
        case ("gaussian")
            do iy = g+1, Bs(2)+g
                y = (x0(2) + dble(iy-g-1)*dx(2) - params_nspp%y0source(iscalar))**2
                do ix = g+1, Bs(1)+g
                    x = (x0(1) + dble(ix-g-1)*dx(1) - params_nspp%x0source(iscalar))**2

                    masksource = dexp( -(x + y) / (params_nspp%widthsource(iscalar))**2  )

                    if (masksource > 1.0d-6) then
                        ! for the source term, we use the usual dirichlet C_eta
                        ! to force scalar to 1
                        source(ix,iy,1) = (masksource - phi(ix,iy,1,j)) / C_eta_apply( int(mask(ix,iy,1,5), kind=2) )
                        ! source(ix,iy,1) = -masksource*(phi(ix,iy,1,j)-1.d0) / C_eta_apply( int(mask(ix,iy,1,5), kind=2) )
                    endif
                end do
            end do

        case("in+outflow")
            do iy = g+1, Bs(2)+g
                y = x0(2) + dble(iy-(g+1))*dx(2)
                do ix = g+1, Bs(1)+g
                    x = x0(1) + dble(ix-(g+1))*dx(1)
                    if ( x <= params_nspp%widthsource(iscalar) ) then
                        ! INFLOW
                        source(ix,iy,1) = (1.0_rk - phi(ix,iy,1,j)) / C_eta_apply( int(mask(ix,iy,1,5), kind=2) ) ! for the source term, we use the usual dirichlet C_eta

                    elseif ( x >= params_nspp%domain_size(1)-params_nspp%widthsource(iscalar) ) then
                        ! OUTFLOW
                        source(ix,iy,1) = (0.0_rk - phi(ix,iy,1,j)) / C_eta_apply( int(mask(ix,iy,1,5), kind=2) )
                    endif
                end do
            end do

        case ("mask_color_emission")
            do iy = g+1, Bs(2)+g
                 do ix = g+1, Bs(1)+g
                     if ( abs(mask(ix,iy,1,5) - params_nspp%widthsource(iscalar)) <=1.0e-8 ) then
                        source(ix,iy,1) = -mask(ix,iy,1,5)*(phi(ix,iy,1,j)-1.d0) / C_eta_apply( int(mask(ix,iy,1,5), kind=2) )
                     endif
                 end do
             end do

        case ("none", "empty")
            ! do nothing.

        case default
            call abort(2608191,"scalar source is unkown.")

        end select


        ! sponge layer
        if (params_nspp%absorbing_sponge) then
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! for the source term, we use the usual dirichlet C_eta
                    ! to force scalar to 0
                    source(ix,iy,1) = source(ix,iy,1) - mask(ix,iy,1,6)*phi(ix,iy,1,j) / C_eta_apply( int(mask(ix,iy,1,5), kind=2) )
                end do
            end do
        endif


        if (params_nspp%scalar_BC_type == "neumann") then
            ! 2nd: compute rhs for this scalar using generalized approach
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ux = phi(ix,iy,1,1)
                    uy = phi(ix,iy,1,2)

                    ! ATTENTION you need to sync the mask
                    chi = mask(ix,iy,1,1)
                    usx = mask(ix,iy,1,2)
                    usy = mask(ix,iy,1,3)

                    ! penalized diffusion coefficient at this point
                    D = kappa*(1.0_rk - chi) + params_nspp%scalar_Ceta(iscalar)*chi

                    ! this is the vector in front of the gradient
                    wx = -((1.d0-chi)*ux + chi*usx)
                    wy = -((1.d0-chi)*uy + chi*usy)

                    ! gradient of passive scalar using generalized stencils
                    gx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,1,j)) * dx_inv
                    gy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,1,j)) * dy_inv

                    ! gradient of mask function (we need that for the diffusive term)
                    ! since this guy reads div( (kappa(1-mask) + eps*mask) * grad(phi) )
                    ! so this boils down to d/dx (D*gx) = D_dx*gx + D*gxx
                    ! so we need D_dx and this is kappa*(1-mask_dx)+ eps*mask_dx
                    chidx = sum(FD1_l(FD1_ls:FD1_le) * mask(ix+FD1_ls:ix+FD1_le,iy,1,1)) * dx_inv
                    chidy = sum(FD1_l(FD1_ls:FD1_le) * mask(ix,iy+FD1_ls:iy+FD1_le,1,1)) * dy_inv

                    D_dx = kappa*(-chidx) + params_nspp%scalar_Ceta(iscalar) * chidx
                    D_dy = kappa*(-chidy) + params_nspp%scalar_Ceta(iscalar) * chidy

                    ! second derivatives of passive scalar using generalized stencils
                    gxx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,1,j)) * dx2_inv
                    gyy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,1,j)) * dy2_inv

                    ! assemble everything
                    rhs(ix,iy,1,j) = wx*gx + wy*gy & ! penalized convection term
                                   + D_dx*gx + D*gxx + D_dy*gy + D*gyy & ! penalized laplacian
                                   + source(ix, iy, 1)
                end do
            end do
        elseif (params_nspp%scalar_BC_type == "dirichlet") then

            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ux = phi(ix,iy,1,1)
                    uy = phi(ix,iy,1,2)

                    chi = mask(ix,iy,1,1)
                    usx = mask(ix,iy,1,2)
                    usy = mask(ix,iy,1,3)

                    ! gradient using generalized stencils
                    phi_dx = sum(FD1_l(FD1_ls:FD1_le) * phi(ix+FD1_ls:ix+FD1_le,iy,1,j)) * dx_inv
                    phi_dy = sum(FD1_l(FD1_ls:FD1_le) * phi(ix,iy+FD1_ls:iy+FD1_le,1,j)) * dy_inv

                    ! laplace using generalized stencils
                    phi_dxdx = sum(FD2(FD2_s:FD2_e) * phi(ix+FD2_s:ix+FD2_e,iy,1,j)) * dx2_inv
                    phi_dydy = sum(FD2(FD2_s:FD2_e) * phi(ix,iy+FD2_s:iy+FD2_e,1,j)) * dy2_inv

                    ! easy RHS for dirichlet BC
                    rhs(ix,iy,1,j) = -ux*phi_dx -uy*phi_dy &
                                   + kappa*(phi_dxdx + phi_dydy) &
                                   - chi*(phi(ix,iy,1,j) - 0.0_rk) / C_eta_apply( int(mask(ix,iy,1,5), kind=2) ) & ! Dirichlet penalization for obstacle mask
                                   + source(ix, iy, 1) ! source term is actually a dirichlet penalization term as well
                end do
            end do

        endif
    end do ! loop over scalars

end subroutine RHS_2D_scalar
