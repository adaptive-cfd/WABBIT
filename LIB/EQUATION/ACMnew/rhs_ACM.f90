!-----------------------------------------------------------------------------
! main level wrapper to set the right hand side on a block. Note this is completely
! independent of the grid and any MPI formalism, neighboring relations and the like.
! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
! from that. Ghost nodes are assumed to be sync'ed.
!-----------------------------------------------------------------------------
subroutine RHS_ACM( time, u, g, x0, dx, rhs, mask, stage, n_domain )
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
    ! NOTE: ACM only supports symmetry BC for the moment (which is handled by wabbit and not ACM)
    integer(kind=2), intent(in) :: n_domain(3)

    real(kind=rk), allocatable, save :: vor(:,:,:,:)
    integer(kind=ik) :: mpierr, i, dim, ix, iy, iz
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: tmp(1:3), tmp2, dV, dV2, penal(1:3), C_eta_inv, x, y, z, f_block(1:3)
    integer(kind=2) :: color

    if (.not. params_acm%initialized) write(*,*) "WARNING: RHS_ACM called but ACM not initialized"

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    dim = params_acm%dim

    ! EXPERIMENTAL. Soft penalization: gently turn on penalization at the beginning (t=0) of a simulation. This makes the initial wave 
    ! travelling through the domain smoother. Its thickness is usually C_eta*C_0 (dimension L). Slowly decreasing C_eta from C_eta_start=1.0
    ! to its final value (small, say C_eta=1e-4) thus drastically increases the shock width. The shock intensity decays as R^-1 (2D) or 
    ! R^-2 (3D), but even in 3D it can still require resolution of the shock wave that can be expensive.
    if (params_acm%soft_penalization_startup) then
        if (time < params_acm%penalization_startup_tau) then
            ! the formulation with EXP seems to be more efficient in damping the initial shock wave; 2D tests showed that nicely.
            ! The new defaults are C_eta_start = 1.0  and  penalization_startup_tau=0.20. Note relatively long before this time the penalization
            ! is (almost) completely activated because of the exponential decay
            params_acm%C_eta = params_acm%C_eta_const + params_acm%C_eta_start*exp(-20.0*time/params_acm%penalization_startup_tau)
        else
            ! constant value of C_eta after the startup phase
            params_acm%C_eta = params_acm%C_eta_const    
        endif
    else
        ! constant value of C_eta
        params_acm%C_eta = params_acm%C_eta_const
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
        if (params_acm%HIT_linear_forcing) then
            params_acm%e_kin = 0.0_rk
            params_acm%enstrophy = 0.0_rk
            params_acm%mean_flow = 0.0_rk
        endif

        if (params_acm%geometry == "Insect") call Update_Insect(time, Insect)

        if (params_acm%use_free_flight_solver) then
            ! reset forces, we compute their value now
            params_acm%force_insect_g = 0.0_rk
            params_acm%moment_insect_g = 0.0_rk
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
            if (maxval(abs(u(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,i))) > 1.0e4_rk) then
                write(*,'("maxval in u(",i2,") = ", es15.8)') i, maxval(abs(u(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,i)))

                ! done by all ranks but well I hope the cluster can take one for the team.
                ! This (empty) file is for scripting purposes on the supercomputers.
                open (77, file='ACM_diverged', status='replace')
                close(77)
                
                call abort(0409201933,"ACM fail: very very large values in state vector.")
            endif
        enddo

        ! Linear Forcing for HIT (Lundgren) requires us to know kinetic energy and dissipation
        ! rate at all times, so compute that, if we use the forcing.
        if (params_acm%HIT_linear_forcing) then
            ! vorticity work array
            if (.not. allocated(vor) ) allocate(vor(1:size(u,1), 1:size(u,2), 1:size(u,3), 1:3 ))
            ! to compute the current dissipation rate
            call compute_vorticity(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, vor(:,:,:,:))

            dV = product(dx(1:params_acm%dim))

            params_acm%mean_flow(1) = dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1))
            params_acm%mean_flow(2) = dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 2))
            if (params_acm%dim==2) then
                params_acm%e_kin = params_acm%e_kin + 0.5_rk*dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g,    :       , 1:params_acm%dim)**2)
                params_acm%enstrophy = params_acm%enstrophy + 0.5_rk*dv*sum(vor(g+1:Bs(1)+g, g+1:Bs(2)+g,    :       , 1)**2)
            else
                params_acm%e_kin = params_acm%e_kin + 0.5_rk*dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:params_acm%dim)**2)
                params_acm%enstrophy = params_acm%enstrophy + 0.5_rk*dv*sum(vor(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:3)**2)
                params_acm%mean_flow(3) = dv*sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 3))
            endif

            ! NOTE: MPI_SUM is perfomed in the post_stage.
        endif


        ! if (params_acm%geometry == "Insect".and. params_acm%use_free_flight_solver) then
        if (params_acm%use_free_flight_solver) then
            if (dim == 2) then
                dV = dx(1)*dx(2)
                C_eta_inv = 1.0_rk / params_acm%C_eta
                f_block = 0.0_rk


                do iy = g+1, Bs(2)+g
                    y = x0(2) + dble(iy-(g+1)) * dx(2) - Insect%xc_body_g(2)
                    do ix = g+1, Bs(1)+g
                        x = x0(1) + dble(ix-(g+1)) * dx(1) - Insect%xc_body_g(1)

                        ! get this points color
                        color = int( mask(ix, iy, 1, 5), kind=2 )

                        ! exclude walls, trees, etc... (they have color 0)
                        ! if (color>0_2 .and. color < 6_2) then
                        ! penalization term
                        penal = -mask(ix,iy,1,1) * (u(ix,iy,1,1:3) - mask(ix,iy,1,2:4)) * C_eta_inv

                        f_block = f_block - penal
                        f_block(3) = 0.0_rk

                        ! x_lev = periodize_coordinate(x_lev, (/xl,yl,zl/))

                        ! moments. For insects, we compute the total moment wrt to the body center
                        ! params_acm%moment_insect_g = params_acm%moment_insect_g - cross((/x, y, z/), penal)*dV
                        ! endif
                    enddo
                enddo

                params_acm%force_insect_g = params_acm%force_insect_g + f_block*dV
            else
                dV = dx(1)*dx(2)*dx(3)
                C_eta_inv = 1.0_rk / params_acm%C_eta
                f_block = 0.0_rk

                do iz = g+1, Bs(3)+g
                    z = x0(3) + dble(iz-(g+1)) * dx(3) - Insect%xc_body_g(3) ! note: x-xc insect
                    do iy = g+1, Bs(2)+g
                        y = x0(2) + dble(iy-(g+1)) * dx(2) - Insect%xc_body_g(2)
                        do ix = g+1, Bs(1)+g
                            x = x0(1) + dble(ix-(g+1)) * dx(1) - Insect%xc_body_g(1)

                            ! get this points color
                            color = int( mask(ix, iy, iz, 5), kind=2 )

                            ! exclude walls, trees, etc... (they have color 0)
                            if (color>0_2 .and. color < 6_2) then
                                ! penalization term
                                penal = -mask(ix,iy,iz,1) * (u(ix,iy,iz,1:3) - mask(ix,iy,iz,2:4)) * C_eta_inv

                                f_block = f_block - penal

                                ! x_lev = periodize_coordinate(x_lev, (/xl,yl,zl/))

                                ! moments. For insects, we compute the total moment wrt to the body center
                                params_acm%moment_insect_g = params_acm%moment_insect_g - cross((/x, y, z/), penal)*dV
                            endif
                        enddo
                    enddo
                enddo
                params_acm%force_insect_g = params_acm%force_insect_g + f_block*dV
            endif ! NOTE: MPI_SUM is perfomed in the post_stage.

        endif

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.

        ! Linear Forcing for HIT (Lundgren) requires us to know kinetic energy and dissipation
        ! rate at all times, so compute that, if we use the forcing.
        if (params_acm%HIT_linear_forcing) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            if (params_acm%dim == 2) then
                params_acm%mean_flow = params_acm%mean_flow / (params_acm%domain_size(1)*params_acm%domain_size(2))
            else
                params_acm%mean_flow = params_acm%mean_flow / (params_acm%domain_size(1)*params_acm%domain_size(2)*params_acm%domain_size(3))
            endif
            params_acm%dissipation = params_acm%enstrophy * params_acm%nu
        endif

        if (params_acm%use_free_flight_solver) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%force_insect_g, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%moment_insect_g, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
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
    character(len=cshort), intent(in)       :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    !> forcing term
    real(kind=rk), dimension(3) :: forcing
    !>
    real(kind=rk) :: dx_inv, dy_inv, dx2_inv, dy2_inv, c_0, nu, C_eta, C_eta_inv, gamma
    real(kind=rk) :: div_U, u_dx, u_dy, u_dxdx, u_dydy, v_dx, v_dy, v_dxdx, &
                     v_dydy, p_dx, p_dy, penalx, penaly, x, y, term_2, spo, p_dxdx, p_dydy, nu_p, &
                     u_dx4, v_dx4, u_dy4, v_dy4, &
                     uu_dx, uv_dy, uw_dz, vu_dx, vv_dy, vw_dz, wu_dx, wv_dy, ww_dz, &
                     C_sponge_inv, &
                     up_dx, vp_dy
    ! loop variables
    integer(kind=ik) :: ix, iy, idir
    ! coefficients for Tam&Webb (4th order 1st derivative)
    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b_FD4(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! 6th order FD scheme
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/) ! 1st derivative
    real(kind=rk), parameter :: b_FD6(-3:3) = (/ 1.0_rk/90.0_rk, -3.0_rk/20.0_rk, 3.0_rk/2.0_rk, -49.0_rk/18.0_rk, 3.0_rk/2.0_rk, -3.0_rk/20.0_rk, 1.0_rk/90.0_rk/) ! 2nd derivative


    ! set parameters for readability
    c_0         = params_acm%c_0
    nu          = params_acm%nu
    nu_p        = params_acm%nu_p
    C_eta       = params_acm%C_eta
    gamma       = params_acm%gamma_p

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    C_eta_inv = 1.0_rk / C_eta

    if (size(phi,1)/=Bs(1)+2*g .or. size(phi,2)/=Bs(2)+2*g .or. size(phi,3)/=params_acm%dim+1+params_acm%N_scalars) then
        call abort(66233,"wrong size, I go for a walk instead.")
    endif


    select case(order_discretization)
    case("FD_2nd_central")
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        if (params_acm%skew_symmetry) then
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    u_dx   = (phi(ix+1,iy,1) - phi(ix-1,iy,1))*dx_inv*0.5_rk
                    v_dx   = (phi(ix+1,iy,2) - phi(ix-1,iy,2))*dx_inv*0.5_rk
                    
                    u_dy   = (phi(ix,iy+1,1) - phi(ix,iy-1,1))*dy_inv*0.5_rk
                    v_dy   = (phi(ix,iy+1,2) - phi(ix,iy-1,2))*dy_inv*0.5_rk

                    uu_dx = (phi(ix+1,iy,1)*phi(ix+1,iy,1) - phi(ix-1,iy,1)*phi(ix-1,iy,1))*dx_inv*0.5_rk
                    uv_dy = (phi(ix,iy+1,1)*phi(ix,iy+1,2) - phi(ix,iy-1,1)*phi(ix,iy-1,2))*dy_inv*0.5_rk

                    vu_dx = (phi(ix+1,iy,2)*phi(ix+1,iy,1) - phi(ix-1,iy,2)*phi(ix-1,iy,1))*dx_inv*0.5_rk
                    vv_dy = (phi(ix,iy+1,2)*phi(ix,iy+1,2) - phi(ix,iy-1,2)*phi(ix,iy-1,2))*dy_inv*0.5_rk

                    u_dxdx = (phi(ix-1,iy,1) -2.0_rk*phi(ix,iy,1) +phi(ix+1,iy,1))*dx2_inv
                    v_dxdx = (phi(ix-1,iy,2) -2.0_rk*phi(ix,iy,2) +phi(ix+1,iy,2))*dx2_inv
                    
                    u_dydy = (phi(ix,iy-1,1) -2.0_rk*phi(ix,iy,1) +phi(ix,iy+1,1))*dy2_inv
                    v_dydy = (phi(ix,iy-1,2) -2.0_rk*phi(ix,iy,2) +phi(ix,iy+1,2))*dy2_inv

                    p_dx   = (phi(ix+1,iy,3) -phi(ix-1,iy,3))*dx_inv*0.5_rk
                    p_dy   = (phi(ix,iy+1,3) -phi(ix,iy-1,3))*dy_inv*0.5_rk

                    div_U = u_dx + v_dy

                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                    ! Rhs in skew-symmetric formulation (of nonlinear term)
                    ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                    rhs(ix,iy,1) = -0.5_rk*(uu_dx + uv_dy   + phi(ix,iy,1)*u_dx + phi(ix,iy,2)*u_dy ) -p_dx + nu*(u_dxdx + u_dydy ) + penalx
                    rhs(ix,iy,2) = -0.5_rk*(vu_dx + vv_dy   + phi(ix,iy,1)*v_dx + phi(ix,iy,2)*v_dy ) -p_dy + nu*(v_dxdx + v_dydy ) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        else
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    u_dx   = (phi(ix+1,iy  ,1) - phi(ix-1,iy  ,1))*dx_inv*0.5_rk
                    u_dy   = (phi(ix  ,iy+1,1) - phi(ix  ,iy-1,1))*dy_inv*0.5_rk

                    v_dx   = (phi(ix+1,iy  ,2) - phi(ix-1,iy  ,2))*dx_inv*0.5_rk
                    v_dy   = (phi(ix  ,iy+1,2) - phi(ix  ,iy-1,2))*dy_inv*0.5_rk

                    u_dxdx = (phi(ix-1,iy  ,1) -2.0_rk*phi(ix,iy,1) +phi(ix+1,iy  ,1))*dx2_inv
                    u_dydy = (phi(ix  ,iy-1,1) -2.0_rk*phi(ix,iy,1) +phi(ix  ,iy+1,1))*dy2_inv

                    v_dxdx = (phi(ix-1,iy  ,2) -2.0_rk*phi(ix,iy,2) +phi(ix+1,iy  ,2))*dx2_inv
                    v_dydy = (phi(ix  ,iy-1,2) -2.0_rk*phi(ix,iy,2) +phi(ix  ,iy+1,2))*dy2_inv

                    p_dx   = (phi(ix+1,iy  ,3) -phi(ix-1,iy  ,3))*dx_inv*0.5_rk
                    p_dy   = (phi(ix  ,iy+1,3) -phi(ix  ,iy-1,3))*dy_inv*0.5_rk

                    div_U = u_dx + v_dy

                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                    ! actual RHS. note mean flow forcing is just a constant and added at the end of the routine
                    rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                    rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        endif

    case("FD_4th_central_optimized")
        !-----------------------------------------------------------------------
        ! 4th order (Tam&web optimized scheme)
        !-----------------------------------------------------------------------
        if (params_acm%skew_symmetry) then
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    ! Note: a(0) does NOT appear (it is zero...)
                    u_dx = (a_TW4(-3)*phi(ix-3,iy,1) +a_TW4(-2)*phi(ix-2,iy,1) +a_TW4(-1)*phi(ix-1,iy,1) +a_TW4(+1)*phi(ix+1,iy,1) +a_TW4(+2)*phi(ix+2,iy,1) +a_TW4(+3)*phi(ix+3,iy,1))*dx_inv
                    p_dx = (a_TW4(-3)*phi(ix-3,iy,3) +a_TW4(-2)*phi(ix-2,iy,3) +a_TW4(-1)*phi(ix-1,iy,3) +a_TW4(+1)*phi(ix+1,iy,3) +a_TW4(+2)*phi(ix+2,iy,3) +a_TW4(+3)*phi(ix+3,iy,3))*dx_inv
                    v_dx = (a_TW4(-3)*phi(ix-3,iy,2) +a_TW4(-2)*phi(ix-2,iy,2) +a_TW4(-1)*phi(ix-1,iy,2) +a_TW4(+1)*phi(ix+1,iy,2) +a_TW4(+2)*phi(ix+2,iy,2) +a_TW4(+3)*phi(ix+3,iy,2))*dx_inv

                    u_dy = (a_TW4(-3)*phi(ix,iy-3,1) +a_TW4(-2)*phi(ix,iy-2,1) +a_TW4(-1)*phi(ix,iy-1,1) +a_TW4(+1)*phi(ix,iy+1,1) +a_TW4(+2)*phi(ix,iy+2,1) +a_TW4(+3)*phi(ix,iy+3,1))*dy_inv
                    v_dy = (a_TW4(-3)*phi(ix,iy-3,2) +a_TW4(-2)*phi(ix,iy-2,2) +a_TW4(-1)*phi(ix,iy-1,2) +a_TW4(+1)*phi(ix,iy+1,2) +a_TW4(+2)*phi(ix,iy+2,2) +a_TW4(+3)*phi(ix,iy+3,2))*dy_inv
                    p_dy = (a_TW4(-3)*phi(ix,iy-3,3) +a_TW4(-2)*phi(ix,iy-2,3) +a_TW4(-1)*phi(ix,iy-1,3) +a_TW4(+1)*phi(ix,iy+1,3) +a_TW4(+2)*phi(ix,iy+2,3) +a_TW4(+3)*phi(ix,iy+3,3))*dy_inv

                    uu_dx = (a_TW4(-3)*phi(ix-3,iy,1)*phi(ix-3,iy,1) +a_TW4(-2)*phi(ix-2,iy,1)*phi(ix-2,iy,1) +a_TW4(-1)*phi(ix-1,iy,1)*phi(ix-1,iy,1) +a_TW4(+1)*phi(ix+1,iy,1)*phi(ix+1,iy,1) +a_TW4(+2)*phi(ix+2,iy,1)*phi(ix+2,iy,1) +a_TW4(+3)*phi(ix+3,iy,1)*phi(ix+3,iy,1))*dx_inv
                    vu_dx = (a_TW4(-3)*phi(ix-3,iy,2)*phi(ix-3,iy,1) +a_TW4(-2)*phi(ix-2,iy,2)*phi(ix-2,iy,1) +a_TW4(-1)*phi(ix-1,iy,2)*phi(ix-1,iy,1) +a_TW4(+1)*phi(ix+1,iy,2)*phi(ix+1,iy,1) +a_TW4(+2)*phi(ix+2,iy,2)*phi(ix+2,iy,1) +a_TW4(+3)*phi(ix+3,iy,2)*phi(ix+3,iy,1))*dx_inv

                    uv_dy = (a_TW4(-3)*phi(ix,iy-3,1)*phi(ix,iy-3,2) +a_TW4(-2)*phi(ix,iy-2,1)*phi(ix,iy-2,2) +a_TW4(-1)*phi(ix,iy-1,1)*phi(ix,iy-1,2) +a_TW4(+1)*phi(ix,iy+1,1)*phi(ix,iy+1,2) +a_TW4(+2)*phi(ix,iy+2,1)*phi(ix,iy+2,2) +a_TW4(+3)*phi(ix,iy+3,1)*phi(ix,iy+3,2))*dy_inv
                    vv_dy = (a_TW4(-3)*phi(ix,iy-3,2)*phi(ix,iy-3,2) +a_TW4(-2)*phi(ix,iy-2,2)*phi(ix,iy-2,2) +a_TW4(-1)*phi(ix,iy-1,2)*phi(ix,iy-1,2) +a_TW4(+1)*phi(ix,iy+1,2)*phi(ix,iy+1,2) +a_TW4(+2)*phi(ix,iy+2,2)*phi(ix,iy+2,2) +a_TW4(+3)*phi(ix,iy+3,2)*phi(ix,iy+3,2))*dy_inv

                    ! second derivatives of u and v
                    u_dxdx = (b_FD4(-2)*phi(ix-2,iy,1) + b_FD4(-1)*phi(ix-1,iy,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix+1,iy,1) + b_FD4(+2)*phi(ix+2,iy,1))*dx2_inv
                    v_dxdx = (b_FD4(-2)*phi(ix-2,iy,2) + b_FD4(-1)*phi(ix-1,iy,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix+1,iy,2) + b_FD4(+2)*phi(ix+2,iy,2))*dx2_inv

                    u_dydy = (b_FD4(-2)*phi(ix,iy-2,1) + b_FD4(-1)*phi(ix,iy-1,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix,iy+1,1) + b_FD4(+2)*phi(ix,iy+2,1))*dy2_inv
                    v_dydy = (b_FD4(-2)*phi(ix,iy-2,2) + b_FD4(-1)*phi(ix,iy-1,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix,iy+1,2) + b_FD4(+2)*phi(ix,iy+2,2))*dy2_inv
                    div_U = u_dx + v_dy

                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                    rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                    rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)

                    ! Rhs in skew-symmetric formulation (of nonlinear term)
                    ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                    rhs(ix,iy,1) = -0.5_rk*(uu_dx + uv_dy   + phi(ix,iy,1)*u_dx + phi(ix,iy,2)*u_dy ) -p_dx + nu*(u_dxdx + u_dydy ) + penalx
                    rhs(ix,iy,2) = -0.5_rk*(vu_dx + vv_dy   + phi(ix,iy,1)*v_dx + phi(ix,iy,2)*v_dy ) -p_dy + nu*(v_dxdx + v_dydy ) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        else
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    ! Note: a(0) does NOT appear (it is zero...)
                    u_dx = (a_TW4(-3)*phi(ix-3,iy,1) + a_TW4(-2)*phi(ix-2,iy,1) + a_TW4(-1)*phi(ix-1,iy,1) + a_TW4(+1)*phi(ix+1,iy,1) + a_TW4(+2)*phi(ix+2,iy,1) + a_TW4(+3)*phi(ix+3,iy,1))*dx_inv
                    p_dx = (a_TW4(-3)*phi(ix-3,iy,3) + a_TW4(-2)*phi(ix-2,iy,3) + a_TW4(-1)*phi(ix-1,iy,3) + a_TW4(+1)*phi(ix+1,iy,3) + a_TW4(+2)*phi(ix+2,iy,3) + a_TW4(+3)*phi(ix+3,iy,3))*dx_inv
                    v_dx = (a_TW4(-3)*phi(ix-3,iy,2) + a_TW4(-2)*phi(ix-2,iy,2) + a_TW4(-1)*phi(ix-1,iy,2) + a_TW4(+1)*phi(ix+1,iy,2) + a_TW4(+2)*phi(ix+2,iy,2) + a_TW4(+3)*phi(ix+3,iy,2))*dx_inv

                    u_dy = (a_TW4(-3)*phi(ix,iy-3,1) + a_TW4(-2)*phi(ix,iy-2,1) + a_TW4(-1)*phi(ix,iy-1,1) + a_TW4(+1)*phi(ix,iy+1,1) + a_TW4(+2)*phi(ix,iy+2,1) + a_TW4(+3)*phi(ix,iy+3,1))*dy_inv
                    v_dy = (a_TW4(-3)*phi(ix,iy-3,2) + a_TW4(-2)*phi(ix,iy-2,2) + a_TW4(-1)*phi(ix,iy-1,2) + a_TW4(+1)*phi(ix,iy+1,2) + a_TW4(+2)*phi(ix,iy+2,2) + a_TW4(+3)*phi(ix,iy+3,2))*dy_inv
                    p_dy = (a_TW4(-3)*phi(ix,iy-3,3) + a_TW4(-2)*phi(ix,iy-2,3) + a_TW4(-1)*phi(ix,iy-1,3) + a_TW4(+1)*phi(ix,iy+1,3) + a_TW4(+2)*phi(ix,iy+2,3) + a_TW4(+3)*phi(ix,iy+3,3))*dy_inv

                    ! second derivatives of u and v
                    u_dxdx = (b_FD4(-2)*phi(ix-2,iy,1) + b_FD4(-1)*phi(ix-1,iy,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix+1,iy,1) + b_FD4(+2)*phi(ix+2,iy,1))*dx2_inv
                    v_dxdx = (b_FD4(-2)*phi(ix-2,iy,2) + b_FD4(-1)*phi(ix-1,iy,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix+1,iy,2) + b_FD4(+2)*phi(ix+2,iy,2))*dx2_inv

                    u_dydy = (b_FD4(-2)*phi(ix,iy-2,1) + b_FD4(-1)*phi(ix,iy-1,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix,iy+1,1) + b_FD4(+2)*phi(ix,iy+2,1))*dy2_inv
                    v_dydy = (b_FD4(-2)*phi(ix,iy-2,2) + b_FD4(-1)*phi(ix,iy-1,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix,iy+1,2) + b_FD4(+2)*phi(ix,iy+2,2))*dy2_inv
                    div_U = u_dx + v_dy

                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                    rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                    rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        endif

    case("FD_4th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard)
        !-----------------------------------------------------------------------
        if (params_acm%skew_symmetry) then
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    ! Note: a(0) does NOT appear (it is zero...)
                    u_dx = (a_FD4(-2)*phi(ix-2,iy,1) + a_FD4(-1)*phi(ix-1,iy,1) + a_FD4(+1)*phi(ix+1,iy,1) + a_FD4(+2)*phi(ix+2,iy,1))*dx_inv
                    v_dx = (a_FD4(-2)*phi(ix-2,iy,2) + a_FD4(-1)*phi(ix-1,iy,2) + a_FD4(+1)*phi(ix+1,iy,2) + a_FD4(+2)*phi(ix+2,iy,2))*dx_inv
                    p_dx = (a_FD4(-2)*phi(ix-2,iy,3) + a_FD4(-1)*phi(ix-1,iy,3) + a_FD4(+1)*phi(ix+1,iy,3) + a_FD4(+2)*phi(ix+2,iy,3))*dx_inv

                    u_dy = (a_FD4(-2)*phi(ix,iy-2,1) + a_FD4(-1)*phi(ix,iy-1,1) + a_FD4(+1)*phi(ix,iy+1,1) + a_FD4(+2)*phi(ix,iy+2,1))*dy_inv
                    v_dy = (a_FD4(-2)*phi(ix,iy-2,2) + a_FD4(-1)*phi(ix,iy-1,2) + a_FD4(+1)*phi(ix,iy+1,2) + a_FD4(+2)*phi(ix,iy+2,2))*dy_inv
                    p_dy = (a_FD4(-2)*phi(ix,iy-2,3) + a_FD4(-1)*phi(ix,iy-1,3) + a_FD4(+1)*phi(ix,iy+1,3) + a_FD4(+2)*phi(ix,iy+2,3))*dy_inv

                    uu_dx = (a_FD4(-2)*phi(ix-2,iy,1)*phi(ix-2,iy,1) + a_FD4(-1)*phi(ix-1,iy,1)*phi(ix-1,iy,1) + a_FD4(+1)*phi(ix+1,iy,1)*phi(ix+1,iy,1) + a_FD4(+2)*phi(ix+2,iy,1)*phi(ix+2,iy,1))*dx_inv
                    uv_dy = (a_FD4(-2)*phi(ix,iy-2,1)*phi(ix,iy-2,2) + a_FD4(-1)*phi(ix,iy-1,1)*phi(ix,iy-1,2) + a_FD4(+1)*phi(ix,iy+1,1)*phi(ix,iy+1,2) + a_FD4(+2)*phi(ix,iy+2,1)*phi(ix,iy+2,2))*dy_inv

                    vu_dx = (a_FD4(-2)*phi(ix-2,iy,2)*phi(ix-2,iy,1) + a_FD4(-1)*phi(ix-1,iy,2)*phi(ix-1,iy,1) + a_FD4(+1)*phi(ix+1,iy,2)*phi(ix+1,iy,1) + a_FD4(+2)*phi(ix+2,iy,2)*phi(ix+2,iy,1))*dx_inv
                    vv_dy = (a_FD4(-2)*phi(ix,iy-2,2)*phi(ix,iy-2,2) + a_FD4(-1)*phi(ix,iy-1,2)*phi(ix,iy-1,2) + a_FD4(+1)*phi(ix,iy+1,2)*phi(ix,iy+1,2) + a_FD4(+2)*phi(ix,iy+2,2)*phi(ix,iy+2,2))*dy_inv

                    ! second derivatives of u and v
                    u_dxdx = (b_FD4(-2)*phi(ix-2,iy,1) + b_FD4(-1)*phi(ix-1,iy,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix+1,iy,1) + b_FD4(+2)*phi(ix+2,iy,1))*dx2_inv
                    v_dxdx = (b_FD4(-2)*phi(ix-2,iy,2) + b_FD4(-1)*phi(ix-1,iy,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix+1,iy,2) + b_FD4(+2)*phi(ix+2,iy,2))*dx2_inv

                    u_dydy = (b_FD4(-2)*phi(ix,iy-2,1) + b_FD4(-1)*phi(ix,iy-1,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix,iy+1,1) + b_FD4(+2)*phi(ix,iy+2,1))*dy2_inv
                    v_dydy = (b_FD4(-2)*phi(ix,iy-2,2) + b_FD4(-1)*phi(ix,iy-1,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix,iy+1,2) + b_FD4(+2)*phi(ix,iy+2,2))*dy2_inv

                    div_U = u_dx + v_dy
                    ! mask (chi) (:,:,1)
                    ! usx (:,:,2)
                    ! usy (:,:,3)
                    ! usz (:,:,4)
                    ! color (.;.5)
                    ! sponge (:,:,6)
                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                    ! Rhs in skew-symmetric formulation (of nonlinear term)
                    ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                    rhs(ix,iy,1) = -0.5_rk*(uu_dx + uv_dy   + phi(ix,iy,1)*u_dx + phi(ix,iy,2)*u_dy ) -p_dx + nu*(u_dxdx + u_dydy ) + penalx
                    rhs(ix,iy,2) = -0.5_rk*(vu_dx + vv_dy   + phi(ix,iy,1)*v_dx + phi(ix,iy,2)*v_dy ) -p_dy + nu*(v_dxdx + v_dydy ) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        else
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    ! Note: a(0) does NOT appear (it is zero...)
                    u_dx = (a_FD4(-2)*phi(ix-2,iy,1) + a_FD4(-1)*phi(ix-1,iy,1) + a_FD4(+1)*phi(ix+1,iy,1) + a_FD4(+2)*phi(ix+2,iy,1))*dx_inv    
                    v_dx = (a_FD4(-2)*phi(ix-2,iy,2) + a_FD4(-1)*phi(ix-1,iy,2) + a_FD4(+1)*phi(ix+1,iy,2) + a_FD4(+2)*phi(ix+2,iy,2))*dx_inv    
                    p_dx = (a_FD4(-2)*phi(ix-2,iy,3) + a_FD4(-1)*phi(ix-1,iy,3) + a_FD4(+1)*phi(ix+1,iy,3) + a_FD4(+2)*phi(ix+2,iy,3))*dx_inv   

                    u_dy = (a_FD4(-2)*phi(ix,iy-2,1) + a_FD4(-1)*phi(ix,iy-1,1) + a_FD4(+1)*phi(ix,iy+1,1) + a_FD4(+2)*phi(ix,iy+2,1))*dy_inv    
                    v_dy = (a_FD4(-2)*phi(ix,iy-2,2) + a_FD4(-1)*phi(ix,iy-1,2) + a_FD4(+1)*phi(ix,iy+1,2) + a_FD4(+2)*phi(ix,iy+2,2))*dy_inv    
                    p_dy = (a_FD4(-2)*phi(ix,iy-2,3) + a_FD4(-1)*phi(ix,iy-1,3) + a_FD4(+1)*phi(ix,iy+1,3) + a_FD4(+2)*phi(ix,iy+2,3))*dy_inv    

                    ! second derivatives of u and v
                    u_dxdx = (b_FD4(-2)*phi(ix-2,iy,1) + b_FD4(-1)*phi(ix-1,iy,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix+1,iy,1) + b_FD4(+2)*phi(ix+2,iy,1))*dx2_inv    
                    v_dxdx = (b_FD4(-2)*phi(ix-2,iy,2) + b_FD4(-1)*phi(ix-1,iy,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix+1,iy,2) + b_FD4(+2)*phi(ix+2,iy,2))*dx2_inv   
                    
                    u_dydy = (b_FD4(-2)*phi(ix,iy-2,1) + b_FD4(-1)*phi(ix,iy-1,1) + b_FD4(0)*phi(ix,iy,1) + b_FD4(+1)*phi(ix,iy+1,1) + b_FD4(+2)*phi(ix,iy+2,1))*dy2_inv    
                    v_dydy = (b_FD4(-2)*phi(ix,iy-2,2) + b_FD4(-1)*phi(ix,iy-1,2) + b_FD4(0)*phi(ix,iy,2) + b_FD4(+1)*phi(ix,iy+1,2) + b_FD4(+2)*phi(ix,iy+2,2))*dy2_inv    
                    
                    div_U = u_dx + v_dy

                    ! mask (chi) (:,:,1)
                    ! usx (:,:,2)
                    ! usy (:,:,3)
                    ! usz (:,:,4)
                    ! color (.;.5)
                    ! sponge (:,:,6)
                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))
    
                    rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                    rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        endif

        select case(params_acm%p_eqn_model)
        case ('acm')
            ! do nothing, is the eqn computed above without additional terms
        case ('diffusive')
            ! add p diffusion []
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    p_dxdx = (b_FD4(-2)*phi(ix-2,iy,3) + b_FD4(-1)*phi(ix-1,iy,3) + b_FD4(0)*phi(ix,iy,3) + b_FD4(+1)*phi(ix+1,iy,3) + b_FD4(+2)*phi(ix+2,iy,3))*dx2_inv
                    p_dydy = (b_FD4(-2)*phi(ix,iy-2,3) + b_FD4(-1)*phi(ix,iy-1,3) + b_FD4(0)*phi(ix,iy,3) + b_FD4(+1)*phi(ix,iy+1,3) + b_FD4(+2)*phi(ix,iy+2,3))*dy2_inv
    
                    rhs(ix,iy,3) = rhs(ix,iy,3) + params_acm%nu_p*(p_dxdx+p_dydy)
                enddo
            enddo   
        case ('EDAC')
            ! add EDAC term AND diffusion [Clausen, Entropically damped form of artificial compressibility for explicit simulation of incompressible flow, 2013]
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    p_dxdx = (b_FD4(-2)*phi(ix-2,iy,3) + b_FD4(-1)*phi(ix-1,iy,3) + b_FD4(0)*phi(ix,iy,3) + b_FD4(+1)*phi(ix+1,iy,3) + b_FD4(+2)*phi(ix+2,iy,3))*dx2_inv
                    p_dydy = (b_FD4(-2)*phi(ix,iy-2,3) + b_FD4(-1)*phi(ix,iy-1,3) + b_FD4(0)*phi(ix,iy,3) + b_FD4(+1)*phi(ix,iy+1,3) + b_FD4(+2)*phi(ix,iy+2,3))*dy2_inv
    
                    up_dx = (a_FD4(-2)*phi(ix-2,iy,1)*phi(ix-2,iy,3) + a_FD4(-1)*phi(ix-1,iy,1)*phi(ix-1,iy,3) + a_FD4(+1)*phi(ix+1,iy,1)*phi(ix+1,iy,3) + a_FD4(+2)*phi(ix+2,iy,1)*phi(ix+2,iy,3))*dx_inv
                    vp_dy = (a_FD4(-2)*phi(ix,iy-2,2)*phi(ix,iy-2,3) + a_FD4(-1)*phi(ix,iy-1,2)*phi(ix,iy-1,3) + a_FD4(+1)*phi(ix,iy+1,2)*phi(ix,iy+1,3) + a_FD4(+2)*phi(ix,iy+2,2)*phi(ix,iy+2,3))*dy_inv

                    rhs(ix,iy,3) = rhs(ix,iy,3) + params_acm%nu_p*(p_dxdx+p_dydy) - up_dx - vp_dy
                enddo
            enddo   
        case ('convective')
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    up_dx = (a_FD4(-2)*phi(ix-2,iy,1)*phi(ix-2,iy,3) + a_FD4(-1)*phi(ix-1,iy,1)*phi(ix-1,iy,3) + a_FD4(+1)*phi(ix+1,iy,1)*phi(ix+1,iy,3) + a_FD4(+2)*phi(ix+2,iy,1)*phi(ix+2,iy,3))*dx_inv
                    vp_dy = (a_FD4(-2)*phi(ix,iy-2,2)*phi(ix,iy-2,3) + a_FD4(-1)*phi(ix,iy-1,2)*phi(ix,iy-1,3) + a_FD4(+1)*phi(ix,iy+1,2)*phi(ix,iy+1,3) + a_FD4(+2)*phi(ix,iy+2,2)*phi(ix,iy+2,3))*dy_inv

                    rhs(ix,iy,3) = rhs(ix,iy,3) - up_dx - vp_dy
                enddo
            enddo    
        case default
            call abort(2501041, "pressure equation model is unkown: "//trim(adjustl(params_acm%p_eqn_model)))
        end select    

    case("FD_6th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard)
        !-----------------------------------------------------------------------
        if (params_acm%skew_symmetry) then
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    ! Note: a(0) does NOT appear (it is zero...)
                    u_dx = (a_FD6(-3)*phi(ix-3,iy,1)+a_FD6(-2)*phi(ix-2,iy,1) +a_FD6(-1)*phi(ix-1,iy,1) +a_FD6(+1)*phi(ix+1,iy,1) +a_FD6(+2)*phi(ix+2,iy,1) +a_FD6(+3)*phi(ix+3,iy,1))*dx_inv
                    v_dx = (a_FD6(-3)*phi(ix-3,iy,2)+a_FD6(-2)*phi(ix-2,iy,2) +a_FD6(-1)*phi(ix-1,iy,2) +a_FD6(+1)*phi(ix+1,iy,2) +a_FD6(+2)*phi(ix+2,iy,2) +a_FD6(+3)*phi(ix+3,iy,2))*dx_inv
                    p_dx = (a_FD6(-3)*phi(ix-3,iy,3)+a_FD6(-2)*phi(ix-2,iy,3) +a_FD6(-1)*phi(ix-1,iy,3) +a_FD6(+1)*phi(ix+1,iy,3) +a_FD6(+2)*phi(ix+2,iy,3) +a_FD6(+3)*phi(ix+3,iy,3))*dx_inv

                    u_dy = (a_FD6(-3)*phi(ix,iy-3,1)+a_FD6(-2)*phi(ix,iy-2,1) +a_FD6(-1)*phi(ix,iy-1,1) +a_FD6(+1)*phi(ix,iy+1,1) +a_FD6(+2)*phi(ix,iy+2,1) +a_FD6(+3)*phi(ix,iy+3,1))*dy_inv
                    v_dy = (a_FD6(-3)*phi(ix,iy-3,2)+a_FD6(-2)*phi(ix,iy-2,2) +a_FD6(-1)*phi(ix,iy-1,2) +a_FD6(+1)*phi(ix,iy+1,2) +a_FD6(+2)*phi(ix,iy+2,2) +a_FD6(+3)*phi(ix,iy+3,2))*dy_inv
                    p_dy = (a_FD6(-3)*phi(ix,iy-3,3)+a_FD6(-2)*phi(ix,iy-2,3) +a_FD6(-1)*phi(ix,iy-1,3) +a_FD6(+1)*phi(ix,iy+1,3) +a_FD6(+2)*phi(ix,iy+2,3) +a_FD6(+3)*phi(ix,iy+3,3))*dy_inv

                    uu_dx = (a_FD6(-3)*phi(ix-3,iy,1)*phi(ix-3,iy,1) +a_FD6(-2)*phi(ix-2,iy,1)*phi(ix-2,iy,1) +a_FD6(-1)*phi(ix-1,iy,1)*phi(ix-1,iy,1) +a_FD6(+1)*phi(ix+1,iy,1)*phi(ix+1,iy,1) +a_FD6(+2)*phi(ix+2,iy,1)*phi(ix+2,iy,1) +a_FD6(+3)*phi(ix+3,iy,1)*phi(ix+3,iy,1))*dx_inv
                    vu_dx = (a_FD6(-3)*phi(ix-3,iy,2)*phi(ix-3,iy,1) +a_FD6(-2)*phi(ix-2,iy,2)*phi(ix-2,iy,1) +a_FD6(-1)*phi(ix-1,iy,2)*phi(ix-1,iy,1) +a_FD6(+1)*phi(ix+1,iy,2)*phi(ix+1,iy,1) +a_FD6(+2)*phi(ix+2,iy,2)*phi(ix+2,iy,1) +a_FD6(+3)*phi(ix+3,iy,2)*phi(ix+3,iy,1))*dx_inv

                    uv_dy = (a_FD6(-3)*phi(ix,iy-3,1)*phi(ix,iy-3,2) +a_FD6(-2)*phi(ix,iy-2,1)*phi(ix,iy-2,2) +a_FD6(-1)*phi(ix,iy-1,1)*phi(ix,iy-1,2) +a_FD6(+1)*phi(ix,iy+1,1)*phi(ix,iy+1,2) +a_FD6(+2)*phi(ix,iy+2,1)*phi(ix,iy+2,2) +a_FD6(+3)*phi(ix,iy+3,1)*phi(ix,iy+3,2))*dy_inv
                    vv_dy = (a_FD6(-3)*phi(ix,iy-3,2)*phi(ix,iy-3,2) +a_FD6(-2)*phi(ix,iy-2,2)*phi(ix,iy-2,2) +a_FD6(-1)*phi(ix,iy-1,2)*phi(ix,iy-1,2) +a_FD6(+1)*phi(ix,iy+1,2)*phi(ix,iy+1,2) +a_FD6(+2)*phi(ix,iy+2,2)*phi(ix,iy+2,2) +a_FD6(+3)*phi(ix,iy+3,2)*phi(ix,iy+3,2))*dy_inv

                    ! second derivatives of u and v
                    u_dxdx = (b_FD6(-3)*phi(ix-3,iy,1) +b_FD6(-2)*phi(ix-2,iy,1) +b_FD6(-1)*phi(ix-1,iy,1) +b_FD6( 0)*phi(ix,iy,1) +b_FD6(+1)*phi(ix+1,iy,1) +b_FD6(+2)*phi(ix+2,iy,1) +b_FD6(+3)*phi(ix+3,iy,1))*dx2_inv
                    v_dxdx = (b_FD6(-3)*phi(ix-3,iy,2) +b_FD6(-2)*phi(ix-2,iy,2) +b_FD6(-1)*phi(ix-1,iy,2) +b_FD6( 0)*phi(ix,iy,2) +b_FD6(+1)*phi(ix+1,iy,2) +b_FD6(+2)*phi(ix+2,iy,2) +b_FD6(+3)*phi(ix+3,iy,2))*dx2_inv
                    
                    u_dydy = (b_FD6(-3)*phi(ix,iy-3,1) +b_FD6(-2)*phi(ix,iy-2,1) +b_FD6(-1)*phi(ix,iy-1,1) +b_FD6( 0)*phi(ix,iy,1) +b_FD6(+1)*phi(ix,iy+1,1) +b_FD6(+2)*phi(ix,iy+2,1) +b_FD6(+3)*phi(ix,iy+3,1))*dy2_inv
                    v_dydy = (b_FD6(-3)*phi(ix,iy-3,2) +b_FD6(-2)*phi(ix,iy-2,2) +b_FD6(-1)*phi(ix,iy-1,2) +b_FD6( 0)*phi(ix,iy,2) +b_FD6(+1)*phi(ix,iy+1,2) +b_FD6(+2)*phi(ix,iy+2,2) +b_FD6(+3)*phi(ix,iy+3,2))*dy2_inv
                    div_U = u_dx + v_dy

                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                    ! Rhs in skew-symmetric formulation (of nonlinear term)
                    ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                    rhs(ix,iy,1) = -0.5_rk*(uu_dx + uv_dy   + phi(ix,iy,1)*u_dx + phi(ix,iy,2)*u_dy ) -p_dx + nu*(u_dxdx + u_dydy ) + penalx
                    rhs(ix,iy,2) = -0.5_rk*(vu_dx + vv_dy   + phi(ix,iy,1)*v_dx + phi(ix,iy,2)*v_dy ) -p_dy + nu*(v_dxdx + v_dydy ) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        else
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    ! Note: a(0) does NOT appear (it is zero...)
                    u_dx = (a_FD6(-3)*phi(ix-3,iy,1)+a_FD6(-2)*phi(ix-2,iy,1) +a_FD6(-1)*phi(ix-1,iy,1) +a_FD6(+1)*phi(ix+1,iy,1) +a_FD6(+2)*phi(ix+2,iy,1) +a_FD6(+3)*phi(ix+3,iy,1))*dx_inv
                    v_dx = (a_FD6(-3)*phi(ix-3,iy,2)+a_FD6(-2)*phi(ix-2,iy,2) +a_FD6(-1)*phi(ix-1,iy,2) +a_FD6(+1)*phi(ix+1,iy,2) +a_FD6(+2)*phi(ix+2,iy,2) +a_FD6(+3)*phi(ix+3,iy,2))*dx_inv
                    p_dx = (a_FD6(-3)*phi(ix-3,iy,3)+a_FD6(-2)*phi(ix-2,iy,3) +a_FD6(-1)*phi(ix-1,iy,3) +a_FD6(+1)*phi(ix+1,iy,3) +a_FD6(+2)*phi(ix+2,iy,3) +a_FD6(+3)*phi(ix+3,iy,3))*dx_inv

                    u_dy = (a_FD6(-3)*phi(ix,iy-3,1)+a_FD6(-2)*phi(ix,iy-2,1) +a_FD6(-1)*phi(ix,iy-1,1) +a_FD6(+1)*phi(ix,iy+1,1) +a_FD6(+2)*phi(ix,iy+2,1) +a_FD6(+3)*phi(ix,iy+3,1))*dy_inv
                    v_dy = (a_FD6(-3)*phi(ix,iy-3,2)+a_FD6(-2)*phi(ix,iy-2,2) +a_FD6(-1)*phi(ix,iy-1,2) +a_FD6(+1)*phi(ix,iy+1,2) +a_FD6(+2)*phi(ix,iy+2,2) +a_FD6(+3)*phi(ix,iy+3,2))*dy_inv
                    p_dy = (a_FD6(-3)*phi(ix,iy-3,3)+a_FD6(-2)*phi(ix,iy-2,3) +a_FD6(-1)*phi(ix,iy-1,3) +a_FD6(+1)*phi(ix,iy+1,3) +a_FD6(+2)*phi(ix,iy+2,3) +a_FD6(+3)*phi(ix,iy+3,3))*dy_inv

                    ! second derivatives of u and v
                    u_dxdx = (b_FD6(-3)*phi(ix-3,iy,1) +b_FD6(-2)*phi(ix-2,iy,1) +b_FD6(-1)*phi(ix-1,iy,1) +b_FD6( 0)*phi(ix,iy,1) +b_FD6(+1)*phi(ix+1,iy,1) +b_FD6(+2)*phi(ix+2,iy,1) +b_FD6(+3)*phi(ix+3,iy,1))*dx2_inv
                    v_dxdx = (b_FD6(-3)*phi(ix-3,iy,2) +b_FD6(-2)*phi(ix-2,iy,2) +b_FD6(-1)*phi(ix-1,iy,2) +b_FD6( 0)*phi(ix,iy,2) +b_FD6(+1)*phi(ix+1,iy,2) +b_FD6(+2)*phi(ix+2,iy,2) +b_FD6(+3)*phi(ix+3,iy,2))*dx2_inv
                    
                    u_dydy = (b_FD6(-3)*phi(ix,iy-3,1) +b_FD6(-2)*phi(ix,iy-2,1) +b_FD6(-1)*phi(ix,iy-1,1) +b_FD6( 0)*phi(ix,iy,1) +b_FD6(+1)*phi(ix,iy+1,1) +b_FD6(+2)*phi(ix,iy+2,1) +b_FD6(+3)*phi(ix,iy+3,1))*dy2_inv
                    v_dydy = (b_FD6(-3)*phi(ix,iy-3,2) +b_FD6(-2)*phi(ix,iy-2,2) +b_FD6(-1)*phi(ix,iy-1,2) +b_FD6( 0)*phi(ix,iy,2) +b_FD6(+1)*phi(ix,iy+1,2) +b_FD6(+2)*phi(ix,iy+2,2) +b_FD6(+3)*phi(ix,iy+3,2))*dy2_inv
                    div_U = u_dx + v_dy

                    penalx = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                    penaly = -mask(ix,iy,1) * C_eta_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                    rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                    rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                    rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
                end do
            end do
        endif

    case default
        call abort(441166, "Discretization unkown "//trim(adjustl(order_discretization))//", you should go play outside. Its nice." )

    end select

    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    ! HACK
    if (.not. params_acm%geometry == "lamballais") then
        if (params_acm%use_sponge) then
            ! avoid division by multiplying with inverse
            C_sponge_inv = 1.0_rk / params_acm%C_sponge

            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! NOTE: the sponge term acts, if active, on ALL components, ux,uy,p
                    ! which is different from the penalization term, which acts only on ux,uy and not p
                    spo = mask(ix,iy,6) * C_sponge_inv

                    rhs(ix,iy,1) = rhs(ix,iy,1) - (phi(ix,iy,1)-params_acm%u_mean_set(1)) * spo
                    rhs(ix,iy,2) = rhs(ix,iy,2) - (phi(ix,iy,2)-params_acm%u_mean_set(2)) * spo
                    rhs(ix,iy,3) = rhs(ix,iy,3) - phi(ix,iy,3)*spo
                enddo
            enddo
        end if
    else
        ! special treatment lamballais, also for p to p_ref in the ring
        ! ring stored in sponge mask array
        ! Gautier, R., Biau, D., Lamballais, E.: A reference solution of the flow over a circular cylinder at Re = 40 , Computers & Fluids 75, 103–111, 2013 

        ! use same C_eta as object, not the sponge value
        C_sponge_inv = 1.0_rk / params_acm%C_eta
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! mask(:,:,4) contains p_ref forcing in the ring
                ! rhs_p      = rhsp         - chi_ring     *( p          - p_lamballais) / C_eta
                rhs(ix,iy,3) = rhs(ix,iy,3) - mask(ix,iy,6)*(phi(ix,iy,3)-mask(ix,iy,4))*C_sponge_inv
            enddo
        enddo
    endif




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
    character(len=cshort), intent(in)       :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time


    !> forcing term
    real(kind=rk), dimension(3) :: forcing

    !> inverse dx, physics/acm parameters
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, c_0, &
                     nu, C_eta, C_eta_inv, gamma, spo, A_forcing, G_gain, t_l_inf, e_kin_set
    !> derivatives
    real(kind=rk) :: u_dx, u_dy, u_dz, u_dxdx, u_dydy, u_dzdz, &
                     v_dx, v_dy, v_dz, v_dxdx, v_dydy, v_dzdz, &
                     w_dx, w_dy, w_dz, w_dxdx, w_dydy, w_dzdz, &
                     p_dx, p_dy, p_dz, penalx, penaly, penalz, u, v, w, p, chi, &
                     uu_dx, uv_dy, uw_dz, vu_dx, vv_dy, vw_dz, wu_dx, wv_dy, ww_dz, &
                     C_sponge_inv, p_dxdx, p_dydy, p_dzdz, pu_dx, pv_dy, pw_dz
    !> loop variables
    integer(kind=ik) :: ix, iy, iz

    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b_FD4(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! 6th order FD scheme
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/) ! 1st derivative
    real(kind=rk), parameter :: b_FD6(-3:3) = (/ 1.0_rk/90.0_rk, -3.0_rk/20.0_rk, 3.0_rk/2.0_rk, -49.0_rk/18.0_rk, 3.0_rk/2.0_rk, -3.0_rk/20.0_rk, 1.0_rk/90.0_rk/) ! 2nd derivative


    ! set parameters for readability
    c_0         = params_acm%c_0
    nu          = params_acm%nu
    C_eta       = params_acm%C_eta
    gamma       = params_acm%gamma_p

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dz_inv = 1.0_rk / dx(3)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    dz2_inv = 1.0_rk / (dx(3)**2)

    C_eta_inv = 1.0_rk / C_eta

    select case(order_discretization)
    case("FD_2nd_central")
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
        if (params_acm%skew_symmetry) then
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

                        uu_dx = (phi(ix+1,iy,iz,1)*phi(ix+1,iy,iz,1) -phi(ix-1,iy,iz,1)*phi(ix-1,iy,iz,1) )*dx_inv*0.5_rk
                        uv_dy = (phi(ix,iy+1,iz,1)*phi(ix,iy+1,iz,2) -phi(ix,iy-1,iz,1)*phi(ix,iy-1,iz,2) )*dy_inv*0.5_rk
                        uw_dz = (phi(ix,iy,iz+1,1)*phi(ix,iy,iz+1,3) -phi(ix,iy,iz-1,1)*phi(ix,iy,iz-1,3) )*dz_inv*0.5_rk

                        vu_dx = (phi(ix+1,iy,iz,2)*phi(ix+1,iy,iz,1) -phi(ix-1,iy,iz,2)*phi(ix-1,iy,iz,1) )*dx_inv*0.5_rk
                        vv_dy = (phi(ix,iy+1,iz,2)*phi(ix,iy+1,iz,2) -phi(ix,iy-1,iz,2)*phi(ix,iy-1,iz,2) )*dy_inv*0.5_rk
                        vw_dz = (phi(ix,iy,iz+1,2)*phi(ix,iy,iz+1,3) -phi(ix,iy,iz-1,2)*phi(ix,iy,iz-1,3) )*dz_inv*0.5_rk

                        wu_dx = (phi(ix+1,iy,iz,3)*phi(ix+1,iy,iz,1) -phi(ix-1,iy,iz,3)*phi(ix-1,iy,iz,1) )*dx_inv*0.5_rk
                        wv_dy = (phi(ix,iy+1,iz,3)*phi(ix,iy+1,iz,2) -phi(ix,iy-1,iz,3)*phi(ix,iy-1,iz,2) )*dy_inv*0.5_rk
                        ww_dz = (phi(ix,iy,iz+1,3)*phi(ix,iy,iz+1,3) -phi(ix,iy,iz-1,3)*phi(ix,iy,iz-1,3) )*dz_inv*0.5_rk

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

                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)

                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))

                        ! Rhs in skew-symmetric formulation (of nonlinear term)
                        ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                        rhs(ix,iy,iz,1) = -0.5_rk*(uu_dx + uv_dy + uw_dz   + u*u_dx + v*u_dy + w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = -0.5_rk*(vu_dx + vv_dy + vw_dz   + u*v_dx + v*v_dy + w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = -0.5_rk*(wu_dx + wv_dy + ww_dz   + u*w_dx + v*w_dy + w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        else
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

                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)

                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))

                        rhs(ix,iy,iz,1) = (-u*u_dx - v*u_dy - w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = (-u*v_dx - v*v_dy - w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = (-u*w_dx - v*w_dy - w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        endif

    case("FD_4th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard scheme)
        !-----------------------------------------------------------------------
        ! Note: a(0) does NOT appear (it is zero...)
        if (params_acm%skew_symmetry) then
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! first derivatives of u, v, p
                        u_dx = (a_FD4(-2)*phi(ix-2,iy,iz,1) +a_FD4(-1)*phi(ix-1,iy,iz,1) +a_FD4(+1)*phi(ix+1,iy,iz,1) +a_FD4(+2)*phi(ix+2,iy,iz,1))*dx_inv
                        v_dx = (a_FD4(-2)*phi(ix-2,iy,iz,2) +a_FD4(-1)*phi(ix-1,iy,iz,2) +a_FD4(+1)*phi(ix+1,iy,iz,2) +a_FD4(+2)*phi(ix+2,iy,iz,2))*dx_inv
                        w_dx = (a_FD4(-2)*phi(ix-2,iy,iz,3) +a_FD4(-1)*phi(ix-1,iy,iz,3) +a_FD4(+1)*phi(ix+1,iy,iz,3) +a_FD4(+2)*phi(ix+2,iy,iz,3))*dx_inv
                        p_dx = (a_FD4(-2)*phi(ix-2,iy,iz,4) +a_FD4(-1)*phi(ix-1,iy,iz,4) +a_FD4(+1)*phi(ix+1,iy,iz,4) +a_FD4(+2)*phi(ix+2,iy,iz,4))*dx_inv

                        u_dy = (a_FD4(-2)*phi(ix,iy-2,iz,1) +a_FD4(-1)*phi(ix,iy-1,iz,1) +a_FD4(+1)*phi(ix,iy+1,iz,1) +a_FD4(+2)*phi(ix,iy+2,iz,1))*dy_inv
                        v_dy = (a_FD4(-2)*phi(ix,iy-2,iz,2) +a_FD4(-1)*phi(ix,iy-1,iz,2) +a_FD4(+1)*phi(ix,iy+1,iz,2) +a_FD4(+2)*phi(ix,iy+2,iz,2))*dy_inv
                        w_dy = (a_FD4(-2)*phi(ix,iy-2,iz,3) +a_FD4(-1)*phi(ix,iy-1,iz,3) +a_FD4(+1)*phi(ix,iy+1,iz,3) +a_FD4(+2)*phi(ix,iy+2,iz,3))*dy_inv
                        p_dy = (a_FD4(-2)*phi(ix,iy-2,iz,4) +a_FD4(-1)*phi(ix,iy-1,iz,4) +a_FD4(+1)*phi(ix,iy+1,iz,4) +a_FD4(+2)*phi(ix,iy+2,iz,4))*dy_inv
                        
                        u_dz = (a_FD4(-2)*phi(ix,iy,iz-2,1) +a_FD4(-1)*phi(ix,iy,iz-1,1) +a_FD4(+1)*phi(ix,iy,iz+1,1) +a_FD4(+2)*phi(ix,iy,iz+2,1))*dz_inv
                        v_dz = (a_FD4(-2)*phi(ix,iy,iz-2,2) +a_FD4(-1)*phi(ix,iy,iz-1,2) +a_FD4(+1)*phi(ix,iy,iz+1,2) +a_FD4(+2)*phi(ix,iy,iz+2,2))*dz_inv
                        w_dz = (a_FD4(-2)*phi(ix,iy,iz-2,3) +a_FD4(-1)*phi(ix,iy,iz-1,3) +a_FD4(+1)*phi(ix,iy,iz+1,3) +a_FD4(+2)*phi(ix,iy,iz+2,3))*dz_inv
                        p_dz = (a_FD4(-2)*phi(ix,iy,iz-2,4) +a_FD4(-1)*phi(ix,iy,iz-1,4) +a_FD4(+1)*phi(ix,iy,iz+1,4) +a_FD4(+2)*phi(ix,iy,iz+2,4))*dz_inv

                        uu_dx = (a_FD4(-2)*phi(ix-2,iy,iz,1)*phi(ix-2,iy,iz,1) +a_FD4(-1)*phi(ix-1,iy,iz,1)*phi(ix-1,iy,iz,1) +a_FD4(+1)*phi(ix+1,iy,iz,1)*phi(ix+1,iy,iz,1) +a_FD4(+2)*phi(ix+2,iy,iz,1)*phi(ix+2,iy,iz,1))*dx_inv
                        uv_dy = (a_FD4(-2)*phi(ix,iy-2,iz,1)*phi(ix,iy-2,iz,2) +a_FD4(-1)*phi(ix,iy-1,iz,1)*phi(ix,iy-1,iz,2) +a_FD4(+1)*phi(ix,iy+1,iz,1)*phi(ix,iy+1,iz,2) +a_FD4(+2)*phi(ix,iy+2,iz,1)*phi(ix,iy+2,iz,2))*dy_inv
                        uw_dz = (a_FD4(-2)*phi(ix,iy,iz-2,1)*phi(ix,iy,iz-2,3) +a_FD4(-1)*phi(ix,iy,iz-1,1)*phi(ix,iy,iz-1,3) +a_FD4(+1)*phi(ix,iy,iz+1,1)*phi(ix,iy,iz+1,3) +a_FD4(+2)*phi(ix,iy,iz+2,1)*phi(ix,iy,iz+2,3))*dz_inv

                        vu_dx = (a_FD4(-2)*phi(ix-2,iy,iz,2)*phi(ix-2,iy,iz,1) +a_FD4(-1)*phi(ix-1,iy,iz,2)*phi(ix-1,iy,iz,1) +a_FD4(+1)*phi(ix+1,iy,iz,2)*phi(ix+1,iy,iz,1) +a_FD4(+2)*phi(ix+2,iy,iz,2)*phi(ix+2,iy,iz,1))*dx_inv
                        vv_dy = (a_FD4(-2)*phi(ix,iy-2,iz,2)*phi(ix,iy-2,iz,2) +a_FD4(-1)*phi(ix,iy-1,iz,2)*phi(ix,iy-1,iz,2) +a_FD4(+1)*phi(ix,iy+1,iz,2)*phi(ix,iy+1,iz,2) +a_FD4(+2)*phi(ix,iy+2,iz,2)*phi(ix,iy+2,iz,2))*dy_inv
                        vw_dz = (a_FD4(-2)*phi(ix,iy,iz-2,2)*phi(ix,iy,iz-2,3) +a_FD4(-1)*phi(ix,iy,iz-1,2)*phi(ix,iy,iz-1,3) +a_FD4(+1)*phi(ix,iy,iz+1,2)*phi(ix,iy,iz+1,3) +a_FD4(+2)*phi(ix,iy,iz+2,2)*phi(ix,iy,iz+2,3))*dz_inv

                        wu_dx = (a_FD4(-2)*phi(ix-2,iy,iz,3)*phi(ix-2,iy,iz,1) +a_FD4(-1)*phi(ix-1,iy,iz,3)*phi(ix-1,iy,iz,1) +a_FD4(+1)*phi(ix+1,iy,iz,3)*phi(ix+1,iy,iz,1) +a_FD4(+2)*phi(ix+2,iy,iz,3)*phi(ix+2,iy,iz,1))*dx_inv
                        wv_dy = (a_FD4(-2)*phi(ix,iy-2,iz,3)*phi(ix,iy-2,iz,2) +a_FD4(-1)*phi(ix,iy-1,iz,3)*phi(ix,iy-1,iz,2) +a_FD4(+1)*phi(ix,iy+1,iz,3)*phi(ix,iy+1,iz,2) +a_FD4(+2)*phi(ix,iy+2,iz,3)*phi(ix,iy+2,iz,2))*dy_inv
                        ww_dz = (a_FD4(-2)*phi(ix,iy,iz-2,3)*phi(ix,iy,iz-2,3) +a_FD4(-1)*phi(ix,iy,iz-1,3)*phi(ix,iy,iz-1,3) +a_FD4(+1)*phi(ix,iy,iz+1,3)*phi(ix,iy,iz+1,3) +a_FD4(+2)*phi(ix,iy,iz+2,3)*phi(ix,iy,iz+2,3))*dz_inv

                        ! second derivatives of u, v and w
                        u_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,1) +b_FD4(-1)*phi(ix-1,iy,iz,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix+1,iy,iz,1) +b_FD4(+2)*phi(ix+2,iy,iz,1))*dx2_inv
                        v_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,2) +b_FD4(-1)*phi(ix-1,iy,iz,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix+1,iy,iz,2) +b_FD4(+2)*phi(ix+2,iy,iz,2))*dx2_inv
                        w_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,3) +b_FD4(-1)*phi(ix-1,iy,iz,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix+1,iy,iz,3) +b_FD4(+2)*phi(ix+2,iy,iz,3))*dx2_inv
                        
                        u_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,1) +b_FD4(-1)*phi(ix,iy-1,iz,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix,iy+1,iz,1) +b_FD4(+2)*phi(ix,iy+2,iz,1))*dy2_inv
                        v_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,2) +b_FD4(-1)*phi(ix,iy-1,iz,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix,iy+1,iz,2) +b_FD4(+2)*phi(ix,iy+2,iz,2))*dy2_inv
                        w_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,3) +b_FD4(-1)*phi(ix,iy-1,iz,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix,iy+1,iz,3) +b_FD4(+2)*phi(ix,iy+2,iz,3))*dy2_inv

                        u_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,1) +b_FD4(-1)*phi(ix,iy,iz-1,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix,iy,iz+1,1) +b_FD4(+2)*phi(ix,iy,iz+2,1))*dz2_inv
                        v_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,2) +b_FD4(-1)*phi(ix,iy,iz-1,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix,iy,iz+1,2) +b_FD4(+2)*phi(ix,iy,iz+2,2))*dz2_inv
                        w_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,3) +b_FD4(-1)*phi(ix,iy,iz-1,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix,iy,iz+1,3) +b_FD4(+2)*phi(ix,iy,iz+2,3))*dz2_inv

                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)

                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))

                        ! Rhs in skew-symmetric formulation (of nonlinear term)
                        ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                        rhs(ix,iy,iz,1) = -0.5_rk*(uu_dx + uv_dy + uw_dz   + u*u_dx + v*u_dy + w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = -0.5_rk*(vu_dx + vv_dy + vw_dz   + u*v_dx + v*v_dy + w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = -0.5_rk*(wu_dx + wv_dy + ww_dz   + u*w_dx + v*w_dy + w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        else
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! first derivatives of u, v, p
                        u_dx = (a_FD4(-2)*phi(ix-2,iy,iz,1) +a_FD4(-1)*phi(ix-1,iy,iz,1) +a_FD4(+1)*phi(ix+1,iy,iz,1) +a_FD4(+2)*phi(ix+2,iy,iz,1))*dx_inv    
                        v_dx = (a_FD4(-2)*phi(ix-2,iy,iz,2) +a_FD4(-1)*phi(ix-1,iy,iz,2) +a_FD4(+1)*phi(ix+1,iy,iz,2) +a_FD4(+2)*phi(ix+2,iy,iz,2))*dx_inv    
                        w_dx = (a_FD4(-2)*phi(ix-2,iy,iz,3) +a_FD4(-1)*phi(ix-1,iy,iz,3) +a_FD4(+1)*phi(ix+1,iy,iz,3) +a_FD4(+2)*phi(ix+2,iy,iz,3))*dx_inv    
                        p_dx = (a_FD4(-2)*phi(ix-2,iy,iz,4) +a_FD4(-1)*phi(ix-1,iy,iz,4) +a_FD4(+1)*phi(ix+1,iy,iz,4) +a_FD4(+2)*phi(ix+2,iy,iz,4))*dx_inv    

                        u_dy = (a_FD4(-2)*phi(ix,iy-2,iz,1) +a_FD4(-1)*phi(ix,iy-1,iz,1) +a_FD4(+1)*phi(ix,iy+1,iz,1) +a_FD4(+2)*phi(ix,iy+2,iz,1))*dy_inv    
                        v_dy = (a_FD4(-2)*phi(ix,iy-2,iz,2) +a_FD4(-1)*phi(ix,iy-1,iz,2) +a_FD4(+1)*phi(ix,iy+1,iz,2) +a_FD4(+2)*phi(ix,iy+2,iz,2))*dy_inv    
                        w_dy = (a_FD4(-2)*phi(ix,iy-2,iz,3) +a_FD4(-1)*phi(ix,iy-1,iz,3) +a_FD4(+1)*phi(ix,iy+1,iz,3) +a_FD4(+2)*phi(ix,iy+2,iz,3))*dy_inv    
                        p_dy = (a_FD4(-2)*phi(ix,iy-2,iz,4) +a_FD4(-1)*phi(ix,iy-1,iz,4) +a_FD4(+1)*phi(ix,iy+1,iz,4) +a_FD4(+2)*phi(ix,iy+2,iz,4))*dy_inv    

                        u_dz = (a_FD4(-2)*phi(ix,iy,iz-2,1) +a_FD4(-1)*phi(ix,iy,iz-1,1) +a_FD4(+1)*phi(ix,iy,iz+1,1) +a_FD4(+2)*phi(ix,iy,iz+2,1))*dz_inv    
                        v_dz = (a_FD4(-2)*phi(ix,iy,iz-2,2) +a_FD4(-1)*phi(ix,iy,iz-1,2) +a_FD4(+1)*phi(ix,iy,iz+1,2) +a_FD4(+2)*phi(ix,iy,iz+2,2))*dz_inv    
                        w_dz = (a_FD4(-2)*phi(ix,iy,iz-2,3) +a_FD4(-1)*phi(ix,iy,iz-1,3) +a_FD4(+1)*phi(ix,iy,iz+1,3) +a_FD4(+2)*phi(ix,iy,iz+2,3))*dz_inv    
                        p_dz = (a_FD4(-2)*phi(ix,iy,iz-2,4) +a_FD4(-1)*phi(ix,iy,iz-1,4) +a_FD4(+1)*phi(ix,iy,iz+1,4) +a_FD4(+2)*phi(ix,iy,iz+2,4))*dz_inv    

                        ! second derivatives of u, v and w
                        u_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,1)  + b_FD4(-1)*phi(ix-1,iy,iz,1)  + b_FD4(0)*phi(ix,iy,iz,1)  + b_FD4(+1)*phi(ix+1,iy,iz,1)  + b_FD4(+2)*phi(ix+2,iy,iz,1))*dx2_inv    
                        v_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,2)  + b_FD4(-1)*phi(ix-1,iy,iz,2)  + b_FD4(0)*phi(ix,iy,iz,2)  + b_FD4(+1)*phi(ix+1,iy,iz,2)  + b_FD4(+2)*phi(ix+2,iy,iz,2))*dx2_inv    
                        w_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,3)  + b_FD4(-1)*phi(ix-1,iy,iz,3)  + b_FD4(0)*phi(ix,iy,iz,3)  + b_FD4(+1)*phi(ix+1,iy,iz,3)  + b_FD4(+2)*phi(ix+2,iy,iz,3))*dx2_inv    

                        u_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,1)  + b_FD4(-1)*phi(ix,iy-1,iz,1)  + b_FD4(0)*phi(ix,iy,iz,1)  + b_FD4(+1)*phi(ix,iy+1,iz,1)  + b_FD4(+2)*phi(ix,iy+2,iz,1))*dy2_inv    
                        v_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,2)  + b_FD4(-1)*phi(ix,iy-1,iz,2)  + b_FD4(0)*phi(ix,iy,iz,2)  + b_FD4(+1)*phi(ix,iy+1,iz,2)  + b_FD4(+2)*phi(ix,iy+2,iz,2))*dy2_inv    
                        w_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,3)  + b_FD4(-1)*phi(ix,iy-1,iz,3)  + b_FD4(0)*phi(ix,iy,iz,3)  + b_FD4(+1)*phi(ix,iy+1,iz,3)  + b_FD4(+2)*phi(ix,iy+2,iz,3))*dy2_inv    

                        u_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,1)  + b_FD4(-1)*phi(ix,iy,iz-1,1)  + b_FD4(0)*phi(ix,iy,iz,1)  + b_FD4(+1)*phi(ix,iy,iz+1,1)  + b_FD4(+2)*phi(ix,iy,iz+2,1))*dz2_inv    
                        v_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,2)  + b_FD4(-1)*phi(ix,iy,iz-1,2)  + b_FD4(0)*phi(ix,iy,iz,2)  + b_FD4(+1)*phi(ix,iy,iz+1,2)  + b_FD4(+2)*phi(ix,iy,iz+2,2))*dz2_inv    
                        w_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,3)  + b_FD4(-1)*phi(ix,iy,iz-1,3)  + b_FD4(0)*phi(ix,iy,iz,3)  + b_FD4(+1)*phi(ix,iy,iz+1,3)  + b_FD4(+2)*phi(ix,iy,iz+2,3))*dz2_inv    
    
                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)
    
                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))
    
                        rhs(ix,iy,iz,1) = (-u*u_dx - v*u_dy - w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = (-u*v_dx - v*v_dy - w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = (-u*w_dx - v*w_dy - w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        endif

        select case(params_acm%p_eqn_model)
        case ('acm')
            ! do nothing, is the eqn computed above without additional terms
        case ('diffusive')
            ! add p diffusion []
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g

                        p_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,4) +b_FD4(-1)*phi(ix-1,iy,iz,4) +b_FD4(0)*phi(ix,iy,iz,4) +b_FD4(+1)*phi(ix+1,iy,iz,4) +b_FD4(+2)*phi(ix+2,iy,iz,4))*dx2_inv
                        p_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,4) +b_FD4(-1)*phi(ix,iy-1,iz,4) +b_FD4(0)*phi(ix,iy,iz,4) +b_FD4(+1)*phi(ix,iy+1,iz,4) +b_FD4(+2)*phi(ix,iy+2,iz,4))*dy2_inv
                        p_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,4) +b_FD4(-1)*phi(ix,iy,iz-1,4) +b_FD4(0)*phi(ix,iy,iz,4) +b_FD4(+1)*phi(ix,iy,iz+1,4) +b_FD4(+2)*phi(ix,iy,iz+2,4))*dz2_inv

                        rhs(ix,iy,iz,4) = rhs(ix,iy,iz,4) + params_acm%nu_p*(p_dxdx + p_dydy + p_dzdz)
                    enddo
                enddo   
            enddo   
        case ('EDAC')
            ! add EDAC term AND diffusion [Clausen, Entropically damped form of artificial compressibility for explicit simulation of incompressible flow, 2013]
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        p_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,4) +b_FD4(-1)*phi(ix-1,iy,iz,4) +b_FD4(0)*phi(ix,iy,iz,4) +b_FD4(+1)*phi(ix+1,iy,iz,4) +b_FD4(+2)*phi(ix+2,iy,iz,4))*dx2_inv
                        p_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,4) +b_FD4(-1)*phi(ix,iy-1,iz,4) +b_FD4(0)*phi(ix,iy,iz,4) +b_FD4(+1)*phi(ix,iy+1,iz,4) +b_FD4(+2)*phi(ix,iy+2,iz,4))*dy2_inv
                        p_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,4) +b_FD4(-1)*phi(ix,iy,iz-1,4) +b_FD4(0)*phi(ix,iy,iz,4) +b_FD4(+1)*phi(ix,iy,iz+1,4) +b_FD4(+2)*phi(ix,iy,iz+2,4))*dz2_inv

                        pu_dx = (a_FD4(-2)*phi(ix-2,iy,iz,4)*phi(ix-2,iy,iz,1) +a_FD4(-1)*phi(ix-1,iy,iz,4)*phi(ix-1,iy,iz,1) +a_FD4(+1)*phi(ix+1,iy,iz,4)*phi(ix+1,iy,iz,1) +a_FD4(+2)*phi(ix+2,iy,iz,4)*phi(ix+2,iy,iz,1))*dx_inv
                        pv_dy = (a_FD4(-2)*phi(ix,iy-2,iz,4)*phi(ix,iy-2,iz,2) +a_FD4(-1)*phi(ix,iy-1,iz,4)*phi(ix,iy-1,iz,2) +a_FD4(+1)*phi(ix,iy+1,iz,4)*phi(ix,iy+1,iz,2) +a_FD4(+2)*phi(ix,iy+2,iz,4)*phi(ix,iy+2,iz,2))*dy_inv
                        pw_dz = (a_FD4(-2)*phi(ix,iy,iz-2,4)*phi(ix,iy,iz-2,3) +a_FD4(-1)*phi(ix,iy,iz-1,4)*phi(ix,iy,iz-1,3) +a_FD4(+1)*phi(ix,iy,iz+1,4)*phi(ix,iy,iz+1,3) +a_FD4(+2)*phi(ix,iy,iz+2,4)*phi(ix,iy,iz+2,3))*dz_inv

                        rhs(ix,iy,iz,4) = rhs(ix,iy,iz,4) + params_acm%nu_p*(p_dxdx + p_dydy + p_dzdz) - (pu_dx+pv_dy+pw_dz)
                    enddo
                enddo   
            enddo   
        case ('convective')
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        pu_dx = (a_FD4(-2)*phi(ix-2,iy,iz,4)*phi(ix-2,iy,iz,1) +a_FD4(-1)*phi(ix-1,iy,iz,4)*phi(ix-1,iy,iz,1) +a_FD4(+1)*phi(ix+1,iy,iz,4)*phi(ix+1,iy,iz,1) +a_FD4(+2)*phi(ix+2,iy,iz,4)*phi(ix+2,iy,iz,1))*dx_inv
                        pv_dy = (a_FD4(-2)*phi(ix,iy-2,iz,4)*phi(ix,iy-2,iz,2) +a_FD4(-1)*phi(ix,iy-1,iz,4)*phi(ix,iy-1,iz,2) +a_FD4(+1)*phi(ix,iy+1,iz,4)*phi(ix,iy+1,iz,2) +a_FD4(+2)*phi(ix,iy+2,iz,4)*phi(ix,iy+2,iz,2))*dy_inv
                        pw_dz = (a_FD4(-2)*phi(ix,iy,iz-2,4)*phi(ix,iy,iz-2,3) +a_FD4(-1)*phi(ix,iy,iz-1,4)*phi(ix,iy,iz-1,3) +a_FD4(+1)*phi(ix,iy,iz+1,4)*phi(ix,iy,iz+1,3) +a_FD4(+2)*phi(ix,iy,iz+2,4)*phi(ix,iy,iz+2,3))*dz_inv

                        rhs(ix,iy,iz,4) = rhs(ix,iy,iz,4) - (pu_dx+pv_dy+pw_dz)
                    enddo
                enddo    
            enddo    
        case default
            call abort(2501041, "pressure equation model is unkown: "//trim(adjustl(params_acm%p_eqn_model)))
        end select   

    case("FD_6th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard scheme)
        !-----------------------------------------------------------------------
        ! Note: a(0) does NOT appear (it is zero...)
        if (params_acm%skew_symmetry) then
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! first derivatives of u, v, p
                        u_dx = (a_FD6(-3)*phi(ix-3,iy,iz,1) +a_FD6(-2)*phi(ix-2,iy,iz,1) +a_FD6(-1)*phi(ix-1,iy,iz,1) +a_FD6(+1)*phi(ix+1,iy,iz,1) +a_FD6(+2)*phi(ix+2,iy,iz,1) +a_FD6(+3)*phi(ix+3,iy,iz,1))*dx_inv
                        v_dx = (a_FD6(-3)*phi(ix-3,iy,iz,2) +a_FD6(-2)*phi(ix-2,iy,iz,2) +a_FD6(-1)*phi(ix-1,iy,iz,2) +a_FD6(+1)*phi(ix+1,iy,iz,2) +a_FD6(+2)*phi(ix+2,iy,iz,2) +a_FD6(+3)*phi(ix+3,iy,iz,2))*dx_inv
                        w_dx = (a_FD6(-3)*phi(ix-3,iy,iz,3) +a_FD6(-2)*phi(ix-2,iy,iz,3) +a_FD6(-1)*phi(ix-1,iy,iz,3) +a_FD6(+1)*phi(ix+1,iy,iz,3) +a_FD6(+2)*phi(ix+2,iy,iz,3) +a_FD6(+3)*phi(ix+3,iy,iz,3))*dx_inv
                        p_dx = (a_FD6(-3)*phi(ix-3,iy,iz,4) +a_FD6(-2)*phi(ix-2,iy,iz,4) +a_FD6(-1)*phi(ix-1,iy,iz,4) +a_FD6(+1)*phi(ix+1,iy,iz,4) +a_FD6(+2)*phi(ix+2,iy,iz,4) +a_FD6(+3)*phi(ix+3,iy,iz,4))*dx_inv

                        u_dy = (a_FD6(-3)*phi(ix,iy-3,iz,1) +a_FD6(-2)*phi(ix,iy-2,iz,1) +a_FD6(-1)*phi(ix,iy-1,iz,1) +a_FD6(+1)*phi(ix,iy+1,iz,1) +a_FD6(+2)*phi(ix,iy+2,iz,1) +a_FD6(+3)*phi(ix,iy+3,iz,1))*dy_inv
                        v_dy = (a_FD6(-3)*phi(ix,iy-3,iz,2) +a_FD6(-2)*phi(ix,iy-2,iz,2) +a_FD6(-1)*phi(ix,iy-1,iz,2) +a_FD6(+1)*phi(ix,iy+1,iz,2) +a_FD6(+2)*phi(ix,iy+2,iz,2) +a_FD6(+3)*phi(ix,iy+3,iz,2))*dy_inv
                        w_dy = (a_FD6(-3)*phi(ix,iy-3,iz,3) +a_FD6(-2)*phi(ix,iy-2,iz,3) +a_FD6(-1)*phi(ix,iy-1,iz,3) +a_FD6(+1)*phi(ix,iy+1,iz,3) +a_FD6(+2)*phi(ix,iy+2,iz,3) +a_FD6(+3)*phi(ix,iy+3,iz,3))*dy_inv
                        p_dy = (a_FD6(-3)*phi(ix,iy-3,iz,4) +a_FD6(-2)*phi(ix,iy-2,iz,4) +a_FD6(-1)*phi(ix,iy-1,iz,4) +a_FD6(+1)*phi(ix,iy+1,iz,4) +a_FD6(+2)*phi(ix,iy+2,iz,4) +a_FD6(+3)*phi(ix,iy+3,iz,4))*dy_inv
                        
                        u_dz = (a_FD6(-3)*phi(ix,iy,iz-3,1) +a_FD6(-2)*phi(ix,iy,iz-2,1) +a_FD6(-1)*phi(ix,iy,iz-1,1) +a_FD6(+1)*phi(ix,iy,iz+1,1) +a_FD6(+2)*phi(ix,iy,iz+2,1) +a_FD6(+3)*phi(ix,iy,iz+3,1))*dz_inv
                        v_dz = (a_FD6(-3)*phi(ix,iy,iz-3,2) +a_FD6(-2)*phi(ix,iy,iz-2,2) +a_FD6(-1)*phi(ix,iy,iz-1,2) +a_FD6(+1)*phi(ix,iy,iz+1,2) +a_FD6(+2)*phi(ix,iy,iz+2,2) +a_FD6(+3)*phi(ix,iy,iz+3,2))*dz_inv
                        w_dz = (a_FD6(-3)*phi(ix,iy,iz-3,3) +a_FD6(-2)*phi(ix,iy,iz-2,3) +a_FD6(-1)*phi(ix,iy,iz-1,3) +a_FD6(+1)*phi(ix,iy,iz+1,3) +a_FD6(+2)*phi(ix,iy,iz+2,3) +a_FD6(+3)*phi(ix,iy,iz+3,3))*dz_inv
                        p_dz = (a_FD6(-3)*phi(ix,iy,iz-3,4) +a_FD6(-2)*phi(ix,iy,iz-2,4) +a_FD6(-1)*phi(ix,iy,iz-1,4) +a_FD6(+1)*phi(ix,iy,iz+1,4) +a_FD6(+2)*phi(ix,iy,iz+2,4) +a_FD6(+3)*phi(ix,iy,iz+3,4))*dz_inv

                        uu_dx = (a_FD6(-3)*phi(ix-3,iy,iz,1)*phi(ix-3,iy,iz,1) +a_FD6(-2)*phi(ix-2,iy,iz,1)*phi(ix-2,iy,iz,1) +a_FD6(-1)*phi(ix-1,iy,iz,1)*phi(ix-1,iy,iz,1) +a_FD6(+1)*phi(ix+1,iy,iz,1)*phi(ix+1,iy,iz,1) +a_FD6(+2)*phi(ix+2,iy,iz,1)*phi(ix+2,iy,iz,1) +a_FD6(+3)*phi(ix+3,iy,iz,1)*phi(ix+3,iy,iz,1))*dx_inv
                        uv_dy = (a_FD6(-3)*phi(ix,iy-3,iz,1)*phi(ix,iy-3,iz,2) +a_FD6(-2)*phi(ix,iy-2,iz,1)*phi(ix,iy-2,iz,2) +a_FD6(-1)*phi(ix,iy-1,iz,1)*phi(ix,iy-1,iz,2) +a_FD6(+1)*phi(ix,iy+1,iz,1)*phi(ix,iy+1,iz,2) +a_FD6(+2)*phi(ix,iy+2,iz,1)*phi(ix,iy+2,iz,2) +a_FD6(+3)*phi(ix,iy+3,iz,1)*phi(ix,iy+3,iz,2))*dy_inv
                        uw_dz = (a_FD6(-3)*phi(ix,iy,iz-3,1)*phi(ix,iy,iz-3,3) +a_FD6(-2)*phi(ix,iy,iz-2,1)*phi(ix,iy,iz-2,3) +a_FD6(-1)*phi(ix,iy,iz-1,1)*phi(ix,iy,iz-1,3) +a_FD6(+1)*phi(ix,iy,iz+1,1)*phi(ix,iy,iz+1,3) +a_FD6(+2)*phi(ix,iy,iz+2,1)*phi(ix,iy,iz+2,3) +a_FD6(+3)*phi(ix,iy,iz+3,1)*phi(ix,iy,iz+3,3))*dz_inv

                        vu_dx = (a_FD6(-3)*phi(ix-3,iy,iz,2)*phi(ix-3,iy,iz,1) +a_FD6(-2)*phi(ix-2,iy,iz,2)*phi(ix-2,iy,iz,1) +a_FD6(-1)*phi(ix-1,iy,iz,2)*phi(ix-1,iy,iz,1) +a_FD6(+1)*phi(ix+1,iy,iz,2)*phi(ix+1,iy,iz,1) +a_FD6(+2)*phi(ix+2,iy,iz,2)*phi(ix+2,iy,iz,1) +a_FD6(+3)*phi(ix+3,iy,iz,2)*phi(ix+3,iy,iz,1))*dx_inv
                        vv_dy = (a_FD6(-3)*phi(ix,iy-3,iz,2)*phi(ix,iy-3,iz,2) +a_FD6(-2)*phi(ix,iy-2,iz,2)*phi(ix,iy-2,iz,2) +a_FD6(-1)*phi(ix,iy-1,iz,2)*phi(ix,iy-1,iz,2) +a_FD6(+1)*phi(ix,iy+1,iz,2)*phi(ix,iy+1,iz,2) +a_FD6(+2)*phi(ix,iy+2,iz,2)*phi(ix,iy+2,iz,2) +a_FD6(+3)*phi(ix,iy+3,iz,2)*phi(ix,iy+3,iz,2))*dy_inv
                        vw_dz = (a_FD6(-3)*phi(ix,iy,iz-3,2)*phi(ix,iy,iz-3,3) +a_FD6(-2)*phi(ix,iy,iz-2,2)*phi(ix,iy,iz-2,3) +a_FD6(-1)*phi(ix,iy,iz-1,2)*phi(ix,iy,iz-1,3) +a_FD6(+1)*phi(ix,iy,iz+1,2)*phi(ix,iy,iz+1,3) +a_FD6(+2)*phi(ix,iy,iz+2,2)*phi(ix,iy,iz+2,3) +a_FD6(+3)*phi(ix,iy,iz+3,2)*phi(ix,iy,iz+3,3))*dz_inv

                        wu_dx = (a_FD6(-3)*phi(ix-3,iy,iz,3)*phi(ix-3,iy,iz,1) +a_FD6(-2)*phi(ix-2,iy,iz,3)*phi(ix-2,iy,iz,1) +a_FD6(-1)*phi(ix-1,iy,iz,3)*phi(ix-1,iy,iz,1) +a_FD6(+1)*phi(ix+1,iy,iz,3)*phi(ix+1,iy,iz,1) +a_FD6(+2)*phi(ix+2,iy,iz,3)*phi(ix+2,iy,iz,1) +a_FD6(+3)*phi(ix+3,iy,iz,3)*phi(ix+3,iy,iz,1))*dx_inv
                        wv_dy = (a_FD6(-3)*phi(ix,iy-3,iz,3)*phi(ix,iy-3,iz,2) +a_FD6(-2)*phi(ix,iy-2,iz,3)*phi(ix,iy-2,iz,2) +a_FD6(-1)*phi(ix,iy-1,iz,3)*phi(ix,iy-1,iz,2) +a_FD6(+1)*phi(ix,iy+1,iz,3)*phi(ix,iy+1,iz,2) +a_FD6(+2)*phi(ix,iy+2,iz,3)*phi(ix,iy+2,iz,2) +a_FD6(+3)*phi(ix,iy+3,iz,3)*phi(ix,iy+3,iz,2))*dy_inv
                        ww_dz = (a_FD6(-3)*phi(ix,iy,iz-3,3)*phi(ix,iy,iz-3,3) +a_FD6(-2)*phi(ix,iy,iz-2,3)*phi(ix,iy,iz-2,3) +a_FD6(-1)*phi(ix,iy,iz-1,3)*phi(ix,iy,iz-1,3) +a_FD6(+1)*phi(ix,iy,iz+1,3)*phi(ix,iy,iz+1,3) +a_FD6(+2)*phi(ix,iy,iz+2,3)*phi(ix,iy,iz+2,3) +a_FD6(+3)*phi(ix,iy,iz+3,3)*phi(ix,iy,iz+3,3))*dz_inv
                        
                        ! second derivatives of u, v and w
                        u_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,1) +b_FD6(-2)*phi(ix-2,iy,iz,1) +b_FD6(-1)*phi(ix-1,iy,iz,1) +b_FD6(0)*phi(ix,iy,iz,1) +b_FD6(+1)*phi(ix+1,iy,iz,1) +b_FD6(+2)*phi(ix+2,iy,iz,1) +b_FD6(+3)*phi(ix+3,iy,iz,1))*dx2_inv
                        u_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,1) +b_FD6(-2)*phi(ix,iy-2,iz,1) +b_FD6(-1)*phi(ix,iy-1,iz,1) +b_FD6(0)*phi(ix,iy,iz,1) +b_FD6(+1)*phi(ix,iy+1,iz,1) +b_FD6(+2)*phi(ix,iy+2,iz,1) +b_FD6(+3)*phi(ix,iy+3,iz,1))*dy2_inv
                        u_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,1) +b_FD6(-2)*phi(ix,iy,iz-2,1) +b_FD6(-1)*phi(ix,iy,iz-1,1) +b_FD6(0)*phi(ix,iy,iz,1) +b_FD6(+1)*phi(ix,iy,iz+1,1) +b_FD6(+2)*phi(ix,iy,iz+2,1) +b_FD6(+3)*phi(ix,iy,iz+3,1))*dz2_inv
                        v_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,2) +b_FD6(-2)*phi(ix-2,iy,iz,2) +b_FD6(-1)*phi(ix-1,iy,iz,2) +b_FD6(0)*phi(ix,iy,iz,2) +b_FD6(+1)*phi(ix+1,iy,iz,2) +b_FD6(+2)*phi(ix+2,iy,iz,2) +b_FD6(+3)*phi(ix+3,iy,iz,2))*dx2_inv
                        v_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,2) +b_FD6(-2)*phi(ix,iy-2,iz,2) +b_FD6(-1)*phi(ix,iy-1,iz,2) +b_FD6(0)*phi(ix,iy,iz,2) +b_FD6(+1)*phi(ix,iy+1,iz,2) +b_FD6(+2)*phi(ix,iy+2,iz,2) +b_FD6(+3)*phi(ix,iy+3,iz,2))*dy2_inv
                        v_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,2) +b_FD6(-2)*phi(ix,iy,iz-2,2) +b_FD6(-1)*phi(ix,iy,iz-1,2) +b_FD6(0)*phi(ix,iy,iz,2) +b_FD6(+1)*phi(ix,iy,iz+1,2) +b_FD6(+2)*phi(ix,iy,iz+2,2) +b_FD6(+3)*phi(ix,iy,iz+3,2))*dz2_inv
                        w_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,3) +b_FD6(-2)*phi(ix-2,iy,iz,3) +b_FD6(-1)*phi(ix-1,iy,iz,3) +b_FD6(0)*phi(ix,iy,iz,3) +b_FD6(+1)*phi(ix+1,iy,iz,3) +b_FD6(+2)*phi(ix+2,iy,iz,3) +b_FD6(+3)*phi(ix+3,iy,iz,3))*dx2_inv
                        w_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,3) +b_FD6(-2)*phi(ix,iy-2,iz,3) +b_FD6(-1)*phi(ix,iy-1,iz,3) +b_FD6(0)*phi(ix,iy,iz,3) +b_FD6(+1)*phi(ix,iy+1,iz,3) +b_FD6(+2)*phi(ix,iy+2,iz,3) +b_FD6(+3)*phi(ix,iy+3,iz,3))*dy2_inv
                        w_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,3) +b_FD6(-2)*phi(ix,iy,iz-2,3) +b_FD6(-1)*phi(ix,iy,iz-1,3) +b_FD6(0)*phi(ix,iy,iz,3) +b_FD6(+1)*phi(ix,iy,iz+1,3) +b_FD6(+2)*phi(ix,iy,iz+2,3) +b_FD6(+3)*phi(ix,iy,iz+3,3))*dz2_inv
                        
                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)

                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))

                        ! Rhs in skew-symmetric formulation (of nonlinear term)
                        ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                        rhs(ix,iy,iz,1) = -0.5_rk*(uu_dx + uv_dy + uw_dz   + u*u_dx + v*u_dy + w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = -0.5_rk*(vu_dx + vv_dy + vw_dz   + u*v_dx + v*v_dy + w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = -0.5_rk*(wu_dx + wv_dy + ww_dz   + u*w_dx + v*w_dy + w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        else
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! first derivatives of u, v, p
                        u_dx = (a_FD6(-3)*phi(ix-3,iy,iz,1) +a_FD6(-2)*phi(ix-2,iy,iz,1) +a_FD6(-1)*phi(ix-1,iy,iz,1) +a_FD6(+1)*phi(ix+1,iy,iz,1) +a_FD6(+2)*phi(ix+2,iy,iz,1) +a_FD6(+3)*phi(ix+3,iy,iz,1))*dx_inv
                        v_dx = (a_FD6(-3)*phi(ix-3,iy,iz,2) +a_FD6(-2)*phi(ix-2,iy,iz,2) +a_FD6(-1)*phi(ix-1,iy,iz,2) +a_FD6(+1)*phi(ix+1,iy,iz,2) +a_FD6(+2)*phi(ix+2,iy,iz,2) +a_FD6(+3)*phi(ix+3,iy,iz,2))*dx_inv
                        w_dx = (a_FD6(-3)*phi(ix-3,iy,iz,3) +a_FD6(-2)*phi(ix-2,iy,iz,3) +a_FD6(-1)*phi(ix-1,iy,iz,3) +a_FD6(+1)*phi(ix+1,iy,iz,3) +a_FD6(+2)*phi(ix+2,iy,iz,3) +a_FD6(+3)*phi(ix+3,iy,iz,3))*dx_inv
                        p_dx = (a_FD6(-3)*phi(ix-3,iy,iz,4) +a_FD6(-2)*phi(ix-2,iy,iz,4) +a_FD6(-1)*phi(ix-1,iy,iz,4) +a_FD6(+1)*phi(ix+1,iy,iz,4) +a_FD6(+2)*phi(ix+2,iy,iz,4) +a_FD6(+3)*phi(ix+3,iy,iz,4))*dx_inv

                        u_dy = (a_FD6(-3)*phi(ix,iy-3,iz,1) +a_FD6(-2)*phi(ix,iy-2,iz,1) +a_FD6(-1)*phi(ix,iy-1,iz,1) +a_FD6(+1)*phi(ix,iy+1,iz,1) +a_FD6(+2)*phi(ix,iy+2,iz,1) +a_FD6(+3)*phi(ix,iy+3,iz,1))*dy_inv
                        v_dy = (a_FD6(-3)*phi(ix,iy-3,iz,2) +a_FD6(-2)*phi(ix,iy-2,iz,2) +a_FD6(-1)*phi(ix,iy-1,iz,2) +a_FD6(+1)*phi(ix,iy+1,iz,2) +a_FD6(+2)*phi(ix,iy+2,iz,2) +a_FD6(+3)*phi(ix,iy+3,iz,2))*dy_inv
                        w_dy = (a_FD6(-3)*phi(ix,iy-3,iz,3) +a_FD6(-2)*phi(ix,iy-2,iz,3) +a_FD6(-1)*phi(ix,iy-1,iz,3) +a_FD6(+1)*phi(ix,iy+1,iz,3) +a_FD6(+2)*phi(ix,iy+2,iz,3) +a_FD6(+3)*phi(ix,iy+3,iz,3))*dy_inv
                        p_dy = (a_FD6(-3)*phi(ix,iy-3,iz,4) +a_FD6(-2)*phi(ix,iy-2,iz,4) +a_FD6(-1)*phi(ix,iy-1,iz,4) +a_FD6(+1)*phi(ix,iy+1,iz,4) +a_FD6(+2)*phi(ix,iy+2,iz,4) +a_FD6(+3)*phi(ix,iy+3,iz,4))*dy_inv
                        
                        u_dz = (a_FD6(-3)*phi(ix,iy,iz-3,1) +a_FD6(-2)*phi(ix,iy,iz-2,1) +a_FD6(-1)*phi(ix,iy,iz-1,1) +a_FD6(+1)*phi(ix,iy,iz+1,1) +a_FD6(+2)*phi(ix,iy,iz+2,1) +a_FD6(+3)*phi(ix,iy,iz+3,1))*dz_inv
                        v_dz = (a_FD6(-3)*phi(ix,iy,iz-3,2) +a_FD6(-2)*phi(ix,iy,iz-2,2) +a_FD6(-1)*phi(ix,iy,iz-1,2) +a_FD6(+1)*phi(ix,iy,iz+1,2) +a_FD6(+2)*phi(ix,iy,iz+2,2) +a_FD6(+3)*phi(ix,iy,iz+3,2))*dz_inv
                        w_dz = (a_FD6(-3)*phi(ix,iy,iz-3,3) +a_FD6(-2)*phi(ix,iy,iz-2,3) +a_FD6(-1)*phi(ix,iy,iz-1,3) +a_FD6(+1)*phi(ix,iy,iz+1,3) +a_FD6(+2)*phi(ix,iy,iz+2,3) +a_FD6(+3)*phi(ix,iy,iz+3,3))*dz_inv
                        p_dz = (a_FD6(-3)*phi(ix,iy,iz-3,4) +a_FD6(-2)*phi(ix,iy,iz-2,4) +a_FD6(-1)*phi(ix,iy,iz-1,4) +a_FD6(+1)*phi(ix,iy,iz+1,4) +a_FD6(+2)*phi(ix,iy,iz+2,4) +a_FD6(+3)*phi(ix,iy,iz+3,4))*dz_inv
                        
                        ! second derivatives of u, v and w
                        u_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,1) +b_FD6(-2)*phi(ix-2,iy,iz,1) +b_FD6(-1)*phi(ix-1,iy,iz,1) +b_FD6(0)*phi(ix,iy,iz,1) +b_FD6(+1)*phi(ix+1,iy,iz,1) +b_FD6(+2)*phi(ix+2,iy,iz,1) +b_FD6(+3)*phi(ix+3,iy,iz,1))*dx2_inv
                        u_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,1) +b_FD6(-2)*phi(ix,iy-2,iz,1) +b_FD6(-1)*phi(ix,iy-1,iz,1) +b_FD6(0)*phi(ix,iy,iz,1) +b_FD6(+1)*phi(ix,iy+1,iz,1) +b_FD6(+2)*phi(ix,iy+2,iz,1) +b_FD6(+3)*phi(ix,iy+3,iz,1))*dy2_inv
                        u_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,1) +b_FD6(-2)*phi(ix,iy,iz-2,1) +b_FD6(-1)*phi(ix,iy,iz-1,1) +b_FD6(0)*phi(ix,iy,iz,1) +b_FD6(+1)*phi(ix,iy,iz+1,1) +b_FD6(+2)*phi(ix,iy,iz+2,1) +b_FD6(+3)*phi(ix,iy,iz+3,1))*dz2_inv
                        v_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,2) +b_FD6(-2)*phi(ix-2,iy,iz,2) +b_FD6(-1)*phi(ix-1,iy,iz,2) +b_FD6(0)*phi(ix,iy,iz,2) +b_FD6(+1)*phi(ix+1,iy,iz,2) +b_FD6(+2)*phi(ix+2,iy,iz,2) +b_FD6(+3)*phi(ix+3,iy,iz,2))*dx2_inv
                        v_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,2) +b_FD6(-2)*phi(ix,iy-2,iz,2) +b_FD6(-1)*phi(ix,iy-1,iz,2) +b_FD6(0)*phi(ix,iy,iz,2) +b_FD6(+1)*phi(ix,iy+1,iz,2) +b_FD6(+2)*phi(ix,iy+2,iz,2) +b_FD6(+3)*phi(ix,iy+3,iz,2))*dy2_inv
                        v_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,2) +b_FD6(-2)*phi(ix,iy,iz-2,2) +b_FD6(-1)*phi(ix,iy,iz-1,2) +b_FD6(0)*phi(ix,iy,iz,2) +b_FD6(+1)*phi(ix,iy,iz+1,2) +b_FD6(+2)*phi(ix,iy,iz+2,2) +b_FD6(+3)*phi(ix,iy,iz+3,2))*dz2_inv
                        w_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,3) +b_FD6(-2)*phi(ix-2,iy,iz,3) +b_FD6(-1)*phi(ix-1,iy,iz,3) +b_FD6(0)*phi(ix,iy,iz,3) +b_FD6(+1)*phi(ix+1,iy,iz,3) +b_FD6(+2)*phi(ix+2,iy,iz,3) +b_FD6(+3)*phi(ix+3,iy,iz,3))*dx2_inv
                        w_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,3) +b_FD6(-2)*phi(ix,iy-2,iz,3) +b_FD6(-1)*phi(ix,iy-1,iz,3) +b_FD6(0)*phi(ix,iy,iz,3) +b_FD6(+1)*phi(ix,iy+1,iz,3) +b_FD6(+2)*phi(ix,iy+2,iz,3) +b_FD6(+3)*phi(ix,iy+3,iz,3))*dy2_inv
                        w_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,3) +b_FD6(-2)*phi(ix,iy,iz-2,3) +b_FD6(-1)*phi(ix,iy,iz-1,3) +b_FD6(0)*phi(ix,iy,iz,3) +b_FD6(+1)*phi(ix,iy,iz+1,3) +b_FD6(+2)*phi(ix,iy,iz+2,3) +b_FD6(+3)*phi(ix,iy,iz+3,3))*dz2_inv
                        
                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)

                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))

                        rhs(ix,iy,iz,1) = (-u*u_dx - v*u_dy - w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = (-u*v_dx - v*v_dy - w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = (-u*w_dx - v*w_dy - w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        endif

    case("FD_4th_central_optimized")
        !-----------------------------------------------------------------------
        ! 4th order (tam&web optimized scheme)
        !-----------------------------------------------------------------------
        ! Note: a(0) does NOT appear (it is zero...)
        if (params_acm%skew_symmetry) then
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! first derivatives of u, v, p
                        u_dx = (a_TW4(-3)*phi(ix-3,iy,iz,1) +a_TW4(-2)*phi(ix-2,iy,iz,1) +a_TW4(-1)*phi(ix-1,iy,iz,1) +a_TW4(+1)*phi(ix+1,iy,iz,1) +a_TW4(+2)*phi(ix+2,iy,iz,1) +a_TW4(+3)*phi(ix+3,iy,iz,1))*dx_inv
                        u_dy = (a_TW4(-3)*phi(ix,iy-3,iz,1) +a_TW4(-2)*phi(ix,iy-2,iz,1) +a_TW4(-1)*phi(ix,iy-1,iz,1) +a_TW4(+1)*phi(ix,iy+1,iz,1) +a_TW4(+2)*phi(ix,iy+2,iz,1) +a_TW4(+3)*phi(ix,iy+3,iz,1))*dy_inv
                        u_dz = (a_TW4(-3)*phi(ix,iy,iz-3,1) +a_TW4(-2)*phi(ix,iy,iz-2,1) +a_TW4(-1)*phi(ix,iy,iz-1,1) +a_TW4(+1)*phi(ix,iy,iz+1,1) +a_TW4(+2)*phi(ix,iy,iz+2,1) +a_TW4(+3)*phi(ix,iy,iz+3,1))*dz_inv
                        v_dx = (a_TW4(-3)*phi(ix-3,iy,iz,2) +a_TW4(-2)*phi(ix-2,iy,iz,2) +a_TW4(-1)*phi(ix-1,iy,iz,2) +a_TW4(+1)*phi(ix+1,iy,iz,2) +a_TW4(+2)*phi(ix+2,iy,iz,2) +a_TW4(+3)*phi(ix+3,iy,iz,2))*dx_inv
                        v_dy = (a_TW4(-3)*phi(ix,iy-3,iz,2) +a_TW4(-2)*phi(ix,iy-2,iz,2) +a_TW4(-1)*phi(ix,iy-1,iz,2) +a_TW4(+1)*phi(ix,iy+1,iz,2) +a_TW4(+2)*phi(ix,iy+2,iz,2) +a_TW4(+3)*phi(ix,iy+3,iz,2))*dy_inv
                        v_dz = (a_TW4(-3)*phi(ix,iy,iz-3,2) +a_TW4(-2)*phi(ix,iy,iz-2,2) +a_TW4(-1)*phi(ix,iy,iz-1,2) +a_TW4(+1)*phi(ix,iy,iz+1,2) +a_TW4(+2)*phi(ix,iy,iz+2,2) +a_TW4(+3)*phi(ix,iy,iz+3,2))*dz_inv
                        w_dx = (a_TW4(-3)*phi(ix-3,iy,iz,3) +a_TW4(-2)*phi(ix-2,iy,iz,3) +a_TW4(-1)*phi(ix-1,iy,iz,3) +a_TW4(+1)*phi(ix+1,iy,iz,3) +a_TW4(+2)*phi(ix+2,iy,iz,3) +a_TW4(+3)*phi(ix+3,iy,iz,3))*dx_inv
                        w_dy = (a_TW4(-3)*phi(ix,iy-3,iz,3) +a_TW4(-2)*phi(ix,iy-2,iz,3) +a_TW4(-1)*phi(ix,iy-1,iz,3) +a_TW4(+1)*phi(ix,iy+1,iz,3) +a_TW4(+2)*phi(ix,iy+2,iz,3) +a_TW4(+3)*phi(ix,iy+3,iz,3))*dy_inv
                        w_dz = (a_TW4(-3)*phi(ix,iy,iz-3,3) +a_TW4(-2)*phi(ix,iy,iz-2,3) +a_TW4(-1)*phi(ix,iy,iz-1,3) +a_TW4(+1)*phi(ix,iy,iz+1,3) +a_TW4(+2)*phi(ix,iy,iz+2,3) +a_TW4(+3)*phi(ix,iy,iz+3,3))*dz_inv
                        p_dx = (a_TW4(-3)*phi(ix-3,iy,iz,4) +a_TW4(-2)*phi(ix-2,iy,iz,4) +a_TW4(-1)*phi(ix-1,iy,iz,4) +a_TW4(+1)*phi(ix+1,iy,iz,4) +a_TW4(+2)*phi(ix+2,iy,iz,4) +a_TW4(+3)*phi(ix+3,iy,iz,4))*dx_inv
                        p_dy = (a_TW4(-3)*phi(ix,iy-3,iz,4) +a_TW4(-2)*phi(ix,iy-2,iz,4) +a_TW4(-1)*phi(ix,iy-1,iz,4) +a_TW4(+1)*phi(ix,iy+1,iz,4) +a_TW4(+2)*phi(ix,iy+2,iz,4) +a_TW4(+3)*phi(ix,iy+3,iz,4))*dy_inv
                        p_dz = (a_TW4(-3)*phi(ix,iy,iz-3,4) +a_TW4(-2)*phi(ix,iy,iz-2,4) +a_TW4(-1)*phi(ix,iy,iz-1,4) +a_TW4(+1)*phi(ix,iy,iz+1,4) +a_TW4(+2)*phi(ix,iy,iz+2,4) +a_TW4(+3)*phi(ix,iy,iz+3,4))*dz_inv
    
                        uu_dx = (a_TW4(-3)*phi(ix-3,iy,iz,1)*phi(ix-3,iy,iz,1) +a_TW4(-2)*phi(ix-2,iy,iz,1)*phi(ix-2,iy,iz,1) +a_TW4(-1)*phi(ix-1,iy,iz,1)*phi(ix-1,iy,iz,1) +a_TW4(+1)*phi(ix+1,iy,iz,1)*phi(ix+1,iy,iz,1) +a_TW4(+2)*phi(ix+2,iy,iz,1)*phi(ix+2,iy,iz,1) +a_TW4(+3)*phi(ix+3,iy,iz,1)*phi(ix+3,iy,iz,1))*dx_inv
                        uv_dy = (a_TW4(-3)*phi(ix,iy-3,iz,1)*phi(ix,iy-3,iz,2) +a_TW4(-2)*phi(ix,iy-2,iz,1)*phi(ix,iy-2,iz,2) +a_TW4(-1)*phi(ix,iy-1,iz,1)*phi(ix,iy-1,iz,2) +a_TW4(+1)*phi(ix,iy+1,iz,1)*phi(ix,iy+1,iz,2) +a_TW4(+2)*phi(ix,iy+2,iz,1)*phi(ix,iy+2,iz,2) +a_TW4(+3)*phi(ix,iy+3,iz,1)*phi(ix,iy+3,iz,2))*dy_inv
                        uw_dz = (a_TW4(-3)*phi(ix,iy,iz-3,1)*phi(ix,iy,iz-3,3) +a_TW4(-2)*phi(ix,iy,iz-2,1)*phi(ix,iy,iz-2,3) +a_TW4(-1)*phi(ix,iy,iz-1,1)*phi(ix,iy,iz-1,3) +a_TW4(+1)*phi(ix,iy,iz+1,1)*phi(ix,iy,iz+1,3) +a_TW4(+2)*phi(ix,iy,iz+2,1)*phi(ix,iy,iz+2,3) +a_TW4(+3)*phi(ix,iy,iz+3,1)*phi(ix,iy,iz+3,3))*dz_inv

                        vu_dx = (a_TW4(-3)*phi(ix-3,iy,iz,2)*phi(ix-3,iy,iz,1) +a_TW4(-2)*phi(ix-2,iy,iz,2)*phi(ix-2,iy,iz,1) +a_TW4(-1)*phi(ix-1,iy,iz,2)*phi(ix-1,iy,iz,1) +a_TW4(+1)*phi(ix+1,iy,iz,2)*phi(ix+1,iy,iz,1) +a_TW4(+2)*phi(ix+2,iy,iz,2)*phi(ix+2,iy,iz,1) +a_TW4(+3)*phi(ix+3,iy,iz,2)*phi(ix+3,iy,iz,1))*dx_inv
                        vv_dy = (a_TW4(-3)*phi(ix,iy-3,iz,2)*phi(ix,iy-3,iz,2) +a_TW4(-2)*phi(ix,iy-2,iz,2)*phi(ix,iy-2,iz,2) +a_TW4(-1)*phi(ix,iy-1,iz,2)*phi(ix,iy-1,iz,2) +a_TW4(+1)*phi(ix,iy+1,iz,2)*phi(ix,iy+1,iz,2) +a_TW4(+2)*phi(ix,iy+2,iz,2)*phi(ix,iy+2,iz,2) +a_TW4(+3)*phi(ix,iy+3,iz,2)*phi(ix,iy+3,iz,2))*dy_inv
                        vw_dz = (a_TW4(-3)*phi(ix,iy,iz-3,2)*phi(ix,iy,iz-3,3) +a_TW4(-2)*phi(ix,iy,iz-2,2)*phi(ix,iy,iz-2,3) +a_TW4(-1)*phi(ix,iy,iz-1,2)*phi(ix,iy,iz-1,3) +a_TW4(+1)*phi(ix,iy,iz+1,2)*phi(ix,iy,iz+1,3) +a_TW4(+2)*phi(ix,iy,iz+2,2)*phi(ix,iy,iz+2,3) +a_TW4(+3)*phi(ix,iy,iz+3,2)*phi(ix,iy,iz+3,3))*dz_inv

                        wu_dx = (a_TW4(-3)*phi(ix-3,iy,iz,3)*phi(ix-3,iy,iz,1) +a_TW4(-2)*phi(ix-2,iy,iz,3)*phi(ix-2,iy,iz,1) +a_TW4(-1)*phi(ix-1,iy,iz,3)*phi(ix-1,iy,iz,1) +a_TW4(+1)*phi(ix+1,iy,iz,3)*phi(ix+1,iy,iz,1) +a_TW4(+2)*phi(ix+2,iy,iz,3)*phi(ix+2,iy,iz,1) +a_TW4(+3)*phi(ix+3,iy,iz,3)*phi(ix+3,iy,iz,1))*dx_inv
                        wv_dy = (a_TW4(-3)*phi(ix,iy-3,iz,3)*phi(ix,iy-3,iz,2) +a_TW4(-2)*phi(ix,iy-2,iz,3)*phi(ix,iy-2,iz,2) +a_TW4(-1)*phi(ix,iy-1,iz,3)*phi(ix,iy-1,iz,2) +a_TW4(+1)*phi(ix,iy+1,iz,3)*phi(ix,iy+1,iz,2) +a_TW4(+2)*phi(ix,iy+2,iz,3)*phi(ix,iy+2,iz,2) +a_TW4(+3)*phi(ix,iy+3,iz,3)*phi(ix,iy+3,iz,2))*dy_inv
                        ww_dz = (a_TW4(-3)*phi(ix,iy,iz-3,3)*phi(ix,iy,iz-3,3) +a_TW4(-2)*phi(ix,iy,iz-2,3)*phi(ix,iy,iz-2,3) +a_TW4(-1)*phi(ix,iy,iz-1,3)*phi(ix,iy,iz-1,3) +a_TW4(+1)*phi(ix,iy,iz+1,3)*phi(ix,iy,iz+1,3) +a_TW4(+2)*phi(ix,iy,iz+2,3)*phi(ix,iy,iz+2,3) +a_TW4(+3)*phi(ix,iy,iz+3,3)*phi(ix,iy,iz+3,3))*dz_inv
                       
                        ! second derivatives of u, v and w
                        u_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,1) +b_FD4(-1)*phi(ix-1,iy,iz,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix+1,iy,iz,1) +b_FD4(+2)*phi(ix+2,iy,iz,1))*dx2_inv
                        u_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,1) +b_FD4(-1)*phi(ix,iy-1,iz,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix,iy+1,iz,1) +b_FD4(+2)*phi(ix,iy+2,iz,1))*dy2_inv
                        u_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,1) +b_FD4(-1)*phi(ix,iy,iz-1,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix,iy,iz+1,1) +b_FD4(+2)*phi(ix,iy,iz+2,1))*dz2_inv
                        v_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,2) +b_FD4(-1)*phi(ix-1,iy,iz,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix+1,iy,iz,2) +b_FD4(+2)*phi(ix+2,iy,iz,2))*dx2_inv
                        v_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,2) +b_FD4(-1)*phi(ix,iy-1,iz,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix,iy+1,iz,2) +b_FD4(+2)*phi(ix,iy+2,iz,2))*dy2_inv
                        v_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,2) +b_FD4(-1)*phi(ix,iy,iz-1,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix,iy,iz+1,2) +b_FD4(+2)*phi(ix,iy,iz+2,2))*dz2_inv
                        w_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,3) +b_FD4(-1)*phi(ix-1,iy,iz,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix+1,iy,iz,3) +b_FD4(+2)*phi(ix+2,iy,iz,3))*dx2_inv
                        w_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,3) +b_FD4(-1)*phi(ix,iy-1,iz,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix,iy+1,iz,3) +b_FD4(+2)*phi(ix,iy+2,iz,3))*dy2_inv
                        w_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,3) +b_FD4(-1)*phi(ix,iy,iz-1,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix,iy,iz+1,3) +b_FD4(+2)*phi(ix,iy,iz+2,3))*dz2_inv
    
                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)
    
                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))
    
                        ! Rhs in skew-symmetric formulation (of nonlinear term)
                        ! see Reiss, J. A Family of Energy Stable, Skew-Symmetric Finite Difference Schemes on Collocated Grids. J Sci Comput 65, 821–838 (2015).
                        rhs(ix,iy,iz,1) = -0.5_rk*(uu_dx + uv_dy + uw_dz   + u*u_dx + v*u_dy + w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = -0.5_rk*(vu_dx + vv_dy + vw_dz   + u*v_dx + v*v_dy + w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = -0.5_rk*(wu_dx + wv_dy + ww_dz   + u*w_dx + v*w_dy + w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        else
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! first derivatives of u, v, p
                        u_dx = (a_TW4(-3)*phi(ix-3,iy,iz,1) +a_TW4(-2)*phi(ix-2,iy,iz,1) +a_TW4(-1)*phi(ix-1,iy,iz,1) +a_TW4(+1)*phi(ix+1,iy,iz,1) +a_TW4(+2)*phi(ix+2,iy,iz,1) +a_TW4(+3)*phi(ix+3,iy,iz,1))*dx_inv
                        u_dy = (a_TW4(-3)*phi(ix,iy-3,iz,1) +a_TW4(-2)*phi(ix,iy-2,iz,1) +a_TW4(-1)*phi(ix,iy-1,iz,1) +a_TW4(+1)*phi(ix,iy+1,iz,1) +a_TW4(+2)*phi(ix,iy+2,iz,1) +a_TW4(+3)*phi(ix,iy+3,iz,1))*dy_inv
                        u_dz = (a_TW4(-3)*phi(ix,iy,iz-3,1) +a_TW4(-2)*phi(ix,iy,iz-2,1) +a_TW4(-1)*phi(ix,iy,iz-1,1) +a_TW4(+1)*phi(ix,iy,iz+1,1) +a_TW4(+2)*phi(ix,iy,iz+2,1) +a_TW4(+3)*phi(ix,iy,iz+3,1))*dz_inv
                        v_dx = (a_TW4(-3)*phi(ix-3,iy,iz,2) +a_TW4(-2)*phi(ix-2,iy,iz,2) +a_TW4(-1)*phi(ix-1,iy,iz,2) +a_TW4(+1)*phi(ix+1,iy,iz,2) +a_TW4(+2)*phi(ix+2,iy,iz,2) +a_TW4(+3)*phi(ix+3,iy,iz,2))*dx_inv
                        v_dy = (a_TW4(-3)*phi(ix,iy-3,iz,2) +a_TW4(-2)*phi(ix,iy-2,iz,2) +a_TW4(-1)*phi(ix,iy-1,iz,2) +a_TW4(+1)*phi(ix,iy+1,iz,2) +a_TW4(+2)*phi(ix,iy+2,iz,2) +a_TW4(+3)*phi(ix,iy+3,iz,2))*dy_inv
                        v_dz = (a_TW4(-3)*phi(ix,iy,iz-3,2) +a_TW4(-2)*phi(ix,iy,iz-2,2) +a_TW4(-1)*phi(ix,iy,iz-1,2) +a_TW4(+1)*phi(ix,iy,iz+1,2) +a_TW4(+2)*phi(ix,iy,iz+2,2) +a_TW4(+3)*phi(ix,iy,iz+3,2))*dz_inv
                        w_dx = (a_TW4(-3)*phi(ix-3,iy,iz,3) +a_TW4(-2)*phi(ix-2,iy,iz,3) +a_TW4(-1)*phi(ix-1,iy,iz,3) +a_TW4(+1)*phi(ix+1,iy,iz,3) +a_TW4(+2)*phi(ix+2,iy,iz,3) +a_TW4(+3)*phi(ix+3,iy,iz,3))*dx_inv
                        w_dy = (a_TW4(-3)*phi(ix,iy-3,iz,3) +a_TW4(-2)*phi(ix,iy-2,iz,3) +a_TW4(-1)*phi(ix,iy-1,iz,3) +a_TW4(+1)*phi(ix,iy+1,iz,3) +a_TW4(+2)*phi(ix,iy+2,iz,3) +a_TW4(+3)*phi(ix,iy+3,iz,3))*dy_inv
                        w_dz = (a_TW4(-3)*phi(ix,iy,iz-3,3) +a_TW4(-2)*phi(ix,iy,iz-2,3) +a_TW4(-1)*phi(ix,iy,iz-1,3) +a_TW4(+1)*phi(ix,iy,iz+1,3) +a_TW4(+2)*phi(ix,iy,iz+2,3) +a_TW4(+3)*phi(ix,iy,iz+3,3))*dz_inv
                        p_dx = (a_TW4(-3)*phi(ix-3,iy,iz,4) +a_TW4(-2)*phi(ix-2,iy,iz,4) +a_TW4(-1)*phi(ix-1,iy,iz,4) +a_TW4(+1)*phi(ix+1,iy,iz,4) +a_TW4(+2)*phi(ix+2,iy,iz,4) +a_TW4(+3)*phi(ix+3,iy,iz,4))*dx_inv
                        p_dy = (a_TW4(-3)*phi(ix,iy-3,iz,4) +a_TW4(-2)*phi(ix,iy-2,iz,4) +a_TW4(-1)*phi(ix,iy-1,iz,4) +a_TW4(+1)*phi(ix,iy+1,iz,4) +a_TW4(+2)*phi(ix,iy+2,iz,4) +a_TW4(+3)*phi(ix,iy+3,iz,4))*dy_inv
                        p_dz = (a_TW4(-3)*phi(ix,iy,iz-3,4) +a_TW4(-2)*phi(ix,iy,iz-2,4) +a_TW4(-1)*phi(ix,iy,iz-1,4) +a_TW4(+1)*phi(ix,iy,iz+1,4) +a_TW4(+2)*phi(ix,iy,iz+2,4) +a_TW4(+3)*phi(ix,iy,iz+3,4))*dz_inv

                        ! second derivatives of u, v and w
                        u_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,1) +b_FD4(-1)*phi(ix-1,iy,iz,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix+1,iy,iz,1) +b_FD4(+2)*phi(ix+2,iy,iz,1))*dx2_inv
                        u_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,1) +b_FD4(-1)*phi(ix,iy-1,iz,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix,iy+1,iz,1) +b_FD4(+2)*phi(ix,iy+2,iz,1))*dy2_inv
                        u_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,1) +b_FD4(-1)*phi(ix,iy,iz-1,1) +b_FD4(0)*phi(ix,iy,iz,1) +b_FD4(+1)*phi(ix,iy,iz+1,1) +b_FD4(+2)*phi(ix,iy,iz+2,1))*dz2_inv
                        v_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,2) +b_FD4(-1)*phi(ix-1,iy,iz,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix+1,iy,iz,2) +b_FD4(+2)*phi(ix+2,iy,iz,2))*dx2_inv
                        v_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,2) +b_FD4(-1)*phi(ix,iy-1,iz,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix,iy+1,iz,2) +b_FD4(+2)*phi(ix,iy+2,iz,2))*dy2_inv
                        v_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,2) +b_FD4(-1)*phi(ix,iy,iz-1,2) +b_FD4(0)*phi(ix,iy,iz,2) +b_FD4(+1)*phi(ix,iy,iz+1,2) +b_FD4(+2)*phi(ix,iy,iz+2,2))*dz2_inv
                        w_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,3) +b_FD4(-1)*phi(ix-1,iy,iz,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix+1,iy,iz,3) +b_FD4(+2)*phi(ix+2,iy,iz,3))*dx2_inv
                        w_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,3) +b_FD4(-1)*phi(ix,iy-1,iz,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix,iy+1,iz,3) +b_FD4(+2)*phi(ix,iy+2,iz,3))*dy2_inv
                        w_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,3) +b_FD4(-1)*phi(ix,iy,iz-1,3) +b_FD4(0)*phi(ix,iy,iz,3) +b_FD4(+1)*phi(ix,iy,iz+1,3) +b_FD4(+2)*phi(ix,iy,iz+2,3))*dz2_inv

                        u = phi(ix, iy, iz, 1)
                        v = phi(ix, iy, iz, 2)
                        w = phi(ix, iy, iz, 3)
                        p = phi(ix, iy, iz, 4)

                        chi = mask(ix,iy,iz,1) * C_eta_inv
                        penalx = -chi * (u - mask(ix,iy,iz,2))
                        penaly = -chi * (v - mask(ix,iy,iz,3))
                        penalz = -chi * (w - mask(ix,iy,iz,4))

                        rhs(ix,iy,iz,1) = (-u*u_dx - v*u_dy - w*u_dz) -p_dx + nu*(u_dxdx + u_dydy + u_dzdz) + penalx
                        rhs(ix,iy,iz,2) = (-u*v_dx - v*v_dy - w*v_dz) -p_dy + nu*(v_dxdx + v_dydy + v_dzdz) + penaly
                        rhs(ix,iy,iz,3) = (-u*w_dx - v*w_dy - w*w_dz) -p_dz + nu*(w_dxdx + w_dydy + w_dzdz) + penalz
                        rhs(ix,iy,iz,4) = -(c_0**2)*(u_dx + v_dy + w_dz) - gamma*p
                    end do
                end do
            end do
        endif

    case default
        call abort(441167, "3d Discretization unkown "//order_discretization//", I ll walk into the light now." )

    end select

    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    if (params_acm%use_sponge) then
        ! avoid division by multiplying with inverse
        C_sponge_inv = 1.0_rk / params_acm%C_sponge

        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! NOTE: the sponge term acts, if active, on ALL components, ux,uy,p
                    ! which is different from the penalization term, which acts only on ux,uy and not p
                    ! NOTE: sponge mask set in hvy_mask
                    spo = mask(ix,iy,iz,6) * C_sponge_inv

                    rhs(ix,iy,iz,1) = rhs(ix,iy,iz,1) - (phi(ix,iy,iz,1)-params_acm%u_mean_set(1)) * spo
                    rhs(ix,iy,iz,2) = rhs(ix,iy,iz,2) - (phi(ix,iy,iz,2)-params_acm%u_mean_set(2)) * spo
                    rhs(ix,iy,iz,3) = rhs(ix,iy,iz,3) - (phi(ix,iy,iz,3)-params_acm%u_mean_set(3)) * spo
                    rhs(ix,iy,iz,4) = rhs(ix,iy,iz,4) - (phi(ix,iy,iz,4))*spo
                end do
            end do
        end do
    end if


    ! --------------------------------------------------------------------------
    ! HIT linear forcing
    ! ATTENTION! This is at last position because I modify phi to avoid a 3-nested do-loop to subtract mean-flow
    ! --------------------------------------------------------------------------
    if (params_acm%HIT_linear_forcing) then
        G_gain = params_acm%HIT_gain
        e_kin_set = params_acm%HIT_energy * product(params_acm%domain_size(1:params_acm%dim))
        t_l_inf = 1.0_rk ! sqrt(nu / epsilon), should be adapted to by setting gain
        ! forcing after Bassene konstant energy (2016)
        A_forcing = (params_acm%dissipation - G_gain * (params_acm%e_kin - e_kin_set) / t_l_inf) / (2.0*params_acm%e_kin)
        
        ! Forcing should not be applied onto the mean-flow, so we subtract it out
        ! ATTENTION! This modifies phi so it should be the last statement with phi
        phi(:,:,:,1) = phi(:,:,:,1) - params_acm%mean_flow(1)
        phi(:,:,:,2) = phi(:,:,:,2) - params_acm%mean_flow(2)
        phi(:,:,:,3) = phi(:,:,:,3) - params_acm%mean_flow(3)
        ! cancel out mean_flow, this is quite brutal but let's see what it does
        rhs(:,:,:,1) = rhs(:,:,:,1) - 1e4*params_acm%mean_flow(1)
        rhs(:,:,:,2) = rhs(:,:,:,2) - 1e4*params_acm%mean_flow(2)
        rhs(:,:,:,3) = rhs(:,:,:,3) - 1e4*params_acm%mean_flow(3)
        ! apply forcing
        rhs(:,:,:,1:3) = rhs(:,:,:,1:3) + A_forcing*phi(:,:,:,1:3)
    endif


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
    character(len=cshort), intent(in)       :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    integer(kind=ik) :: ix, iy, iz, iscalar, j

    ! coefficients for Tam&Webb (4th order 1st derivative)
    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b_FD4(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! 6th order FD scheme
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/) ! 1st derivative
    real(kind=rk), parameter :: b_FD6(-3:3) = (/ 1.0_rk/90.0_rk, -3.0_rk/20.0_rk, 3.0_rk/2.0_rk, -49.0_rk/18.0_rk, 3.0_rk/2.0_rk, -3.0_rk/20.0_rk, 1.0_rk/90.0_rk/) ! 2nd derivative

    real(kind=rk) :: kappa, x, y, z, masksource, nu, R, R0sq
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

    nu = params_acm%nu

    select case(order_discretization)
    case("FD_2nd_central")
        call abort(2208191, "passive scalar implemented only with 4th order....sorry, I am lazy")

    case("FD_4th_central_optimized")
        !-----------------------------------------------------------------------
        ! 4th order
        !-----------------------------------------------------------------------
        do iscalar = 1, params_acm%N_scalars
            ! actual index of this scalar in the array
            j = iscalar + (params_acm%dim + 1)

            ! compute diffusivity from schmidt number (and fluid viscosity)
            kappa = nu / params_acm%schmidt_numbers(iscalar)

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

            case ("inflow-x")
                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            x = x0(1) + dble(ix-(g+1))*dx(1)
                            if ( x <= params_acm%widthsource(iscalar) ) then
                                ! INFLOW
                                source(ix,iy,iz) = (1.0_rk - phi(ix,iy,iz,j)) / params_acm%C_eta ! for the source term, we use the usual dirichlet C_eta
                            endif
                        end do
                    end do
                end do

            case ("in+outflow-x")
                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            x = x0(1) + dble(ix-(g+1))*dx(1)
                            if ( x <= params_acm%widthsource(iscalar) ) then
                                ! INFLOW
                                source(ix,iy,iz) = (1.0_rk - phi(ix,iy,iz,j)) / params_acm%C_eta ! for the source term, we use the usual dirichlet C_eta
                            endif
                            if ( x >= params_acm%domain_size(1)-params_acm%widthsource(iscalar) ) then
                                ! OUTFLOW
                                source(ix,iy,iz) = (0.0_rk - phi(ix,iy,iz,j)) / params_acm%C_eta ! for the source term, we use the usual dirichlet C_eta
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

            if (params_acm%scalar_BC_type == "neumann") then
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
                            gx = (a_TW4(-3)*phi(ix-3,iy,iz,j)&
                                 +a_TW4(-2)*phi(ix-2,iy,iz,j)&
                                 +a_TW4(-1)*phi(ix-1,iy,iz,j)&
                                 +a_TW4(+2)*phi(ix+2,iy,iz,j)&
                                 +a_TW4(+3)*phi(ix+3,iy,iz,j)&
                                 +a_TW4(+1)*phi(ix+1,iy,iz,j))*dx_inv

                            gy = (a_TW4(-3)*phi(ix,iy-3,iz,j)&
                                 +a_TW4(-2)*phi(ix,iy-2,iz,j)&
                                 +a_TW4(-1)*phi(ix,iy-1,iz,j)&
                                 +a_TW4(+2)*phi(ix,iy+2,iz,j)&
                                 +a_TW4(+3)*phi(ix,iy+3,iz,j)&
                                 +a_TW4(+1)*phi(ix,iy+1,iz,j))*dy_inv

                            gz = (a_TW4(-3)*phi(ix,iy,iz-3,j)&
                                 +a_TW4(-2)*phi(ix,iy,iz-2,j)&
                                 +a_TW4(-1)*phi(ix,iy,iz-1,j)&
                                 +a_TW4(+2)*phi(ix,iy,iz+2,j)&
                                 +a_TW4(+3)*phi(ix,iy,iz+3,j)&
                                 +a_TW4(+1)*phi(ix,iy,iz+1,j))*dz_inv

                            ! gradient of mask function ( we need that for the diffusive term)
                            ! since this guy reads div( (kappa(1-mask) + eps*mask) * grad(phi) )
                            ! so this boils down to d/dx (D*gx) = D_dx*gx + D*gxx
                            ! so we need D_dx and this is kappa*(1-mask_dx)+ eps*mask_dx
                            chidx = (a_TW4(-3)*mask(ix-3,iy,iz, 1)&
                                    +a_TW4(-2)*mask(ix-2,iy,iz, 1)&
                                    +a_TW4(-1)*mask(ix-1,iy,iz, 1)&
                                    +a_TW4(+3)*mask(ix+3,iy,iz, 1)&
                                    +a_TW4(+2)*mask(ix+2,iy,iz, 1)&
                                    +a_TW4(+1)*mask(ix+1,iy,iz, 1))*dx_inv

                            chidy = (a_TW4(-3)*mask(ix,iy-3,iz, 1)&
                                    +a_TW4(-2)*mask(ix,iy-2,iz, 1)&
                                    +a_TW4(-1)*mask(ix,iy-1,iz, 1)&
                                    +a_TW4(+3)*mask(ix,iy+3,iz, 1)&
                                    +a_TW4(+2)*mask(ix,iy+2,iz, 1)&
                                    +a_TW4(+1)*mask(ix,iy+1,iz, 1))*dy_inv

                            chidz = (a_TW4(-3)*mask(ix,iy,iz-3, 1)&
                                    +a_TW4(-2)*mask(ix,iy,iz-2, 1)&
                                    +a_TW4(-1)*mask(ix,iy,iz-1, 1)&
                                    +a_TW4(+3)*mask(ix,iy,iz+3, 1)&
                                    +a_TW4(+2)*mask(ix,iy,iz+2, 1)&
                                    +a_TW4(+1)*mask(ix,iy,iz+1, 1))*dz_inv

                            D_dx = kappa*(1.0_rk-chidx) + params_acm%scalar_Ceta(iscalar) * chidx
                            D_dy = kappa*(1.0_rk-chidy) + params_acm%scalar_Ceta(iscalar) * chidy
                            D_dz = kappa*(1.0_rk-chidz) + params_acm%scalar_Ceta(iscalar) * chidz

                            ! second derivatives of passive scalar
                            gxx = (b_FD4(-2)*phi(ix-2,iy,iz ,j)&
                                  +b_FD4(-1)*phi(ix-1,iy,iz ,j)&
                                  +b_FD4( 0)*phi(ix  ,iy,iz ,j)&
                                  +b_FD4(+1)*phi(ix+1,iy,iz ,j)&
                                  +b_FD4(+2)*phi(ix+2,iy,iz ,j))*dx2_inv
                            gyy = (b_FD4(-2)*phi(ix,iy-2,iz ,j)&
                                  +b_FD4(-1)*phi(ix,iy-1,iz ,j)&
                                  +b_FD4( 0)*phi(ix,iy  ,iz ,j)&
                                  +b_FD4(+1)*phi(ix,iy+1,iz ,j)&
                                  +b_FD4(+2)*phi(ix,iy+2,iz ,j))*dy2_inv
                            gzz = (b_FD4(-2)*phi(ix,iy,iz-2 ,j)&
                                  +b_FD4(-1)*phi(ix,iy,iz-1 ,j)&
                                  +b_FD4( 0)*phi(ix,iy,iz   ,j)&
                                  +b_FD4(+1)*phi(ix,iy,iz+1 ,j)&
                                  +b_FD4(+2)*phi(ix,iy,iz+2 ,j))*dz2_inv

                            ! assemble everything
                            rhs(ix,iy,iz,j) = wx*gx + wy*gy + wz*gz & ! penalized convection term
                            + D_dx*gx + D*gxx & ! penalized laplacian
                            + D_dy*gy + D*gyy &
                            + D_dz*gz + D*gzz &
                            + source(ix,iy,iz)
                        end do
                    end do
                end do
            elseif (params_acm%scalar_BC_type == "dirichlet") then
                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            ux = phi(ix,iy,iz,1)
                            uy = phi(ix,iy,iz,2)
                            uz = phi(ix,iy,iz,3)

                            chi = mask(ix,iy,iz,1)

                            ! gradient
                            phi_dx = (a_TW4(-3)*phi(ix-3,iy,iz,j) &
                                    + a_TW4(-2)*phi(ix-2,iy,iz,j) &
                                    + a_TW4(-1)*phi(ix-1,iy,iz,j) &
                                    + a_TW4(+1)*phi(ix+1,iy,iz,j) &
                                    + a_TW4(+2)*phi(ix+2,iy,iz,j) &
                                    + a_TW4(+3)*phi(ix+3,iy,iz,j))*dx_inv

                            phi_dy = (a_TW4(-3)*phi(ix,iy-3,iz,j) &
                                    + a_TW4(-2)*phi(ix,iy-2,iz,j) &
                                    + a_TW4(-1)*phi(ix,iy-1,iz,j) &
                                    + a_TW4(+1)*phi(ix,iy+1,iz,j) &
                                    + a_TW4(+2)*phi(ix,iy+2,iz,j) &
                                    + a_TW4(+3)*phi(ix,iy+3,iz,j))*dy_inv

                            phi_dz = (a_TW4(-3)*phi(ix,iy,iz-3,j) &
                                    + a_TW4(-2)*phi(ix,iy,iz-2,j) &
                                    + a_TW4(-1)*phi(ix,iy,iz-1,j) &
                                    + a_TW4(+1)*phi(ix,iy,iz+1,j) &
                                    + a_TW4(+2)*phi(ix,iy,iz+2,j) &
                                    + a_TW4(+3)*phi(ix,iy,iz+3,j))*dz_inv

                            ! laplace
                            phi_dxdx = (  b_FD4(-2)*phi(ix-2,iy,iz,j) &
                                        + b_FD4(-1)*phi(ix-1,iy,iz,j) &
                                        + b_FD4( 0)*phi(ix  ,iy,iz,j) &
                                        + b_FD4(+1)*phi(ix+1,iy,iz,j) &
                                        + b_FD4(+2)*phi(ix+2,iy,iz,j))*dx2_inv

                            phi_dydy = (  b_FD4(-2)*phi(ix,iy-2,iz,j) &
                                        + b_FD4(-1)*phi(ix,iy-1,iz,j) &
                                        + b_FD4( 0)*phi(ix,iy  ,iz,j) &
                                        + b_FD4(+1)*phi(ix,iy+1,iz,j) &
                                        + b_FD4(+2)*phi(ix,iy+2,iz,j))*dy2_inv

                            phi_dzdz = (  b_FD4(-2)*phi(ix,iy,iz-2,j) &
                                        + b_FD4(-1)*phi(ix,iy,iz-1,j) &
                                        + b_FD4( 0)*phi(ix,iy,iz  ,j) &
                                        + b_FD4(+1)*phi(ix,iy,iz+1,j) &
                                        + b_FD4(+2)*phi(ix,iy,iz+2,j))*dz2_inv

                            ! easy RHS for dirichlet BC
                            rhs(ix,iy,iz,j) = -ux*phi_dx -uy*phi_dy -uz*phi_dz &
                            + kappa*(phi_dxdx + phi_dydy + phi_dzdz) &
                            - chi*(phi(ix,iy,iz,j) - 0.0_rk) / params_acm%C_eta & ! Dirichlet penalization for obstacle mask (instead of Neumann penalization)
                            + source(ix, iy, iz) ! source term is actually a dirichlet penalization term as well
                        enddo
                    enddo
                enddo
            else
                call abort(22092234, "scalar_BC_type unkown"//trim(adjustl(params_acm%scalar_BC_type))//". Time to fake a smile :)" )
            endif
        end do ! loop over scalars


    case default
        call abort(441167, "3d Discretization unkown "//trim(adjustl(order_discretization))//", I ll walk into the light now." )
    end select
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
    character(len=cshort), intent(in)       :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    integer(kind=ik) :: ix, iy, iscalar, j

    real(kind=rk) :: kappa, x, y, masksource, nu, R
    real(kind=rk) :: dx_inv, dy_inv, dx2_inv, dy2_inv
    real(kind=rk) :: ux, uy, usx, usy, wx, wy, gx, gy, D, chi, chidx, chidy, D_dx, D_dy, gxx, gyy
    real(kind=rk) :: phi_dx, phi_dy, phi_dxdx, phi_dydy
    ! coefficients for Tam&Webb (4th order 1st derivative)
    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b_FD4(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! 6th order FD scheme
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/) ! 1st derivative
    real(kind=rk), parameter :: b_FD6(-3:3) = (/ 1.0_rk/90.0_rk, -3.0_rk/20.0_rk, 3.0_rk/2.0_rk, -49.0_rk/18.0_rk, 3.0_rk/2.0_rk, -3.0_rk/20.0_rk, 1.0_rk/90.0_rk/) ! 2nd derivative

    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're negligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: source(:,:,:)

    if (.not. allocated(source)) allocate(source(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
    source = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    nu = params_acm%nu


    if (order_discretization == "FD_4th_central_optimized") then
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

            case("in+outflow")
                do iy = g+1, Bs(2)+g
                    y = x0(2) + dble(iy-(g+1))*dx(2)
                    do ix = g+1, Bs(1)+g
                        x = x0(1) + dble(ix-(g+1))*dx(1)
                        if ( x <= params_acm%widthsource(iscalar) ) then
                            ! INFLOW
                            source(ix,iy,1) = (1.0_rk - phi(ix,iy,1,j)) / params_acm%C_eta ! for the source term, we use the usual dirichlet C_eta

                        elseif ( x >= params_acm%domain_size(1)-params_acm%widthsource(iscalar) ) then
                            ! OUTFLOW
                            source(ix,iy,1) = (0.0_rk - phi(ix,iy,1,j)) / params_acm%C_eta
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


            if (params_acm%scalar_BC_type == "neumann") then
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
                        gx = (a_TW4(-3)*phi(ix-3,iy,1,j) &
                             +a_TW4(-2)*phi(ix-2,iy,1,j) &
                             +a_TW4(-1)*phi(ix-1,iy,1,j) &
                             +a_TW4( 0)*phi(ix  ,iy,1,j) &
                             +a_TW4(+1)*phi(ix+1,iy,1,j) &
                             +a_TW4(+2)*phi(ix+2,iy,1,j) &
                             +a_TW4(+3)*phi(ix+3,iy,1,j))*dx_inv

                        gy = (a_TW4(-3)*phi(ix,iy-3,1,j) &
                             +a_TW4(-2)*phi(ix,iy-2,1,j) &
                             +a_TW4(-1)*phi(ix,iy-1,1,j) &
                             +a_TW4( 0)*phi(ix,iy  ,1,j) &
                             +a_TW4(+1)*phi(ix,iy+1,1,j) &
                             +a_TW4(+2)*phi(ix,iy+2,1,j) &
                             +a_TW4(+3)*phi(ix,iy+3,1,j))*dy_inv

                        ! gradient of mask function ( we need that for the diffusive term)
                        ! since this guy reads div( (kappa(1-mask) + eps*mask) * grad(phi) )
                        ! so this boils down to d/dx (D*gx) = D_dx*gx + D*gxx
                        ! so we need D_dx and this is kappa*(1-mask_dx)+ eps*mask_dx
                        chidx = (a_TW4(-3)*mask(ix-3,iy,1, 1) &
                                +a_TW4(-2)*mask(ix-2,iy,1, 1) &
                                +a_TW4(-1)*mask(ix-1,iy,1, 1) &
                                +a_TW4( 0)*mask(ix  ,iy,1, 1) &
                                +a_TW4(+1)*mask(ix+1,iy,1, 1) &
                                +a_TW4(+2)*mask(ix+2,iy,1, 1) &
                                +a_TW4(+3)*mask(ix+3,iy,1, 1))*dx_inv

                        chidy = (a_TW4(-3)*mask(ix,iy-3,1, 1) &
                                +a_TW4(-2)*mask(ix,iy-2,1, 1) &
                                +a_TW4(-1)*mask(ix,iy-1,1, 1) &
                                +a_TW4( 0)*mask(ix,iy  ,1, 1) &
                                +a_TW4(+1)*mask(ix,iy+1,1, 1) &
                                +a_TW4(+2)*mask(ix,iy+2,1, 1) &
                                +a_TW4(+3)*mask(ix,iy+3,1, 1))*dy_inv

                        D_dx = kappa*(-chidx) + params_acm%scalar_Ceta(iscalar) * chidx
                        D_dy = kappa*(-chidy) + params_acm%scalar_Ceta(iscalar) * chidy

                        ! second derivatives of passive scalar
                        gxx = (b_FD4(-2)*phi(ix-2,iy,1 ,j) &
                              +b_FD4(-1)*phi(ix-1,iy,1 ,j) &
                              +b_FD4( 0)*phi(ix  ,iy,1 ,j) &
                              +b_FD4(+1)*phi(ix+1,iy,1 ,j) &
                              +b_FD4(+2)*phi(ix+2,iy,1 ,j))*dx2_inv
                        gyy = (b_FD4(-2)*phi(ix,iy-2,1 ,j) &
                              +b_FD4(-1)*phi(ix,iy-1,1 ,j) &
                              +b_FD4( 0)*phi(ix,iy  ,1 ,j) &
                              +b_FD4(+1)*phi(ix,iy+1,1 ,j) &
                              +b_FD4(+2)*phi(ix,iy+2,1 ,j))*dy2_inv

                        ! assemble everything
                        rhs(ix,iy,1,j) = wx*gx + wy*gy & ! penalized convection term
                                       + D_dx*gx + D*gxx + D_dy*gy + D*gyy & ! penalized laplacian
                                       + source(ix, iy, 1)
                    end do
                end do
            elseif (params_acm%scalar_BC_type == "dirichlet") then

                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ux = phi(ix,iy,1,1)
                        uy = phi(ix,iy,1,2)

                        chi = mask(ix,iy,1,1)
                        usx = mask(ix,iy,1,2)
                        usy = mask(ix,iy,1,3)

                        ! gradient
                        phi_dx = (a_TW4(-3)*phi(ix-3,iy,1,j) &
                                + a_TW4(-2)*phi(ix-2,iy,1,j) &
                                + a_TW4(-1)*phi(ix-1,iy,1,j) &
                                + a_TW4( 0)*phi(ix  ,iy,1,j) &
                                + a_TW4(+1)*phi(ix+1,iy,1,j) &
                                + a_TW4(+2)*phi(ix+2,iy,1,j) &
                                + a_TW4(+3)*phi(ix+3,iy,1,j))*dx_inv
                        phi_dy = (a_TW4(-3)*phi(ix,iy-3,1,j) &
                                + a_TW4(-2)*phi(ix,iy-2,1,j) &
                                + a_TW4(-1)*phi(ix,iy-1,1,j) &
                                + a_TW4( 0)*phi(ix,iy  ,1,j) &
                                + a_TW4(+1)*phi(ix,iy+1,1,j) &
                                + a_TW4(+2)*phi(ix,iy+2,1,j) &
                                + a_TW4(+3)*phi(ix,iy+3,1,j))*dy_inv

                        ! laplace
                        phi_dxdx = (  b_FD4(-2)*phi(ix-2,iy,1,j) &
                                    + b_FD4(-1)*phi(ix-1,iy,1,j) &
                                    + b_FD4( 0)*phi(ix  ,iy,1,j) &
                                    + b_FD4(+1)*phi(ix+1,iy,1,j) &
                                    + b_FD4(+2)*phi(ix+2,iy,1,j))*dx2_inv
                        phi_dydy = (  b_FD4(-2)*phi(ix,iy-2,1,j) &
                                    + b_FD4(-1)*phi(ix,iy-1,1,j) &
                                    + b_FD4( 0)*phi(ix,iy  ,1,j) &
                                    + b_FD4(+1)*phi(ix,iy+1,1,j) &
                                    + b_FD4(+2)*phi(ix,iy+2,1,j))*dy2_inv

                        ! easy RHS for dirichlet BC
                        rhs(ix,iy,1,j) = -ux*phi_dx -uy*phi_dy &
                                       + kappa*(phi_dxdx + phi_dydy) &
                                       - chi*(phi(ix,iy,1,j) - 0.0_rk) / params_acm%C_eta & ! Dirichlet penalization for obstacle mask
                                       + source(ix, iy, 1) ! source term is actually a dirichlet penalization term as well
                    end do
                end do

            endif
        end do ! loop over scalars


    else
        call abort(2109231, "Passive scalar discretization unkown "//trim(adjustl(order_discretization))//". Just as well. I was tired anyways, sleep tight!" )
    end if
end subroutine RHS_2D_scalar



! on a block compute the dissipation rate,
! i.e. \varepsilon = 2 nu (du_i/dx_j) (du_i/dx_j) 
! and thus do not make use of the vorticity
subroutine dissipation_ACM_block(Bs, g, dx, u, dissipation_rate)
    implicit none

    !> grid parameter
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> spacing of the block
    real(kind=rk), dimension(3), intent(in) :: dx
    !> datafields
    real(kind=rk), intent(inout) :: u(:,:,:,:)
    real(kind=rk), intent(inout) :: dissipation_rate

    !> derivatives
    real(kind=rk) :: ux_dx, ux_dy, ux_dz
    real(kind=rk) :: uy_dx, uy_dy, uy_dz
    real(kind=rk) :: uz_dx, uz_dy, uz_dz

    !> inverse of dx, dy, dz
    real(kind=rk) :: dx_inv, dy_inv, dz_inv
    ! loop variables
    integer(kind=ik) :: ix, iy, iz
    ! coefficients for Tam&Webb (4th order 1st derivative)
    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    ! 4th order coefficients for second derivative
    real(kind=rk), parameter :: b_FD4(-2:2) = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)
    ! 6th order FD scheme
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/) ! 1st derivative
    real(kind=rk), parameter :: b_FD6(-3:3) = (/ 1.0_rk/90.0_rk, -3.0_rk/20.0_rk, 3.0_rk/2.0_rk, -49.0_rk/18.0_rk, 3.0_rk/2.0_rk, -3.0_rk/20.0_rk, 1.0_rk/90.0_rk/) ! 2nd derivative


    if (.not. params_acm%initialized) write(*,*) "WARNING: vorticity_ACM_block called but ACM not initialized"

    dissipation_rate = 0.0_rk

    if ( params_acm%dim == 2) then
        dx_inv = 1.0_rk / dx(1)
        dy_inv = 1.0_rk / dx(2)

        iz = 1

        select case(params_acm%discretization)
        case("FD_2nd_central")
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ux_dx = (u(ix+1,iy,iz,1) - u(ix-1,iy,iz,1))*dx_inv*0.5_rk
                    ux_dy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dx_inv*0.5_rk

                    uy_dx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dy_inv*0.5_rk
                    uy_dy = (u(ix,iy+1,iz,2) - u(ix,iy-1,iz,2))*dy_inv*0.5_rk

                    dissipation_rate = dissipation_rate + ux_dx**2 + ux_dy**2 + uy_dx**2 + uy_dy**2
                end do
            end do

        case("FD_4th_central")
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ux_dx = ( a_FD4(-2)*u(ix-2,iy,iz,1) &
                            + a_FD4(-1)*u(ix-1,iy,iz,1) &
                            + a_FD4(0) *u(ix,iy,iz,1) &
                            + a_FD4(+1)*u(ix+1,iy,iz,1) &
                            + a_FD4(+2)*u(ix+2,iy,iz,1) )*dx_inv
                            
                    ux_dy = ( a_FD4(-2)*u(ix,iy-2,iz,1) &
                            + a_FD4(-1)*u(ix,iy-1,iz,1) &
                            + a_FD4(0) *u(ix,iy,iz,1) &
                            + a_FD4(+1)*u(ix,iy+1,iz,1) &
                            + a_FD4(+2)*u(ix,iy+2,iz,1) )*dy_inv

                    uy_dx = ( a_FD4(-2)*u(ix-2,iy,iz,2) &
                            + a_FD4(-1)*u(ix-1,iy,iz,2) &
                            + a_FD4(0) *u(ix,iy,iz,2) &
                            + a_FD4(+1)*u(ix+1,iy,iz,2) &
                            + a_FD4(+2)*u(ix+2,iy,iz,2) )*dx_inv
                            
                    uy_dy = ( a_FD4(-2)*u(ix,iy-2,iz,2) &
                            + a_FD4(-1)*u(ix,iy-1,iz,2) &
                            + a_FD4(0) *u(ix,iy,iz,2) &
                            + a_FD4(+1)*u(ix,iy+1,iz,2) &
                            + a_FD4(+2)*u(ix,iy+2,iz,2) )*dy_inv

                    dissipation_rate = dissipation_rate + ux_dx**2 + ux_dy**2 + uy_dx**2 + uy_dy**2
                end do
            end do

        case("FD_6th_central")
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ux_dx = (a_FD6(-3)*u(ix-3,iy,iz,1) &
                            + a_FD6(-2)*u(ix-2,iy,iz,1) &
                            + a_FD6(-1)*u(ix-1,iy,iz,1) &
                            + a_FD6(0) *u(ix,iy,iz,1) &
                            + a_FD6(+1)*u(ix+1,iy,iz,1) &
                            + a_FD6(+2)*u(ix+2,iy,iz,1) &
                            + a_FD6(+3)*u(ix+3,iy,iz,1))*dx_inv

                    ux_dy = (a_FD6(-3)*u(ix,iy-3,iz,1) &
                            + a_FD6(-2)*u(ix,iy-2,iz,1) &
                            + a_FD6(-1)*u(ix,iy-1,iz,1) &
                            + a_FD6(0) *u(ix,iy,iz,1) &
                            + a_FD6(+1)*u(ix,iy+1,iz,1) &
                            + a_FD6(+2)*u(ix,iy+2,iz,1) &
                            + a_FD6(+3)*u(ix,iy+3,iz,1))*dy_inv

                    uy_dx = (a_FD6(-3)*u(ix-3,iy,iz,2) &
                            + a_FD6(-2)*u(ix-2,iy,iz,2) &
                            + a_FD6(-1)*u(ix-1,iy,iz,2) &
                            + a_FD6(0) *u(ix,iy,iz,2) &
                            + a_FD6(+1)*u(ix+1,iy,iz,2) &
                            + a_FD6(+2)*u(ix+2,iy,iz,2) &
                            + a_FD6(+3)*u(ix+3,iy,iz,2))*dx_inv

                    uy_dy = (a_FD6(-3)*u(ix,iy-3,iz,2) &
                            + a_FD6(-2)*u(ix,iy-2,iz,2) &
                            + a_FD6(-1)*u(ix,iy-1,iz,2) &
                            + a_FD6(0) *u(ix,iy,iz,2) &
                            + a_FD6(+1)*u(ix,iy+1,iz,2) &
                            + a_FD6(+2)*u(ix,iy+2,iz,2) &
                            + a_FD6(+3)*u(ix,iy+3,iz,2))*dy_inv

                    dissipation_rate = dissipation_rate + ux_dx**2 + ux_dy**2 + uy_dx**2 + uy_dy**2
                end do
            end do

        case("FD_4th_central_optimized")
            call abort(2362976, "dissipation rate with FD_4th_central_optimized not yet done.")

        case default
            call abort(1902201, "unknown order_discretization in ACM vorticity")
        end select

    else
        dx_inv = 1.0_rk / dx(1)
        dy_inv = 1.0_rk / dx(2)
        dz_inv = 1.0_rk / dx(3)

        select case(params_acm%discretization)
        case("FD_2nd_central")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ux_dx = (u(ix+1,iy,iz,1) - u(ix-1,iy,iz,1))*dx_inv*0.5_rk
                        ux_dy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dx_inv*0.5_rk
                        ux_dz = (u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1))*dx_inv*0.5_rk

                        uy_dx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dy_inv*0.5_rk
                        uy_dy = (u(ix,iy+1,iz,2) - u(ix,iy-1,iz,2))*dy_inv*0.5_rk
                        uy_dz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dy_inv*0.5_rk

                        uz_dx = (u(ix+1,iy,iz,3) - u(ix-1,iy,iz,3))*dz_inv*0.5_rk
                        uz_dy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dz_inv*0.5_rk
                        uz_dz = (u(ix,iy,iz+1,3) - u(ix,iy,iz-1,3))*dz_inv*0.5_rk

                        dissipation_rate = dissipation_rate + &
                        ux_dx**2 + ux_dy**2 + ux_dz**2 + uy_dx**2 + uy_dy**2 + uy_dz**2 + uz_dx**2 + uz_dy**2 + uz_dz**2
                    end do
                end do
            end do

        case("FD_4th_central")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ux_dx = ( a_FD4(-2)*u(ix-2,iy,iz,1) &
                                + a_FD4(-1)*u(ix-1,iy,iz,1) &
                                + a_FD4(0) *u(ix,iy,iz,1) &
                                + a_FD4(+1)*u(ix+1,iy,iz,1) &
                                + a_FD4(+2)*u(ix+2,iy,iz,1) )*dx_inv
                                
                        ux_dy = ( a_FD4(-2)*u(ix,iy-2,iz,1) &
                                + a_FD4(-1)*u(ix,iy-1,iz,1) &
                                + a_FD4(0) *u(ix,iy,iz,1) &
                                + a_FD4(+1)*u(ix,iy+1,iz,1) &
                                + a_FD4(+2)*u(ix,iy+2,iz,1) )*dy_inv

                        ux_dz = ( a_FD4(-2)*u(ix,iy,iz-2,1) &
                                + a_FD4(-1)*u(ix,iy,iz-1,1) &
                                + a_FD4(0) *u(ix,iy,iz,1) &
                                + a_FD4(+1)*u(ix,iy,iz+1,1) &
                                + a_FD4(+2)*u(ix,iy,iz+2,1) )*dz_inv

                        uy_dx = ( a_FD4(-2)*u(ix-2,iy,iz,2) &
                                + a_FD4(-1)*u(ix-1,iy,iz,2) &
                                + a_FD4(0) *u(ix,iy,iz,2) &
                                + a_FD4(+1)*u(ix+1,iy,iz,2) &
                                + a_FD4(+2)*u(ix+2,iy,iz,2) )*dx_inv
                                
                        uy_dy = ( a_FD4(-2)*u(ix,iy-2,iz,2) &
                                + a_FD4(-1)*u(ix,iy-1,iz,2) &
                                + a_FD4(0) *u(ix,iy,iz,2) &
                                + a_FD4(+1)*u(ix,iy+1,iz,2) &
                                + a_FD4(+2)*u(ix,iy+2,iz,2) )*dy_inv

                        uy_dz = ( a_FD4(-2)*u(ix,iy,iz-2,2) &
                                + a_FD4(-1)*u(ix,iy,iz-1,2) &
                                + a_FD4(0) *u(ix,iy,iz,2) &
                                + a_FD4(+1)*u(ix,iy,iz+1,2) &
                                + a_FD4(+2)*u(ix,iy,iz+2,2) )*dz_inv
                        
                        uz_dx = ( a_FD4(-2)*u(ix-2,iy,iz,3) &
                                + a_FD4(-1)*u(ix-1,iy,iz,3) &
                                + a_FD4(0) *u(ix,iy,iz,3) &
                                + a_FD4(+1)*u(ix+1,iy,iz,3) &
                                + a_FD4(+2)*u(ix+2,iy,iz,3) )*dx_inv
                                
                        uz_dy = ( a_FD4(-2)*u(ix,iy-2,iz,3) &
                                + a_FD4(-1)*u(ix,iy-1,iz,3) &
                                + a_FD4(0) *u(ix,iy,iz,3) &
                                + a_FD4(+1)*u(ix,iy+1,iz,3) &
                                + a_FD4(+2)*u(ix,iy+2,iz,3) )*dy_inv

                        uz_dz = ( a_FD4(-2)*u(ix,iy,iz-2,3) &
                                + a_FD4(-1)*u(ix,iy,iz-1,3) &
                                + a_FD4(0) *u(ix,iy,iz,3) &
                                + a_FD4(+1)*u(ix,iy,iz+1,3) &
                                + a_FD4(+2)*u(ix,iy,iz+2,3) )*dz_inv

                        dissipation_rate = dissipation_rate + &
                        ux_dx**2 + ux_dy**2 + ux_dz**2 + uy_dx**2 + uy_dy**2 + uy_dz**2 + uz_dx**2 + uz_dy**2 + uz_dz**2
                    end do
                end do
            end do

        case("FD_6th_central")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ux_dx = (a_FD6(-3)*u(ix-3,iy,iz,1) &
                              + a_FD6(-2)*u(ix-2,iy,iz,1) &
                              + a_FD6(-1)*u(ix-1,iy,iz,1) &
                              + a_FD6(0) *u(ix,iy,iz,1) &
                              + a_FD6(+1)*u(ix+1,iy,iz,1) &
                              + a_FD6(+2)*u(ix+2,iy,iz,1) &
                              + a_FD6(+3)*u(ix+3,iy,iz,1))*dx_inv

                        ux_dy = (a_FD6(-3)*u(ix,iy-3,iz,1) &
                              + a_FD6(-2)*u(ix,iy-2,iz,1) &
                              + a_FD6(-1)*u(ix,iy-1,iz,1) &
                              + a_FD6(0) *u(ix,iy,iz,1) &
                              + a_FD6(+1)*u(ix,iy+1,iz,1) &
                              + a_FD6(+2)*u(ix,iy+2,iz,1) &
                              + a_FD6(+3)*u(ix,iy+3,iz,1))*dy_inv

                        ux_dz = (a_FD6(-3)*u(ix,iy,iz-3,1) &
                              + a_FD6(-2)*u(ix,iy,iz-2,1) &
                              + a_FD6(-1)*u(ix,iy,iz-1,1) &
                              + a_FD6(0) *u(ix,iy,iz,1) &
                              + a_FD6(+1)*u(ix,iy,iz+1,1) &
                              + a_FD6(+2)*u(ix,iy,iz+2,1) &
                              + a_FD6(+3)*u(ix,iy,iz+3,1))*dz_inv

                        uy_dx = (a_FD6(-3)*u(ix-3,iy,iz,2) &
                              + a_FD6(-2)*u(ix-2,iy,iz,2) &
                              + a_FD6(-1)*u(ix-1,iy,iz,2) &
                              + a_FD6(0) *u(ix,iy,iz,2) &
                              + a_FD6(+1)*u(ix+1,iy,iz,2) &
                              + a_FD6(+2)*u(ix+2,iy,iz,2) &
                              + a_FD6(+3)*u(ix+3,iy,iz,2))*dx_inv

                        uy_dy = (a_FD6(-3)*u(ix,iy-3,iz,2) &
                              + a_FD6(-2)*u(ix,iy-2,iz,2) &
                              + a_FD6(-1)*u(ix,iy-1,iz,2) &
                              + a_FD6(0) *u(ix,iy,iz,2) &
                              + a_FD6(+1)*u(ix,iy+1,iz,2) &
                              + a_FD6(+2)*u(ix,iy+2,iz,2) &
                              + a_FD6(+3)*u(ix,iy+3,iz,2))*dy_inv

                        uy_dz = (a_FD6(-3)*u(ix,iy,iz-3,2) &
                              + a_FD6(-2)*u(ix,iy,iz-2,2) &
                              + a_FD6(-1)*u(ix,iy,iz-1,2) &
                              + a_FD6(0) *u(ix,iy,iz,2) &
                              + a_FD6(+1)*u(ix,iy,iz+1,2) &
                              + a_FD6(+2)*u(ix,iy,iz+2,2) &
                              + a_FD6(+3)*u(ix,iy,iz+3,2))*dz_inv

                        uz_dx = (a_FD6(-3)*u(ix-3,iy,iz,3) &
                              + a_FD6(-2)*u(ix-2,iy,iz,3) &
                              + a_FD6(-1)*u(ix-1,iy,iz,3) &
                              + a_FD6(0) *u(ix,iy,iz,3) &
                              + a_FD6(+1)*u(ix+1,iy,iz,3) &
                              + a_FD6(+2)*u(ix+2,iy,iz,3) &
                              + a_FD6(+3)*u(ix+3,iy,iz,3))*dx_inv

                        uz_dy = (a_FD6(-3)*u(ix,iy-3,iz,3) &
                              + a_FD6(-2)*u(ix,iy-2,iz,3) &
                              + a_FD6(-1)*u(ix,iy-1,iz,3) &
                              + a_FD6(0) *u(ix,iy,iz,3) &
                              + a_FD6(+1)*u(ix,iy+1,iz,3) &
                              + a_FD6(+2)*u(ix,iy+2,iz,3) &
                              + a_FD6(+3)*u(ix,iy+3,iz,3))*dy_inv

                        uz_dz = (a_FD6(-3)*u(ix,iy,iz-3,3) &
                              + a_FD6(-2)*u(ix,iy,iz-2,3) &
                              + a_FD6(-1)*u(ix,iy,iz-1,3) &
                              + a_FD6(0) *u(ix,iy,iz,3) &
                              + a_FD6(+1)*u(ix,iy,iz+1,3) &
                              + a_FD6(+2)*u(ix,iy,iz+2,3) &
                              + a_FD6(+3)*u(ix,iy,iz+3,3))*dz_inv


                        dissipation_rate = dissipation_rate + &
                            ux_dx**2 + ux_dy**2 + ux_dz**2 + uy_dx**2 + uy_dy**2 + uy_dz**2 + uz_dx**2 + uz_dy**2 + uz_dz**2
                    end do
                end do
            end do

        case("FD_4th_central_optimized")
            call abort(2362976, "dissipation rate with FD_4th_central_optimized not yet done.")

        case default
            call abort(1902201, "unknown order_discretization in ACM vorticity")
        end select

    endif

    dissipation_rate = dissipation_rate * product(dx(1:params_acm%dim))


end subroutine
