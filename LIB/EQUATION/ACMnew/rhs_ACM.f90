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

    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.
        ! performs initializations in the RHS module, such as resetting integrals

        ! Linear Forcing for HIT (Lundgren) requires us to know kinetic energy and dissipation
        ! rate at all times, so compute that, if we use the forcing.
        if (params_acm%use_HIT_linear_forcing) then
            params_acm%e_kin = 0.0_rk
            params_acm%enstrophy = 0.0_rk
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
                call abort(0409201933,"ACM fail: very very large values in state vector.")
            endif
        enddo

        ! Linear Forcing for HIT (Lundgren) requires us to know kinetic energy and dissipation
        ! rate at all times, so compute that, if we use the forcing.
        if (params_acm%use_HIT_linear_forcing) then
            ! vorticity work array
            if (.not. allocated(vor) ) allocate(vor(1:size(u,1), 1:size(u,2), 1:size(u,3), 1:3 ))
            ! to compute the current dissipation rate
            call vorticity_ACM_block(Bs, g, dx, u, vor)

            if (dim == 2) then
                dV = dx(1)*dx(2)
                dV2 = dV * 0.5_rk

                iz = 1
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        ! note dV2 contains the 0.5 from energy as well as the spacing
                        params_acm%e_kin = params_acm%e_kin + dv2*( u(ix,iy,iz,1)*u(ix,iy,iz,1) + u(ix,iy,iz,2)*u(ix,iy,iz,2) )

                        params_acm%enstrophy = params_acm%enstrophy + dv * vor(ix,iy,iz,1)*vor(ix,iy,iz,1)
                    enddo
                enddo
            else
                dV = dx(1)*dx(2)*dx(3)
                dV2 = dV * 0.5_rk

                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            ! note dV2 contains the 0.5 from energy as well as the spacing
                            params_acm%e_kin = params_acm%e_kin + dv2*( u(ix,iy,iz,1)*u(ix,iy,iz,1) + u(ix,iy,iz,2)*u(ix,iy,iz,2) &
                            + u(ix,iy,iz,3)*u(ix,iy,iz,3) )

                            params_acm%enstrophy = params_acm%enstrophy + dv*( vor(ix,iy,iz,1)*vor(ix,iy,iz,1) &
                            + vor(ix,iy,iz,2)*vor(ix,iy,iz,2) + vor(ix,iy,iz,3)*vor(ix,iy,iz,3) )
                        enddo
                    enddo
                enddo

            endif ! NOTE: MPI_SUM is perfomed in the post_stage.
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
        if (params_acm%use_HIT_linear_forcing) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
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
    character(len=cshort), intent(in)           :: order_discretization
    !> time
    real(kind=rk), intent(in)               :: time

    !> forcing term
    real(kind=rk), dimension(3) :: forcing
    !>
    real(kind=rk) :: dx_inv, dy_inv, dx2_inv, dy2_inv, c_0, nu, eps, eps_inv, gamma
    real(kind=rk) :: div_U, u_dx, u_dy, u_dxdx, u_dydy, v_dx, v_dy, v_dxdx, &
                     v_dydy, p_dx, p_dy, penalx, penaly, x, y, term_2, spo, p_dxdx, p_dydy, nu_p, &
                     u_dx4, v_dx4, u_dy4, v_dy4
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

    select case(order_discretization)
    case ("FD_2nd_central")
        !-----------------------------------------------------------------------
        ! 2nd order
        !-----------------------------------------------------------------------
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

                penalx = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                penaly = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                ! actual RHS. note mean flow forcing is just a constant and added at the end of the routine
                rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
            end do
        end do

    case("FD_4th_central_optimized")
        !-----------------------------------------------------------------------
        ! 4th order (Tam&web optimized scheme)
        !-----------------------------------------------------------------------
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! first derivatives of u, v, p
                ! Note: a(0) does NOT appear (it is zero...)
                u_dx = (a_TW4(-3)*phi(ix-3,iy,1) &
                      + a_TW4(-2)*phi(ix-2,iy,1) &
                      + a_TW4(-1)*phi(ix-1,iy,1) &
                      + a_TW4(+1)*phi(ix+1,iy,1) &
                      + a_TW4(+2)*phi(ix+2,iy,1) &
                      + a_TW4(+3)*phi(ix+3,iy,1))*dx_inv

                u_dy = (a_TW4(-3)*phi(ix,iy-3,1) &
                      + a_TW4(-2)*phi(ix,iy-2,1) &
                      + a_TW4(-1)*phi(ix,iy-1,1) &
                      + a_TW4(+1)*phi(ix,iy+1,1) &
                      + a_TW4(+2)*phi(ix,iy+2,1) &
                      + a_TW4(+3)*phi(ix,iy+3,1))*dy_inv

                v_dx = (a_TW4(-3)*phi(ix-3,iy,2) &
                      + a_TW4(-2)*phi(ix-2,iy,2) &
                      + a_TW4(-1)*phi(ix-1,iy,2) &
                      + a_TW4(+1)*phi(ix+1,iy,2) &
                      + a_TW4(+2)*phi(ix+2,iy,2) &
                      + a_TW4(+3)*phi(ix+3,iy,2))*dx_inv

                v_dy = (a_TW4(-3)*phi(ix,iy-3,2) &
                      + a_TW4(-2)*phi(ix,iy-2,2) &
                      + a_TW4(-1)*phi(ix,iy-1,2) &
                      + a_TW4(+1)*phi(ix,iy+1,2) &
                      + a_TW4(+2)*phi(ix,iy+2,2) &
                      + a_TW4(+3)*phi(ix,iy+3,2))*dy_inv

                p_dx = (a_TW4(-3)*phi(ix-3,iy,3) &
                      + a_TW4(-2)*phi(ix-2,iy,3) &
                      + a_TW4(-1)*phi(ix-1,iy,3) &
                      + a_TW4(+1)*phi(ix+1,iy,3) &
                      + a_TW4(+2)*phi(ix+2,iy,3) &
                      + a_TW4(+3)*phi(ix+3,iy,3))*dx_inv

                p_dy = (a_TW4(-3)*phi(ix,iy-3,3) &
                      + a_TW4(-2)*phi(ix,iy-2,3) &
                      + a_TW4(-1)*phi(ix,iy-1,3) &
                      + a_TW4(+1)*phi(ix,iy+1,3) &
                      + a_TW4(+2)*phi(ix,iy+2,3) &
                      + a_TW4(+3)*phi(ix,iy+3,3))*dy_inv

                ! second derivatives of u and v
                u_dxdx = (b_FD4(-2)*phi(ix-2,iy,1) &
                        + b_FD4(-1)*phi(ix-1,iy,1) &
                        + b_FD4(0)*phi(ix,iy,1) &
                        + b_FD4(+1)*phi(ix+1,iy,1) &
                        + b_FD4(+2)*phi(ix+2,iy,1))*dx2_inv

                u_dydy = (b_FD4(-2)*phi(ix,iy-2,1) &
                        + b_FD4(-1)*phi(ix,iy-1,1) &
                        + b_FD4(0)*phi(ix,iy,1) &
                        + b_FD4(+1)*phi(ix,iy+1,1) &
                        + b_FD4(+2)*phi(ix,iy+2,1))*dy2_inv

                v_dxdx = (b_FD4(-2)*phi(ix-2,iy,2) &
                        + b_FD4(-1)*phi(ix-1,iy,2) &
                        + b_FD4(0)*phi(ix,iy,2) &
                        + b_FD4(+1)*phi(ix+1,iy,2) &
                        + b_FD4(+2)*phi(ix+2,iy,2))*dx2_inv

                v_dydy = (b_FD4(-2)*phi(ix,iy-2,2) &
                        + b_FD4(-1)*phi(ix,iy-1,2) &
                        + b_FD4(0)*phi(ix,iy,2) &
                        + b_FD4(+1)*phi(ix,iy+1,2) &
                        + b_FD4(+2)*phi(ix,iy+2,2))*dy2_inv

                div_U = u_dx + v_dy

                penalx = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                penaly = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
            end do
        end do

    case("FD_4th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard)
        !-----------------------------------------------------------------------
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! first derivatives of u, v, p
                ! Note: a(0) does NOT appear (it is zero...)
                u_dx = (a_FD4(-2)*phi(ix-2,iy,1) &
                      + a_FD4(-1)*phi(ix-1,iy,1) &
                      + a_FD4(+1)*phi(ix+1,iy,1) &
                      + a_FD4(+2)*phi(ix+2,iy,1))*dx_inv

                u_dy = (a_FD4(-2)*phi(ix,iy-2,1) &
                      + a_FD4(-1)*phi(ix,iy-1,1) &
                      + a_FD4(+1)*phi(ix,iy+1,1) &
                      + a_FD4(+2)*phi(ix,iy+2,1))*dy_inv

                v_dx = (a_FD4(-2)*phi(ix-2,iy,2) &
                      + a_FD4(-1)*phi(ix-1,iy,2) &
                      + a_FD4(+1)*phi(ix+1,iy,2) &
                      + a_FD4(+2)*phi(ix+2,iy,2))*dx_inv

                v_dy = (a_FD4(-2)*phi(ix,iy-2,2) &
                      + a_FD4(-1)*phi(ix,iy-1,2) &
                      + a_FD4(+1)*phi(ix,iy+1,2) &
                      + a_FD4(+2)*phi(ix,iy+2,2))*dy_inv

                p_dx = (a_FD4(-2)*phi(ix-2,iy,3) &
                      + a_FD4(-1)*phi(ix-1,iy,3) &
                      + a_FD4(+1)*phi(ix+1,iy,3) &
                      + a_FD4(+2)*phi(ix+2,iy,3))*dx_inv

                p_dy = (a_FD4(-2)*phi(ix,iy-2,3) &
                      + a_FD4(-1)*phi(ix,iy-1,3) &
                      + a_FD4(+1)*phi(ix,iy+1,3) &
                      + a_FD4(+2)*phi(ix,iy+2,3))*dy_inv

              ! second derivatives of u and v
              u_dxdx = (b_FD4(-2)*phi(ix-2,iy,1) &
                      + b_FD4(-1)*phi(ix-1,iy,1) &
                      + b_FD4(0)*phi(ix,iy,1) &
                      + b_FD4(+1)*phi(ix+1,iy,1) &
                      + b_FD4(+2)*phi(ix+2,iy,1))*dx2_inv

              u_dydy = (b_FD4(-2)*phi(ix,iy-2,1) &
                      + b_FD4(-1)*phi(ix,iy-1,1) &
                      + b_FD4(0)*phi(ix,iy,1) &
                      + b_FD4(+1)*phi(ix,iy+1,1) &
                      + b_FD4(+2)*phi(ix,iy+2,1))*dy2_inv

              v_dxdx = (b_FD4(-2)*phi(ix-2,iy,2) &
                      + b_FD4(-1)*phi(ix-1,iy,2) &
                      + b_FD4(0)*phi(ix,iy,2) &
                      + b_FD4(+1)*phi(ix+1,iy,2) &
                      + b_FD4(+2)*phi(ix+2,iy,2))*dx2_inv

              v_dydy = (b_FD4(-2)*phi(ix,iy-2,2) &
                      + b_FD4(-1)*phi(ix,iy-1,2) &
                      + b_FD4(0)*phi(ix,iy,2) &
                      + b_FD4(+1)*phi(ix,iy+1,2) &
                      + b_FD4(+2)*phi(ix,iy+2,2))*dy2_inv

                div_U = u_dx + v_dy
! mask (chi) (:,:,1)
                ! usx (:,:,2)
                ! usy (:,:,3)
                ! usz (:,:,4)
                ! color (.;.5)
                ! sponge (:,:,6)
                penalx = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                penaly = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
            end do
        end do

    case("FD_6th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard)
        !-----------------------------------------------------------------------
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! first derivatives of u, v, p
                ! Note: a(0) does NOT appear (it is zero...)
                u_dx = (a_FD6(-3)*phi(ix-3,iy,1)&
                       +a_FD6(-2)*phi(ix-2,iy,1) &
                       +a_FD6(-1)*phi(ix-1,iy,1) &
                       +a_FD6(+1)*phi(ix+1,iy,1) &
                       +a_FD6(+2)*phi(ix+2,iy,1) &
                       +a_FD6(+3)*phi(ix+3,iy,1))*dx_inv

                u_dy = (a_FD6(-3)*phi(ix,iy-3,1)&
                       +a_FD6(-2)*phi(ix,iy-2,1) &
                       +a_FD6(-1)*phi(ix,iy-1,1) &
                       +a_FD6(+1)*phi(ix,iy+1,1) &
                       +a_FD6(+2)*phi(ix,iy+2,1) &
                       +a_FD6(+3)*phi(ix,iy+3,1))*dy_inv

                v_dx = (a_FD6(-3)*phi(ix-3,iy,2)&
                       +a_FD6(-2)*phi(ix-2,iy,2) &
                       +a_FD6(-1)*phi(ix-1,iy,2) &
                       +a_FD6(+1)*phi(ix+1,iy,2) &
                       +a_FD6(+2)*phi(ix+2,iy,2) &
                       +a_FD6(+3)*phi(ix+3,iy,2))*dx_inv

                v_dy = (a_FD6(-3)*phi(ix,iy-3,2)&
                       +a_FD6(-2)*phi(ix,iy-2,2) &
                       +a_FD6(-1)*phi(ix,iy-1,2) &
                       +a_FD6(+1)*phi(ix,iy+1,2) &
                       +a_FD6(+2)*phi(ix,iy+2,2) &
                       +a_FD6(+3)*phi(ix,iy+3,2))*dy_inv

                p_dx = (a_FD6(-3)*phi(ix-3,iy,3)&
                       +a_FD6(-2)*phi(ix-2,iy,3) &
                       +a_FD6(-1)*phi(ix-1,iy,3) &
                       +a_FD6(+1)*phi(ix+1,iy,3) &
                       +a_FD6(+2)*phi(ix+2,iy,3) &
                       +a_FD6(+3)*phi(ix+3,iy,3))*dx_inv

                p_dy = (a_FD6(-3)*phi(ix,iy-3,3)&
                       +a_FD6(-2)*phi(ix,iy-2,3) &
                       +a_FD6(-1)*phi(ix,iy-1,3) &
                       +a_FD6(+1)*phi(ix,iy+1,3) &
                       +a_FD6(+2)*phi(ix,iy+2,3) &
                       +a_FD6(+3)*phi(ix,iy+3,3))*dy_inv

                       ! second derivatives of u and v
                u_dxdx = (b_FD6(-3)*phi(ix-3,iy,1) &
                        + b_FD6(-2)*phi(ix-2,iy,1) &
                        + b_FD6(-1)*phi(ix-1,iy,1) &
                        + b_FD6( 0)*phi(ix,iy,1) &
                        + b_FD6(+1)*phi(ix+1,iy,1) &
                        + b_FD6(+2)*phi(ix+2,iy,1) &
                        + b_FD6(+3)*phi(ix+3,iy,1))*dx2_inv

                u_dydy = (b_FD6(-3)*phi(ix,iy-3,1) &
                        + b_FD6(-2)*phi(ix,iy-2,1) &
                        + b_FD6(-1)*phi(ix,iy-1,1) &
                        + b_FD6( 0)*phi(ix,iy,1) &
                        + b_FD6(+1)*phi(ix,iy+1,1) &
                        + b_FD6(+2)*phi(ix,iy+2,1) &
                        + b_FD6(+3)*phi(ix,iy+3,1))*dy2_inv

                v_dxdx = (b_FD6(-3)*phi(ix-3,iy,2) &
                        + b_FD6(-2)*phi(ix-2,iy,2) &
                        + b_FD6(-1)*phi(ix-1,iy,2) &
                        + b_FD6( 0)*phi(ix,iy,2) &
                        + b_FD6(+1)*phi(ix+1,iy,2) &
                        + b_FD6(+2)*phi(ix+2,iy,2) &
                        + b_FD6(+3)*phi(ix+3,iy,2))*dx2_inv

                v_dydy = (b_FD6(-3)*phi(ix,iy-3,2) &
                        + b_FD6(-2)*phi(ix,iy-2,2) &
                        + b_FD6(-1)*phi(ix,iy-1,2) &
                        + b_FD6( 0)*phi(ix,iy,2) &
                        + b_FD6(+1)*phi(ix,iy+1,2) &
                        + b_FD6(+2)*phi(ix,iy+2,2) &
                        + b_FD6(+3)*phi(ix,iy+3,2))*dy2_inv

                div_U = u_dx + v_dy

                penalx = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,1) -mask(ix,iy,2))
                penaly = -mask(ix,iy,1) * eps_inv * (phi(ix,iy,2) -mask(ix,iy,3))

                rhs(ix,iy,1) = -phi(ix,iy,1)*u_dx - phi(ix,iy,2)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
                rhs(ix,iy,2) = -phi(ix,iy,1)*v_dx - phi(ix,iy,2)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
                rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*phi(ix,iy,3)
            end do
        end do

    case default
        call abort(441166, "Discretization unkown "//trim(adjustl(order_discretization))//", you should go play outside. Its nice." )

    end select

    ! ---------------------------------------------------------------------------
    ! EXPERIMENTAL
    ! hyperviscosity, with 2nd order 4th derivative
    ! ---------------------------------------------------------------------------
    ! if (abs(nu_p) > 1.0e-13_rk) then
    !     do iy = g+1, Bs(2)+g
    !         do ix = g+1, Bs(1)+g
    !             u_dx4  = ( 1.0_rk*phi(ix-2,iy,1) + (-4.0_rk)*phi(ix-1,iy,1) + 6.0_rk*phi(ix,iy,1) &
    !             + (-4.0_rk)*phi(ix+1,iy,1) + 1.0_rk*phi(ix+2,iy,1))*dx2_inv*dx2_inv
    !             u_dy4  = ( 1.0_rk*phi(ix,iy-2,1) + (-4.0_rk)*phi(ix,iy-1,1) + 6.0_rk*phi(ix,iy,1) &
    !             + (-4.0_rk)*phi(ix,iy+1,1) + 1.0_rk*phi(ix,iy+2,1))*dy2_inv*dy2_inv
    !
    !             v_dx4  = ( 1.0_rk*phi(ix-2,iy,2) + (-4.0_rk)*phi(ix-1,iy,2) + 6.0_rk*phi(ix,iy,2) &
    !             + (-4.0_rk)*phi(ix+1,iy,2) + 1.0_rk*phi(ix+2,iy,2))*dx2_inv*dx2_inv
    !             v_dy4  = ( 1.0_rk*phi(ix,iy-2,2) + (-4.0_rk)*phi(ix,iy-1,2) + 6.0_rk*phi(ix,iy,2) &
    !             + (-4.0_rk)*phi(ix,iy+1,2) + 1.0_rk*phi(ix,iy+2,2))*dy2_inv*dy2_inv
    !
    !
    !             rhs(ix,iy,1) = rhs(ix,iy,1) - nu_p*(u_dx4+u_dy4)
    !             rhs(ix,iy,2) = rhs(ix,iy,2) - nu_p*(v_dx4+v_dy4)
    !         end do
    !     end do
    ! endif

    ! ---------------------------------------------------------------------------
    ! EXPERIMENTAL
    ! Pressure diffusion term, experimental. (Hence not integrated in the above loops,
    ! performance does not yet matter, only if this term turns out to be super useful)
    ! ---------------------------------------------------------------------------
    ! if (nu_p > 1.0e-13_rk) then
    !     select case(order_discretization)
    !     case ("FD_2nd_central")
    !         ! 2nd order
    !         do iy = g+1, Bs(2)+g
    !             do ix = g+1, Bs(1)+g
    !                 p_dxdx = (phi(ix-1,iy  ,3) -2.0_rk*phi(ix,iy,3) +phi(ix+1,iy  ,3))*dx2_inv
    !                 p_dydy = (phi(ix  ,iy-1,3) -2.0_rk*phi(ix,iy,3) +phi(ix  ,iy+1,3))*dy2_inv
    !
    !                 rhs(ix,iy,3) = rhs(ix,iy,3) + nu_p*(p_dxdx + p_dydy)
    !             end do
    !         end do
    !
    !     case("FD_4th_central_optimized","FD_4th_central")
    !         do iy = g+1, Bs(2)+g
    !             do ix = g+1, Bs(1)+g
    !                 p_dxdx = (b(-2)*phi(ix-2,iy,3) + b(-1)*phi(ix-1,iy,3) + b(0)*phi(ix,iy,3) &
    !                        +  b(+1)*phi(ix+1,iy,3) + b(+2)*phi(ix+2,iy,3))*dx2_inv
    !                 p_dydy = (b(-2)*phi(ix,iy-2,3) + b(-1)*phi(ix,iy-1,3) + b(0)*phi(ix,iy,3) &
    !                        +  b(+1)*phi(ix,iy+1,3) + b(+2)*phi(ix,iy+2,3))*dy2_inv
    !
    !                 rhs(ix,iy,3) = rhs(ix,iy,3) + nu_p*(p_dxdx + p_dydy)
    !             end do
    !         end do
    !
    !     case default
    !         call abort(2204041, "Discretization unkown "//order_discretization//", I ll walk into the light now." )
    !
    !     end select
    ! endif


    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    ! HACK
    if (.not. params_acm%geometry == "lamballais") then
        if (params_acm%use_sponge) then
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
    else
        ! special treatment lamballais, also for p to p_ref in the ring
        ! ring stored in sponge mask array
        ! Gautier, R., Biau, D., Lamballais, E.: A reference solution of the flow over a circular cylinder at Re = 40 , Computers & Fluids 75, 103–111, 2013 

        ! use same C_eta as object, not the sponge value
        eps_inv = 1.0_rk / params_acm%C_eta
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! mask(:,:,4) contains p_ref forcing in the ring
                rhs(ix,iy,3) = rhs(ix,iy,3) - mask(ix,iy,6)*(phi(ix,iy,3)-mask(ix,iy,4))*eps_inv
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
                     nu, eps, eps_inv, gamma, spo, A_forcing, G_gain, t_l_inf, e_kin_set
    !> derivatives
    real(kind=rk) :: div_U, u_dx, u_dy, u_dz, u_dxdx, u_dydy, u_dzdz, &
                     v_dx, v_dy, v_dz, v_dxdx, v_dydy, v_dzdz, &
                     w_dx, w_dy, w_dz, w_dxdx, w_dydy, w_dzdz, &
                     p_dx, p_dy, p_dz, penalx, penaly, penalz, u, v, w, p, chi
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
    eps         = params_acm%C_eta
    gamma       = params_acm%gamma_p

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dz_inv = 1.0_rk / dx(3)

    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)
    dz2_inv = 1.0_rk / (dx(3)**2)

    eps_inv = 1.0_rk / eps

    select case(order_discretization)
    case("FD_2nd_central")
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

    case("FD_4th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard scheme)
        !-----------------------------------------------------------------------
        ! Note: a(0) does NOT appear (it is zero...)
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    u_dx = (a_FD4(-2)*phi(ix-2,iy,iz,1) &
                           +a_FD4(-1)*phi(ix-1,iy,iz,1) &
                           +a_FD4(+1)*phi(ix+1,iy,iz,1) &
                           +a_FD4(+2)*phi(ix+2,iy,iz,1))*dx_inv

                    u_dy = (a_FD4(-2)*phi(ix,iy-2,iz,1) &
                           +a_FD4(-1)*phi(ix,iy-1,iz,1) &
                           +a_FD4(+1)*phi(ix,iy+1,iz,1) &
                           +a_FD4(+2)*phi(ix,iy+2,iz,1))*dy_inv

                    u_dz = (a_FD4(-2)*phi(ix,iy,iz-2,1) &
                           +a_FD4(-1)*phi(ix,iy,iz-1,1) &
                           +a_FD4(+1)*phi(ix,iy,iz+1,1) &
                           +a_FD4(+2)*phi(ix,iy,iz+2,1))*dz_inv

                    v_dx = (a_FD4(-2)*phi(ix-2,iy,iz,2) &
                           +a_FD4(-1)*phi(ix-1,iy,iz,2) &
                           +a_FD4(+1)*phi(ix+1,iy,iz,2) &
                           +a_FD4(+2)*phi(ix+2,iy,iz,2))*dx_inv

                    v_dy = (a_FD4(-2)*phi(ix,iy-2,iz,2) &
                           +a_FD4(-1)*phi(ix,iy-1,iz,2) &
                           +a_FD4(+1)*phi(ix,iy+1,iz,2) &
                           +a_FD4(+2)*phi(ix,iy+2,iz,2))*dy_inv

                    v_dz = (a_FD4(-2)*phi(ix,iy,iz-2,2) &
                           +a_FD4(-1)*phi(ix,iy,iz-1,2) &
                           +a_FD4(+1)*phi(ix,iy,iz+1,2) &
                           +a_FD4(+2)*phi(ix,iy,iz+2,2))*dz_inv

                    w_dx = (a_FD4(-2)*phi(ix-2,iy,iz,3) &
                           +a_FD4(-1)*phi(ix-1,iy,iz,3) &
                           +a_FD4(+1)*phi(ix+1,iy,iz,3) &
                           +a_FD4(+2)*phi(ix+2,iy,iz,3))*dx_inv

                    w_dy = (a_FD4(-2)*phi(ix,iy-2,iz,3) &
                           +a_FD4(-1)*phi(ix,iy-1,iz,3) &
                           +a_FD4(+1)*phi(ix,iy+1,iz,3) &
                           +a_FD4(+2)*phi(ix,iy+2,iz,3))*dy_inv

                    w_dz = (a_FD4(-2)*phi(ix,iy,iz-2,3) &
                           +a_FD4(-1)*phi(ix,iy,iz-1,3) &
                           +a_FD4(+1)*phi(ix,iy,iz+1,3) &
                           +a_FD4(+2)*phi(ix,iy,iz+2,3))*dz_inv

                    p_dx = (a_FD4(-2)*phi(ix-2,iy,iz,4) &
                           +a_FD4(-1)*phi(ix-1,iy,iz,4) &
                           +a_FD4(+1)*phi(ix+1,iy,iz,4) &
                           +a_FD4(+2)*phi(ix+2,iy,iz,4))*dx_inv

                    p_dy = (a_FD4(-2)*phi(ix,iy-2,iz,4) &
                           +a_FD4(-1)*phi(ix,iy-1,iz,4) &
                           +a_FD4(+1)*phi(ix,iy+1,iz,4) &
                           +a_FD4(+2)*phi(ix,iy+2,iz,4))*dy_inv

                    p_dz = (a_FD4(-2)*phi(ix,iy,iz-2,4) &
                           +a_FD4(-1)*phi(ix,iy,iz-1,4) &
                           +a_FD4(+1)*phi(ix,iy,iz+1,4) &
                           +a_FD4(+2)*phi(ix,iy,iz+2,4))*dz_inv

                    ! second derivatives of u, v and w
                    u_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,1) &
                            + b_FD4(-1)*phi(ix-1,iy,iz,1) &
                            + b_FD4(0)*phi(ix,iy,iz,1) &
                            + b_FD4(+1)*phi(ix+1,iy,iz,1) &
                            + b_FD4(+2)*phi(ix+2,iy,iz,1))*dx2_inv

                    u_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,1) &
                            + b_FD4(-1)*phi(ix,iy-1,iz,1) &
                            + b_FD4(0)*phi(ix,iy,iz,1) &
                            + b_FD4(+1)*phi(ix,iy+1,iz,1) &
                            + b_FD4(+2)*phi(ix,iy+2,iz,1))*dy2_inv

                    u_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,1) &
                            + b_FD4(-1)*phi(ix,iy,iz-1,1) &
                            + b_FD4(0)*phi(ix,iy,iz,1) &
                            + b_FD4(+1)*phi(ix,iy,iz+1,1) &
                            + b_FD4(+2)*phi(ix,iy,iz+2,1))*dz2_inv

                    v_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,2) &
                            + b_FD4(-1)*phi(ix-1,iy,iz,2) &
                            + b_FD4(0)*phi(ix,iy,iz,2) &
                            + b_FD4(+1)*phi(ix+1,iy,iz,2) &
                            + b_FD4(+2)*phi(ix+2,iy,iz,2))*dx2_inv

                    v_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,2) &
                            + b_FD4(-1)*phi(ix,iy-1,iz,2) &
                            + b_FD4(0)*phi(ix,iy,iz,2) &
                            + b_FD4(+1)*phi(ix,iy+1,iz,2) &
                            + b_FD4(+2)*phi(ix,iy+2,iz,2))*dy2_inv

                    v_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,2) &
                            + b_FD4(-1)*phi(ix,iy,iz-1,2) &
                            + b_FD4(0)*phi(ix,iy,iz,2) &
                            + b_FD4(+1)*phi(ix,iy,iz+1,2) &
                            + b_FD4(+2)*phi(ix,iy,iz+2,2))*dz2_inv

                    w_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,3) &
                            + b_FD4(-1)*phi(ix-1,iy,iz,3) &
                            + b_FD4(0)*phi(ix,iy,iz,3) &
                            + b_FD4(+1)*phi(ix+1,iy,iz,3) &
                            + b_FD4(+2)*phi(ix+2,iy,iz,3))*dx2_inv

                    w_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,3) &
                            + b_FD4(-1)*phi(ix,iy-1,iz,3) &
                            + b_FD4(0)*phi(ix,iy,iz,3) &
                            + b_FD4(+1)*phi(ix,iy+1,iz,3) &
                            + b_FD4(+2)*phi(ix,iy+2,iz,3))*dy2_inv

                    w_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,3) &
                            + b_FD4(-1)*phi(ix,iy,iz-1,3) &
                            + b_FD4(0)*phi(ix,iy,iz,3) &
                            + b_FD4(+1)*phi(ix,iy,iz+1,3) &
                            + b_FD4(+2)*phi(ix,iy,iz+2,3))*dz2_inv

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

    case("FD_6th_central")
        !-----------------------------------------------------------------------
        ! 4th order (standard scheme)
        !-----------------------------------------------------------------------
        ! Note: a(0) does NOT appear (it is zero...)
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    u_dx = (a_FD6(-3)*phi(ix-3,iy,iz,1) &
                           +a_FD6(-2)*phi(ix-2,iy,iz,1) &
                           +a_FD6(-1)*phi(ix-1,iy,iz,1) &
                           +a_FD6(+1)*phi(ix+1,iy,iz,1) &
                           +a_FD6(+2)*phi(ix+2,iy,iz,1) &
                           +a_FD6(+3)*phi(ix+3,iy,iz,1))*dx_inv

                    u_dy = (a_FD6(-3)*phi(ix,iy-3,iz,1) &
                           +a_FD6(-2)*phi(ix,iy-2,iz,1) &
                           +a_FD6(-1)*phi(ix,iy-1,iz,1) &
                           +a_FD6(+1)*phi(ix,iy+1,iz,1) &
                           +a_FD6(+2)*phi(ix,iy+2,iz,1) &
                           +a_FD6(+3)*phi(ix,iy+3,iz,1))*dy_inv

                    u_dz = (a_FD6(-3)*phi(ix,iy,iz-3,1) &
                           +a_FD6(-2)*phi(ix,iy,iz-2,1) &
                           +a_FD6(-1)*phi(ix,iy,iz-1,1) &
                           +a_FD6(+1)*phi(ix,iy,iz+1,1) &
                           +a_FD6(+2)*phi(ix,iy,iz+2,1) &
                           +a_FD6(+3)*phi(ix,iy,iz+3,1))*dz_inv

                    v_dx = (a_FD6(-3)*phi(ix-3,iy,iz,2) &
                           +a_FD6(-2)*phi(ix-2,iy,iz,2) &
                           +a_FD6(-1)*phi(ix-1,iy,iz,2) &
                           +a_FD6(+1)*phi(ix+1,iy,iz,2) &
                           +a_FD6(+2)*phi(ix+2,iy,iz,2) &
                           +a_FD6(+3)*phi(ix+3,iy,iz,2))*dx_inv

                    v_dy = (a_FD6(-3)*phi(ix,iy-3,iz,2) &
                           +a_FD6(-2)*phi(ix,iy-2,iz,2) &
                           +a_FD6(-1)*phi(ix,iy-1,iz,2) &
                           +a_FD6(+1)*phi(ix,iy+1,iz,2) &
                           +a_FD6(+2)*phi(ix,iy+2,iz,2) &
                           +a_FD6(+3)*phi(ix,iy+3,iz,2))*dy_inv

                    v_dz = (a_FD6(-3)*phi(ix,iy,iz-3,2) &
                           +a_FD6(-2)*phi(ix,iy,iz-2,2) &
                           +a_FD6(-1)*phi(ix,iy,iz-1,2) &
                           +a_FD6(+1)*phi(ix,iy,iz+1,2) &
                           +a_FD6(+2)*phi(ix,iy,iz+2,2) &
                           +a_FD6(+3)*phi(ix,iy,iz+3,2))*dz_inv

                    w_dx = (a_FD6(-3)*phi(ix-3,iy,iz,3) &
                           +a_FD6(-2)*phi(ix-2,iy,iz,3) &
                           +a_FD6(-1)*phi(ix-1,iy,iz,3) &
                           +a_FD6(+1)*phi(ix+1,iy,iz,3) &
                           +a_FD6(+2)*phi(ix+2,iy,iz,3) &
                           +a_FD6(+3)*phi(ix+3,iy,iz,3))*dx_inv

                    w_dy = (a_FD6(-3)*phi(ix,iy-3,iz,3) &
                           +a_FD6(-2)*phi(ix,iy-2,iz,3) &
                           +a_FD6(-1)*phi(ix,iy-1,iz,3) &
                           +a_FD6(+1)*phi(ix,iy+1,iz,3) &
                           +a_FD6(+2)*phi(ix,iy+2,iz,3) &
                           +a_FD6(+3)*phi(ix,iy+3,iz,3))*dy_inv

                    w_dz = (a_FD6(-3)*phi(ix,iy,iz-3,3) &
                           +a_FD6(-2)*phi(ix,iy,iz-2,3) &
                           +a_FD6(-1)*phi(ix,iy,iz-1,3) &
                           +a_FD6(+1)*phi(ix,iy,iz+1,3) &
                           +a_FD6(+2)*phi(ix,iy,iz+2,3) &
                           +a_FD6(+3)*phi(ix,iy,iz+3,3))*dz_inv

                    p_dx = (a_FD6(-3)*phi(ix-3,iy,iz,4) &
                           +a_FD6(-2)*phi(ix-2,iy,iz,4) &
                           +a_FD6(-1)*phi(ix-1,iy,iz,4) &
                           +a_FD6(+1)*phi(ix+1,iy,iz,4) &
                           +a_FD6(+2)*phi(ix+2,iy,iz,4) &
                           +a_FD6(+3)*phi(ix+3,iy,iz,4))*dx_inv

                    p_dy = (a_FD6(-3)*phi(ix,iy-3,iz,4) &
                           +a_FD6(-2)*phi(ix,iy-2,iz,4) &
                           +a_FD6(-1)*phi(ix,iy-1,iz,4) &
                           +a_FD6(+1)*phi(ix,iy+1,iz,4) &
                           +a_FD6(+2)*phi(ix,iy+2,iz,4) &
                           +a_FD6(+3)*phi(ix,iy+3,iz,4))*dy_inv

                    p_dz = (a_FD6(-3)*phi(ix,iy,iz-3,4) &
                           +a_FD6(-2)*phi(ix,iy,iz-2,4) &
                           +a_FD6(-1)*phi(ix,iy,iz-1,4) &
                           +a_FD6(+1)*phi(ix,iy,iz+1,4) &
                           +a_FD6(+2)*phi(ix,iy,iz+2,4) &
                           +a_FD6(+3)*phi(ix,iy,iz+3,4))*dz_inv

                    ! second derivatives of u, v and w
                    u_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,1) &
                            + b_FD6(-2)*phi(ix-2,iy,iz,1) &
                            + b_FD6(-1)*phi(ix-1,iy,iz,1) &
                            + b_FD6(0)*phi(ix,iy,iz,1) &
                            + b_FD6(+1)*phi(ix+1,iy,iz,1) &
                            + b_FD6(+2)*phi(ix+2,iy,iz,1) &
                            + b_FD6(+3)*phi(ix+3,iy,iz,1))*dx2_inv

                    u_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,1) &
                            + b_FD6(-2)*phi(ix,iy-2,iz,1) &
                            + b_FD6(-1)*phi(ix,iy-1,iz,1) &
                            + b_FD6(0)*phi(ix,iy,iz,1) &
                            + b_FD6(+1)*phi(ix,iy+1,iz,1) &
                            + b_FD6(+2)*phi(ix,iy+2,iz,1) &
                            + b_FD6(+3)*phi(ix,iy+3,iz,1))*dy2_inv

                    u_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,1) &
                            + b_FD6(-2)*phi(ix,iy,iz-2,1) &
                            + b_FD6(-1)*phi(ix,iy,iz-1,1) &
                            + b_FD6(0)*phi(ix,iy,iz,1) &
                            + b_FD6(+1)*phi(ix,iy,iz+1,1) &
                            + b_FD6(+2)*phi(ix,iy,iz+2,1) &
                            + b_FD6(+3)*phi(ix,iy,iz+3,1))*dz2_inv

                    v_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,2) &
                            + b_FD6(-2)*phi(ix-2,iy,iz,2) &
                            + b_FD6(-1)*phi(ix-1,iy,iz,2) &
                            + b_FD6(0)*phi(ix,iy,iz,2) &
                            + b_FD6(+1)*phi(ix+1,iy,iz,2) &
                            + b_FD6(+2)*phi(ix+2,iy,iz,2) &
                            + b_FD6(+3)*phi(ix+3,iy,iz,2))*dx2_inv

                    v_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,2) &
                            + b_FD6(-2)*phi(ix,iy-2,iz,2) &
                            + b_FD6(-1)*phi(ix,iy-1,iz,2) &
                            + b_FD6(0)*phi(ix,iy,iz,2) &
                            + b_FD6(+1)*phi(ix,iy+1,iz,2) &
                            + b_FD6(+2)*phi(ix,iy+2,iz,2) &
                            + b_FD6(+3)*phi(ix,iy+3,iz,2))*dy2_inv

                    v_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,2) &
                            + b_FD6(-2)*phi(ix,iy,iz-2,2) &
                            + b_FD6(-1)*phi(ix,iy,iz-1,2) &
                            + b_FD6(0)*phi(ix,iy,iz,2) &
                            + b_FD6(+1)*phi(ix,iy,iz+1,2) &
                            + b_FD6(+2)*phi(ix,iy,iz+2,2) &
                            + b_FD6(+3)*phi(ix,iy,iz+3,2))*dz2_inv

                    w_dxdx = (b_FD6(-3)*phi(ix-3,iy,iz,3) &
                            + b_FD6(-2)*phi(ix-2,iy,iz,3) &
                            + b_FD6(-1)*phi(ix-1,iy,iz,3) &
                            + b_FD6(0)*phi(ix,iy,iz,3) &
                            + b_FD6(+1)*phi(ix+1,iy,iz,3) &
                            + b_FD6(+2)*phi(ix+2,iy,iz,3) &
                            + b_FD6(+3)*phi(ix+3,iy,iz,3))*dx2_inv

                    w_dydy = (b_FD6(-3)*phi(ix,iy-3,iz,3) &
                            + b_FD6(-2)*phi(ix,iy-2,iz,3) &
                            + b_FD6(-1)*phi(ix,iy-1,iz,3) &
                            + b_FD6(0)*phi(ix,iy,iz,3) &
                            + b_FD6(+1)*phi(ix,iy+1,iz,3) &
                            + b_FD6(+2)*phi(ix,iy+2,iz,3) &
                            + b_FD6(+3)*phi(ix,iy+3,iz,3))*dy2_inv

                    w_dzdz = (b_FD6(-3)*phi(ix,iy,iz-3,3) &
                            + b_FD6(-2)*phi(ix,iy,iz-2,3) &
                            + b_FD6(-1)*phi(ix,iy,iz-1,3) &
                            + b_FD6(0)*phi(ix,iy,iz,3) &
                            + b_FD6(+1)*phi(ix,iy,iz+1,3) &
                            + b_FD6(+2)*phi(ix,iy,iz+2,3) &
                            + b_FD6(+3)*phi(ix,iy,iz+3,3))*dz2_inv

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

    case("FD_4th_central_optimized")
        !-----------------------------------------------------------------------
        ! 4th order (tam&web optimized scheme)
        !-----------------------------------------------------------------------
        ! Note: a(0) does NOT appear (it is zero...)
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! first derivatives of u, v, p
                    u_dx = (a_TW4(-3)*phi(ix-3,iy,iz,1) &
                           +a_TW4(-2)*phi(ix-2,iy,iz,1) &
                           +a_TW4(-1)*phi(ix-1,iy,iz,1) &
                           +a_TW4(+1)*phi(ix+1,iy,iz,1) &
                           +a_TW4(+2)*phi(ix+2,iy,iz,1) &
                           +a_TW4(+3)*phi(ix+3,iy,iz,1))*dx_inv

                    u_dy = (a_TW4(-3)*phi(ix,iy-3,iz,1) &
                           +a_TW4(-2)*phi(ix,iy-2,iz,1) &
                           +a_TW4(-1)*phi(ix,iy-1,iz,1) &
                           +a_TW4(+1)*phi(ix,iy+1,iz,1) &
                           +a_TW4(+2)*phi(ix,iy+2,iz,1) &
                           +a_TW4(+3)*phi(ix,iy+3,iz,1))*dy_inv

                    u_dz = (a_TW4(-3)*phi(ix,iy,iz-3,1) &
                           +a_TW4(-2)*phi(ix,iy,iz-2,1) &
                           +a_TW4(-1)*phi(ix,iy,iz-1,1) &
                           +a_TW4(+1)*phi(ix,iy,iz+1,1) &
                           +a_TW4(+2)*phi(ix,iy,iz+2,1) &
                           +a_TW4(+3)*phi(ix,iy,iz+3,1))*dz_inv

                    v_dx = (a_TW4(-3)*phi(ix-3,iy,iz,2) &
                           +a_TW4(-2)*phi(ix-2,iy,iz,2) &
                           +a_TW4(-1)*phi(ix-1,iy,iz,2) &
                           +a_TW4(+1)*phi(ix+1,iy,iz,2) &
                           +a_TW4(+2)*phi(ix+2,iy,iz,2) &
                           +a_TW4(+3)*phi(ix+3,iy,iz,2))*dx_inv

                    v_dy = (a_TW4(-3)*phi(ix,iy-3,iz,2) &
                           +a_TW4(-2)*phi(ix,iy-2,iz,2) &
                           +a_TW4(-1)*phi(ix,iy-1,iz,2) &
                           +a_TW4(+1)*phi(ix,iy+1,iz,2) &
                           +a_TW4(+2)*phi(ix,iy+2,iz,2) &
                           +a_TW4(+3)*phi(ix,iy+3,iz,2))*dy_inv

                    v_dz = (a_TW4(-3)*phi(ix,iy,iz-3,2) &
                           +a_TW4(-2)*phi(ix,iy,iz-2,2) &
                           +a_TW4(-1)*phi(ix,iy,iz-1,2) &
                           +a_TW4(+1)*phi(ix,iy,iz+1,2) &
                           +a_TW4(+2)*phi(ix,iy,iz+2,2) &
                           +a_TW4(+3)*phi(ix,iy,iz+3,2))*dz_inv

                    w_dx = (a_TW4(-3)*phi(ix-3,iy,iz,3) &
                           +a_TW4(-2)*phi(ix-2,iy,iz,3) &
                           +a_TW4(-1)*phi(ix-1,iy,iz,3) &
                           +a_TW4(+1)*phi(ix+1,iy,iz,3) &
                           +a_TW4(+2)*phi(ix+2,iy,iz,3) &
                           +a_TW4(+3)*phi(ix+3,iy,iz,3))*dx_inv

                    w_dy = (a_TW4(-3)*phi(ix,iy-3,iz,3) &
                           +a_TW4(-2)*phi(ix,iy-2,iz,3) &
                           +a_TW4(-1)*phi(ix,iy-1,iz,3) &
                           +a_TW4(+1)*phi(ix,iy+1,iz,3) &
                           +a_TW4(+2)*phi(ix,iy+2,iz,3) &
                           +a_TW4(+3)*phi(ix,iy+3,iz,3))*dy_inv

                    w_dz = (a_TW4(-3)*phi(ix,iy,iz-3,3) &
                           +a_TW4(-2)*phi(ix,iy,iz-2,3) &
                           +a_TW4(-1)*phi(ix,iy,iz-1,3) &
                           +a_TW4(+1)*phi(ix,iy,iz+1,3) &
                           +a_TW4(+2)*phi(ix,iy,iz+2,3) &
                           +a_TW4(+3)*phi(ix,iy,iz+3,3))*dz_inv

                    p_dx = (a_TW4(-3)*phi(ix-3,iy,iz,4) &
                           +a_TW4(-2)*phi(ix-2,iy,iz,4) &
                           +a_TW4(-1)*phi(ix-1,iy,iz,4) &
                           +a_TW4(+1)*phi(ix+1,iy,iz,4) &
                           +a_TW4(+2)*phi(ix+2,iy,iz,4) &
                           +a_TW4(+3)*phi(ix+3,iy,iz,4))*dx_inv

                    p_dy = (a_TW4(-3)*phi(ix,iy-3,iz,4) &
                           +a_TW4(-2)*phi(ix,iy-2,iz,4) &
                           +a_TW4(-1)*phi(ix,iy-1,iz,4) &
                           +a_TW4(+1)*phi(ix,iy+1,iz,4) &
                           +a_TW4(+2)*phi(ix,iy+2,iz,4) &
                           +a_TW4(+3)*phi(ix,iy+3,iz,4))*dy_inv

                    p_dz = (a_TW4(-3)*phi(ix,iy,iz-3,4) &
                           +a_TW4(-2)*phi(ix,iy,iz-2,4) &
                           +a_TW4(-1)*phi(ix,iy,iz-1,4) &
                           +a_TW4(+1)*phi(ix,iy,iz+1,4) &
                           +a_TW4(+2)*phi(ix,iy,iz+2,4) &
                           +a_TW4(+3)*phi(ix,iy,iz+3,4))*dz_inv


                   ! second derivatives of u, v and w
                   u_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,1) &
                           + b_FD4(-1)*phi(ix-1,iy,iz,1) &
                           + b_FD4(0)*phi(ix,iy,iz,1) &
                           + b_FD4(+1)*phi(ix+1,iy,iz,1) &
                           + b_FD4(+2)*phi(ix+2,iy,iz,1))*dx2_inv

                   u_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,1) &
                           + b_FD4(-1)*phi(ix,iy-1,iz,1) &
                           + b_FD4(0)*phi(ix,iy,iz,1) &
                           + b_FD4(+1)*phi(ix,iy+1,iz,1) &
                           + b_FD4(+2)*phi(ix,iy+2,iz,1))*dy2_inv

                   u_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,1) &
                           + b_FD4(-1)*phi(ix,iy,iz-1,1) &
                           + b_FD4(0)*phi(ix,iy,iz,1) &
                           + b_FD4(+1)*phi(ix,iy,iz+1,1) &
                           + b_FD4(+2)*phi(ix,iy,iz+2,1))*dz2_inv

                   v_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,2) &
                           + b_FD4(-1)*phi(ix-1,iy,iz,2) &
                           + b_FD4(0)*phi(ix,iy,iz,2) &
                           + b_FD4(+1)*phi(ix+1,iy,iz,2) &
                           + b_FD4(+2)*phi(ix+2,iy,iz,2))*dx2_inv

                   v_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,2) &
                           + b_FD4(-1)*phi(ix,iy-1,iz,2) &
                           + b_FD4(0)*phi(ix,iy,iz,2) &
                           + b_FD4(+1)*phi(ix,iy+1,iz,2) &
                           + b_FD4(+2)*phi(ix,iy+2,iz,2))*dy2_inv

                   v_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,2) &
                           + b_FD4(-1)*phi(ix,iy,iz-1,2) &
                           + b_FD4(0)*phi(ix,iy,iz,2) &
                           + b_FD4(+1)*phi(ix,iy,iz+1,2) &
                           + b_FD4(+2)*phi(ix,iy,iz+2,2))*dz2_inv

                   w_dxdx = (b_FD4(-2)*phi(ix-2,iy,iz,3) &
                           + b_FD4(-1)*phi(ix-1,iy,iz,3) &
                           + b_FD4(0)*phi(ix,iy,iz,3) &
                           + b_FD4(+1)*phi(ix+1,iy,iz,3) &
                           + b_FD4(+2)*phi(ix+2,iy,iz,3))*dx2_inv

                   w_dydy = (b_FD4(-2)*phi(ix,iy-2,iz,3) &
                           + b_FD4(-1)*phi(ix,iy-1,iz,3) &
                           + b_FD4(0)*phi(ix,iy,iz,3) &
                           + b_FD4(+1)*phi(ix,iy+1,iz,3) &
                           + b_FD4(+2)*phi(ix,iy+2,iz,3))*dy2_inv

                   w_dzdz = (b_FD4(-2)*phi(ix,iy,iz-2,3) &
                           + b_FD4(-1)*phi(ix,iy,iz-1,3) &
                           + b_FD4(0)*phi(ix,iy,iz,3) &
                           + b_FD4(+1)*phi(ix,iy,iz+1,3) &
                           + b_FD4(+2)*phi(ix,iy,iz+2,3))*dz2_inv

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

    case default
        call abort(441167, "3d Discretization unkown "//order_discretization//", I ll walk into the light now." )

    end select

    ! --------------------------------------------------------------------------
    ! HIT linear forcing
    ! --------------------------------------------------------------------------
    ! if (params_acm%use_HIT_linear_forcing) then
    !     G_gain = 100.0_rk
    !     e_kin_set = 0.5_rk * (product(params_acm%domain_size))
    !     t_l_inf = 1.0_rk ! (2.0_rk/3.0_rk) * e_kin_set /
    !     A_forcing = (params_acm%dissipation - G_gain * (params_acm%e_kin - e_kin_set) / t_l_inf) / (2.0*params_acm%e_kin)
    !
    !     rhs(:,:,:,1) = rhs(:,:,:,1) + A_forcing*phi(:,:,:,1)
    !     rhs(:,:,:,2) = rhs(:,:,:,2) + A_forcing*phi(:,:,:,2)
    !     rhs(:,:,:,3) = rhs(:,:,:,3) + A_forcing*phi(:,:,:,3)
    ! endif

    ! --------------------------------------------------------------------------
    ! sponge term.
    ! --------------------------------------------------------------------------
    if (params_acm%use_sponge) then
        ! avoid division by multiplying with inverse
        eps_inv = 1.0_rk / params_acm%C_sponge

        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! NOTE: the sponge term acts, if active, on ALL components, ux,uy,p
                    ! which is different from the penalization term, which acts only on ux,uy and not p
                    ! NOTE: sponge mask set in hvy_mask
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




subroutine vorticity_ACM_block(Bs, g, dx, u, vor)
    implicit none

    !> grid parameter
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> spacing of the block
    real(kind=rk), dimension(3), intent(in) :: dx
    !> datafields
    real(kind=rk), intent(inout) :: u(:,:,:,:)
    real(kind=rk), intent(inout) :: vor(:,:,:,:)

    !> derivatives
    real(kind=rk) :: u_dy, u_dz, v_dx, v_dz, w_dx, w_dy
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

    if ( params_acm%dim == 2) then
        !-----------------------------------------------------------------------
        ! 2D, vorticity is a scalar
        !-----------------------------------------------------------------------
        dx_inv = 1.0_rk / dx(1)
        dy_inv = 1.0_rk / dx(2)

        select case (params_acm%discretization)
        case("FD_2nd_central")
            iz = 1
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dy = (u(ix,iy+1,iz,1)-u(ix,iy-1,iz,1))*dy_inv*0.5_rk
                    v_dx = (u(ix+1,iy,iz,2)-u(ix-1,iy,iz,2))*dx_inv*0.5_rk

                    vor(ix,iy,iz,1) = v_dx - u_dy
                end do
            end do

        case("FD_4th_central")
            iz = 1
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dy = (a_FD4(-2)*u(ix,iy-2,iz,1) &
                           +a_FD4(-1)*u(ix,iy-1,iz,1) &
                           +a_FD4(+1)*u(ix,iy+1,iz,1) &
                           +a_FD4(+2)*u(ix,iy+2,iz,1))*dy_inv

                    v_dx = (a_FD4(-2)*u(ix-2,iy,iz,2) &
                           +a_FD4(-1)*u(ix-1,iy,iz,2) &
                           +a_FD4(+1)*u(ix+1,iy,iz,2) &
                           +a_FD4(+2)*u(ix+2,iy,iz,2))*dx_inv

                    vor(ix,iy,iz,1) = v_dx - u_dy
                end do
            end do

        case("FD_6th_central")
            iz = 1
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dy = (a_FD6(-3)*u(ix,iy-3,iz,1) &
                           +a_FD6(-2)*u(ix,iy-2,iz,1) &
                           +a_FD6(-1)*u(ix,iy-1,iz,1) &
                           +a_FD6(+1)*u(ix,iy+1,iz,1) &
                           +a_FD6(+2)*u(ix,iy+2,iz,1) &
                           +a_FD6(+3)*u(ix,iy+3,iz,1))*dy_inv

                    v_dx = (a_FD6(-3)*u(ix-3,iy,iz,2) &
                           +a_FD6(-2)*u(ix-2,iy,iz,2) &
                           +a_FD6(-1)*u(ix-1,iy,iz,2) &
                           +a_FD6(+1)*u(ix+1,iy,iz,2) &
                           +a_FD6(+2)*u(ix+2,iy,iz,2) &
                           +a_FD6(+3)*u(ix+3,iy,iz,2))*dx_inv

                    vor(ix,iy,iz,1) = v_dx - u_dy
                end do
            end do

        case("FD_4th_central_optimized")
            iz = 1
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dy = (a_TW4(-3)*u(ix,iy-3,iz,1) &
                          + a_TW4(-2)*u(ix,iy-2,iz,1) &
                          + a_TW4(-1)*u(ix,iy-1,iz,1) &
                          + a_TW4(0) *u(ix,iy,iz,1) &
                          + a_TW4(+1)*u(ix,iy+1,iz,1) &
                          + a_TW4(+2)*u(ix,iy+2,iz,1) &
                          + a_TW4(+3)*u(ix,iy+3,iz,1))*dy_inv

                    v_dx = (a_TW4(-3)*u(ix-3,iy,iz,2) &
                          + a_TW4(-2)*u(ix-2,iy,iz,2) &
                          + a_TW4(-1)*u(ix-1,iy,iz,2) &
                          + a_TW4(0) *u(ix,iy,iz,2) &
                          + a_TW4(+1)*u(ix+1,iy,iz,2) &
                          + a_TW4(+2)*u(ix+2,iy,iz,2) &
                          + a_TW4(+3)*u(ix+3,iy,iz,2))*dx_inv

                    vor(ix,iy,iz,1) = v_dx - u_dy
                end do
            end do

        case default
            call abort(1902201, "unknown order_discretization in ACM vorticity")
        end select
    else
        !-----------------------------------------------------------------------
        ! 3D, vorticity is a vector
        !-----------------------------------------------------------------------
        dx_inv = 1.0_rk / dx(1)
        dy_inv = 1.0_rk / dx(2)
        dz_inv = 1.0_rk / dx(3)

        select case(params_acm%discretization)
        case("FD_2nd_central")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        u_dy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dy_inv*0.5_rk
                        u_dz = (u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1))*dz_inv*0.5_rk
                        v_dx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dx_inv*0.5_rk
                        v_dz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dz_inv*0.5_rk
                        w_dx = (u(ix+1,iy,iz,3) - u(ix-1,iy,iz,3))*dx_inv*0.5_rk
                        w_dy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dy_inv*0.5_rk

                        vor(ix,iy,iz,1) = w_dy - v_dz
                        vor(ix,iy,iz,2) = u_dz - w_dx
                        vor(ix,iy,iz,3) = v_dx - u_dy
                    end do
                end do
            end do

        case("FD_4th_central")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        u_dy = (a_FD4(-2)*u(ix,iy-2,iz,1) &
                              + a_FD4(-1)*u(ix,iy-1,iz,1) &
                              + a_FD4(0) *u(ix,iy,iz,1) &
                              + a_FD4(+1)*u(ix,iy+1,iz,1) &
                              + a_FD4(+2)*u(ix,iy+2,iz,1) )*dy_inv

                        u_dz = (a_FD4(-2)*u(ix,iy,iz-2,1) &
                              + a_FD4(-1)*u(ix,iy,iz-1,1) &
                              + a_FD4(0) *u(ix,iy,iz,1) &
                              + a_FD4(+1)*u(ix,iy,iz+1,1) &
                              + a_FD4(+2)*u(ix,iy,iz+2,1) )*dz_inv

                        v_dx = (a_FD4(-2)*u(ix-2,iy,iz,2)  &
                              + a_FD4(-1)*u(ix-1,iy,iz,2) &
                              + a_FD4(0) *u(ix,iy,iz,2) &
                              + a_FD4(+1)*u(ix+1,iy,iz,2) &
                              + a_FD4(+2)*u(ix+2,iy,iz,2) )*dx_inv

                        v_dz = (a_FD4(-2)*u(ix,iy,iz-2,2) &
                              + a_FD4(-1)*u(ix,iy,iz-1,2) &
                              + a_FD4(0) *u(ix,iy,iz,2) &
                              + a_FD4(+1)*u(ix,iy,iz+1,2) &
                              + a_FD4(+2)*u(ix,iy,iz+2,2) )*dz_inv

                        w_dx = (a_FD4(-2)*u(ix-2,iy,iz,3) &
                              + a_FD4(-1)*u(ix-1,iy,iz,3) &
                              + a_FD4(0) *u(ix,iy,iz,3) &
                              + a_FD4(+1)*u(ix+1,iy,iz,3) &
                              + a_FD4(+2)*u(ix+2,iy,iz,3) )*dx_inv

                        w_dy = (a_FD4(-2)*u(ix,iy-2,iz,3) &
                              + a_FD4(-1)*u(ix,iy-1,iz,3) &
                              + a_FD4(0) *u(ix,iy,iz,3) &
                              + a_FD4(+1)*u(ix,iy+1,iz,3) &
                              + a_FD4(+2)*u(ix,iy+2,iz,3) )*dy_inv

                        vor(ix,iy,iz,1) = w_dy - v_dz
                        vor(ix,iy,iz,2) = u_dz - w_dx
                        vor(ix,iy,iz,3) = v_dx - u_dy
                    end do
                end do
            end do

        case("FD_6th_central")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        u_dy = (a_FD6(-3)*u(ix,iy-3,iz,1) &
                              + a_FD6(-2)*u(ix,iy-2,iz,1) &
                              + a_FD6(-1)*u(ix,iy-1,iz,1) &
                              + a_FD6(0) *u(ix,iy,iz,1) &
                              + a_FD6(+1)*u(ix,iy+1,iz,1) &
                              + a_FD6(+2)*u(ix,iy+2,iz,1) &
                              + a_FD6(+3)*u(ix,iy+3,iz,1))*dy_inv

                        u_dz = (a_FD6(-3)*u(ix,iy,iz-3,1) &
                              + a_FD6(-2)*u(ix,iy,iz-2,1) &
                              + a_FD6(-1)*u(ix,iy,iz-1,1) &
                              + a_FD6(0) *u(ix,iy,iz,1) &
                              + a_FD6(+1)*u(ix,iy,iz+1,1) &
                              + a_FD6(+2)*u(ix,iy,iz+2,1) &
                              + a_FD6(+3)*u(ix,iy,iz+3,1))*dz_inv

                        v_dx = (a_FD6(-3)*u(ix-3,iy,iz,2) &
                              + a_FD6(-2)*u(ix-2,iy,iz,2)  &
                              + a_FD6(-1)*u(ix-1,iy,iz,2) &
                              + a_FD6(0) *u(ix,iy,iz,2) &
                              + a_FD6(+1)*u(ix+1,iy,iz,2) &
                              + a_FD6(+2)*u(ix+2,iy,iz,2) &
                              + a_FD6(+3)*u(ix+3,iy,iz,2))*dx_inv

                        v_dz = (a_FD6(-3)*u(ix,iy,iz-3,2) &
                              + a_FD6(-2)*u(ix,iy,iz-2,2) &
                              + a_FD6(-1)*u(ix,iy,iz-1,2) &
                              + a_FD6(0) *u(ix,iy,iz,2) &
                              + a_FD6(+1)*u(ix,iy,iz+1,2) &
                              + a_FD6(+2)*u(ix,iy,iz+2,2) &
                              + a_FD6(+3)*u(ix,iy,iz+3,2))*dz_inv

                        w_dx = (a_FD6(-3)*u(ix-3,iy,iz,3) &
                              + a_FD6(-2)*u(ix-2,iy,iz,3) &
                              + a_FD6(-1)*u(ix-1,iy,iz,3) &
                              + a_FD6(0) *u(ix,iy,iz,3) &
                              + a_FD6(+1)*u(ix+1,iy,iz,3) &
                              + a_FD6(+2)*u(ix+2,iy,iz,3) &
                              + a_FD6(+3)*u(ix+3,iy,iz,3))*dx_inv

                        w_dy = (a_FD6(-3)*u(ix,iy-3,iz,3) &
                              + a_FD6(-2)*u(ix,iy-2,iz,3) &
                              + a_FD6(-1)*u(ix,iy-1,iz,3) &
                              + a_FD6(0) *u(ix,iy,iz,3) &
                              + a_FD6(+1)*u(ix,iy+1,iz,3) &
                              + a_FD6(+2)*u(ix,iy+2,iz,3) &
                              + a_FD6(+3)*u(ix,iy+3,iz,3))*dy_inv

                        vor(ix,iy,iz,1) = w_dy - v_dz
                        vor(ix,iy,iz,2) = u_dz - w_dx
                        vor(ix,iy,iz,3) = v_dx - u_dy
                    end do
                end do
            end do

        case("FD_4th_central_optimized")
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        u_dy = (a_TW4(-3)*u(ix,iy-3,iz,1) &
                              + a_TW4(-2)*u(ix,iy-2,iz,1) &
                              + a_TW4(-1)*u(ix,iy-1,iz,1) &
                              + a_TW4(0) *u(ix,iy,iz,1) &
                              + a_TW4(+1)*u(ix,iy+1,iz,1) &
                              + a_TW4(+2)*u(ix,iy+2,iz,1) &
                              + a_TW4(+3)*u(ix,iy+3,iz,1))*dy_inv

                        u_dz = (a_TW4(-3)*u(ix,iy,iz-3,1) &
                              + a_TW4(-2)*u(ix,iy,iz-2,1) &
                              + a_TW4(-1)*u(ix,iy,iz-1,1) &
                              + a_TW4(0) *u(ix,iy,iz,1) &
                              + a_TW4(+1)*u(ix,iy,iz+1,1) &
                              + a_TW4(+2)*u(ix,iy,iz+2,1) &
                              + a_TW4(+3)*u(ix,iy,iz+3,1))*dz_inv

                        v_dx = (a_TW4(-3)*u(ix-3,iy,iz,2) &
                              + a_TW4(-2)*u(ix-2,iy,iz,2)  &
                              + a_TW4(-1)*u(ix-1,iy,iz,2) &
                              + a_TW4(0) *u(ix,iy,iz,2) &
                              + a_TW4(+1)*u(ix+1,iy,iz,2) &
                              + a_TW4(+2)*u(ix+2,iy,iz,2) &
                              + a_TW4(+3)*u(ix+3,iy,iz,2))*dx_inv

                        v_dz = (a_TW4(-3)*u(ix,iy,iz-3,2) &
                              + a_TW4(-2)*u(ix,iy,iz-2,2) &
                              + a_TW4(-1)*u(ix,iy,iz-1,2) &
                              + a_TW4(0) *u(ix,iy,iz,2) &
                              + a_TW4(+1)*u(ix,iy,iz+1,2) &
                              + a_TW4(+2)*u(ix,iy,iz+2,2) &
                              + a_TW4(+3)*u(ix,iy,iz+3,2))*dz_inv

                        w_dx = (a_TW4(-3)*u(ix-3,iy,iz,3) &
                              + a_TW4(-2)*u(ix-2,iy,iz,3) &
                              + a_TW4(-1)*u(ix-1,iy,iz,3) &
                              + a_TW4(0) *u(ix,iy,iz,3) &
                              + a_TW4(+1)*u(ix+1,iy,iz,3) &
                              + a_TW4(+2)*u(ix+2,iy,iz,3) &
                              + a_TW4(+3)*u(ix+3,iy,iz,3))*dx_inv

                        w_dy = (a_TW4(-3)*u(ix,iy-3,iz,3) &
                              + a_TW4(-2)*u(ix,iy-2,iz,3) &
                              + a_TW4(-1)*u(ix,iy-1,iz,3) &
                              + a_TW4(0) *u(ix,iy,iz,3) &
                              + a_TW4(+1)*u(ix,iy+1,iz,3) &
                              + a_TW4(+2)*u(ix,iy+2,iz,3) &
                              + a_TW4(+3)*u(ix,iy+3,iz,3))*dy_inv

                        vor(ix,iy,iz,1) = w_dy - v_dz
                        vor(ix,iy,iz,2) = u_dz - w_dx
                        vor(ix,iy,iz,3) = v_dx - u_dy
                    end do
                end do
            end do

        case default
            call abort(1902201, "unknown order_discretization in ACM vorticity")
        end select
    endif


end subroutine
