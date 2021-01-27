!-----------------------------------------------------------------------------
! main level wrapper to compute statistics (such as mean flow, global energy,
! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
! on individual blocks. This requires one to use the same staging concept as for the RHS.
!-----------------------------------------------------------------------------
subroutine STATISTICS_ACM( time, dt, u, g, x0, dx, stage, work, mask )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time, dt

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! work data, for mask, vorticity etc. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(inout) :: mask(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in) :: stage

    ! local variables
    integer(kind=ik) :: mpierr, ix, iy, iz, k
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: tmp(1:6), meanflow_block(1:3), residual_block(1:3), ekin_block, tmp_volume, tmp_volume2
    real(kind=rk) :: force_block(1:3, 0:6), moment_block(1:3,0:6), x_glob(1:3), x_lev(1:3), penalpower_block
    real(kind=rk) :: x0_moment(1:3,0:6), ipowtotal=0.0_rk, apowtotal=0.0_rk
    real(kind=rk) :: CFL, CFL_eta, CFL_nu
    real(kind=rk) :: C_eta_inv, dV, x, y, z, penal(1:3)
    real(kind=rk), dimension(3) :: dxyz
    real(kind=rk), dimension(1:3,1:5) :: iwmoment
    real(kind=rk), save :: umag, umax, dx_min
    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're ngeligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: div(:,:,:)
    ! Color         Description
    !   0           Boring parts (channel walls, cavity)
    !   1           Interesting parts (e.g. a cylinder), for the insects this is BODY
    !   2           Other parts, for the insects, this is LEFT WING
    !   3           For the insects, this is RIGHT WING
    !   4           Other parts, for the insects, this is 2ND LEFT WING
    !   5           For the insects, this is 2ND RIGHT WING
    integer(kind=2), allocatable, save :: mask_color(:,:,:)
    integer(kind=2) :: color
    logical :: is_insect
    !> coefficients for Tam&Webb
    real(kind=rk), parameter :: a(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    real(kind=rk) :: ux_dx, ux_dy, ux_dz, uy_dx, uy_dy, uy_dz, uz_dx, uz_dy, uz_dz, &
    vorx, vory, vorz, vortex_stretching, dx_inv, dy_inv, dz_inv



    if (.not. params_acm%initialized) write(*,*) "WARNING: STATISTICS_ACM called but ACM not initialized"

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    if (params_acm%dim==3) then
        if (.not. allocated(mask_color)) allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g))
        if (.not. allocated(div)) allocate(div(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g))
        dV = dx(1)*dx(2)*dx(3)
    else
        if (.not. allocated(mask_color)) allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
        if (.not. allocated(div)) allocate(div(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
        dV = dx(1)*dx(2)
    endif

    ! save some computing time by using a logical and not comparing strings every time
    is_insect = .false.
    if (params_acm%geometry == "Insect") is_insect = .true.



    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.
        ! performs initializations in the RHS module, such as resetting integrals
        params_acm%mean_flow = 0.0_rk
        params_acm%force_color = 0.0_rk
        params_acm%moment_color = 0.0_rk
        params_acm%e_kin = 0.0_rk
        params_acm%enstrophy = 0.0_rk
        params_acm%mask_volume = 0.0_rk
        params_acm%sponge_volume = 0.0_rk
        params_acm%vortex_stretching = 0.0_rk
        params_acm%u_residual = 0.0_rk
        params_acm%div_max = 0.0_rk
        params_acm%div_min = 0.0_rk
        params_acm%penalpower = 0.0_rk
        dx_min = 90.0e9_rk

        if (is_insect) then
            call Update_Insect(time, Insect)
        endif

    case ("integral_stage")
        ! minium spacing of current grid, not smallest possible one.
        dx_min = min( dx_min, minval(dx(1:params_acm%dim)) )
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! This stage contains all operations which are running on the blocks
        !
        ! called for each block.

        do k = 1, size(u,4)
            if (maxval(abs(u(:,:,:,k))) > 1.0e4_rk) then
                write(*,'("maxval in u(:,:,:,",i2,") = ", es15.8)') k, maxval(abs(u(:,:,:,k)))
                call abort(0409201934,"ACM fail: very very large values in state vector.")
            endif
        enddo

        ! tmp values for computing the current block only
        meanflow_block = 0.0_rk
        force_block = 0.0_rk
        moment_block = 0.0_rk
        residual_block = 0.0_rk
        ekin_block = 0.0_rk
        tmp_volume = 0.0_rk
        tmp_volume2 = 0.0_rk
        penalpower_block = 0.0_rk

        if (params_acm%dim == 2) then
            ! --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D ---
            C_eta_inv = 1.0_rk / params_acm%C_eta

            ! note in 2D case, uz is ignored, so we pass p just for fun.
            call divergence( u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, div)

            ! mask divergence inside the solid body
            where (mask(:,:,:,1)>0.0_rk)
                div = 0.00_rk
            end where

            do iy = g+1, Bs(2)+g-1 ! Note: loops skip redundant points
            do ix = g+1, Bs(1)+g-1
                ! coloring not implemented for 2D
                color = 0_2

                ! compute mean flow for output in statistics
                meanflow_block(1) = meanflow_block(1) + u(ix,iy,1,1)
                meanflow_block(2) = meanflow_block(2) + u(ix,iy,1,2)

                ! volume of mask (useful to see if it is properly generated)
                tmp_volume = tmp_volume + mask(ix,iy,1,1)
                tmp_volume2 = tmp_volume2 + mask(ix,iy,1,6)

                ! forces acting on body
                force_block(1, color) = force_block(1, color) + (u(ix,iy,1,1)-mask(ix,iy,1,2))*mask(ix,iy,1,1)*C_eta_inv
                force_block(2, color) = force_block(2, color) + (u(ix,iy,1,2)-mask(ix,iy,1,3))*mask(ix,iy,1,1)*C_eta_inv

                ! residual velocity in the solid domain
                residual_block(1) = max( residual_block(1), (u(ix,iy,1,1)-mask(ix,iy,1,2)) * mask(ix,iy,1,1))
                residual_block(2) = max( residual_block(2), (u(ix,iy,1,2)-mask(ix,iy,1,3)) * mask(ix,iy,1,1))

                ! kinetic energy
                ekin_block = ekin_block + 0.5_rk*sum( u(ix,iy,1,1:2)**2 )

                ! maximum of velocity in the field
                params_acm%umag = max( params_acm%umag, u(ix,iy,1,1)*u(ix,iy,1,1) + u(ix,iy,1,2)*u(ix,iy,1,2) )

                ! maximum/min divergence in velocity field
                params_acm%div_max = max( params_acm%div_max, div(ix,iy,1) )
                params_acm%div_min = min( params_acm%div_min, div(ix,iy,1) )
            enddo
            enddo

            ! this is a 2D case and u has only two components: the 3rd one (in fact thats the pressure)
            ! is ignored by the vorticity routine.
            call compute_vorticity(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, work(:,:,:,:))

            ! compute mean enstrophy over the entire domain (including the solid body, in case of penalization)
            params_acm%enstrophy = params_acm%enstrophy + sum(work(g+1:Bs(1)+g-1,g+1:Bs(2)+g-1,1,1)**2)*dx(1)*dx(2)

        else
            ! --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D ---
            C_eta_inv = 1.0_rk / params_acm%C_eta

            ! compute divergence on this block
            call divergence( u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, div)

            ! mask divergence inside the solid body
            where (mask(:,:,:,1)>0.0_rk)
                div = 0.00_rk
            end where


            do iz = g+1, Bs(3)+g-1 ! Note: loops skip redundant points
                z = x0(3) + dble(iz-(g+1)) * dx(3)
                do iy = g+1, Bs(2)+g-1
                    y = x0(2) + dble(iy-(g+1)) * dx(2)
                    do ix = g+1, Bs(1)+g-1
                        x = x0(1) + dble(ix-(g+1)) * dx(1)

                        ! compute mean flow for output in statistics
                        meanflow_block(1:3) = meanflow_block(1:3) + u(ix,iy,iz,1:3)
                        ekin_block = ekin_block + 0.5_rk*sum( u(ix,iy,iz,1:3)**2 )

                        ! maximum of velocity in the field
                        params_acm%umag = max( params_acm%umag, u(ix,iy,iz,1)*u(ix,iy,iz,1) + u(ix,iy,iz,2)*u(ix,iy,iz,2) + u(ix,iy,iz,3)*u(ix,iy,iz,3) )

                        ! maximum/min divergence in velocity field
                        params_acm%div_max = max( params_acm%div_max, div(ix,iy,iz) )
                        params_acm%div_min = min( params_acm%div_min, div(ix,iy,iz) )

                        ! volume of mask (useful to see if it is properly generated)
                        ! NOTE: in wabbit, mask is really the mask: it is not divided by C_eta yet.
                        tmp_volume = tmp_volume + mask(ix, iy, iz, 1)
                        tmp_volume2 = tmp_volume2 + mask(ix, iy, iz, 6)

                        ! get this points color
                        color = int( mask(ix, iy, iz, 5), kind=2 )

                        ! this is a bugfix: if the color is smaller 0 or greater 6, something
                        ! went wrong with mask color (happens when not actually using penalization, for instance)
                        if (color>0_2 .and. color < 6_2) then
                            ! penalization term
                            penal = -mask(ix,iy,iz,1) * (u(ix,iy,iz,1:3) - mask(ix,iy,iz,2:4)) * C_eta_inv

                            ! penalization power: the contribution to the energy eqn from the penalization
                            ! (ie moving boundarys act as energy sink/source)
                            ! is the scalar product of penalization term and u (see Engels et al J. Comp. Phys 2015) eqn.34
                            penalpower_block = penalpower_block + sum( penal*u(ix,iy,iz,1:3) )

                            ! forces acting on this color
                            force_block(1:3, color) = force_block(1:3, color) - penal

                            ! moments. For insects, we compute the total moment wrt to the body center, and
                            ! the wing moments wrt to the hinge points. The latter two are used to compute the
                            ! aerodynamic power. Makes sense only in 3D.
                            if (is_insect) then
                                ! point of reference for the moments
                                x0_moment = 0.0_rk
                                ! body moment
                                x0_moment(1:3, Insect%color_body) = Insect%xc_body_g
                                ! left wing
                                x0_moment(1:3, Insect%color_l) = Insect%x_pivot_l_g
                                ! right wing
                                x0_moment(1:3, Insect%color_r) = Insect%x_pivot_r_g
                                ! second left and second right wings
                                if (Insect%second_wing_pair) then
                                  x0_moment(1:3, Insect%color_l2) = Insect%x_pivot_l2_g
                                  x0_moment(1:3, Insect%color_r2) = Insect%x_pivot_r2_g
                                endif

                                ! exclude walls, trees, etc...
                                if (color > 0_2) then
                                    ! moment with color-dependent lever
                                    x_lev(1:3) = (/x, y, z/) - x0_moment(1:3, color)

                                    ! is the obstacle is near the boundary, parts of it may cross the periodic
                                    ! boundary. therefore, ensure that xlev is periodized:
                                    ! x_lev = periodize_coordinate(x_lev, (/xl,yl,zl/))

                                    ! Compute moments relative to each part
                                    moment_block(:,color) = moment_block(:,color) - cross(x_lev, penal)

                                    ! in the seventh color, we compute the total moment for the whole
                                    ! insect wrt the center point (body+wings)
                                    x_lev(1:3) = (/x, y, z/) - Insect%xc_body_g(1:3)
                                    moment_block(:,6)  = moment_block(:,6) - cross(x_lev, penal)
                                endif

                            endif

                            ! residual velocity in the solid domain
                            residual_block(1) = max( residual_block(1), (u(ix,iy,iz,1)-mask(ix,iy,iz,2))*mask(ix,iy,iz,1) )
                            residual_block(2) = max( residual_block(2), (u(ix,iy,iz,2)-mask(ix,iy,iz,3))*mask(ix,iy,iz,1) )
                            residual_block(3) = max( residual_block(3), (u(ix,iy,iz,3)-mask(ix,iy,iz,4))*mask(ix,iy,iz,1) )
                        endif
                    enddo
                enddo
            enddo

            ! For the enstrophy equation, we require the vorticity:
            call compute_vorticity(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, work(:,:,:,:))

            ! compute mean enstrophy over the entire domain (including the solid body, in case of penalization)
            params_acm%enstrophy = params_acm%enstrophy + sum( work(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 1:3)**2 )*dx(1)*dx(2)*dx(3)

            ! hack only for 4th order scheme
            if (params_acm%discretization == "FD_4th_central_optimized") then
                dx_inv = 1.0_rk / dx(1)
                dy_inv = 1.0_rk / dx(2)
                dz_inv = 1.0_rk / dx(3)

                do iz = g+1, Bs(3)+g-1 ! Note: loops skip redundant points
                    do iy = g+1, Bs(2)+g-1
                        do ix = g+1, Bs(1)+g-1
                            ! vorx = work(ix,iy,iz,1)
                            ! vory = work(ix,iy,iz,2)
                            ! vorz = work(ix,iy,iz,3)
                            !
                            ! ux = u(ix,iy,iz,1)
                            ! uy = u(ix,iy,iz,1)
                            ! uz = u(ix,iy,iz,1)

                            ux_dx = (a(-3)*u(ix-3,iy,iz,1) + a(-2)*u(ix-2,iy,iz,1) + a(-1)*u(ix-1,iy,iz,1) + a(0)*u(ix,iy,iz,1) &
                                  +  a(+1)*u(ix+1,iy,iz,1) + a(+2)*u(ix+2,iy,iz,1) + a(+3)*u(ix+3,iy,iz,1))*dx_inv
                            ux_dy = (a(-3)*u(ix,iy-3,iz,1) + a(-2)*u(ix,iy-2,iz,1) + a(-1)*u(ix,iy-1,iz,1) + a(0)*u(ix,iy,iz,1) &
                                  +  a(+1)*u(ix,iy+1,iz,1) + a(+2)*u(ix,iy+2,iz,1) + a(+3)*u(ix,iy+3,iz,1))*dy_inv
                            ux_dz = (a(-3)*u(ix,iy,iz-3,1) + a(-2)*u(ix,iy,iz-2,1) + a(-1)*u(ix,iy,iz-1,1) + a(0)*u(ix,iy,iz,1) &
                                  +  a(+1)*u(ix,iy,iz+1,1) + a(+2)*u(ix,iy,iz+2,1) + a(+3)*u(ix,iy,iz+3,1))*dz_inv

                            uy_dx = (a(-3)*u(ix-3,iy,iz,2) + a(-2)*u(ix-2,iy,iz,2) + a(-1)*u(ix-1,iy,iz,2) + a(0)*u(ix,iy,iz,2) &
                                  +  a(+1)*u(ix+1,iy,iz,2) + a(+2)*u(ix+2,iy,iz,2) + a(+3)*u(ix+3,iy,iz,2))*dx_inv
                            uy_dy = (a(-3)*u(ix,iy-3,iz,2) + a(-2)*u(ix,iy-2,iz,2) + a(-1)*u(ix,iy-1,iz,2) + a(0)*u(ix,iy,iz,2) &
                                  +  a(+1)*u(ix,iy+1,iz,2) + a(+2)*u(ix,iy+2,iz,2) + a(+3)*u(ix,iy+3,iz,2))*dy_inv
                            uy_dz = (a(-3)*u(ix,iy,iz-3,2) + a(-2)*u(ix,iy,iz-2,2) + a(-1)*u(ix,iy,iz-1,2) + a(0)*u(ix,iy,iz,2) &
                                  +  a(+1)*u(ix,iy,iz+1,2) + a(+2)*u(ix,iy,iz+2,2) + a(+3)*u(ix,iy,iz+3,2))*dz_inv

                            uz_dx = (a(-3)*u(ix-3,iy,iz,3) + a(-2)*u(ix-2,iy,iz,3) + a(-1)*u(ix-1,iy,iz,3) + a(0)*u(ix,iy,iz,3) &
                                  +  a(+1)*u(ix+1,iy,iz,3) + a(+2)*u(ix+2,iy,iz,3) + a(+3)*u(ix+3,iy,iz,3))*dx_inv
                            uz_dy = (a(-3)*u(ix,iy-3,iz,3) + a(-2)*u(ix,iy-2,iz,3) + a(-1)*u(ix,iy-1,iz,3) + a(0)*u(ix,iy,iz,3) &
                                  +  a(+1)*u(ix,iy+1,iz,3) + a(+2)*u(ix,iy+2,iz,3) + a(+3)*u(ix,iy+3,iz,3))*dy_inv
                            uz_dz = (a(-3)*u(ix,iy,iz-3,3) + a(-2)*u(ix,iy,iz-2,3) + a(-1)*u(ix,iy,iz-1,3) + a(0)*u(ix,iy,iz,3) &
                                  +  a(+1)*u(ix,iy,iz+1,3) + a(+2)*u(ix,iy,iz+2,3) + a(+3)*u(ix,iy,iz+3,3))*dz_inv

                            vorx = uz_dy - uy_dz
                            vory = ux_dz - uz_dx
                            vorz = uy_dx - ux_dy

                            vortex_stretching = vorx*vorx*ux_dx + vorx*vory*ux_dy + vorx*vorz*ux_dz + &
                                                vory*vorx*uy_dx + vory*vory*uy_dy + vory*vorz*uy_dz + &
                                                vorz*vorx*uz_dx + vorz*vory*uz_dy + vorz*vorz*uz_dz

                            work(ix,iy,iz,1) = vortex_stretching
                        enddo
                    enddo
                enddo

                ! sum up block values (process-local, MPI_ALLREDUCE will follow later)
                ! Note: skip redundant points (integration)
                params_acm%vortex_stretching = params_acm%vortex_stretching + sum( work(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 1) ) * dV
            endif

        endif ! 2D / 3D case

        ! we just computed the values on the current block, which we now add to the
        ! existing blocks in the variables (recall normalization by dV)
        params_acm%u_residual    = params_acm%u_residual    + residual_block * dV
        params_acm%mean_flow     = params_acm%mean_flow     + meanflow_block * dV
        params_acm%penalpower    = params_acm%penalpower    + penalpower_block * dV
        params_acm%mask_volume   = params_acm%mask_volume   + tmp_volume * dV
        params_acm%sponge_volume = params_acm%sponge_volume + tmp_volume2 * dV
        params_acm%force_color   = params_acm%force_color   + force_block * dV
        params_acm%moment_color  = params_acm%moment_color  + moment_block * dV
        params_acm%e_kin         = params_acm%e_kin         + ekin_block * dV

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.

        !-------------------------------------------------------------------------
        ! mean flow
        tmp(1:3) = params_acm%mean_flow
        call MPI_ALLREDUCE(tmp(1:3), params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        if (params_acm%dim == 2) then
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%domain_size(1)*params_acm%domain_size(2))
        else
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%domain_size(1)*params_acm%domain_size(2)*params_acm%domain_size(3))
        endif

        !-------------------------------------------------------------------------
        ! force & moment
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%force_color, size(params_acm%force_color), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%moment_color, size(params_acm%moment_color), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! residual velocity in solid domain
        tmp(1:3) = params_acm%u_residual
        call MPI_ALLREDUCE(tmp(1:3), params_acm%u_residual, 3, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! volume of mask (useful to see if it is properly generated)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%mask_volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%sponge_volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%penalpower, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic energy
        tmp(1) = params_acm%e_kin
        call MPI_ALLREDUCE(tmp(1), params_acm%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! divergence
        tmp(1) = params_acm%div_min
        call MPI_ALLREDUCE(tmp(1), params_acm%div_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, WABBIT_COMM, mpierr)

        tmp(1) = params_acm%div_max
        call MPI_ALLREDUCE(tmp(1), params_acm%div_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic enstrophy
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%vortex_stretching, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        tmp(1) = dx_min ! minium spacing of current grid, not smallest possible one.
        call MPI_ALLREDUCE(tmp(1), dx_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, WABBIT_COMM, mpierr)

        tmp(1) = params_acm%umag
        call MPI_ALLREDUCE(tmp(1), params_acm%umag, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        umag = params_acm%umag

        ! compute aerodynamic power
        if (is_insect) then
            ! store moments for the insect (so that it can compute the aerodynamic power)
            Insect%PartIntegrals( Insect%color_body )%Torque = params_acm%moment_color(1:3, Insect%color_body )
            Insect%PartIntegrals( Insect%color_l )%Torque = params_acm%moment_color(1:3, Insect%color_l )
            Insect%PartIntegrals( Insect%color_r )%Torque = params_acm%moment_color(1:3, Insect%color_r )
            if (Insect%second_wing_pair) then
                Insect%PartIntegrals( Insect%color_l2 )%Torque = params_acm%moment_color(1:3, Insect%color_l2 )
                Insect%PartIntegrals( Insect%color_r2 )%Torque = params_acm%moment_color(1:3, Insect%color_r2 )
            endif

            call aero_power (Insect, apowtotal)
            call inert_power(Insect, ipowtotal, iwmoment)
        endif

        !-------------------------------------------------------------------------
        ! write statistics to ascii files.
        if (params_acm%mpirank == 0) then
            if (umag > 1.0e-12_rk) then
                call append_t_file( 'umag.t', (/time, sqrt(umag), params_acm%c_0, &
                params_acm%c_0/sqrt(umag), sqrt(umag) + sqrt(params_acm%c_0**2 + umag)/) )
            else
                call append_t_file( 'umag.t', (/time, sqrt(umag), params_acm%c_0, &
                0.0_rk, sqrt(umag) + sqrt(params_acm%c_0**2 + umag) /) )
            endif

            CFL   = dt * (sqrt(umag) + sqrt(params_acm%c_0**2 + umag)) / dx_min
            CFL_nu = dt * params_acm%nu / dx_min**2
            CFL_eta = dt / params_acm%C_eta

            call append_t_file( 'CFL.t', (/time, CFL, CFL_nu, CFL_eta/) )
            call append_t_file( 'meanflow.t', (/time, params_acm%mean_flow/) )
            call append_t_file( 'div.t', (/time, params_acm%div_max, params_acm%div_min/) )
            call append_t_file( 'forces.t', (/time, sum(params_acm%force_color(1,:)), &
            sum(params_acm%force_color(2,:)), sum(params_acm%force_color(3,:)) /) )

            if (is_insect) then
                call append_t_file( 'aero_power.t', (/time, apowtotal, ipowtotal/) )

                color = 6_2
                call append_t_file( 'moments.t', (/time, params_acm%moment_color(:,color)/) )

                ! body
                color = Insect%color_body
                call append_t_file( 'forces_body.t', (/time, params_acm%force_color(:,color)/) )
                call append_t_file( 'moments_body.t', (/time, params_acm%moment_color(:,color)/) )

                ! left wing
                color = Insect%color_l
                call append_t_file( 'forces_leftwing.t', (/time, params_acm%force_color(:,color)/) )
                call append_t_file( 'moments_leftwing.t', (/time, params_acm%moment_color(:,color), iwmoment(:,color)/) )

                ! right wing
                color = Insect%color_r
                call append_t_file( 'forces_rightwing.t', (/time, params_acm%force_color(:,color)/) )
                call append_t_file( 'moments_rightwing.t', (/time, params_acm%moment_color(:,color), iwmoment(:,color)/) )

                ! kinematics data ('kinematics.t')
                if (Insect%second_wing_pair) then
                    ! second left wing
                    color = Insect%color_l2
                    call append_t_file( 'forces_leftwing2.t', (/time, params_acm%force_color(:,color)/) )
                    call append_t_file( 'moments_leftwing2.t', (/time, params_acm%moment_color(:,color), iwmoment(:,color)/) )

                    ! second right wing
                    color = Insect%color_r2
                    call append_t_file( 'forces_rightwing2.t', (/time, params_acm%force_color(:,color)/) )
                    call append_t_file( 'moments_rightwing2.t', (/time, params_acm%moment_color(:,color), iwmoment(:,color)/) )

                    call append_t_file( 'kinematics.t', (/time, Insect%xc_body_g, Insect%psi, Insect%beta, &
                    Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
                    Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
                    Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
                    Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w, &
                    Insect%alpha_l2, Insect%phi_l2, Insect%theta_l2, &
                    Insect%alpha_r2, Insect%phi_r2, Insect%theta_r2, &
                    Insect%rot_rel_wing_l2_w, Insect%rot_rel_wing_r2_w, &
                    Insect%rot_dt_wing_l2_w, Insect%rot_dt_wing_r2_w/) )
                else
                    call append_t_file( 'kinematics.t', (/time, Insect%xc_body_g, Insect%psi, Insect%beta, &
                    Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
                    Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
                    Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
                    Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w/) )
                endif
            endif

            call append_t_file( 'e_kin.t', (/time, params_acm%e_kin/) )
            call append_t_file( 'enstrophy.t', (/time, params_acm%enstrophy, params_acm%vortex_stretching/) )
            call append_t_file( 'penalpower.t', (/time, params_acm%penalpower/) )
            call append_t_file( 'mask_volume.t', (/time, params_acm%mask_volume, params_acm%sponge_volume/) )
            call append_t_file( 'u_residual.t', (/time, params_acm%u_residual/) )
        end if

    case default
        call abort(7772,"the STATISTICS wrapper requests a stage this physics module cannot handle.")

    end select


end subroutine STATISTICS_ACM
