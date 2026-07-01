!-----------------------------------------------------------------------------
! main level wrapper to compute statistics (such as mean flow, global energy,
! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
! on individual blocks. This requires one to use the same staging concept as for the RHS.
!-----------------------------------------------------------------------------
subroutine STATISTICS_NSPP( time, dt, u, g, x0, dx, stage, work, mask )
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
    integer(kind=ik) :: mpierr, ix, iy, iz, k, Bs(1:3), x1,x2,y1,y2,z1,z2
    real(kind=rk) :: meanflow_block(1:3), residual_block(1:3), ekin_block, mask_volume(1:ncolors), sponge_volume, meanflow_channel_block(1:3)
    real(kind=rk) :: force_block(1:3, 1:ncolors), moment_block(1:3, 1:ncolors), x_glob(1:3), x_lev(1:3)
    real(kind=rk) :: x0_moment(1:3, 1:ncolors), ipowtotal(1:n_insects), apowtotal(1:n_insects)
    real(kind=rk) :: CFL, CFL_eta, CFL_nu, penal_power_block(1:3), usx, usy, usz, chi, chi_sponge, dissipation_block
    real(kind=rk) :: C_eta_inv, C_sponge_inv, dV, x, y, z, penal(1:3), V_channel
    real(kind=rk) :: dxyz(1:3), iwmoment(1:n_insects, 1:3,1:5)
    real(kind=rk), save :: umag, umax, dx_min, scalar_removal_block, dissipation, u_RMS
    ! Color defines which objects belong together, default values are:
    ! Color         Description
    !   0           Boring parts (channel walls, cavity)
    !   1           Interesting parts (e.g. a cylinder), for the insects this is BODY
    !   2           Other parts, for the insects, this is LEFT WING
    !   3           For the insects, this is RIGHT WING
    !   4           Other parts, for the insects, this is 2ND LEFT WING
    !   5           For the insects, this is 2ND RIGHT WING
    integer(kind=2) :: color
    logical :: is_insect, has_two_wings
    integer :: i_insect
    character(len=64) :: fname

    if (.not. params_nspp%initialized) write(*,*) "WARNING: STATISTICS_NSPP called but NSPP not initialized"

    ! compute the size of blocks
    Bs = 1
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    if (params_nspp%dim == 3) Bs(3) = size(u,3) - 2*g

    ! block interior slicing excluding ghost points, for easier reading
    x1 = g+1
    x2 = Bs(1)+g
    y1 = g+1
    y2 = Bs(2)+g
    if (params_nspp%dim==3) then
        z1 = g+1
        z2 = Bs(3)+g
    else
        z1 = 1
        z2 = 1
    endif
    dV = product(dx(1:params_nspp%dim))

    C_eta_inv    = 1.0_rk / params_nspp%C_eta
    C_sponge_inv = 1.0_rk / params_nspp%C_sponge

    ! save some computing time by using a logical and not comparing every time
    is_insect = .false.
    if (any(params_nspp%geometries(:) == "Insect").or.any(params_nspp%geometries(:)=="active_grid") &
        .or. any(params_nspp%geometries(:)=="cylinder-free").or.any(params_nspp%geometries(:)=="sphere-free").or.any(params_nspp%geometries(:)=="plate-free")) is_insect = .true.

    has_two_wings = .false.
    do i_insect = 1, n_insects
        has_two_wings = has_two_wings .or. (insects(i_insect)%second_wing_pair)
    enddo



    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.
        ! performs initializations in the RHS module, such as resetting integrals
        params_nspp%mean_flow = 0.0_rk
        params_nspp%force_color = 0.0_rk
        params_nspp%meanflow_channel = 0.0_rk
        params_nspp%moment_color = 0.0_rk
        params_nspp%e_kin = 0.0_rk
        params_nspp%enstrophy = 0.0_rk
        params_nspp%max_vort = 0.0_rk
        params_nspp%helicity = 0.0_rk
        params_nspp%mask_volume = 0.0_rk
        params_nspp%sponge_volume = 0.0_rk
        params_nspp%u_residual = 0.0_rk
        params_nspp%div_max = 0.0_rk
        params_nspp%div_min = 0.0_rk
        params_nspp%penal_power = 0.0_rk
        params_nspp%scalar_removal = 0.0_rk
        params_nspp%dissipation = 0.0_rk
        params_nspp%umag = 0.0_rk
        if (params_nspp%time_statistics) then
            params_nspp%time_statistics_mean = 0.0_rk
            params_nspp%time_statistics_maxabs = 0.0_rk
        endif

        dx_min = 90.0e9_rk

        ! let's update all insects. If there is none, then this is just an empty loop, so no problemo
        call Update_All_Insects(time)

    case ("integral_stage")
        ! minium spacing of current grid, not smallest possible one.
        dx_min = min( dx_min, minval(dx(1:params_nspp%dim)) )
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! This stage contains all operations which are running on the blocks
        !
        ! called for each block.

        do k = 1, size(u,4)
            if (maxval(abs(u(:,:,:,k))) > LIM_DIVERGED) then
                write(*,'("Statistics: maxval in u(:,:,:,",i2,") = ", es10.3, ", Block with origin", 3(1x,es9.2), " and dx", 3(1x,es9.2))') k, maxval(abs(u(:,:,:,k))), x0, dx

                call dump_block_fancy(u(:,:,:,k:k), "block_NSPP_diverged_statistics.txt", Bs, g)

                ! done by all ranks but well I hope the cluster can take one for the team
                ! This (empty) file is for scripting purposes on the supercomputers.
                open (77, file='NSPP_diverged', status='replace')
                close(77)

                call abort(0409201934,"NSPP fail: very very large values in state vector.")
            endif
        enddo

        ! tmp values for computing the current block only
        meanflow_channel_block = 0.0_rk
        meanflow_block = 0.0_rk
        force_block = 0.0_rk
        moment_block = 0.0_rk
        residual_block = 0.0_rk
        ekin_block = 0.0_rk
        mask_volume = 0.0_rk
        sponge_volume = 0.0_rk
        penal_power_block = 0.0_rk
        scalar_removal_block = 0.0_rk

        ! default for moment computation is the centre point of the object (=lever)
        ! In general, this is not block-dependent, but it's cheap to compute anyways, so I won't make it persistent
        x0_moment(1,:) = params_nspp%x_cntr(1)
        x0_moment(2,:) = params_nspp%x_cntr(2)
        x0_moment(3,:) = params_nspp%x_cntr(3)
        ! in insects, other levers are used
        ! some preparations that we do not want to compute each time within the loop
        do i_insect = 1, n_insects
            ! body moment
            x0_moment(1:3, Insects(i_insect)%color_body) = Insects(i_insect)%xc_body_g
            ! left wing
            x0_moment(1:3, Insects(i_insect)%color_l) = Insects(i_insect)%x_pivot_l_g
            ! right wing
            x0_moment(1:3, Insects(i_insect)%color_r) = Insects(i_insect)%x_pivot_r_g
            ! second left and second right wings
            if (Insects(i_insect)%second_wing_pair) then
                x0_moment(1:3, Insects(i_insect)%color_l2) = Insects(i_insect)%x_pivot_l2_g
                x0_moment(1:3, Insects(i_insect)%color_r2) = Insects(i_insect)%x_pivot_r2_g
            endif
        enddo


        ! most integral quantities can be computed with block-wise function calls
        ! compute mean flow for output in statistics
        do k = 1, params_nspp%dim
            meanflow_block(k) = sum(u(x1:x2, y1:y2, z1:z2, k))
        enddo

        ! kinetic energy
        ekin_block = 0.5_rk*sum( u(x1:x2, y1:y2, z1:z2,1:params_nspp%dim)**2 )

        ! square of maximum of velocity in the field
        params_nspp%umag = max( params_nspp%umag, maxval(sum(u(x1:x2, y1:y2, z1:z2,1:params_nspp%dim)**2, dim=4)))

        ! maximum/min divergence in velocity field
        call compute_divergence( u(:,:,:,1:params_nspp%dim), dx, Bs, g, params_nspp%discretization, work(:,:,:,1))

        ! mask divergence inside the solid body
        where (mask(:,:,:,1)>0.0_rk)
            work(:,:,:,1) = 0.00_rk
        end where
        params_nspp%div_max = max( params_nspp%div_max, maxval(work(x1:x2, y1:y2, z1:z2,1)))
        params_nspp%div_min = min( params_nspp%div_min, minval(work(x1:x2, y1:y2, z1:z2,1)))

        ! penalization part
        if (params_nspp%penalization .or. params_nspp%use_sponge) then
            ! volume of sponge
            if (params_nspp%use_sponge) sponge_volume = sum(mask(x1:x2, y1:y2, z1:z2, 6))

            ! penalization power (input from solid motion), see Engels et al. J. Comput. Phys. 2015
            ! energy input from solid : (usx * (u-uxs) + usy * (v-usy) + usz * (w-usz)) * chi / C_eta
            penal_power_block(1) = sum( &
                 (mask(x1:x2, y1:y2, z1:z2, 2)*(u(x1:x2, y1:y2, z1:z2,1)-mask(x1:x2, y1:y2, z1:z2, 2)) &
                + mask(x1:x2, y1:y2, z1:z2, 3)*(u(x1:x2, y1:y2, z1:z2,2)-mask(x1:x2, y1:y2, z1:z2, 3))) &
                * mask(x1:x2, y1:y2, z1:z2, 1)*C_eta_inv)
            if (params_nspp%dim == 3) then
                penal_power_block(1) = penal_power_block(1) + sum( &
                      mask(x1:x2, y1:y2, z1:z2, 4)*(u(x1:x2, y1:y2, z1:z2,3)-mask(x1:x2, y1:y2, z1:z2, 4)) &
                    * mask(x1:x2, y1:y2, z1:z2, 1)*C_eta_inv )
            endif
            ! dissipation inside solid : ((u-uxs)**2 + (v-usy)**2 + (w-usz)**2) * chi / C_eta
            penal_power_block(2) = sum( &
                 ((u(x1:x2, y1:y2, z1:z2,1)-mask(x1:x2, y1:y2, z1:z2, 2))**2 &
                + (u(x1:x2, y1:y2, z1:z2,2)-mask(x1:x2, y1:y2, z1:z2, 3))**2) &
                * mask(x1:x2, y1:y2, z1:z2, 1)*C_eta_inv )
            if (params_nspp%dim == 3) then
                penal_power_block(2) = penal_power_block(2) + sum(((u(x1:x2, y1:y2, z1:z2,3)-mask(x1:x2, y1:y2, z1:z2, 4))**2) * mask(x1:x2, y1:y2, z1:z2, 1)*C_eta_inv )
            endif
            ! sponge input : ( u*(u-ubar) + v*(v-vbar) + w*(w-wbar) ) * chi_sponge / C_sponge
            if (params_nspp%use_sponge) then
                penal_power_block(3) = sum( &
                     (u(x1:x2, y1:y2, z1:z2,1)*(u(x1:x2, y1:y2, z1:z2,1)-params_nspp%u_mean_set(1)) &
                    + u(x1:x2, y1:y2, z1:z2,2)*(u(x1:x2, y1:y2, z1:z2,2)-params_nspp%u_mean_set(2))) &
                    * mask(x1:x2, y1:y2, z1:z2, 6) * C_sponge_inv )
                if (params_nspp%dim == 3) then
                    penal_power_block(3) = penal_power_block(3) + sum( &
                        u(x1:x2, y1:y2, z1:z2,3)*(u(x1:x2, y1:y2, z1:z2,3)-params_nspp%u_mean_set(3)) &
                        * mask(x1:x2, y1:y2, z1:z2, 6) * C_sponge_inv )
                endif
            endif

            ! residual velocity in the solid domain
            do k = 1, params_nspp%dim
                residual_block(k) = maxval( abs(u(x1:x2, y1:y2, z1:z2,k)-mask(x1:x2, y1:y2, z1:z2, k+1)) * mask(x1:x2, y1:y2, z1:z2, 1) )
            enddo

            ! some computations need point-wise loops, as we are picking the color of the points
            do iz = z1,z2
                if (params_nspp%dim == 2) then
                    z = 0.0_rk
                else
                    z = x0(3) + dble(iz-(g+1)) * dx(3)
                endif
                do iy = y1,y2
                    y = x0(2) + dble(iy-(g+1)) * dx(2) 
                    do ix = x1,x2
                        x = x0(1) + dble(ix-(g+1)) * dx(1)

                        color = int( mask(ix, iy, iz, 5), kind=2 )

                        ! mask volume
                        mask_volume(color) = mask_volume(color) + mask(ix,iy,iz,1)

                        ! penalization term
                        penal = -mask(ix,iy,iz,1) * (u(ix,iy,iz,1:3) - mask(ix,iy,iz,2:4)) * C_eta_inv

                        ! forces acting on this color
                        force_block(1:params_nspp%dim, color) = force_block(1:params_nspp%dim, color) - penal(1:params_nspp%dim)

                        ! moment with color-dependent lever
                        x_lev(1:3) = (/x, y, z/) - x0_moment(1:3, color)

                        ! is the obstacle is near the boundary, parts of it may cross the periodic
                        ! boundary. therefore, ensure that xlev is periodized:
                        ! x_lev = periodize_coordinate(x_lev, (/xl,yl,zl/))
                        
                        ! moment acting on this color
                        moment_block(1:3,color) = moment_block(1:3,color) - cross(x_lev, penal)


                        ! moments. For insects, we compute the total moment wrt to the body center, and
                        ! the wing moments wrt to the hinge points. The latter two are used to compute the
                        ! aerodynamic power. Makes sense only in 3D.
                        if (is_insect .and. params_nspp%dim == 3) then
                            do i_insect = 1, n_insects
                                if (any(color == (/insects(i_insect)%color_body, insects(i_insect)%color_l, insects(i_insect)%color_r, insects(i_insect)%color_l2, insects(i_insect)%color_r2/))) then
                                    ! in the geometry color of the insect, we compute the total mask volume, force and moment for the whole
                                    ! insect wrt the center point (body+wings)
                                    mask_volume(insects(i_insect)%color_geometry) = mask_volume(insects(i_insect)%color_geometry) + mask(ix,iy,iz,1)
                                    force_block(1:params_nspp%dim, insects(i_insect)%color_geometry) = force_block(1:params_nspp%dim, insects(i_insect)%color_geometry) - penal(1:params_nspp%dim)
                                    x_lev(1:3) = (/x, y, z/) - Insects(1)%xc_body_g(1:3)
                                    moment_block(:,insects(i_insect)%color_geometry)  = moment_block(:,insects(i_insect)%color_geometry) - cross(x_lev, penal)
                                endif
                            enddo
                        endif
                    enddo
                enddo
            enddo

            ! if the scalar BC is Dirichlet, then the solid absorbs some scalar, and it makes
            ! sense to keep track of this. however, note that with Neumann BC, that makes no
            ! sense (zero flux is imposed and the solution for the scalar inside the solid is arbitrary)
            if (params_nspp%use_passive_scalar .and. params_nspp%scalar_BC_type == "dirichlet") then
                ! should we ever use more than 1 scalar seriously, this has to be adopted
                ! because it uses only the 1st one (:,:,:,5)
                ! assumes *homogeneous* dirichlet condition
                scalar_removal_block = sum(mask(x1:x2, y1:y2, z1:z2, 1) * u(x1:x2, y1:y2, z1:z2, params_nspp%dim+2)) * C_eta_inv
            endif
        endif

        ! we just computed the values on the current block, which we now add to the
        ! existing blocks in the variables (recall normalization by dV)
        params_nspp%u_residual     = params_nspp%u_residual     + residual_block * dV
        params_nspp%mean_flow      = params_nspp%mean_flow      + meanflow_block * dV
        params_nspp%meanflow_channel = params_nspp%meanflow_channel + meanflow_channel_block * dV
        params_nspp%mask_volume    = params_nspp%mask_volume    + mask_volume * dV
        params_nspp%sponge_volume  = params_nspp%sponge_volume  + sponge_volume * dV
        params_nspp%force_color    = params_nspp%force_color    + force_block * dV
        params_nspp%moment_color   = params_nspp%moment_color   + moment_block * dV
        params_nspp%e_kin          = params_nspp%e_kin          + ekin_block * dV
        params_nspp%penal_power    = params_nspp%penal_power    + penal_power_block * dV
        params_nspp%scalar_removal = params_nspp%scalar_removal + scalar_removal_block * dV

        !-------------------------------------------------------------------------
        ! compute enstrophy in the whole domain (including penalized regions)
        ! note in 2D case, uz is ignored, so we pass p=u(:,:,:,3) just for fun.
        call compute_vorticity(u(:,:,:,1:params_nspp%dim), dx, Bs, g, params_nspp%discretization, work(:,:,:,:))

        if (params_nspp%dim == 2) then
            params_nspp%enstrophy = params_nspp%enstrophy + 0.5_rk*sum(work(x1:x2, y1:y2, z1:z2, 1)**2)*dV
            params_nspp%max_vort = max( params_nspp%max_vort, maxval( abs(work(x1:x2, y1:y2, z1:z2, 1)) ) )
        else
            params_nspp%enstrophy = params_nspp%enstrophy + 0.5_rk*sum(work(x1:x2, y1:y2, z1:z2, 1:3)**2)*dV
            params_nspp%max_vort = max( params_nspp%max_vort, maxval( sqrt(sum( work(x1:x2, y1:y2, z1:z2, 1:3)**2, dim=4 )) ) )
            params_nspp%helicity = params_nspp%helicity + 0.5_rk*sum(work(x1:x2, y1:y2, z1:z2, 1:3) * u(x1:x2, y1:y2, z1:z2, 1:3))*dV
        end if

        if (params_nspp%nu > 0.0_rk) then
            call compute_dissipation(u(:,:,:,1:params_nspp%dim), dx, Bs, g, params_nspp%discretization, work(:,:,:,1))
            params_nspp%dissipation = params_nspp%dissipation - params_nspp%nu * sum(work(x1:x2, y1:y2, z1:z2, 1)) * dV
        endif

        ! we want to output some measurements for the timestatistics, just to see if everything is in order
        if (params_nspp%time_statistics) then
            do k = 1, params_nspp%n_time_statistics
                params_nspp%time_statistics_mean(k) = params_nspp%time_statistics_mean(k) + dV * sum( u(x1:x2, y1:y2, z1:z2,params_nspp%dim+1+params_nspp%N_scalars+k) )
                params_nspp%time_statistics_maxabs(k) = max( params_nspp%time_statistics_maxabs(k), maxval(abs(u(x1:x2, y1:y2, z1:z2,params_nspp%dim+1+params_nspp%N_scalars+k))) )
            enddo
        endif

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.

        ! mean flow (in entire domain)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        params_nspp%mean_flow = params_nspp%mean_flow / product(params_nspp%domain_size(1:params_nspp%dim))

        if (params_nspp%use_channel_forcing) then
            ! mean flow but only in fluid domain
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%meanflow_channel, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

            V_channel = params_nspp%domain_size(1)*params_nspp%domain_size(3)*(params_nspp%domain_size(2)-2.0_rk*params_nspp%h_channel)
            params_nspp%meanflow_channel = params_nspp%meanflow_channel / V_channel
        endif

        if (params_nspp%penalization .or. params_nspp%use_sponge) then
            !-------------------------------------------------------------------------
            ! volume of mask (useful to see if it is properly generated)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%mask_volume, size(params_nspp%mask_volume), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            ! volume of sponge
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%sponge_volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

            ! force & moment
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%force_color, size(params_nspp%force_color), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%moment_color, size(params_nspp%moment_color), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            ! residual velocity in solid domain
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%u_residual, 3, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%penal_power, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        endif

        !-------------------------------------------------------------------------
        ! scalar removal
        if (params_nspp%use_passive_scalar .and. params_nspp%scalar_BC_type == "dirichlet") then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%scalar_removal, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        endif

        !-------------------------------------------------------------------------
        ! kinetic energy
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! divergence
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%div_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%div_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic enstrophy
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%max_vort, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%helicity, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        if (params_nspp%nu > 0.0_rk) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%dissipation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        endif

        ! minium spacing of current grid, not smallest possible one.
        call MPI_ALLREDUCE(MPI_IN_PLACE, dx_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%umag, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        ! time statistics
        if (params_nspp%time_statistics) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%time_statistics_mean, params_nspp%n_time_statistics, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
            params_nspp%time_statistics_mean = params_nspp%time_statistics_mean / product(params_nspp%domain_size(1:params_nspp%dim))
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_nspp%time_statistics_maxabs, params_nspp%n_time_statistics, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        endif

        umag = params_nspp%umag

        ! compute aerodynamic power
        if (is_insect) then
            do i_insect = 1, n_insects
                ! store moments for the insect (so that it can compute the aerodynamic power) - body, wing_l, wing_r, (wing_l2, wing_r2)
                Insects(i_insect)%PartIntegrals( 1 )%Torque = params_nspp%moment_color(1:3, Insects(i_insect)%color_body )
                Insects(i_insect)%PartIntegrals( 2 )%Torque = params_nspp%moment_color(1:3, Insects(i_insect)%color_l )
                Insects(i_insect)%PartIntegrals( 3 )%Torque = params_nspp%moment_color(1:3, Insects(i_insect)%color_r )
                if (Insects(i_insect)%second_wing_pair) then
                    Insects(i_insect)%PartIntegrals( 4 )%Torque = params_nspp%moment_color(1:3, Insects(i_insect)%color_l2 )
                    Insects(i_insect)%PartIntegrals( 5 )%Torque = params_nspp%moment_color(1:3, Insects(i_insect)%color_r2 )
                endif

                call aero_power (Insects(i_insect), apowtotal(i_insect))
                call inert_power(Insects(i_insect), ipowtotal(i_insect), iwmoment(i_insect, :, :))

                ! store simplified insect state (which is NOT the same as the state vector that we use in free_flight simulations)
                call append_t_file( 'insect_state.t', (/time, Insects(i_insect)%xc_body_g(1:3), Insects(i_insect)%vc_body_g(1:3), &
                    Insects(i_insect)%psi, Insects(i_insect)%beta, Insects(i_insect)%gamma /) )
            enddo
        endif

        !-------------------------------------------------------------------------
        ! write statistics to ascii files.
        if (params_nspp%mpirank == 0) then
            call append_t_file( 'umag.t', (/time, sqrt(umag)/) )

            CFL   = dt * (sqrt(umag)) / dx_min
            CFL_nu = dt * params_nspp%nu / dx_min**2
            CFL_eta = dt / params_nspp%C_eta

            call append_t_file( 'CFL.t', (/time, CFL, CFL_nu, CFL_eta/) )
            call append_t_file( 'meanflow.t', (/time, params_nspp%mean_flow/) )

            if (params_nspp%use_channel_forcing) then
                call append_t_file( 'meanflow_channel.t', (/time, params_nspp%meanflow_channel   /) )
            endif

            call append_t_file( 'div.t', (/time, params_nspp%div_max, params_nspp%div_min/) )
            if (params_nspp%penalization .or. params_nspp%use_sponge) then
                ! total force (excluding insect parts, they are considered by the full geometry, otherwise we count them twice)
                call append_t_file( 'forces.t', (/time, sum(params_nspp%force_color(1,1:params_nspp%n_geometries)), &
                sum(params_nspp%force_color(2,1:params_nspp%n_geometries)), sum(params_nspp%force_color(3,1:params_nspp%n_geometries)) /) )

                ! forces/moment for individual colors.
                ! This is what we should have done in the first place. For the insects below,
                ! some colors are (also) stored to different names. That's unfortunate but not a big deal.
                ! Reshape flattens the array, so that we'd have for a(2,2): (/a(1,1), a(2,1), a(1,2), a(2,2)/)
                call append_t_file( "forces_color.t", (/time, reshape(params_nspp%force_color(:,:), (/ 3*ncolors/))/) )

                ! save moment for each color in one file
                call append_t_file( "moments_color.t", (/time, reshape(params_nspp%moment_color(:,:), (/ 3*ncolors/))/) )

                if (is_insect) then
                    call append_t_file( 'aero_power.t', (/time, apowtotal(:), ipowtotal(:)/) )

                    ! kinematics
                    call Write_insect_data( time )

                    ! information for each insect
                    do i_insect = 1, n_insects
                        ! total moment w.r.t body center is computed in zeroth color slot:
                        color = insects(i_insect)%color_geometry
                        call append_t_file( 'moments.t', (/time, params_nspp%moment_color(:,color), dble(i_insect)/) )

                        ! body
                        color = insects(i_insect)%color_body
                        call append_t_file( 'forces_body.t', (/time, params_nspp%force_color(:,color), dble(i_insect)/) )
                        call append_t_file( 'moments_body.t', (/time, params_nspp%moment_color(:,color), dble(i_insect)/) )

                        ! left wing
                        color = insects(i_insect)%color_l
                        call append_t_file( 'forces_leftwing.t', (/time, params_nspp%force_color(:,color), dble(i_insect)/) )
                        call append_t_file( 'moments_leftwing.t', (/time, params_nspp%moment_color(:,color), iwmoment(i_insect, :, color), dble(i_insect)/) )

                        ! right wing
                        color = insects(i_insect)%color_r
                        call append_t_file( 'forces_rightwing.t', (/time, params_nspp%force_color(:,color), dble(i_insect)/) )
                        call append_t_file( 'moments_rightwing.t', (/time, params_nspp%moment_color(:,color), iwmoment(i_insect, :, color), dble(i_insect)/) )

                        ! kinematics data ('kinematics.t')
                        if (Insects(i_insect)%second_wing_pair) then
                            ! second left wing
                            color = insects(i_insect)%color_l2
                            call append_t_file( 'forces_leftwing2.t', (/time, params_nspp%force_color(:,color), dble(i_insect)/) )
                            call append_t_file( 'moments_leftwing2.t', (/time, params_nspp%moment_color(:,color), iwmoment(i_insect, :, color), dble(i_insect)/) )

                            ! second right wing
                            color = insects(i_insect)%color_r2
                            call append_t_file( 'forces_rightwing2.t', (/time, params_nspp%force_color(:,color), dble(i_insect)/) )
                            call append_t_file( 'moments_rightwing2.t', (/time, params_nspp%moment_color(:,color), iwmoment(i_insect, :, color), dble(i_insect)/) )

                        endif
                    enddo
                endif

                call append_t_file( 'mask_volume.t', (/time, params_nspp%mask_volume(1:ncolors), params_nspp%sponge_volume/) )
                call append_t_file( 'penal_power.t', (/time, params_nspp%penal_power/) )
                call append_t_file( 'u_residual.t', (/time, params_nspp%u_residual/) )    
            endif

            if (params_nspp%use_passive_scalar .and. params_nspp%scalar_BC_type == "dirichlet") then
                call append_t_file( 'scalar_removal.t', (/time, params_nspp%scalar_removal/) )
            endif
            call append_t_file( 'e_kin.t', (/time, params_nspp%e_kin/) )
            call append_t_file( 'enstrophy.t', (/time, params_nspp%enstrophy, params_nspp%max_vort/) )
            if (params_nspp%dim == 3) then
                call append_t_file( 'helicity.t', (/time, params_nspp%helicity/) )
            endif
            if (params_nspp%nu > 0.0_rk)  then
                call append_t_file( 'dissipation.t', (/time, params_nspp%dissipation/) )
            endif

            ! turbulent statistics - these are normed by the volume!
            if (params_nspp%nu*params_nspp%enstrophy > 0.0_rk .and. params_nspp%HIT_linear_forcing) then
                ! dissipation = 2*params_nspp%nu*params_nspp%enstrophy/product(params_nspp%domain_size(1:params_nspp%dim))
                dissipation = params_nspp%dissipation/product(params_nspp%domain_size(1:params_nspp%dim))
                u_RMS = sqrt(2*params_nspp%e_kin/product(params_nspp%domain_size(1:params_nspp%dim))/3)
                call append_t_file( 'turbulent_statistics.t', (/time, dissipation, params_nspp%e_kin/product(params_nspp%domain_size(1:params_nspp%dim)), u_RMS, &
                    (params_nspp%nu**3.0_rk / dissipation)**0.25_rk, sqrt(params_nspp%nu/dissipation), (params_nspp%nu*dissipation)**0.25_rk, &
                    sqrt(15.0_rk*params_nspp%nu*u_RMS**2/dissipation), sqrt(15.0_rk*params_nspp%nu*u_RMS**2/dissipation)*u_RMS/params_nspp%nu/))
            endif

            ! time statistics
            if (params_nspp%time_statistics) then
                call append_t_file( 'time_statistics_mean.t', (/time, params_nspp%time_statistics_mean/) )
                call append_t_file( 'time_statistics_maxabs.t', (/time, params_nspp%time_statistics_maxabs/) )
            endif

            ! this file is to simply keep track of simulations, should they be restarted with different parameters.
            ! We just keep track of the most important and most frequently changed parameters.
            call append_t_file( 'parameters.t', (/time, params_nspp%C_eta, params_nspp%C_sponge, params_nspp%nu/) )
        end if

    case default
        call abort(7772,"the STATISTICS wrapper requests a stage this physics module cannot handle.")

    end select


end subroutine STATISTICS_NSPP
