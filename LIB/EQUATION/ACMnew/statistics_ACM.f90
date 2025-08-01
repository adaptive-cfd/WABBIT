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
    integer(kind=ik) :: mpierr, ix, iy, iz, k, Bs(1:3)
    real(kind=rk) :: meanflow_block(1:3), residual_block(1:3), ekin_block, tmp_volume, tmp_volume2, meanflow_channel_block(1:3)
    real(kind=rk) :: force_block(1:3, 0:ncolors), moment_block(1:3, 0:ncolors), x_glob(1:3), x_lev(1:3)
    real(kind=rk) :: x0_moment(1:3, 0:ncolors), ipowtotal=0.0_rk, apowtotal=0.0_rk
    real(kind=rk) :: CFL, CFL_eta, CFL_nu, penal_power_block(1:3), usx, usy, usz, chi, chi_sponge, dissipation_block
    real(kind=rk) :: C_eta_inv, C_sponge_inv, dV, x, y, z, penal(1:3), ACM_energy_block, C_0, V_channel
    real(kind=rk), dimension(3) :: dxyz
    real(kind=rk), dimension(1:3,1:5) :: iwmoment
    real(kind=rk), save :: umag, umax, dx_min, scalar_removal_block, dissipation, u_RMS
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
    logical :: is_insect, use_color
    character(len=64) :: fname

    if (.not. params_acm%initialized) write(*,*) "WARNING: STATISTICS_ACM called but ACM not initialized"

    ! compute the size of blocks
    Bs = 1
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    if (params_acm%dim == 3) Bs(3) = size(u,3) - 2*g

    C_eta_inv    = 1.0_rk / params_acm%C_eta
    C_sponge_inv = 1.0_rk / params_acm%C_sponge
    C_0          = params_acm%C_0

    if (params_acm%dim==3) then
        if (.not. allocated(mask_color)) allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g))
        if (.not. allocated(div)) allocate(div(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g))
        dV = dx(1)*dx(2)*dx(3)
    else
        if (.not. allocated(mask_color)) allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
        if (.not. allocated(div)) allocate(div(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
        dV = dx(1)*dx(2)
    endif

    ! save some computing time by using a logical and not comparing every time
    is_insect = .false.
    use_color = .false.
    if (params_acm%geometry == "Insect") is_insect = .true.
    if (params_acm%geometry == "two-moving-cylinders") use_color = .true.



    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.
        ! performs initializations in the RHS module, such as resetting integrals
        params_acm%mean_flow = 0.0_rk
        params_acm%force_color = 0.0_rk
        params_acm%meanflow_channel = 0.0_rk
        params_acm%moment_color = 0.0_rk
        params_acm%e_kin = 0.0_rk
        params_acm%enstrophy = 0.0_rk
        params_acm%helicity = 0.0_rk
        params_acm%mask_volume = 0.0_rk
        params_acm%sponge_volume = 0.0_rk
        params_acm%u_residual = 0.0_rk
        params_acm%div_max = 0.0_rk
        params_acm%div_min = 0.0_rk
        params_acm%penal_power = 0.0_rk
        params_acm%scalar_removal = 0.0_rk
        params_acm%dissipation = 0.0_rk
        params_acm%umag = 0.0_rk
        params_acm%ACM_energy = 0.0_rk

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
                write(*,'("Statistics: maxval in u(:,:,:,",i2,") = ", es10.4, ", Block with origin", 3(1x,es9.2), " and dx", 3(1x,es8.2))') k, maxval(abs(u(:,:,:,k))), x0, dx

                call dump_block_fancy(u(:,:,:,k:k), "block_ACM_diverged_statistics.txt", Bs, g)

                ! done by all ranks but well I hope the cluster can take one for the team
                ! This (empty) file is for scripting purposes on the supercomputers.
                open (77, file='ACM_diverged', status='replace')
                close(77)

                call abort(0409201934,"ACM fail: very very large values in state vector.")
            endif
        enddo

        ! tmp values for computing the current block only
        meanflow_channel_block = 0.0_rk
        meanflow_block = 0.0_rk
        force_block = 0.0_rk
        moment_block = 0.0_rk
        residual_block = 0.0_rk
        ekin_block = 0.0_rk
        tmp_volume = 0.0_rk
        tmp_volume2 = 0.0_rk
        penal_power_block = 0.0_rk
        scalar_removal_block = 0.0_rk
        ACM_energy_block = 0.0_rk

        ! default for moment computation is the centre point of the object (=lever)
        x0_moment(1,:) = params_acm%x_cntr(1)
        x0_moment(2,:) = params_acm%x_cntr(2)
        x0_moment(3,:) = params_acm%x_cntr(3)
        ! in insects, other levers are used
        ! some preparations that we do not want to compute each time within the loop
        if (is_insect) then
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
        endif

        if (params_acm%dim == 2) then
            ! --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D ---
            ! note in 2D case, uz is ignored, so we pass p=u(:,:,:,3) just for fun.
            call divergence( u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, div)

            ! mask divergence inside the solid body
            where (mask(:,:,:,1)>0.0_rk)
                div = 0.00_rk
            end where

            ! most integral quantities can be computed with block-wise function calls
            ! compute mean flow for output in statistics
            meanflow_block(1) = sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1))
            meanflow_block(2) = sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 2))

            ! kinetic energy
            ekin_block = 0.5_rk*sum( u(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1:2)**2 )

            ! acm energy (kinetic energy + artificial pressure energy, 0.5*p**2/C0**2). As C0->\infty, 
            ! the latter vanishes, just as in the incompressible limit, where we have only kinetic energy.
            ACM_energy_block = 0.5_rk*sum( u(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 3)**2 )/C_0**2 + ekin_block

            ! square of maximum of velocity in the field
            params_acm%umag = max( params_acm%umag, maxval(sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g,:,1:2)**2, dim=4)))

            ! maximum/min divergence in velocity field
            params_acm%div_max = max( params_acm%div_max, maxval(div(g+1:Bs(1)+g, g+1:Bs(2)+g, 1)))
            params_acm%div_min = min( params_acm%div_min, minval(div(g+1:Bs(1)+g, g+1:Bs(2)+g, 1)))

            ! volume of mask (useful to see if it is properly generated)
            tmp_volume = sum(mask(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1))
            ! volume of sponge
            tmp_volume2 = sum(mask(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 6))

            ! some computations need point-wise loops
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    color = int( mask(ix, iy, 1, 5), kind=2 )

                    chi = mask(ix,iy,1,1) * C_eta_inv
                    usx = mask(ix,iy,1,2)
                    usy = mask(ix,iy,1,3)

                    ! forces acting on body
                    force_block(1, color) = force_block(1, color) + (u(ix,iy,1,1)-usx)*chi
                    force_block(2, color) = force_block(2, color) + (u(ix,iy,1,2)-usy)*chi

                    ! penalization power (input from solid motion), see Engels et al. J. Comput. Phys. 2015
                    ! energy input from solid:
                    penal_power_block(1) = penal_power_block(1) + (usx*(u(ix,iy,1,1)-usx) + usy*(u(ix,iy,1,2)-usy))*chi
                    ! dissipation inside solid:
                    penal_power_block(2) = penal_power_block(2) + ( (u(ix,iy,1,1)-usx)**2 + (u(ix,iy,1,2)-usy)**2) *chi
                    ! sponge input:
                    penal_power_block(3) = penal_power_block(3) + ( u(ix,iy,1,1)*(u(ix,iy,1,1)-params_acm%u_mean_set(1)) + u(ix,iy,1,2)*(u(ix,iy,1,2)-params_acm%u_mean_set(2))) *mask(ix,iy,1,6)


                    ! residual velocity in the solid domain
                    residual_block(1) = max( residual_block(1), (u(ix,iy,1,1)-usx) * mask(ix,iy,1,1))
                    residual_block(2) = max( residual_block(2), (u(ix,iy,1,2)-usy) * mask(ix,iy,1,1))
                enddo
            enddo

            ! if the scalar BC is Dirichlet, then the solid absorbs some scalar, and it makes
            ! sense to keep track of this. however, note that with Neumann BC, that makes no
            ! sense (zero flux is imposed and the solution for the scalar inside the solid is arbitrary)
            if (params_acm%use_passive_scalar .and. params_acm%scalar_BC_type == "dirichlet") then
                ! should we ever use more than 1 scalar seriously, this has to be adopted
                ! because it uses only the 1st one (:,:,1,4)
                ! assumes *homogeneous* dirichlet condition
                scalar_removal_block = scalar_removal_block + sum(mask(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1) * u(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 4)) * C_eta_inv
            endif
        else
            ! --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D ---
            ! compute divergence on this block
            call divergence( u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, div)

            ! mask divergence inside the solid body
            where (mask(:,:,:,1)>0.0_rk)
                div = 0.00_rk
            end where

            ! block-wise function calls for quantities
            ! compute mean flow for output in statistics
            meanflow_block(1) = sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1))
            meanflow_block(2) = sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 2))
            meanflow_block(3) = sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 3))

            ! relevant in 3D only
            if (params_acm%use_channel_forcing) then
                do iy = g+1, Bs(2)+g
                    y = x0(2) + dble(iy-(g+1)) * dx(2)

                    ! exclude channel walls
                    if (( y>params_acm%h_channel).and.(y<params_acm%domain_size(2)-params_acm%h_channel)) then
                        meanflow_channel_block(1) = meanflow_channel_block(1) + sum(u(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 1))
                        meanflow_channel_block(2) = meanflow_channel_block(2) + sum(u(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 2))
                        meanflow_channel_block(3) = meanflow_channel_block(3) + sum(u(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 3))
                    endif
                enddo
            endif

            ! kinetic energy
            ekin_block = 0.5_rk*sum( u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:3)**2 )

            ! acm energy (kinetic energy + artificial pressure energy, 0.5*p**2/C0**2). As C0->\infty, 
            ! the latter vanishes, just as in the incompressible limit, where we have only kinetic energy.
            ACM_energy_block = 0.5_rk*sum( u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 4)**2 )/C_0**2 + ekin_block

            ! square of maximum of velocity in the field
            params_acm%umag = max( params_acm%umag, maxval(sum(u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:3)**2,dim=4)))

            ! maximum/min divergence in velocity field
            params_acm%div_max = max( params_acm%div_max, maxval(div(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g)))
            params_acm%div_min = min( params_acm%div_min, minval(div(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g)))

            ! volume of mask (useful to see if it is properly generated)
            ! NOTE: in wabbit, mask is really the mask: it is not divided by C_eta yet.
            tmp_volume = sum(mask(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1))
            ! volume of sponge
            tmp_volume2 = sum(mask(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 6))

            ! point-wise loop for quantities
            do iz = g+1, Bs(3)+g
                z = x0(3) + dble(iz-(g+1)) * dx(3)
                do iy = g+1, Bs(2)+g
                    y = x0(2) + dble(iy-(g+1)) * dx(2)
                    do ix = g+1, Bs(1)+g
                        x = x0(1) + dble(ix-(g+1)) * dx(1)

                        chi        = mask(ix,iy,iz,1) * C_eta_inv
                        usx        = mask(ix,iy,iz,2)
                        usy        = mask(ix,iy,iz,3)
                        usz        = mask(ix,iy,iz,4)
                        color      = int( mask(ix, iy, iz, 5), kind=2 )
                        chi_sponge = mask(ix,iy,iz,6) * C_sponge_inv

                        ! penalization term
                        penal = -mask(ix,iy,iz,1) * (u(ix,iy,iz,1:3) - mask(ix,iy,iz,2:4)) * C_eta_inv

                        ! residual velocity in the solid domain
                        residual_block(1) = max( residual_block(1), (u(ix,iy,iz,1)-mask(ix,iy,iz,2))*mask(ix,iy,iz,1) )
                        residual_block(2) = max( residual_block(2), (u(ix,iy,iz,2)-mask(ix,iy,iz,3))*mask(ix,iy,iz,1) )
                        residual_block(3) = max( residual_block(3), (u(ix,iy,iz,3)-mask(ix,iy,iz,4))*mask(ix,iy,iz,1) )

                        ! penalization power (input from solid motion), see Engels et al. J. Comput. Phys. 2015
                        ! NOTE: as this is qty used to show the global energy budget, it includes ALL mask colors
                        ! 1/ energy input from solid:
                        penal_power_block(1) = penal_power_block(1) + (usx*(u(ix,iy,iz,1)-usx) + usy*(u(ix,iy,iz,2)-usy) + usz*(u(ix,iy,iz,3)-usz))*chi
                        ! 2/ dissipation inside solid:
                        penal_power_block(2) = penal_power_block(2) + ((u(ix,iy,iz,1)-usx)**2 + (u(ix,iy,iz,2)-usy)**2 + (u(ix,iy,iz,3)-usz)**2)*chi
                        ! 3/ sponge input:
                        ! attention division by C_sponge C0**2 because of the way the energy eqn is derived (it is divided by C_0**2
                        ! such that dp/dt 1/C0**2 is on the right side)
                        penal_power_block(3) = penal_power_block(3) + ( u(ix,iy,iz,1)*(u(ix,iy,iz,1)-params_acm%u_mean_set(1)) &
                                                                      + u(ix,iy,iz,2)*(u(ix,iy,iz,2)-params_acm%u_mean_set(2)) &
                                                                      + u(ix,iy,iz,3)*(u(ix,iy,iz,3)-params_acm%u_mean_set(3)) &
                                                                      + u(ix,iy,iz,4)*(u(ix,iy,iz,4)-0.0_rk)/C_0**2 ) *chi_sponge 

                        ! forces acting on this color
                        force_block(1:3, color) = force_block(1:3, color) - penal

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
                        if ((is_insect).and.(color>0_2 .and. color < 6_2)) then ! Insect has colors (1,2,3,4,5) at most
                            ! in the zeroth color, we compute the total moment for the whole
                            ! insect wrt the center point (body+wings)
                            x_lev(1:3) = (/x, y, z/) - Insect%xc_body_g(1:3)
                            moment_block(:,0)  = moment_block(:,0) - cross(x_lev, penal)
                        endif

                    enddo
                enddo
            enddo

            ! if the scalar BC is Dirichlet, then the solid absorbs some scalar, and it makes
            ! sense to keep track of this. however, note that with Neumann BC, that makes no
            ! sense (zero flux is imposed and the solution for the scalar inside the solid is arbitrary)
            if (params_acm%use_passive_scalar .and. params_acm%scalar_BC_type == "dirichlet") then
                ! should we ever use more than 1 scalar seriously, this has to be adopted
                ! because it uses only the 1st one (:,:,:,5)
                ! assumes *homogeneous* dirichlet condition
                scalar_removal_block = sum(mask(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1) * u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 5)) * C_eta_inv
            endif

        endif

        ! we just computed the values on the current block, which we now add to the
        ! existing blocks in the variables (recall normalization by dV)
        params_acm%u_residual     = params_acm%u_residual     + residual_block * dV
        params_acm%mean_flow      = params_acm%mean_flow      + meanflow_block * dV
        params_acm%meanflow_channel = params_acm%meanflow_channel + meanflow_channel_block * dV
        params_acm%mask_volume    = params_acm%mask_volume    + tmp_volume * dV
        params_acm%sponge_volume  = params_acm%sponge_volume  + tmp_volume2 * dV
        params_acm%force_color    = params_acm%force_color    + force_block * dV
        params_acm%moment_color   = params_acm%moment_color   + moment_block * dV
        params_acm%e_kin          = params_acm%e_kin          + ekin_block * dV
        params_acm%penal_power    = params_acm%penal_power    + penal_power_block * dV
        params_acm%scalar_removal = params_acm%scalar_removal + scalar_removal_block * dV
        params_acm%ACM_energy     = params_acm%ACM_energy     + ACM_energy_block * dV

        !-------------------------------------------------------------------------
        ! compute enstrophy in the whole domain (including penalized regions)
        ! note in 2D case, uz is ignored, so we pass p=u(:,:,:,3) just for fun.
        call compute_vorticity(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, work(:,:,:,:))

        if (params_acm%dim == 2) then
            params_acm%enstrophy = params_acm%enstrophy + 0.5_rk*sum(work(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1)**2)*dV
        else
            params_acm%enstrophy = params_acm%enstrophy + 0.5_rk*sum(work(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:3)**2)*dV
            params_acm%helicity = params_acm%helicity + 0.5_rk*sum(work(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:3) * u(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1:3))*dV
        end if

        call dissipation_ACM_block(Bs, g, dx, u(:,:,:,1:3), dissipation_block)
        params_acm%dissipation = params_acm%dissipation + dissipation_block

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.

        ! volume of mask (useful to see if it is properly generated)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%mask_volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        ! volume of sponge
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%sponge_volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        ! mean flow (in entire domain)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        params_acm%mean_flow = params_acm%mean_flow / product(params_acm%domain_size(1:params_acm%dim))

        if (params_acm%use_channel_forcing) then
            ! mean flow but only in fluid domain
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%meanflow_channel, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

            V_channel = params_acm%domain_size(1)*params_acm%domain_size(3)*(params_acm%domain_size(2)-2.0_rk*params_acm%h_channel)
            params_acm%meanflow_channel = params_acm%meanflow_channel / V_channel
        endif

        ! force & moment
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%force_color, size(params_acm%force_color), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%moment_color, size(params_acm%moment_color), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        ! residual velocity in solid domain
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%u_residual, 3, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%penal_power, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! scalar removal
        if (params_acm%use_passive_scalar .and. params_acm%scalar_BC_type == "dirichlet") then
            call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%scalar_removal, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        endif

        !-------------------------------------------------------------------------
        ! kinetic energy
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%ACM_energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! divergence
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%div_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, WABBIT_COMM, mpierr)

        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%div_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic enstrophy
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%helicity, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%dissipation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        ! minium spacing of current grid, not smallest possible one.
        call MPI_ALLREDUCE(MPI_IN_PLACE, dx_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_acm%umag, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

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

            if (params_acm%use_channel_forcing) then
                call append_t_file( 'meanflow_channel.t', (/time, params_acm%meanflow_channel   /) )
            endif

            call append_t_file( 'div.t', (/time, params_acm%div_max, params_acm%div_min/) )
            if (params_acm%penalization .or. params_acm%use_sponge) then
                ! total force (excluding zero color)
                call append_t_file( 'forces.t', (/time, sum(params_acm%force_color(1,1:ncolors)), &
                sum(params_acm%force_color(2,1:ncolors)), sum(params_acm%force_color(3,1:ncolors)) /) )

                ! forces/moment for individual colors.
                ! This is what we should have done in the first place. For the insects below,
                ! some colors are (also) stored to different names. That's unfortunate but not a big deal.
                do color = 0, ncolors
                    ! save force for each color
                    write(fname,'("forces_color",i1,".t")') color
                    call append_t_file( fname, (/time, params_acm%force_color(:,color)/) )

                    ! save moment for each color
                    write(fname,'("moments_color",i1,".t")') color
                    call append_t_file( fname, (/time, params_acm%moment_color(:,color)/) )
                enddo

                if (is_insect) then
                    call append_t_file( 'aero_power.t', (/time, apowtotal, ipowtotal/) )

                    ! total moment w.r.t body center is computed in zeroth color slot:
                    color = 0_2
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

                call append_t_file( 'mask_volume.t', (/time, params_acm%mask_volume, params_acm%sponge_volume/) )
                call append_t_file( 'penal_power.t', (/time, params_acm%penal_power/) )
                call append_t_file( 'u_residual.t', (/time, params_acm%u_residual/) )    
            endif

            if (params_acm%use_passive_scalar .and. params_acm%scalar_BC_type == "dirichlet") then
                call append_t_file( 'scalar_removal.t', (/time, params_acm%scalar_removal/) )
            endif
            call append_t_file( 'e_kin.t', (/time, params_acm%e_kin, params_acm%ACM_energy/) )
            call append_t_file( 'enstrophy.t', (/time, params_acm%enstrophy/) )
            if (params_acm%dim == 3) then
                call append_t_file( 'helicity.t', (/time, params_acm%helicity/) )
            endif
            call append_t_file( 'dissipation.t', (/time, params_acm%dissipation/) )

            ! turbulent statistics - these are normed by the volume!
            if (params_acm%nu*params_acm%enstrophy > 0.0_rk) then
                dissipation = 2*params_acm%nu*params_acm%enstrophy/product(params_acm%domain_size(1:params_acm%dim))
                u_RMS = sqrt(2*params_acm%e_kin/product(params_acm%domain_size(1:params_acm%dim))/3)
                call append_t_file( 'turbulent_statistics.t', (/time, dissipation, params_acm%e_kin/product(params_acm%domain_size(1:params_acm%dim)), u_RMS, &
                    (params_acm%nu**3.0_rk / dissipation)**0.25_rk, sqrt(params_acm%nu/dissipation), (params_acm%nu*dissipation)**0.25_rk, &
                    sqrt(15.0_rk*params_acm%nu*u_RMS**2/dissipation), sqrt(15.0_rk*params_acm%nu*u_RMS**2/dissipation)*u_RMS/params_acm%nu/))
            endif

            ! this file is to simply keep track of simulations, should they be restarted with different parameters.
            ! We just keep track of the most important and most frequently changed parameters.
            call append_t_file( 'parameters.t', (/time, params_acm%C_0, params_acm%C_eta, params_acm%C_sponge, params_acm%nu/) )
        end if

    case default
        call abort(7772,"the STATISTICS wrapper requests a stage this physics module cannot handle.")

    end select


end subroutine STATISTICS_ACM
