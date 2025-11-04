!-----------------------------------------------------------------------------
! main level wrapper to compute local statistics over time, like average or minmax of specific quantities
! these are then also adapted over time
!-----------------------------------------------------------------------------
subroutine TIME_STATISTICS_ACM( time, dt, time_start, u, g, x0, dx, work, mask )
    use module_operators

    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time, dt, time_start

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

    integer(kind=ik) :: iB, ix, iy, iz, Bs(3), i_ts, N_offset, mean_idx1, mean_idx2, j
    real(kind=rk) :: time_diff

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    N_offset = params_acm%dim + 1 + params_acm%N_scalars
    time_diff = time-time_start

    do i_ts = 1, params_acm%N_time_statistics
        ! if time is right at time_start, we do nothing, we already set the values to zero or load them in earlier
        if ( time - time_start < 1.0e-12_rk ) then
            ! u(:,:,:,N_offset + i_ts) = 0.0_rk
            cycle
        end if

        select case (trim(params_acm%time_statistics_names(i_ts)))
        case ("ux-avg")
            ! compute the average of ux over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,1)
        case ("ux-var")
            ! compute the variance of ux over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "ux-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'ux-avg' somewhere after 'ux-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,1))
        case ("uy-avg")
            ! compute the average of uy over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,2)
        case ("uy-var")
            ! compute the variance of uy over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "uy-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'uy-avg' somewhere after 'uy-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,2) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,2) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,2))
        case ("uz-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uz not available in 2D simulations.")
            ! compute the average of uz over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,3)
        case ("uz-var")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uz not available in 2D simulations.")
            ! compute the variance of uz over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "uz-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'uz-avg' somewhere after 'uz-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,3) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,3) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,3))
        case ("p-avg")
            ! compute the average of p over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,params_acm%dim+1)
        case ("p-var")
            ! compute the variance of p over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "p-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'p-avg' or 'p-mean' somewhere after 'p-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,params_acm%dim+1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,params_acm%dim+1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,params_acm%dim+1))
        
        case ("umag-avg")
            ! compute the average of velocity magnitude over time
            work (:,:,:,1) = u(:,:,:,1)**2 + u(:,:,:,2)**2
            if (params_acm%dim == 3) work(:,:,:,1) = work(:,:,:,1) + u(:,:,:,3)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * sqrt( work(:,:,:,1) )
        case ("umag-var")
            ! compute the variance of velocity magnitude over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "umag-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'umag-avg' or 'umag-mean' somewhere after 'umag-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            work (:,:,:,1) = u(:,:,:,1)**2 + u(:,:,:,2)**2
            if (params_acm%dim == 3) work(:,:,:,1) = work(:,:,:,1) + u(:,:,:,3)**2
            work(:,:,:,1) = sqrt( work(:,:,:,1) )  ! current umag value
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        case ("divu-avg", "div-avg")
            ! compute the average of div(u) over time
            call compute_divergence(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("divu-var", "div-var")
            ! compute the variance of div(u) over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "divu-avg", mean_idx1)
            call find_single_mean_index(i_ts, "div-avg", mean_idx2)
            mean_idx1 = max(mean_idx1, mean_idx2)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'divu-avg' or 'div-avg' somewhere after 'divu-var' or 'div-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            call compute_divergence(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))
        case ("divu-minmax", "div-minmax")
            ! compute the minmax of div(u) over time
            call compute_divergence(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            do iz = merge(1,g+1,params_acm%dim==2), merge(1,Bs(3)+g,params_acm%dim==2)
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        if (abs(work(ix,iy,iz,1)) > u(ix,iy,iz,N_offset + i_ts)) u(ix,iy,iz,N_offset + i_ts) = abs(work(ix,iy,iz,1))
                    end do
                end do
            end do
        
        case ("vort-avg", "vor-avg")
            if (params_acm%dim == 3) call abort(250916, "[TIME_STATISTICS_ACM]: vorticity is not a scalar in 3D simulations - use vorx, vory, vorz or vorabs instead.")
            ! compute the average of vorticity magnitude over time
            call compute_vorticity(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,:))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("vort-var", "vor-var")
            if (params_acm%dim == 3) call abort(251016, "[TIME_STATISTICS_ACM]: vorticity is not a scalar in 3D simulations - use vorx, vory, vorz or vorabs instead.")
            ! compute the variance of vorticity magnitude over time using Welford's method
            ! This needs the computation of the average in the variable afterwards to work
            ! NOTE: Using manual search loop instead of find_single_mean_index helper because
            ! this case has 4 possible mean names (vort-avg, vor-avg, vort-mean, vor-mean)
            ! and the helper function only supports 2 name options
            call find_single_mean_index(i_ts, "vort-avg", mean_idx1)
            call find_single_mean_index(i_ts, "vor-avg", mean_idx2)
            mean_idx1 = max(mean_idx1, mean_idx2)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'vort-avg', 'vor-avg' somewhere after 'vort-var' or 'vor-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            call compute_vorticity(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,:))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))
        
        case ("vorabs-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: vorticity is a scalar in 2D simulations - use vor instead.")
            ! compute the average of vorticity magnitude over time
            call compute_vorticity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("vorabs-var")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: vorticity is a scalar in 2D simulations - use vor instead.")
            ! compute the variance of vorticity magnitude over time using Welford's method
            ! This needs the computation of the average in the variable afterwards to work
            call find_single_mean_index(i_ts, "vorabs-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'vorabs-avg' somewhere after 'vorabs-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            call compute_vorticity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        ! Helicity magnitude cases
        case ("helabs-avg")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: helabs not available in 2D simulations.")
            ! compute the average of helicity magnitude over time
            call compute_helicity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("helabs-var")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: helabs not available in 2D simulations.")
            ! compute the variance of helicity magnitude over time using Welford's method
            call find_single_mean_index(i_ts, "helabs-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251021_ik, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'helabs-avg' somewhere after 'helabs-var' in time_statistics_names.")
            call compute_helicity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        case ("dissipation-avg")
            ! compute the average of dissipation over time
            call compute_dissipation(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("dissipation-var")
            ! compute the variance of dissipation over time using Welford's method
            ! This needs the computation of the average in the variable afterwards to work
            call find_single_mean_index(i_ts, "dissipation-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'dissipation-avg' somewhere after 'dissipation-var' in time_statistics_names.")
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            call compute_dissipation(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        ! Velocity derivative cases
        case ("uxdx-avg")
            ! compute the average of ∂ux/∂x over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdy-avg")
            ! compute the average of ∂ux/∂y over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdz-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uxdz not available in 2D simulations.")
            ! compute the average of ∂ux/∂z over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydx-avg")
            ! compute the average of ∂uy/∂x over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydy-avg")
            ! compute the average of ∂uy/∂y over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydz-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uydz not available in 2D simulations.")
            ! compute the average of ∂uy/∂z over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdx-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdx not available in 2D simulations.")
            ! compute the average of ∂uz/∂x over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdy-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdy not available in 2D simulations.")
            ! compute the average of ∂uz/∂y over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdz-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdz not available in 2D simulations.")
            ! compute the average of ∂uz/∂z over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        
        ! Squared velocity derivative cases
        case ("uxdx2-avg")
            ! compute the average of (∂ux/∂x)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdy2-avg")
            ! compute the average of (∂ux/∂y)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdz2-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uxdz2 not available in 2D simulations.")
            ! compute the average of (∂ux/∂z)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydx2-avg")
            ! compute the average of (∂uy/∂x)² over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydy2-avg")
            ! compute the average of (∂uy/∂y)² over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydz2-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uydz2 not available in 2D simulations.")
            ! compute the average of (∂uy/∂z)² over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdx2-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdx2 not available in 2D simulations.")
            ! compute the average of (∂uz/∂x)² over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdy2-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdy2 not available in 2D simulations.")
            ! compute the average of (∂uz/∂y)² over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdz2-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdz2 not available in 2D simulations.")
            ! compute the average of (∂uz/∂z)² over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        
        ! Squared velocity component cases
        case ("ux-2-avg")
            ! compute the average of ux² over time
            work(:,:,:,1) = u(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uy-2-avg")
            ! compute the average of uy² over time
            work(:,:,:,1) = u(:,:,:,2)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uz-2-avg")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uz-2 not available in 2D simulations.")
            ! compute the average of uz² over time
            work(:,:,:,1) = u(:,:,:,3)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)

        ! Cross-product velocity component averages
        case ("uxuy-avg")
            ! compute the average of ux*uy over time
            work(:,:,:,1) = u(:,:,:,1) * u(:,:,:,2)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxuz-avg")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uxuz-avg not available in 2D simulations.")
            ! compute the average of ux*uz over time
            work(:,:,:,1) = u(:,:,:,1) * u(:,:,:,3)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uyuz-avg")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uyuz-avg not available in 2D simulations.")
            ! compute the average of uy*uz over time
            work(:,:,:,1) = u(:,:,:,2) * u(:,:,:,3)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)

        ! Velocity-pressure cross-product averages
        case ("uxp-avg")
            ! compute the average of ux*p over time
            work(:,:,:,1) = u(:,:,:,1) * u(:,:,:,params_acm%dim+1)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uyp-avg")
            ! compute the average of uy*p over time
            work(:,:,:,1) = u(:,:,:,2) * u(:,:,:,params_acm%dim+1)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzp-avg")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uzp-avg not available in 2D simulations.")
            ! compute the average of uz*p over time
            work(:,:,:,1) = u(:,:,:,3) * u(:,:,:,params_acm%dim+1)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)

        ! Covariance cases using Welford's online algorithm
        case ("uxuy-cov")
            ! compute the covariance of ux and uy over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "ux-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'ux-avg' somewhere after 'uxuy-cov' in time_statistics_names.")
            call find_single_mean_index(i_ts, "uy-avg", mean_idx2)
            if (mean_idx2 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uy-avg' somewhere after 'uxuy-cov' in time_statistics_names.")
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            ! mean_x_new = (time_diff-dt)/time_diff*mean_x_old + dt/time_diff*x
            ! mean_y_new = (time_diff-dt)/time_diff*mean_y_old + dt/time_diff*y
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,2) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,2))
        case ("uxuz-cov")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uxuz-cov not available in 2D simulations.")
            ! compute the covariance of ux and uz over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "ux-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'ux-avg' somewhere after 'uxuz-cov' in time_statistics_names.")
            call find_single_mean_index(i_ts, "uz-avg", mean_idx2)
            if (mean_idx2 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uz-avg' somewhere after 'uxuz-cov' in time_statistics_names.")
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,3) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,3))
        case ("uyuz-cov")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uyuz-cov not available in 2D simulations.")
            ! compute the covariance of uy and uz over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "uy-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uy-avg' somewhere after 'uyuz-cov' in time_statistics_names.")
            call find_single_mean_index(i_ts, "uz-avg", mean_idx2)
            if (mean_idx2 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uz-avg' somewhere after 'uyuz-cov' in time_statistics_names.")
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,2) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,3) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,3))

        ! Velocity-pressure covariances
        case ("uxp-cov")
            ! compute the covariance of ux and p over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "ux-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'ux-avg' somewhere after 'uxp-cov' in time_statistics_names.")
            call find_single_mean_index(i_ts, "p-avg", mean_idx2)
            if (mean_idx2 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'p-avg' somewhere after 'uxp-cov' in time_statistics_names.")
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,params_acm%dim+1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,params_acm%dim+1))
        case ("uyp-cov")
            ! compute the covariance of uy and p over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "uy-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uy-avg' somewhere after 'uyp-cov' in time_statistics_names.")
            call find_single_mean_index(i_ts, "p-avg", mean_idx2)
            if (mean_idx2 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'p-avg' somewhere after 'uyp-cov' in time_statistics_names.")
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,2) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,params_acm%dim+1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,params_acm%dim+1))
        case ("uzp-cov")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uzp-cov not available in 2D simulations.")
            ! compute the covariance of uz and p over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "uz-avg", mean_idx1)
            if (mean_idx1 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uz-avg' somewhere after 'uzp-cov' in time_statistics_names.")
            call find_single_mean_index(i_ts, "p-avg", mean_idx2)
            if (mean_idx2 == -1) call abort(251016_ik, "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'p-avg' somewhere after 'uzp-cov' in time_statistics_names.")
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,3) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,params_acm%dim+1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,params_acm%dim+1))

        case default
            call abort(2153000, "[TIME_STATISTICS_ACM]: time_statistics_name "// &
                 trim(params_acm%time_statistics_names(i_ts))//" not recognized.")
        end select
    end do

end subroutine TIME_STATISTICS_ACM

!-----------------------------------------------------------------------------
! Helper subroutine to find a single mean value index after a given position
!
! This helper function was created to eliminate code duplication across variance 
! and covariance calculations. Previously, each variance/covariance case contained
! identical search loops to find the corresponding mean values. This helper:
!
! 1. Reduces code duplication by centralizing the search logic
! 2. Provides consistent validation behavior across all cases
! 3. Enables flexible positioning of mean values anywhere somewhere after the variance/covariance
! 5. Provides consistent error handling with custom messages and error codes
!
!-----------------------------------------------------------------------------
subroutine find_single_mean_index(current_idx, mean_name, mean_idx)
    implicit none
    
    integer(kind=ik), intent(in) :: current_idx
    character(len=*), intent(in) :: mean_name  ! One possible name (e.g., "ux-avg")
    integer(kind=ik), intent(out) :: mean_idx
    
    integer(kind=ik) :: j
    
    mean_idx = -1  ! Initialize to not found
    
    ! Search for mean name after current index
    do j = current_idx + 1, params_acm%N_time_statistics
        if (trim(params_acm%time_statistics_names(j)) == trim(mean_name)) then
            mean_idx = j
            exit
        end if
    end do
    
end subroutine find_single_mean_index