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
        case ("ux-int")
            ! compute the integral of ux over time
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* u(:,:,:,1)
        case ("ux-avg", "ux-mean")
            ! compute the average of ux over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,1)
        case ("ux-var")
            ! compute the variance of ux over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "ux-avg", "ux-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'ux-avg' or 'ux-mean' somewhere after 'ux-var' in time_statistics_names.", 251016_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,1))
        
        case ("uy-int")
            ! compute the integral of uy over time
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* u(:,:,:,2)
        case ("uy-avg", "uy-mean")
            ! compute the average of uy over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,2)
        case ("uy-var")
            ! compute the variance of uy over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "uy-avg", "uy-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'uy-avg' or 'uy-mean' somewhere after 'uy-var' in time_statistics_names.", 251016_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,2) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,2) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,2))
        
        case ("uz-int")
            if (params_acm%dim == 2) call abort(250915, "[TIME_STATISTICS_ACM]: uz not available in 2D simulations.")
            ! compute the integral of uz over time
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* u(:,:,:,3)
        case ("uz-avg", "uz-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uz not available in 2D simulations.")
            ! compute the average of uz over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,3)
        case ("uz-var")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uz not available in 2D simulations.")
            ! compute the variance of uz over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "uz-avg", "uz-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'uz-avg' or 'uz-mean' somewhere after 'uz-var' in time_statistics_names.", 251016_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,3) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,3) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,3))
        
        case ("p-int")
            ! compute the integral of p over time
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* u(:,:,:,params_acm%dim+1)
        case ("p-avg", "p-mean")
            ! compute the average of p over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * u(:,:,:,params_acm%dim+1)
        case ("p-var")
            ! compute the variance of p over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "p-avg", "p-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'p-avg' or 'p-mean' somewhere after 'p-var' in time_statistics_names.", 251016_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,params_acm%dim+1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,params_acm%dim+1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,params_acm%dim+1))
        
        case ("umag-int")
            ! compute the integral of velocity magnitude over time
            work(:,:,:,1) = u(:,:,:,1)**2 + u(:,:,:,2)**2
            if (params_acm%dim == 3) work(:,:,:,1) = work(:,:,:,1) + u(:,:,:,3)**2
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* sqrt( work(:,:,:,1) )
        case ("umag-avg", "umag-mean")
            ! compute the average of velocity magnitude over time
            work (:,:,:,1) = u(:,:,:,1)**2 + u(:,:,:,2)**2
            if (params_acm%dim == 3) work(:,:,:,1) = work(:,:,:,1) + u(:,:,:,3)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * sqrt( work(:,:,:,1) )
        case ("umag-var")
            ! compute the variance of velocity magnitude over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, "umag-avg", "umag-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'umag-avg' or 'umag-mean' somewhere after 'umag-var' in time_statistics_names.", 251016_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            work (:,:,:,1) = u(:,:,:,1)**2 + u(:,:,:,2)**2
            if (params_acm%dim == 3) work(:,:,:,1) = work(:,:,:,1) + u(:,:,:,3)**2
            work(:,:,:,1) = sqrt( work(:,:,:,1) )  ! current umag value
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        case ("divu-int", "div-int")
            ! compute the integral, average or minmax of div(u) over time
            call compute_divergence(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("divu-avg", "div-avg", "divu-mean", "div-mean")
            ! compute the average of div(u) over time
            call compute_divergence(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("divu-var", "div-var")
            ! compute the variance of div(u) over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            ! NOTE: Using manual search loop instead of find_single_mean_index helper because
            ! this case has 4 possible mean names (divu-avg, div-avg, divu-mean, div-mean)
            ! and the helper function only supports 2 name options
            mean_idx1 = -1
            do j = i_ts + 1, params_acm%N_time_statistics
                if (trim(params_acm%time_statistics_names(j)) == "divu-avg" .or. &
                    trim(params_acm%time_statistics_names(j)) == "div-avg" .or. &
                    trim(params_acm%time_statistics_names(j)) == "divu-mean" .or. &
                    trim(params_acm%time_statistics_names(j)) == "div-mean") then
                    mean_idx1 = j
                    exit
                end if
            end do
            if (mean_idx1 == -1) then
                call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'divu-avg', 'div-avg', 'divu-mean', or 'div-mean' somewhere after 'divu-var' or 'div-var' in time_statistics_names.")
            end if
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
        
        case ("vort-int", "vor-int")
            if (params_acm%dim == 3) call abort(250916, "[TIME_STATISTICS_ACM]: vorticity is not a scalar in 3D simulations - use vorx, vory, vorz or vorabs instead.")
            ! compute the integral of vorticity magnitude over time
            call compute_vorticity(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,:))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("vort-avg", "vor-avg", "vort-mean", "vor-mean")
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
            mean_idx1 = -1
            do j = i_ts + 1, params_acm%N_time_statistics
                if (trim(params_acm%time_statistics_names(j)) == "vort-avg" .or. &
                    trim(params_acm%time_statistics_names(j)) == "vor-avg" .or. &
                    trim(params_acm%time_statistics_names(j)) == "vort-mean" .or. &
                    trim(params_acm%time_statistics_names(j)) == "vor-mean") then
                    mean_idx1 = j
                    exit
                end if
            end do
            if (mean_idx1 == -1) then
                call abort(251016, "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'vort-avg', 'vor-avg', 'vort-mean', or 'vor-mean' AFTER 'vort-var' or 'vor-var' in time_statistics_names.")
            end if
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            call compute_vorticity(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,:))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))
        
        case ("vorabs-int")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: vorticity is a scalar in 2D simulations - use vor instead.")
            ! compute the integral of vorticity magnitude over time
            call compute_vorticity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("vorabs-avg", "vorabs-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: vorticity is a scalar in 2D simulations - use vor instead.")
            ! compute the average of vorticity magnitude over time
            call compute_vorticity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("vorabs-var")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: vorticity is a scalar in 2D simulations - use vor instead.")
            ! compute the variance of vorticity magnitude over time using Welford's method
            ! This needs the computation of the average in the variable afterwards to work
            call find_single_mean_index(i_ts, "vorabs-avg", "vorabs-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'vorabs-avg' or 'vorabs-mean' AFTER 'vorabs-var' in time_statistics_names.", 251016_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            call compute_vorticity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        ! Helicity magnitude cases
        case ("helabs-int")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: helabs not available in 2D simulations.")
            ! compute the integral of helicity magnitude over time
            call compute_helicity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("helabs-avg", "helabs-mean")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: helabs not available in 2D simulations.")
            ! compute the average of helicity magnitude over time
            call compute_helicity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("helabs-var")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: helabs not available in 2D simulations.")
            ! compute the variance of helicity magnitude over time using Welford's method
            call find_single_mean_index(i_ts, "helabs-avg", "helabs-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'helabs-avg' or 'helabs-mean' AFTER 'helabs-var' in time_statistics_names.", 251021_ik)
            call compute_helicity_abs(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        case ("dissipation-int")
            ! compute the integral of dissipation over time
            call compute_dissipation(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("dissipation-avg", "dissipation-mean")
            ! compute the average of dissipation over time
            call compute_dissipation(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("dissipation-var")
            ! compute the variance of dissipation over time using Welford's method
            ! This needs the computation of the average in the variable afterwards to work
            call find_single_mean_index(i_ts, "dissipation-avg", "dissipation-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the variance together with the mean value. Insert 'dissipation-avg' or 'dissipation-mean' AFTER 'dissipation-var' in time_statistics_names.", 251016_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            call compute_dissipation(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (work(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (work(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*work(:,:,:,1))

        ! Velocity derivative cases
        case ("uxdx-avg", "uxdx-mean")
            ! compute the average of ∂ux/∂x over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdy-avg", "uxdy-mean")
            ! compute the average of ∂ux/∂y over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdz-avg", "uxdz-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uxdz not available in 2D simulations.")
            ! compute the average of ∂ux/∂z over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydx-avg", "uydx-mean")
            ! compute the average of ∂uy/∂x over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydy-avg", "uydy-mean")
            ! compute the average of ∂uy/∂y over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydz-avg", "uydz-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uydz not available in 2D simulations.")
            ! compute the average of ∂uy/∂z over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdx-avg", "uzdx-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdx not available in 2D simulations.")
            ! compute the average of ∂uz/∂x over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdy-avg", "uzdy-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdy not available in 2D simulations.")
            ! compute the average of ∂uz/∂y over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdz-avg", "uzdz-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdz not available in 2D simulations.")
            ! compute the average of ∂uz/∂z over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        
        ! Squared velocity derivative cases
        case ("uxdx2-avg", "uxdx2-mean")
            ! compute the average of (∂ux/∂x)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdy2-avg", "uxdy2-mean")
            ! compute the average of (∂ux/∂y)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxdz2-avg", "uxdz2-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uxdz2 not available in 2D simulations.")
            ! compute the average of (∂ux/∂z)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydx2-avg", "uydx2-mean")
            ! compute the average of (∂uy/∂x)² over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydy2-avg", "uydy2-mean")
            ! compute the average of (∂uy/∂y)² over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uydz2-avg", "uydz2-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uydz2 not available in 2D simulations.")
            ! compute the average of (∂uy/∂z)² over time
            call compute_derivative(u(:,:,:,2), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdx2-avg", "uzdx2-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdx2 not available in 2D simulations.")
            ! compute the average of (∂uz/∂x)² over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdy2-avg", "uzdy2-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdy2 not available in 2D simulations.")
            ! compute the average of (∂uz/∂y)² over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzdz2-avg", "uzdz2-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uzdz2 not available in 2D simulations.")
            ! compute the average of (∂uz/∂z)² over time
            call compute_derivative(u(:,:,:,3), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        
        ! Squared velocity component cases
        case ("ux-2-avg", "ux-2-mean")
            ! compute the average of ux² over time
            work(:,:,:,1) = u(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uy-2-avg", "uy-2-mean")
            ! compute the average of uy² over time
            work(:,:,:,1) = u(:,:,:,2)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uz-2-avg", "uz-2-mean")
            if (params_acm%dim == 2) call abort(250916, "[TIME_STATISTICS_ACM]: uz-2 not available in 2D simulations.")
            ! compute the average of uz² over time
            work(:,:,:,1) = u(:,:,:,3)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)

        ! Cross-product velocity component averages
        case ("uxuy-avg", "uxuy-mean")
            ! compute the average of ux*uy over time
            work(:,:,:,1) = u(:,:,:,1) * u(:,:,:,2)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uxuz-avg", "uxuz-mean")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uxuz-avg not available in 2D simulations.")
            ! compute the average of ux*uz over time
            work(:,:,:,1) = u(:,:,:,1) * u(:,:,:,3)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uyuz-avg", "uyuz-mean")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uyuz-avg not available in 2D simulations.")
            ! compute the average of uy*uz over time
            work(:,:,:,1) = u(:,:,:,2) * u(:,:,:,3)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)

        ! Velocity-pressure cross-product averages
        case ("uxp-avg", "uxp-mean")
            ! compute the average of ux*p over time
            work(:,:,:,1) = u(:,:,:,1) * u(:,:,:,params_acm%dim+1)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uyp-avg", "uyp-mean")
            ! compute the average of uy*p over time
            work(:,:,:,1) = u(:,:,:,2) * u(:,:,:,params_acm%dim+1)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)
        case ("uzp-avg", "uzp-mean")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uzp-avg not available in 2D simulations.")
            ! compute the average of uz*p over time
            work(:,:,:,1) = u(:,:,:,3) * u(:,:,:,params_acm%dim+1)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/time_diff * u(:,:,:,N_offset + i_ts) + dt/time_diff * work(:,:,:,1)

        ! Covariance cases using Welford's online algorithm
        case ("uxuy-cov")
            ! compute the covariance of ux and uy over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "ux-avg", "ux-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'ux-avg' or 'ux-mean' somewhere after 'uxuy-cov' in time_statistics_names.", 251016_ik)
            call find_single_mean_index(i_ts, "uy-avg", "uy-mean", mean_idx2, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uy-avg' or 'uy-mean' somewhere after 'uxuy-cov' in time_statistics_names.", 251016_ik)
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
            call find_single_mean_index(i_ts, "ux-avg", "ux-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'ux-avg' or 'ux-mean' somewhere after 'uxuz-cov' in time_statistics_names.", 251016_ik)
            call find_single_mean_index(i_ts, "uz-avg", "uz-mean", mean_idx2, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uz-avg' or 'uz-mean' somewhere after 'uxuz-cov' in time_statistics_names.", 251016_ik)
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,3) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,3))
        case ("uyuz-cov")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uyuz-cov not available in 2D simulations.")
            ! compute the covariance of uy and uz over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "uy-avg", "uy-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uy-avg' or 'uy-mean' somewhere after 'uyuz-cov' in time_statistics_names.", 251016_ik)
            call find_single_mean_index(i_ts, "uz-avg", "uz-mean", mean_idx2, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uz-avg' or 'uz-mean' somewhere after 'uyuz-cov' in time_statistics_names.", 251016_ik)
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,2) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,3) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,3))

        ! Velocity-pressure covariances
        case ("uxp-cov")
            ! compute the covariance of ux and p over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "ux-avg", "ux-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'ux-avg' or 'ux-mean' somewhere after 'uxp-cov' in time_statistics_names.", 251016_ik)
            call find_single_mean_index(i_ts, "p-avg", "p-mean", mean_idx2, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'p-avg' or 'p-mean' somewhere after 'uxp-cov' in time_statistics_names.", 251016_ik)
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,params_acm%dim+1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,params_acm%dim+1))
        case ("uyp-cov")
            ! compute the covariance of uy and p over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "uy-avg", "uy-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uy-avg' or 'uy-mean' somewhere after 'uyp-cov' in time_statistics_names.", 251016_ik)
            call find_single_mean_index(i_ts, "p-avg", "p-mean", mean_idx2, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'p-avg' or 'p-mean' somewhere after 'uyp-cov' in time_statistics_names.", 251016_ik)
            ! cov_new = (time_diff-dt)/time_diff*cov_old + dt/time_diff*(x - mean_x_old)*(y - mean_y_new)
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* &
                 (u(:,:,:,2) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,params_acm%dim+1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx2) - dt/(time_diff)*u(:,:,:,params_acm%dim+1))
        case ("uzp-cov")
            if (params_acm%dim == 2) call abort(251016, "[TIME_STATISTICS_ACM]: uzp-cov not available in 2D simulations.")
            ! compute the covariance of uz and p over time using Welford's method
            ! This needs the computation of both averages somewhere after this field
            call find_single_mean_index(i_ts, "uz-avg", "uz-mean", mean_idx1, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'uz-avg' or 'uz-mean' somewhere after 'uzp-cov' in time_statistics_names.", 251016_ik)
            call find_single_mean_index(i_ts, "p-avg", "p-mean", mean_idx2, &
                "[TIME_STATISTICS_ACM]: You need to compute the covariance together with both mean values. Insert 'p-avg' or 'p-mean' somewhere after 'uzp-cov' in time_statistics_names.", 251016_ik)
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
! 3. Enables flexible positioning of mean values anywhere AFTER the variance/covariance
! 4. Supports exactly 2 possible mean names (e.g., "ux-avg", "ux-mean")
! 5. Provides consistent error handling with custom messages and error codes
!
! Note: Cases with more than 2 possible mean names (like vort-var with 4 options:
! vort-avg, vor-avg, vort-mean, vor-mean) still use manual search loops since
! this helper is designed for the common case of 2 name variants.
!-----------------------------------------------------------------------------
subroutine find_single_mean_index(current_idx, mean_name1, mean_name2, mean_idx, error_message, error_code)
    implicit none
    
    integer(kind=ik), intent(in) :: current_idx
    character(len=*), intent(in) :: mean_name1, mean_name2  ! Two possible names (e.g., "ux-avg", "ux-mean")
    integer(kind=ik), intent(out) :: mean_idx
    character(len=*), intent(in) :: error_message
    integer(kind=ik), intent(in) :: error_code
    
    integer(kind=ik) :: j
    
    mean_idx = -1  ! Initialize to not found
    
    ! Search for either mean name after current index
    do j = current_idx + 1, params_acm%N_time_statistics
        if (trim(params_acm%time_statistics_names(j)) == trim(mean_name1) .or. &
            trim(params_acm%time_statistics_names(j)) == trim(mean_name2)) then
            mean_idx = j
            exit
        end if
    end do
    
    ! If mean is not found, abort with error
    if (mean_idx == -1) then
        call abort(error_code, error_message)
    end if
    
end subroutine find_single_mean_index