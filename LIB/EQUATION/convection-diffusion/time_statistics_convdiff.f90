!-----------------------------------------------------------------------------
! main level wrapper to compute local statistics over time, like average or minmax of specific quantities
! these are then also adapted over time
!-----------------------------------------------------------------------------
subroutine TIME_STATISTICS_convdiff( time, dt, time_start, u, g, x0, dx, work, mask )
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

    integer(kind=ik) :: iB, ix, iy, iz, Bs(3), i_ts, N_offset
    real(kind=rk) :: time_diff

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    N_offset = params_convdiff%N_scalars
    time_diff = time-time_start

    do i_ts = 1, params_convdiff%N_time_statistics
        ! at the beginning of the simulation, initialize the array to zero
        if ( abs(time - time_start) < 1.0e-12_rk ) then
            u(:,:,:,N_offset + i_ts) = 0.0_rk
            cycle
        end if

        select case (trim(params_convdiff%time_statistics_names(i_ts)))
        case ("phi-int")
            ! compute the integral of phi over time
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* u(:,:,:,1)
        case ("phi-avg", "phi-mean")
            ! compute the average of phi over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* u(:,:,:,1)
        case ("phi-var")
            ! compute the variance of phi over time using Welford's method
            ! This needs the computation of the average in the variable afterwards to work
            if (i_ts == params_convdiff%N_time_statistics .or. (trim(params_convdiff%time_statistics_names(i_ts+1)) /= "phi-avg" .and. trim(params_convdiff%time_statistics_names(i_ts+1)) /= "phi-mean")) then
                call abort(2153001, "[TIME_STATISTICS_CONVDIFF]: You need to compute the variance together with the mean value. Insert 'phi-avg' right after 'phi-var' in time_statistics_names.")
            end if
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff) * &
                 (u(:,:,:,1) - u(:,:,:,N_offset + i_ts + 1)) * &
                 (u(:,:,:,1) - ((time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts + 1) + dt/(time_diff)*u(:,:,:,1)))
        case ("phi-minmax")
            ! compute the minmax of phi over time
            do iz = merge(1,g+1,params_convdiff%dim==2), merge(1,Bs(3)+g,params_convdiff%dim==2)
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        if (abs(u(ix,iy,iz,1)) > u(ix,iy,iz,N_offset + i_ts)) u(ix,iy,iz,N_offset + i_ts) = abs(u(ix,iy,iz,1))
                    end do
                end do
            end do

        case ("phi2-int")
            ! compute the integral of phi² over time
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* u(:,:,:,1)**2
        case ("phi2-avg", "phi2-mean")
            ! compute the average of phi² over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* u(:,:,:,1)**2
        
        ! Phi derivative cases
        case ("phidx-int")
            ! compute the integral of ∂phi/∂x over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("phidx-avg", "phidx-mean")
            ! compute the average of ∂phi/∂x over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        case ("phidy-int")
            ! compute the integral of ∂phi/∂y over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("phidy-avg", "phidy-mean")
            ! compute the average of ∂phi/∂y over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        case ("phidz-int")
            if (params_convdiff%dim == 2) call abort(250916, "[TIME_STATISTICS_CONVDIFF]: phidz not available in 2D simulations.")
            ! compute the integral of ∂phi/∂z over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("phidz-avg", "phidz-mean")
            if (params_convdiff%dim == 2) call abort(250916, "[TIME_STATISTICS_CONVDIFF]: phidz not available in 2D simulations.")
            ! compute the average of ∂phi/∂z over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        
        ! Squared phi derivative cases
        case ("phidx2-int")
            ! compute the integral of (∂phi/∂x)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("phidx2-avg", "phidx2-mean")
            ! compute the average of (∂phi/∂x)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        case ("phidy2-int")
            ! compute the integral of (∂phi/∂y)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("phidy2-avg", "phidy2-mean")
            ! compute the average of (∂phi/∂y)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        case ("phidz2-int")
            if (params_convdiff%dim == 2) call abort(250916, "[TIME_STATISTICS_CONVDIFF]: phidz2 not available in 2D simulations.")
            ! compute the integral of (∂phi/∂z)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = u(:,:,:,N_offset + i_ts) + dt* work(:,:,:,1)
        case ("phidz2-avg", "phidz2-mean")
            if (params_convdiff%dim == 2) call abort(250916, "[TIME_STATISTICS_CONVDIFF]: phidz2 not available in 2D simulations.")
            ! compute the average of (∂phi/∂z)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)

        case default
            call abort(2153000, "[TIME_STATISTICS_CONVDIFF]: time_statistics_name "// &
                 trim(params_convdiff%time_statistics_names(i_ts))//" not recognized.")
        end select
    end do

end subroutine TIME_STATISTICS_convdiff