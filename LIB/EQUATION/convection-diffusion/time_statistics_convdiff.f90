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

    integer(kind=ik) :: iB, ix, iy, iz, Bs(3), i_ts, N_offset, i_scalar, mean_idx1
    real(kind=rk) :: time_diff
    character(len=cshort) :: name_phi, stat_name

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

        ! store the trimmed statistic name for easier comparison
        stat_name = trim(params_convdiff%time_statistics_names(i_ts))

        ! we need to know which scalar we are working on, it is in the form phiX-var or phidx-avg etc.
        ! so we read the number after "phi" and before "-"
        i_scalar = 0
        if (stat_name(1:3) == "phi") then
            ! check if the 4th character is actually a digit
            if (.not. (stat_name(4:4) >= '0' .and. stat_name(4:4) <= '9')) then
                call abort(250931, "[TIME_STATISTICS_CONVDIFF]: time_statistics_name '"//trim(stat_name)//"' has invalid format. Expected 'phiX-[statistic]' where X is a digit (1-9), but found '"//stat_name(4:4)//"' at position 4.")
            end if
            
            read( stat_name(4:4), * ) i_scalar
            if (i_scalar < 1 .or. i_scalar > params_convdiff%N_scalars) call abort(250929, "[TIME_STATISTICS_CONVDIFF]: In time_statistics_names, the scalar index after 'phi' is out of range.")
        end if
        write( name_phi, '(A,I0)' ) "phi", i_scalar

        if (stat_name == trim(name_phi) // "-avg") then
            ! compute the average of phi over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* u(:,:,:,1)
        elseif (stat_name == trim(name_phi) // "-var") then
            ! compute the variance of phi over time using Welford's method
            ! This needs the computation of the average somewhere after this field
            call find_single_mean_index(i_ts, trim(name_phi)//"-avg", mean_idx1, &
                "[TIME_STATISTICS_CONVDIFF]: You need to compute the variance together with the mean value. Insert '"//trim(name_phi)//"-avg' somewhere after '"//trim(name_phi)//"-var' in time_statistics_names.", 251023_ik)
            ! var_new = (time_diff-dt)/time_diff*var_old + dt/time_diff*(x - mean_old)*(x - mean_new)
            ! mean_new = (time_diff-dt)/time_diff*mean_old + dt/time_diff*x
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff) * &
                 (u(:,:,:,1) - u(:,:,:,N_offset + mean_idx1)) * &
                 (u(:,:,:,1) - (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + mean_idx1) - dt/(time_diff)*u(:,:,:,1))
            elseif (stat_name == trim(name_phi) // "-avg") then
                ! compute the average of phi over time (required for variance)
                u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* u(:,:,:,1)
        elseif (stat_name == trim(name_phi) // "-minmax") then
            ! compute the minmax of phi over time
            do iz = merge(1,g+1,params_convdiff%dim==2), merge(1,Bs(3)+g,params_convdiff%dim==2)
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        if (abs(u(ix,iy,iz,1)) > u(ix,iy,iz,N_offset + i_ts)) u(ix,iy,iz,N_offset + i_ts) = abs(u(ix,iy,iz,1))
                    end do
                end do
            end do

        elseif (stat_name == trim(name_phi) // "-2-avg") then
            ! compute the average of phi² over time
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* u(:,:,:,1)**2
        
        ! Phi derivative cases
        elseif (stat_name == trim(name_phi) // "dx-avg") then
            ! compute the average of ∂phi/∂x over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        elseif (stat_name == trim(name_phi) // "dy-avg") then
            ! compute the average of ∂phi/∂y over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        elseif (stat_name == trim(name_phi) // "dz-avg") then
            if (params_convdiff%dim == 2) call abort(250916, "[TIME_STATISTICS_CONVDIFF]: phidz not available in 2D simulations.")
            ! compute the average of ∂phi/∂z over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_convdiff%discretization, work(:,:,:,1))
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        
        ! Squared phi derivative cases
        elseif (stat_name == trim(name_phi) // "dx2-avg") then
            ! compute the average of (∂phi/∂x)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 1, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        elseif (stat_name == trim(name_phi) // "dy2-avg") then
            ! compute the average of (∂phi/∂y)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 2, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        elseif (stat_name == trim(name_phi) // "dz2-avg") then
            if (params_convdiff%dim == 2) call abort(250916, "[TIME_STATISTICS_CONVDIFF]: phidz2 not available in 2D simulations.")
            ! compute the average of (∂phi/∂z)² over time
            call compute_derivative(u(:,:,:,1), dx, Bs, g, 3, 1, params_convdiff%discretization, work(:,:,:,1))
            work(:,:,:,1) = work(:,:,:,1)**2
            u(:,:,:,N_offset + i_ts) = (time_diff-dt)/(time_diff)*u(:,:,:,N_offset + i_ts) + dt/(time_diff)* work(:,:,:,1)
        
        ! ToDo: Implement cases for velocity

        else
            call abort(2153000, "[TIME_STATISTICS_CONVDIFF]: time_statistics_name "// &
                 stat_name//" not recognized.")
        end if
    end do

end subroutine TIME_STATISTICS_convdiff


!----------------------------------------------------------------------------- 
! Helper subroutine to find a single mean value index _after a given position_
! This is adapted from the ACM version for flexible mean search in variance/covariance cases
!----------------------------------------------------------------------------- 
subroutine find_single_mean_index(current_idx, mean_name, mean_idx, error_message, error_code)
    implicit none
    integer(kind=ik), intent(in) :: current_idx
    character(len=*), intent(in) :: mean_name
    integer(kind=ik), intent(out) :: mean_idx
    character(len=*), intent(in) :: error_message
    integer(kind=ik), intent(in) :: error_code
    integer(kind=ik) :: j
    mean_idx = -1
    do j = current_idx + 1, params_convdiff%N_time_statistics
        if (trim(params_convdiff%time_statistics_names(j)) == trim(mean_name)) then
            mean_idx = j
            exit
        end if
    end do
    if (mean_idx == -1) then
        call abort(error_code, error_message)
    end if
end subroutine find_single_mean_index