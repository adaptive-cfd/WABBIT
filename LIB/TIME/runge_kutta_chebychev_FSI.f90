subroutine RungeKuttaChebychev_FSI(time, dt, iteration, params, hvy_block, hvy_work, &
    hvy_mask, hvy_tmp, tree_ID)
    ! it is unfortunate but this routine breaks the encapsulation concept, as it requires to
    ! have access to the "INSECT" structure for time stepping. hence, but ACM and INSECT modules
    ! have to be loaded here.
    use module_acm
    use module_insects
    implicit none

    real(kind=rk), intent(inout)        :: time, dt
    integer(kind=ik), intent(in)        :: iteration
    type (type_params), intent(in)      :: params                       !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)   !> heavy work data array - block data
    !> hvy_tmp are qty that depend on the grid and not explicitly on time.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    ! in fortran, we work with indices:
    integer                             :: y0=3, y1=4, y2=5, F1=6, tmp(1:3)
    integer, parameter                  :: y00=1, F0=2
    integer                             :: i, k, s, hvy_id
    real(kind=rk)                       :: tau
    logical, save                       :: setup_complete = .false.
    logical, save                       :: informed = .false.

    if (.not. setup_complete) then
        call setup_RKC_coefficients(params)
        setup_complete = .true.
    endif


    ! s is the number of stages
    s = params%s

    if (.not. informed) then
        if (params%rank==0) then
            write(*,'(80("-"))')
            write(*,*) "Runge-Kutta-Chebychev method (FSI version!)"
            write(*,'("Using s=",i2," stages")') s
            write(*,'(80("-"))')
        endif
        informed = .true.
    endif

    !---------------------------------------------------------------------------
    ! check if we have enough RHS slots available
    if (size(hvy_work,6) < 6 ) then
        call abort(06091916, "not enough registers for runge-kutta-chebychev method")
    endif

    if (.not. params%physics_type=="ACM-new") then
        call abort(1202211,"this special FSI time stepper works only with the ACM")
    endif

    if (.not. allocated(Insect%rhs)) then
        allocate(Insect%rhs(1:20, 1:size(hvy_work,6)))
    endif

    if (s<4) then
        call abort(1715929,"runge-kutta-chebychev: s cannot be less than 4")
    endif

    ! synchronize ghost nodes
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID), skipDiagonalNeighbors=.true. )

    ! calculate time step
    call calculate_time_step(params, time, iteration, hvy_block, dt, tree_ID)

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! Y0 = u (in matlab: y00 = u;)
        hvy_work(:,:,:,:,hvy_id, y00 ) = hvy_block(:,:,:,:,hvy_id)
        ! we need two copies (one is an iteration variable, the other (above) is kept constant)
        hvy_work(:,:,:,:,hvy_id, y0  ) = hvy_block(:,:,:,:,hvy_id)
    enddo
    Insect%rhs(:,y00) = Insect%STATE
    Insect%rhs(:,y0 ) = Insect%STATE

    ! F0 (RHS at initial time, old time level)
    ! note: call sync_ghosts on input data before
    call RHS_wrapper( time, params, hvy_block, hvy_work(:,:,:,:,:,F0), hvy_mask, hvy_tmp, tree_ID )

    ! the rhs wrapper has computed params_acm%force_insect_g and moment_insect_g
    call rigid_solid_rhs( time, iteration, Insect%STATE, Insect%rhs(:,F0), &
    params_acm%force_insect_g, params_acm%moment_insect_g, Insect )

    ! euler step
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        ! y1 = y0 + mu_tilde(1) * dt * F0;
        hvy_work(:,:,:,:,hvy_id, y1 ) = hvy_work(:,:,:,:,hvy_id, y0) &
        + mu_tilde(s,1) * dt * hvy_work(:,:,:,:,hvy_id, F0)
    enddo
    Insect%rhs(:,y1) = Insect%rhs(:,y0) + mu_tilde(s,1) * dt * Insect%rhs(:,F0)

    ! runge-kutta-chebychev stages
    do i = 2, s
        ! for explicitly time-dependend RHS [tau = time + c(i-1)*dt;]
        tau = time + c(s,i-1)*dt

        ! F1 = rhs(y1);
        ! note: call sync_ghosts on input data before
        call sync_ghosts( params, lgt_block, hvy_work(:,:,:,:,:,y1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID), skipDiagonalNeighbors=.true. )

        call RHS_wrapper( tau, params, hvy_work(:,:,:,:,:,y1), hvy_work(:,:,:,:,:,F1), hvy_mask, hvy_tmp, tree_ID )
        ! the rhs wrapper has computed params_acm%force_insect_g and moment_insect_g
        call rigid_solid_rhs(tau, iteration, Insect%rhs(:,y1), Insect%rhs(:,F1), &
        params_acm%force_insect_g, params_acm%moment_insect_g, Insect)

        ! main formula
        ! y2 = (1-mu(i)-nu(i)) * y00 + mu(i) * y1 + nu(i) * y0 + mu_tilde(i)*dt*F1 + gamma_tilde(i)*dt*F0;
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k,tree_ID)

            hvy_work(:,:,:,:,hvy_id, y2 ) = (1.0_rk-mu(s,i)-nu(s,i))*hvy_work(:,:,:,:,hvy_id, y00) &
            + mu(s,i) * hvy_work(:,:,:,:,hvy_id, y1 ) &
            + nu(s,i) * hvy_work(:,:,:,:,hvy_id, y0 ) &
            + mu_tilde(s,i) * dt * hvy_work(:,:,:,:,hvy_id, F1 ) &
            + gamma_tilde(s,i) * dt * hvy_work(:,:,:,:,hvy_id, F0 )
        enddo

        Insect%rhs(:,y2) = (1.0_rk-mu(s,i)-nu(s,i))*Insect%rhs(:,y00) &
        + mu(s,i) * Insect%rhs(:,y1) &
        + nu(s,i) * Insect%rhs(:,y0) &
        + mu_tilde(s,i) * dt * Insect%rhs(:,F1) &
        + gamma_tilde(s,i) * dt * Insect%rhs(:,F0)

        ! iteration. note in last stage, do not iterate (we return y0 otherwise)
        if ( i < s) then
            tmp = (/y0, y1, y2/)
            y0 = tmp(2)
            y1 = tmp(3)
            y2 = tmp(1)
        endif
    end do

    ! return result
    !!     u = y2;
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)

        hvy_block(:,:,:,:,hvy_id ) = hvy_work(:,:,:,:,hvy_id, y2 )
    enddo

    Insect%STATE = Insect%rhs(:,y2)

    call append_t_file( 'insect_state_vector.t', (/time+dt, Insect%STATE, params_acm%force_insect_g/) )

    Insect%time = Insect%time+dt

end subroutine RungeKuttaChebychev_FSI
