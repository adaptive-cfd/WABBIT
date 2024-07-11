! Time stepper to advance coupled fluid/rigid solid system with the same order and
! method for both subproblems
subroutine RungeKuttaGeneric_FSI(time, dt, iteration, params, hvy_block, hvy_work, hvy_mask, hvy_tmp, tree_ID)
    ! it is unfortunate but this routine breaks the encapsulation concept, as it requires to
    ! have access to the "INSECT" structure for time stepping. hence, but ACM and INSECT modules
    ! have to be loaded here.
    use module_acm
    use module_insects
    implicit none

    real(kind=rk), intent(inout)        :: time, dt
    integer(kind=ik), intent(in)        :: iteration
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)
    !> hvy_tmp are qty that depend on the grid and not explicitly on time.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: j, k, hvy_id, z1, z2, g, l, grhs
    real(kind=rk) :: t
    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable, save  :: rk_coeffs(:,:)

    Bs   = params%Bs
    g    = params%g
    grhs = params%g_rhs

    if (params%dim==2) then
        z1 = 1
        z2 = 1
    else
        z1 = g+1
        z2 = Bs(3)+g
    endif

    if (.not. params%physics_type=="ACM-new") then
        call abort(1202211,"this special FSI time stepper works only with the ACM")
    endif

    if (.not. allocated(Insect%rhs)) then
        allocate(Insect%rhs(1:20, 1:size(hvy_work,6)))
    endif

    if (.not.allocated(rk_coeffs)) allocate(rk_coeffs(size(params%butcher_tableau,1),size(params%butcher_tableau,2)) )
    dt = 9.0e9_rk
    ! set rk_coeffs
    rk_coeffs = params%butcher_tableau

    ! synchronize ghost nodes
    call sync_ghosts_RHS_tree( params, hvy_block, tree_ID, g_minus=grhs, g_plus=grhs)

    ! calculate time step
    call calculate_time_step(params, time, iteration, hvy_block, dt, tree_ID)

    ! first stage, call to RHS. note the resulting RHS is stored in hvy_work(), first
    ! slot after the copy of the state vector (hence 2)
    call RHS_wrapper(time + dt*rk_coeffs(1,1), params, hvy_block, hvy_work(:,:,:,:,:,2), hvy_mask, hvy_tmp, tree_ID)

    ! the rhs wrapper has computed params_acm%force_insect_g and moment_insect_g
    call rigid_solid_rhs(time + dt*rk_coeffs(1,1), iteration, Insect%STATE, Insect%rhs(:,2), &
    params_acm%force_insect_g, params_acm%moment_insect_g, Insect)

    ! save data at time t to heavy work array
    ! copy state vector content to work array. NOTE: 09/04/2018: moved this after RHS_wrapper
    ! since we can allow the RHS wrapper to modify the state vector (eg for mean flow fixing)
    ! if the copy part is above, the changes in state vector are ignored
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        ! first slot in hvy_work is previous time step (time level at start of time step)
        hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id, 1 ) = &
        hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id )
    end do
    Insect%rhs(:,1) = Insect%STATE(:)



    ! compute k_1, k_2, .... (coefficients for final stage)
    do j = 2, size(rk_coeffs, 1) - 1
        !-----------------------------------------------------------------------
        ! prepare input for the RK substep
        !-----------------------------------------------------------------------
        ! gives back the input for the RHS (from which in the final stage the next
        ! time step is computed).\n
        !
        ! k_j = RHS(t+dt*c_j,  datafield(t) + dt*sum(a_jl*k_l))
        ! (e.g. k3 = RHS(t+dt*c_3, data_field(t) + dt*(a31*k1+a32*k2)) ) \n

        ! first: k_j = RHS(data_field(t) + ...
        ! loop over all active heavy data blocks
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            ! first slot in hvy_work is previous time step
            hvy_block(g+1:Bs(1)+g,g+1:Bs(2)+g,z1:z2,:,hvy_id) = &
            hvy_work(g+1:Bs(1)+g, g+1:Bs(2)+g,z1:z2,:,hvy_id,1)
        end do
        Insect%STATE(:) = Insect%rhs(:,1)

        do l = 2, j
            ! check if coefficient is zero - if so, avoid loop over all components and active blocks
            if (abs(rk_coeffs(j,l)) < 1.0e-8_rk) then
                cycle
            end if

            ! loop over all active heavy data blocks
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                ! new input for computation of k-coefficients
                ! k_j = RHS((t+dt*c_j, data_field(t) + sum(a_jl*k_l))
                hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id) = &
                hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id) &
                + dt * rk_coeffs(j,l) * hvy_work(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id, l)
            end do

            Insect%STATE(:) = Insect%STATE(:) + dt * rk_coeffs(j,l) * Insect%rhs(:,l)
        end do

        !-----------------------------------------------------------------------
        ! RK substep
        !-----------------------------------------------------------------------
        ! synchronize ghost nodes for new input
        call sync_ghosts_RHS_tree( params, hvy_block, tree_ID, g_minus=grhs, g_plus=grhs)

        ! note substeps are at different times, use temporary time "t"
        t = time + dt*rk_coeffs(j,1)

        ! computes mask from Insect%STATE, computes aerodyn forces using hvy_block, both hvy_block and Insect%STATE are updated above
        call RHS_wrapper(t, params, hvy_block, hvy_work(:,:,:,:,:,j+1), hvy_mask, hvy_tmp, tree_ID)

        ! the rhs wrapper has computed params_acm%force_insect_g and moment_insect_g
        call rigid_solid_rhs(t, iteration, Insect%STATE, Insect%rhs(:,j+1), &
        params_acm%force_insect_g, params_acm%moment_insect_g, Insect)
    end do



    ! final stage (actual final update of state vector)
    ! for the RK4 the final stage looks like this:
    ! data_field(t+dt) = data_field(t) + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! u_n = u_n +...
        hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id) = &
        hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id, 1)

        do j = 2, size(rk_coeffs, 2)
            ! check if coefficient is zero - if so, avoid loop over all components and active blocks
            if ( abs(rk_coeffs(size(rk_coeffs, 1),j)) < 1.0e-8_rk) then
                cycle
            endif

            ! ... dt*(b1*k1 + b2*k2+ ..)
            ! rk_coeffs(size(rk_coeffs,1)) , since we want to access last line,
            ! e.g. b1 = butcher(last line,2)
            hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id) = &
            hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id) + &
            dt*rk_coeffs(size(rk_coeffs,1),j) * hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id, j)
        end do
    end do

    Insect%STATE(:) = Insect%rhs(:,1)
    do j = 2, size(rk_coeffs, 2)
        if ( abs(rk_coeffs(size(rk_coeffs, 1),j)) < 1.0e-8_rk) then
            cycle
        endif
        Insect%STATE = Insect%STATE + dt*rk_coeffs(size(rk_coeffs,1),j)*Insect%rhs(:,j)
    enddo

    call append_t_file( 'insect_state_vector.t', (/time+dt, Insect%STATE, params_acm%force_insect_g/) )

    Insect%time = Insect%time+dt
end subroutine RungeKuttaGeneric_FSI
