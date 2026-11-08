! Time stepper to advance coupled fluid/rigid solid system with the same order and
! method for both subproblems
subroutine RungeKuttaGeneric_FSI(time, dt, iteration, params, hvy_block, hvy_work, hvy_mask, hvy_tmp, tree_ID)
    ! it is unfortunate but this routine breaks the encapsulation concept, as it requires to
    ! have access to the "INSECT" structure for time stepping. hence, but ACM and INSECT modules
    ! have to be loaded here.
    use module_acm
    use module_insects
    use module_solid_model

    implicit none

    real(kind=rk), intent(inout)        :: time, dt
    integer(kind=ik), intent(inout)     :: iteration
    !> user defined parameter structure
    type (type_params), intent(inout)   :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)
    !> hvy_tmp are qty that depend on the grid and not explicitly on time.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: j, k, hvy_id, z1, z2, g, l, grhs, Neqn_RHS, i_insect
    real(kind=rk) :: t, t_call, t_stage
    real(kind=rk), allocatable, save    :: state_vector_write(:)
    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable, save  :: rk_coeffs(:,:)

    Bs   = params%Bs
    g    = params%g
    grhs = params%g_rhs
    Neqn_RHS = params%n_eqn_rhs

#ifdef DEV
    if (grhs>g) then
        write(*,*) "grhs=", grhs, "g=", g
        call abort(180225,"ERROR: grhs>g")
    endif
#endif

    if (params%dim==2) then
        z1 = 1
        z2 = 1
    else
        z1 = g+1
        z2 = Bs(3)+g
    endif

    
    if (.not. (params%physics_type=="ACM-new" .or. params%physics_type=="NSPP")) then
        call abort(1202211,"this special FSI time stepper works only with ACM or NSPP for now")
    endif
    if (.not. allocated(state_vector_write)) allocate(state_vector_write(1:26 * n_insects + 1))


    if (.not.(params_acm%use_solid_model .or. params_acm%use_free_flight_solver)) then
        call abort(2611031, "Using the FSI time stepper without params_acm%use_solid_model .or. params_acm%use_free_flight_solver makes no sense.")
    endif

    ! Rigid body dynamics (free-flight) solver. Solves the EOM (Newtons law) together with the fluid.
    ! For both fluid and solid, the same time stepper is used.
    do i_insect = 1,n_insects
        if (.not. allocated(Insects(i_insect)%rhs)) then
            allocate(Insects(i_insect)%rhs(1:20, 1:size(hvy_work,6)))
        endif
    enddo

    if (.not.allocated(rk_coeffs)) allocate(rk_coeffs(size(params%butcher_tableau,1),size(params%butcher_tableau,2)) )
    dt = 9.0e9_rk
    ! set rk_coeffs
    rk_coeffs = params%butcher_tableau

    ! synchronize ghost nodes
    t_call = MPI_wtime()
    call sync_ghosts_RHS_tree( params, hvy_block(:,:,:,1:Neqn_RHS,:), tree_ID, g_minus=grhs, g_plus=grhs)
    call toc( "timestep (sync ghosts)", 20, MPI_wtime()-t_call)


    if (params_acm%use_solid_model) then
        ! --- FSI WITH BEAM-EQN ----
        ! NOTE: the FSI solver uses two different time steppers for fluid and solid.
        !
        ! Pressure interpolation at the old time level (in fact, only used in the first CN2 step).
        ! Note how grhs > N_support must be satisfied if the delta interpolation is used.
        if (grhs < 3) call abort(26092609, "Change grhs to 3 for the moment ")

        ! update top/bottom surfaces - these are the points for the interpolation (xq_solid)
        call plate_geometry( time )

        ! interpolate pressure at the top- and bottom surface (the Beam eqn will only see their difference)
        ! Note: interpolation is a grid level routine and thus needs to be done here, on the WABBIT level
        call interpolatePointCloud_tree(params, hvy_block(:,:,:,4:4,:), 1, xq_solid, pq_solid_old, interpolationMethod, .false. )
    endif

    ! calculate time step
    call calculate_time_step(params, time, iteration, hvy_block, dt, tree_ID)

    !------------------------------------------------------------------------------------------------------------------
    ! save data at time t to heavy work array
    !------------------------------------------------------------------------------------------------------------------
    ! copy state vector content to work array. NOTE: 09/04/2018: moved this after RHS_wrapper
    ! since we can allow the RHS wrapper to modify the state vector (eg for mean flow fixing)
    ! if the copy part is above, the changes in state vector are ignored
    ! NOTE: 02/02/2026 by JB: I moved it back before, because we would like to re-use the space in hvy_block for projection methods
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        ! first slot in hvy_work is previous time step (time level at start of time step)
        hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id, 1 ) = &
            hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id )
    end do

    ! Rigid body dynamics (free-flight) solver. Solves the EOM (Newtons law) together with the fluid.
    ! For both fluid and solid, the same time stepper is used.
    if (params_acm%use_free_flight_solver) then
        do i_insect = 1, n_insects
            Insects(i_insect)%rhs(:,1) = Insects(i_insect)%STATE
        enddo
    endif

    !------------------------------------------------------------------------------------------------------------------
    ! first RK stage
    !------------------------------------------------------------------------------------------------------------------
    ! first stage, call to RHS. note the resulting RHS is stored in hvy_work(), first
    ! slot after the copy of the state vector (hence 2)
    t_call = MPI_wtime()
    call RHS_wrapper(time + dt*rk_coeffs(1,1), params, hvy_block, hvy_work(:,:,:,:,:,2), hvy_mask, hvy_tmp, tree_ID)
    call toc( "timestep (RHS wrapper)", 21, MPI_wtime()-t_call)

    ! Rigid body dynamics (free-flight) solver. Solves the EOM (Newtons law) together with the fluid.
    ! For both fluid and solid, the same time stepper is used.
    if (params_acm%use_free_flight_solver) then
        t_call = MPI_wtime()
        do i_insect = 1, n_insects
            ! the rhs wrapper has computed insects(i_insect)%force_g and moment_g
            call rigid_solid_rhs( time + dt*rk_coeffs(1,1), iteration, Insects(i_insect)%STATE, Insects(i_insect)%rhs(:,2), &
                insects(i_insect)%force_g, insects(i_insect)%moment_g, Insects(i_insect) )
        enddo
        call toc( "timestep (RHS rigid solid)", 24, MPI_wtime()-t_call)
    endif

    !------------------------------------------------------------------------------------------------------------------
    ! other RK stages
    !------------------------------------------------------------------------------------------------------------------
    ! compute k_1, k_2, .... (coefficients for final stage)
    do j = 2, size(rk_coeffs, 1) - 1
        t_stage = MPI_wtime()
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
            hvy_block(g+1:Bs(1)+g,g+1:Bs(2)+g,z1:z2,1:Neqn_RHS,hvy_id) = &
                hvy_work(g+1:Bs(1)+g, g+1:Bs(2)+g,z1:z2,1:Neqn_RHS,hvy_id,1)
        end do

        ! Rigid body dynamics (free-flight) solver. Solves the EOM (Newtons law) together with the fluid.
        ! For both fluid and solid, the same time stepper is used.
        if (params_acm%use_free_flight_solver) then
            do i_insect = 1, n_insects
                Insects(i_insect)%STATE = Insects(i_insect)%rhs(:,1)
            enddo
        endif

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
                hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id) = &
                    hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id) &
                    + dt * rk_coeffs(j,l) * hvy_work(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id, l)
            end do

            ! Rigid body dynamics (free-flight) solver. Solves the EOM (Newtons law) together with the fluid.
            ! For both fluid and solid, the same time stepper is used.
            if (params_acm%use_free_flight_solver) then
                do i_insect = 1, n_insects
                    Insects(i_insect)%STATE(:) = Insects(i_insect)%STATE(:) + dt * rk_coeffs(j,l) * Insects(i_insect)%rhs(:,l)
                enddo
            endif
        end do

        !-----------------------------------------------------------------------
        ! RK substep
        !-----------------------------------------------------------------------
        ! synchronize ghost nodes for new input
        t_call = MPI_wtime()
        call sync_ghosts_RHS_tree( params, hvy_block(:,:,:,1:Neqn_RHS,:), tree_ID, g_minus=grhs, g_plus=grhs)
        call toc( "timestep (sync ghosts)", 20, MPI_wtime()-t_call)

        ! note substeps are at different times, use temporary time "t"
        t = time + dt*rk_coeffs(j,1)

        ! computes mask from Insect%STATE, computes aerodyn forces using hvy_block, both hvy_block and Insect%STATE are updated above
        t_call = MPI_wtime()
        call RHS_wrapper(t, params, hvy_block, hvy_work(:,:,:,:,:,j+1), hvy_mask, hvy_tmp, tree_ID)
        call toc( "timestep (RHS wrapper)", 21, MPI_wtime()-t_call)

        ! Rigid body dynamics (free-flight) solver. Solves the EOM (Newtons law) together with the fluid.
        ! For both fluid and solid, the same time stepper is used.
        if (params_acm%use_free_flight_solver) then
            t_call = MPI_wtime()
            do i_insect = 1, n_insects
                ! the rhs wrapper has computed insects(i_insect)%force_g and moment_g
                call rigid_solid_rhs(t, iteration, Insects(i_insect)%STATE, Insects(i_insect)%rhs(:,j+1), &
                    insects(i_insect)%force_g, insects(i_insect)%moment_g, Insects(i_insect))
            enddo
            call toc( "timestep (RHS rigid solid)", 22, MPI_wtime()-t_call)
        endif

        call toc( "timestep (RK stage)", 23, MPI_wtime()-t_stage)

    end do


    ! final stage (actual final update of state vector)
    ! for the RK4 the final stage looks like this:
    ! data_field(t+dt) = data_field(t) + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! u_n = u_n +...
        hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id) = &
            hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id, 1)

        do j = 2, size(rk_coeffs, 2)
            ! check if coefficient is zero - if so, avoid loop over all components and active blocks
            if ( abs(rk_coeffs(size(rk_coeffs, 1),j)) < 1.0e-8_rk) then
                cycle
            endif

            ! ... dt*(b1*k1 + b2*k2+ ..)
            ! rk_coeffs(size(rk_coeffs,1)) , since we want to access last line,
            ! e.g. b1 = butcher(last line,2)
            hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id) = &
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id) + &
                dt*rk_coeffs(size(rk_coeffs,1),j) * hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn_RHS, hvy_id, j)
        end do
    end do

    ! Rigid body dynamics (free-flight) solver. Solves the EOM (Newtons law) together with the fluid.
    ! For both fluid and solid, the same time stepper is used.
    if (params_acm%use_free_flight_solver) then
        do i_insect = 1, n_insects
            Insects(i_insect)%STATE(:) = Insects(i_insect)%rhs(:,1)
        enddo

        do j = 2, size(rk_coeffs, 2)
            if ( abs(rk_coeffs(size(rk_coeffs, 1),j)) < 1.0e-8_rk) then
                cycle
            endif
            do i_insect = 1, n_insects
                Insects(i_insect)%STATE = Insects(i_insect)%STATE + dt*rk_coeffs(size(rk_coeffs,1),j)*Insects(i_insect)%rhs(:,j)
            enddo
        enddo

        ! We write the insect state vector
        ! in order to avoid writing many files, everything is set into one, we only have to write the exact data
        state_vector_write = 0.0_rk
        state_vector_write(1) = time + dt
        do i_insect = 1, n_insects
            state_vector_write( (i_insect-1)*26 + 2 : i_insect*26 + 1) = (/ Insects(i_insect)%STATE, insects(i_insect)%force_g, insects(i_insect)%moment_g /)
        enddo
        call append_t_file( 'insect_state_vector.t', state_vector_write )

        do i_insect = 1, n_insects
            Insects(i_insect)%time = Insects(i_insect)%time + dt
        enddo
    endif

    
    if (params_acm%use_solid_model) then
        ! The fluid is now advanced in time: t -> t+dt
        ! During the fluid time stepping, the beam did not deform. However, the imposed motion at the 
        ! leading edge is correctly taken into account at different sub-time steps of the fluid time step.
        ! However, the deformation is constant over the fluid time step -- that is an approximation.
        
        
        ! update surfaces
        ! Note: we do not have the beam at the new time+dt yet - the call plate_geometry only updates the imposed motion
        ! at the leading edge. The deformation at time+dt is not yet computed - this is the major error contribution
        ! (approximation) of this time stepper. It's called semi-implicit staggered FSI method.
        call plate_geometry( time+dt )
        
        
        ! Interpolate forces using the advanced pressure field of the fluid (time+dt) at the positions
        ! returned by plate_geometry. There is a time lag here because we cannot easily solve the 
        ! nonlinear implicit system. (In other words, pressure_new is not actually the pressure_new, because 
        ! the deformation is at the old time level).
        ! Note: interpolation is a grid level routine and thus needs to be done here, on the WABBIT level
        call interpolatePointCloud_tree(params, hvy_block(:,:,:,4:4,:), 1, xq_solid, pq_solid_new, interpolationMethod, .true. )
        
        ! Make the time step of the solid, t->t+dt
        ! Input time is the old time level.
        call SolidSolverWrapper( time, dt )
    endif

end subroutine RungeKuttaGeneric_FSI
