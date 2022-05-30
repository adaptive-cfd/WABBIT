subroutine RungeKuttaGeneric(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, &
    hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist, tree_id)
    implicit none

    real(kind=rk), intent(inout)        :: time, dt
    integer(kind=ik), intent(in)        :: iteration
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)
    !> hvy_tmp are qty that depend on the grid and not explicitly on time.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:,:)
    !> is only needed for testing purposes
    integer(kind=ik), intent(in), optional :: tree_id


    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: j, k, hvy_id, g, Neqn, l, tree_id_evolv !z1, z2
    real(kind=rk) :: t
    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable, save  :: rk_coeffs(:,:)

    Neqn = params%n_eqn
    Bs   = params%Bs
    g    = params%n_ghosts

    ! if (params%dim==2) then
    !     z1 = 1
    !     z2 = 1
    ! else
    !     z1 = g+1
    !     z2 = Bs(3)+g
    ! endif


    if (present(tree_id)) then
        tree_id_evolv = tree_id
    else
        tree_id_evolv = tree_ID_flow
    end if

    if (.not.allocated(rk_coeffs)) allocate(rk_coeffs(size(params%butcher_tableau,1),size(params%butcher_tableau,2)) )
    dt = 9.0e9_rk
    ! set rk_coeffs
    rk_coeffs = params%butcher_tableau

    ! synchronize ghost nodes
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_evolv), hvy_n(tree_id_evolv) )

    ! calculate time step
    call calculate_time_step(params, time, iteration, hvy_block, hvy_active(:,tree_id_evolv), hvy_n(tree_id_evolv), lgt_block, &
        lgt_active(:,tree_id_evolv), lgt_n(tree_id_evolv), dt)
    ! first stage, call to RHS. note the resulting RHS is stored in hvy_work(), first
    ! slot after the copy of the state vector (hence 2)
    call RHS_wrapper(time + dt*rk_coeffs(1,1), params, hvy_block, hvy_work(:,:,:,:,:,2), &
    hvy_mask, hvy_tmp, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor, tree_id_evolv )

    ! save data at time t to heavy work array
    ! copy state vector content to work array. NOTE: 09/04/2018: moved this after RHS_wrapper
    ! since we can allow the RHS wrapper to modify the state vector (eg for mean flow fixing)
    ! if the copy part is above, the changes in state vector are ignored
    do k = 1, hvy_n(tree_id_evolv)
        hvy_id = hvy_active(k, tree_id_evolv)
        ! first slot in hvy_work is previous time step (time level at start of time step)
        ! hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id, 1 ) = hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id )
        hvy_work( :, :, :, :, hvy_id, 1 ) = hvy_block( :,:,:, :, hvy_id )
    end do


    ! compute k_1, k_2, .... (coefficients for final stage)
    do j = 2, size(rk_coeffs, 1) - 1
        ! prepare input for the RK substep
        ! gives back the input for the RHS (from which in the final stage the next
        ! time step is computed).\n
        !
        ! k_j = RHS(t+dt*c_j,  datafield(t) + dt*sum(a_jl*k_l))
        ! (e.g. k3 = RHS(t+dt*c_3, data_field(t) + dt*(a31*k1+a32*k2)) ) \n

        ! first: k_j = RHS(data_field(t) + ...
        ! loop over all active heavy data blocks
        do k = 1, hvy_n(tree_id_evolv)
            ! first slot in hvy_work is previous time step
            ! hvy_block(g+1:Bs(1)+g,g+1:Bs(2)+g,z1:z2,:,hvy_active(k,tree_id_evolv)) = &
            ! hvy_work(g+1:Bs(1)+g, g+1:Bs(2)+g,z1:z2,:,hvy_active(k,tree_id_evolv),1)

            hvy_block(:,:,:,:,hvy_active(k,tree_id_evolv)) = &
            hvy_work(:,:,:,:,hvy_active(k,tree_id_evolv),1)
        end do

        do l = 2, j
            ! check if coefficient is zero - if so, avoid loop over all data fields and active blocks
            if (abs(rk_coeffs(j,l)) < 1.0e-8_rk) then
                cycle
            end if

            ! loop over all active heavy data blocks
            do k = 1, hvy_n(tree_id_evolv)
                ! new input for computation of k-coefficients
                ! k_j = RHS((t+dt*c_j, data_field(t) + sum(a_jl*k_l))
                hvy_block(:, :, :, :, hvy_active(k,tree_id_evolv)) = &
                hvy_block(:, :, :, :, hvy_active(k,tree_id_evolv)) &
                + dt * rk_coeffs(j,l) * hvy_work(:,:,:,:,hvy_active(k, tree_id_evolv), l)
                ! hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_active(k,tree_id_evolv)) = &
                ! hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_active(k,tree_id_evolv)) &
                ! + dt * rk_coeffs(j,l) * hvy_work(g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_active(k, tree_id_evolv), l)
            end do
        end do

        ! synchronize ghost nodes for new input
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_evolv), hvy_n(tree_id_evolv) )

        ! note substeps are at different times, use temporary time "t"
        t = time + dt*rk_coeffs(j,1)

        call RHS_wrapper(t, params, hvy_block, hvy_work(:,:,:,:,:,j+1), hvy_mask, hvy_tmp, lgt_block, &
        lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor, tree_id_evolv )
    end do



    ! final stage (actual final update of state vector)
    ! for the RK4 the final stage looks like this:
    ! data_field(t+dt) = data_field(t) + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
    do k = 1, hvy_n(tree_id_evolv)
        ! u_n = u_n +...
        hvy_block( :, :, :, 1:Neqn, hvy_active(k,tree_id_evolv)) = &
        hvy_work( :, :, :, 1:Neqn, hvy_active(k,tree_id_evolv), 1)
        ! hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn, hvy_active(k,tree_id_evolv)) = &
        ! hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn,hvy_active(k,tree_id_evolv), 1)

        do j = 2, size(rk_coeffs, 2)
            ! check if coefficient is zero - if so, avoid loop over all data fields and active blocks
            if ( abs(rk_coeffs(size(rk_coeffs, 1),j)) < 1.0e-8_rk) then
                cycle
            endif

            ! ... dt*(b1*k1 + b2*k2+ ..)
            ! rk_coeffs(size(rk_coeffs,1)) , since we want to access last line,
            ! e.g. b1 = butcher(last line,2)
            hvy_block( :,:,:, 1:Neqn, hvy_active(k,tree_id_evolv)) = hvy_block( :,:,:, 1:Neqn, hvy_active(k,tree_id_evolv)) &
                   + dt*rk_coeffs(size(rk_coeffs,1),j) * hvy_work( :,:,:, 1:Neqn, hvy_active(k,tree_id_evolv), j)
            ! hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn, hvy_active(k,tree_id_evolv)) = hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn, hvy_active(k,tree_id_evolv)) &
            !        + dt*rk_coeffs(size(rk_coeffs,1),j) * hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, 1:Neqn, hvy_active(k,tree_id_evolv), j)
        end do
    end do

end subroutine RungeKuttaGeneric
