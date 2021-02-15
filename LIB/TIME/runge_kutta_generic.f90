subroutine RungeKuttaGeneric(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, &
    hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist)
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

    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: j, k, hvy_id, z1, z2, g, Neqn
    real(kind=rk) :: t
    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable, save  :: rk_coeffs(:,:)

    Neqn = params%n_eqn
    Bs   = params%Bs
    g    = params%n_ghosts

    if (params%dim==2) then
        z1 = 1
        z2 = 1
    else
        z1 = g+1
        z2 = Bs(3)+g
    endif

    if (.not.allocated(rk_coeffs)) allocate(rk_coeffs(size(params%butcher_tableau,1),size(params%butcher_tableau,2)) )
    dt = 9.0e9_rk
    ! set rk_coeffs
    rk_coeffs = params%butcher_tableau

    ! synchronize ghost nodes
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

    ! calculate time step
    call calculate_time_step(params, time, iteration, hvy_block, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), lgt_block, &
        lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), dt)

    ! first stage, call to RHS. note the resulting RHS is stored in hvy_work(), first
    ! slot after the copy of the state vector (hence 2)
    call RHS_wrapper(time + dt*rk_coeffs(1,1), params, hvy_block, hvy_work(:,:,:,:,:,2), &
    hvy_mask, hvy_tmp, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor )

    ! save data at time t to heavy work array
    ! copy state vector content to work array. NOTE: 09/04/2018: moved this after RHS_wrapper
    ! since we can allow the RHS wrapper to modify the state vector (eg for mean flow fixing)
    ! if the copy part is above, the changes in state vector are ignored
    do k = 1, hvy_n(tree_ID_flow)
        hvy_id = hvy_active(k, tree_ID_flow)
        ! first slot in hvy_work is previous time step
        hvy_work( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id, 1 ) = hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, z1:z2, :, hvy_id )
    end do


    ! compute k_1, k_2, .... (coefficients for final stage)
    do j = 2, size(rk_coeffs, 1) - 1
        ! prepare input for the RK substep
        call set_RK_input(dt, params, rk_coeffs(j,:), j, hvy_block, hvy_work, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))

        ! synchronize ghost nodes for new input
        ! further ghost nodes synchronization, fixed grid
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

        ! note substeps are at different times, use temporary time "t"
        t = time + dt*rk_coeffs(j,1)

        call RHS_wrapper(t, params, hvy_block, hvy_work(:,:,:,:,:,j+1), hvy_mask, hvy_tmp, lgt_block, &
        lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor)
    end do

    ! final stage
    call final_stage_RK(params, dt, hvy_work, hvy_block, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), rk_coeffs)

end subroutine RungeKuttaGeneric
