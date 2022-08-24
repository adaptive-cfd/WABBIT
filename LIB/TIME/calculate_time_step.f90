!> \brief Get the time step according to the conditions of the physics module
subroutine calculate_time_step( params, time, iteration, hvy_block, dt, tree_ID )

    use module_physics_metamodule, only : GET_DT_BLOCK_meta

    !--------------------------------------------------------------
    implicit none
    type (type_params), intent(in):: params                    !< user defined parameter structure
    real(kind=rk), intent(in)     :: time                      !< current time of the simulation
    integer(kind=ik), intent(in)  :: iteration
    real(kind=rk), intent(in)     :: hvy_block(:, :, :, :, :)  !< heavy data array contains the block data of the statevector
    integer(kind=ik), intent(in)  :: tree_ID
    real(kind=rk), intent(out)    :: dt                         !< time step dt
    !--------------------------------------------------------------
    ! MPI error variable
    integer(kind=ik) :: ierr, Jmax, k, lgt_id, hvy_id
    reaL(kind=rk ) :: ddx(1:3), xx0(1:3), dt_tmp

    dt = 9.0e9_rk

    ! is there a fixed timestep set?
    if (params%dt_fixed > 0.0) then
       dt = params%dt_fixed
    else
      ! --------------------------------------------------------------------------
      ! physics module restrictions on time step
      ! --------------------------------------------------------------------------
      do k = 1, hvy_n(tree_ID)
          hvy_id = hvy_active(k, tree_ID)
          ! as the local CFL condition depends on the blocks, we give the routine
          ! the block grid (origin/spacing)
          call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
          call get_block_spacing_origin( params, lgt_id, xx0, ddx )

          ! physics modules dictate some restrictions due to CFL conditions, penalization
          ! or other operators. Everything that is physics-dependent goes here. it is
          ! computed for each block, then the minimum is used.
          call GET_DT_BLOCK_meta( params%physics_type, time, iteration, hvy_block(:,:,:,:,hvy_id), &
              params%Bs,params%n_ghosts, xx0, ddx, dt_tmp)

          dt = min( dt, dt_tmp )
      end do

      ! synchronize time steps
      ! store local (per process) time step
      dt_tmp = dt
      ! global minimum time step
      call MPI_Allreduce(dt_tmp, dt, 1, MPI_REAL8, MPI_MIN, WABBIT_COMM, ierr)

      ! --------------------------------------------------------------------------
      ! other constraints (saving, final time, etc.)
      ! --------------------------------------------------------------------------
      ! is there an upper limit for the time step set in parameter file?
      if (params%dt_max > 0.0) dt = min( params%dt_max, dt)
    endif


    if ( params%write_method == 'fixed_time' ) then
        ! time step should also fit in output time step size
        ! criterion: check in time+dt above next output time
        ! if so: truncate time+dt
        if ( time+dt > params%next_write_time .and. time<params%next_write_time ) then
            dt = params%next_write_time - time
        end if
    end if


    if ( abs(params%tsave_stats-9999999.9_rk)>1e-1_rk ) then
        ! time step should also fit in statistics output time step size
        ! criterion: check in time+dt above next output time
        ! if so: truncate time+dt
        if ( time+dt > params%next_stats_time .and. time<params%next_stats_time ) then
            dt = params%next_stats_time - time
        end if
    end if


    ! do not jump past final time
    if (time + dt > params%time_max .and. time<=params%time_max) dt = params%time_max - time

    ! --------------------------------------------------------------------------
    ! log time step to accii file
    ! --------------------------------------------------------------------------
    call append_t_file('dt.t', (/time, dt/))

end subroutine calculate_time_step
