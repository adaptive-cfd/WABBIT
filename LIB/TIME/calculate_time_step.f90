!> \brief Get the time step according to the conditions of the physics module
!> \details
!> \author sm
!> \date 18/04/17 - create
!> \date 03/09/18 - removed unused variables (P.Krah) commit 0b79d945422eedda70fa5f874cd2889a10ef8287
subroutine calculate_time_step( params, time, hvy_block, hvy_active, hvy_n, lgt_block, lgt_active, lgt_n, dt )

    use module_physics_metamodule, only : GET_DT_BLOCK

    !--------------------------------------------------------------
    implicit none
    type (type_params), intent(in):: params                    !< user defined parameter structure
    real(kind=rk), intent(in)     :: time                      !< current time of the simulation
    real(kind=rk), intent(in)     :: hvy_block(:, :, :, :, :)  !< heavy data array contains the block data of the statevector
    integer(kind=ik), intent(in)  :: hvy_active(:),hvy_n       !< list of active blocks (heavy data) and number of active blocks
    integer(kind=ik), intent(in)  :: lgt_block(:, :),lgt_active(:),lgt_n!< light data array,active list, number of active blocks
    real(kind=rk), intent(out)    :: dt                         !< time step dt
    !--------------------------------------------------------------
    ! MPI error variable
    integer(kind=ik) :: ierr, Jmax, k, lgt_id
    reaL(kind=rk ) :: ddx(1:3), xx0(1:3), dt_tmp

    dt = 9.0e9_rk

    ! is there a fixed timestep set?
    if (params%dt_fixed > 0.0) then
       dt = params%dt_fixed
    else
      ! --------------------------------------------------------------------------
      ! physics module restrictions on time step
      ! --------------------------------------------------------------------------
      do k = 1, hvy_n
          ! as the local CFL condition depends on the blocks, we give the routine
          ! the block grid (origin/spacing)
          call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
          call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

          ! physics modules dictate some restrictions due to CFL conditions, penalization
          ! or other operators. Everything that is physics-dependent goes here. it is
          ! computed for each block, then the minimum is used.
          call GET_DT_BLOCK( params%physics_type, time, hvy_block(:,:,:,:,hvy_active(k)), &
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

    if (dt <= 0.0_rk) then
      call abort(12131,"For some reason, we ended up with a negative or zero time step. This is not back to the future!!!")
    endif

    ! --------------------------------------------------------------------------
    ! log time step to accii file
    ! --------------------------------------------------------------------------
    if (params%rank==0) then
        open(14,file='dt.t',status='unknown',position='append')
        write (14,'(2(g15.8,1x))') time, dt
        close(14)
    endif

end subroutine calculate_time_step
