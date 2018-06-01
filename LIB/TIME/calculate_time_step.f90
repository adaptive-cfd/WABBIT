!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name calculate_time_step.f90
!> \version 0.5
!> \author sm
!
!> \brief calculate time step
!
!>
!! input:    -  params \n
!! output:   -  time step dt \n
!!
!!
!! = log ======================================================================================
!! \n
!! 18/04/17 - create
!
! ********************************************************************************************
subroutine calculate_time_step( params, time, hvy_block, hvy_active, hvy_n, lgt_block, lgt_active, lgt_n, dt )

    use module_physics_metamodule, only : GET_DT_BLOCK

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    real(kind=rk), intent(in)           :: time

    !> heavy data array - block data
    real(kind=rk), intent(in)           :: hvy_block(:, :, :, :, :)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)

    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> time step dt
    real(kind=rk), intent(out)          :: dt

    ! loop variables
    integer(kind=ik)                    :: dF, N_dF, UxF, UyF, UzF

    ! norm of vector u
    real(kind=rk)                       :: norm_u, my_norm_u

    ! MPI error variable
    integer(kind=ik)                    :: ierr

    ! maximal mesh level
    integer(kind=ik)                    :: Jmax

    reaL(kind=rk ) :: ddx(1:3), xx0(1:3), dt_tmp
    integer(kind=ik) :: k, lgt_id

!---------------------------------------------------------------------------------------------
! variables initialization

    N_dF  = params%number_data_fields

    UxF  = 0
    UyF  = 0
    UzF  = 0

!---------------------------------------------------------------------------------------------
! main body
  dt = 9.0e9_rk

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
            params%number_block_nodes,params%number_ghost_nodes, xx0, ddx, dt_tmp)

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
    ! is there a fixed timestep set?
    if (params%dt_fixed > 0.0) dt = params%dt_fixed

    if ( params%write_method == 'fixed_time' ) then
        ! time step should also fit in output time step size
        ! criterion: check in time+dt above next output time
        ! if so: truncate time+dt
        if ( time+dt > params%next_write_time .and. time<params%next_write_time ) then
            dt = params%next_write_time - time
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
