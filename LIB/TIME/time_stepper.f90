!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name time_stepper.f90
!> \version 0.5
!> \author msr, sm
!
!> \brief time step main function
!
!>
!! data_field(t) = data_field(t) + sum(b_j*k_j) 
!! with k_j = RHS(t+dt*c_j, datafield(t) + dt*sum(a_ji*k_i))         \n
!!
!! input:    
!!           - time variable
!!           - params
!!           - light and heavy data
!!           - neighbor list 
!!
!! output:   
!!           - time variable
!!           - heavy data array
!!
!!
!! physics: \n
! --------
!> - convection/diffusion: works only for one datafield \n
!!
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!!
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 07/12/16 - now uses heavy work data array and work for different physics \n
!! 31/01/17 - switch to 3D, v0.5 \n
!! 23/05/17 - new structure for time_stepper, now works for any explicit Runge Kutta method (up to RK of order 4)
!
! ********************************************************************************************

subroutine time_stepper( time, params, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> time varible
    real(kind=rk), intent(inout)        :: time

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)

    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)

    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)         :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, lgt_id, d, j

    ! time step, dx
    real(kind=rk)                       :: dt, dx, my_dx
    ! spacing and origin of a block
    real(kind=rk)                       :: xx0(1:3), ddx(1:3)
    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! new time after timestep dt, so current time can be used in RHS
    real(kind=rk)                       :: time_dt

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable          :: rk_coeffs(:,:)

    integer(kind=ik)                    :: allocate_error

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    sub_t0 = MPI_Wtime()

    allocate(rk_coeffs(size(params%butcher_tableau,1),size(params%butcher_tableau,2)), stat=allocate_error)
    call check_allocation(allocate_error)

    ! reset dx
    my_dx = 9.0e9_rk
    dx    = 9.0e9_rk
    ! reset dt
    dt    = 9.0e9_rk

    time_sum = 0.0_rk

    ! set MPI parameter
    rank         = params%rank
    d = 2
    if ( params%threeD_case) d = 3

    ! set rk_coeffs
    rk_coeffs = params%butcher_tableau

!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! calculate time step
    ! loop over all active blocks (heavy data)
    ! FIXME: you could also look over light data, as ddx is available only from that. no mpi
    do k = 1, hvy_n
        ! light id of this block
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
        ! compute blocks' spacing from treecode
        call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )
        ! find smallest dx of active blocks
        my_dx = min(my_dx, minval(ddx(1:d)) )
    end do

    ! find globally smallest dx
    call MPI_Allreduce(my_dx, dx, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! calculate dt
    call calculate_time_step(params, dx, dt)

    ! calculate value for time after one dt
    time_dt = time + dt
    ! last timestep should fit in maximal time
    if ( time_dt >= params%time_max) then
        dt = params%time_max - time
        time_dt = params%time_max
    end if

    ! time measurement without ghost nodes synchronization
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)

    ! synchronize ghost nodes
    ! first ghost nodes synchronization, so grid has changed
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

    ! restart time
    sub_t0 = MPI_Wtime()
    
    ! save data at time t
    call save_data_t(params, hvy_work, hvy_block, hvy_active, hvy_n)

    call RHS_wrapper(time, dt, params, hvy_work, rk_coeffs(1,1), 1, lgt_block, hvy_active, hvy_n, hvy_block)

    ! compute k_1, k_2, .... (coefficients for final stage)
    do j = 2, size(rk_coeffs, 1)-1

        call set_RK_input(dt, params, rk_coeffs(j,:), j, hvy_block, hvy_work, hvy_active, hvy_n)
        
        ! time measurement without ghost nodes synchronization
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        sub_t1   = MPI_Wtime()
        time_sum = time_sum + (sub_t1 - sub_t0)

        ! synchronize ghost nodes for new input
        ! further ghost nodes synchronization, fixed grid
        call synchronize_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .false., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer)

        ! restart time
        sub_t0 = MPI_Wtime()

        call RHS_wrapper(time, dt, params, hvy_work, rk_coeffs(j,1), j, lgt_block, hvy_active, hvy_n, hvy_block)

    end do

    ! final stage
    call final_stage_RK(params, dt, hvy_work, hvy_block, hvy_active, hvy_n, rk_coeffs)
    
    ! increase time variable after all RHS substeps
    time = time_dt

    ! end time
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)
    
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "time_step (w/o ghost synch.)" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "time_step (w/o ghost synch.)"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + time_sum
    end if

    deallocate(rk_coeffs, stat=allocate_error)

end subroutine time_stepper
