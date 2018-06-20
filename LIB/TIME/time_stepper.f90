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
!! Runge-Kutta: \n
!! data_field(t) = data_field(t) + sum(b_j*k_j)
!! with k_j = RHS(t+dt*c_j, datafield(t) + dt*sum(a_ji*k_i))         \n
!!
!! \image html time_step.svg "Time-stepper" width=300
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
!
! = log ======================================================================================
!
!> \date 08/11/16 - switch to v0.4 \n
!! \date 07/12/16 - now uses heavy work data array and work for different physics \n
!! \date 31/01/17 - switch to 3D, v0.5 \n
!! \date 23/05/17 - new structure for time_stepper, now works for any explicit Runge Kutta method
!! (up to RK of order 4)
! ********************************************************************************************

subroutine time_stepper(time, params, lgt_block, hvy_block, hvy_work, &
    hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n)
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
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    ! loop variables
    integer(kind=ik)                    :: k, j, neq
    ! time step, dx
    real(kind=rk)                       :: dt
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, sub_t1, t_sum, t
    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable          :: rk_coeffs(:,:)
    logical::test
!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    t0 = MPI_Wtime()
    t_sum = 0.0_rk
    neq = params%number_data_fields

    allocate(rk_coeffs(size(params%butcher_tableau,1),size(params%butcher_tableau,2)) )
    dt = 9.0e9_rk
    ! set rk_coeffs
    rk_coeffs = params%butcher_tableau

!---------------------------------------------------------------------------------------------
! main body

    ! synchronize ghost nodes
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    ! ----------------------------------------------------------------------------------------
    ! calculate time step
    call calculate_time_step(params, time, hvy_block, hvy_active, hvy_n, lgt_block, &
        lgt_active, lgt_n, dt)

    ! first stage, call to RHS. note the resulting RHS is stored in hvy_work(), first
    ! slots after the copy of the state vector, which is in the first 1:neq slots
    j = 1
    call RHS_wrapper(time + dt*rk_coeffs(1,1), params, hvy_block(:,:,:,1:neq,:),&
        hvy_work(:,:,:,j*neq+1:(j+1)*neq,:), lgt_block, hvy_active, hvy_n)

    ! save data at time t to heavy work array
    ! copy state vector content to work array. NOTE: 09/04/2018: moved this after RHS_wrapper
    ! since we can allow the RHS wrapper to modify the state vector (eg for mean flow fixing)
    ! if the copy part is above, the changes in state vector are ignored
    do k = 1, hvy_n
      hvy_work( :, :, :, 1:neq, hvy_active(k) ) = hvy_block( :, :, :, 1:neq, hvy_active(k) )
    end do


    ! compute k_1, k_2, .... (coefficients for final stage)
    do j = 2, size(rk_coeffs, 1)-1
        ! prepare input for the RK substep
        call set_RK_input(dt, params, rk_coeffs(j,:), j, hvy_block, hvy_work, hvy_active, hvy_n)

        ! synchronize ghost nodes for new input
        ! further ghost nodes synchronization, fixed grid
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        ! note substeps are at different times, use temporary time "t"
        t = time + dt*rk_coeffs(j,1)
        call RHS_wrapper(t, params, hvy_block(:,:,:,1:neq,:), &
            hvy_work(:,:,:,j*neq+1:(j+1)*neq,:), lgt_block, hvy_active, hvy_n)
    end do

    ! final stage
    call final_stage_RK(params, dt, hvy_work, hvy_block, hvy_active, hvy_n, rk_coeffs)


    ! increase time variable after all RHS substeps
    time = time + dt
    deallocate(rk_coeffs )


    call toc( params, "time_step (everything incl ghosts)", MPI_Wtime()-t0)
end subroutine time_stepper
