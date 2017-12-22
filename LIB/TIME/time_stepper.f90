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
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 07/12/16 - now uses heavy work data array and work for different physics \n
!! 31/01/17 - switch to 3D, v0.5 \n
!! 23/05/17 - new structure for time_stepper, now works for any explicit Runge Kutta method (up to RK of order 4)
!
! ********************************************************************************************

subroutine time_stepper( time, params, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, com_lists, com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

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
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)

    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)

    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)         :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, j, neq

    ! time step, dx
    real(kind=rk)                       :: dt

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, sub_t1, t_sum, t

    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable          :: rk_coeffs(:,:)


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

    ! ----------------------------------------------------------------------------------------
    ! calculate time step
    call calculate_time_step(params, time, hvy_block, hvy_active, hvy_n, lgt_block, lgt_active, lgt_n, dt)


    ! synchronize ghost nodes
    ! first ghost nodes synchronization, so grid has changed
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )


call cp_state_vect_to_work(params, hvy_work, hvy_block, hvy_active, hvy_n)

    ! save data at time t to heavy work array
    ! call save_data_t(params, hvy_work, hvy_block, hvy_active, hvy_n)
    j = 1
    call RHS_wrapper(time + dt*rk_coeffs(1,1), params, hvy_block(:,:,:,1:neq,:), hvy_work(:,:,:,j*neq+1:(j+1)*neq,:), lgt_block, hvy_active, hvy_n)


    ! function u1 =RK4(u0,dt,t,v0,mask,mask_dx,rhs)
    !     % Runge Kutta 4. Ordnung
    !
    !     % die vier Aufrufe der rechten Seiten;
    !     % egal, ob rhs zahl, vektor  oder vektoren=matrix zurueckgibt:
    !     % solange Dimension u entspricht geht alles glatt.
    !
    !     k1=rhs(u0          , t     ,v0,mask,mask_dx);
    !     k2=rhs(u0+dt/2* k1 , t+dt/2,v0,mask,mask_dx);
    !     k3=rhs(u0+dt/2* k2 , t+dt/2,v0,mask,mask_dx);
    !     k4=rhs(u0+dt* k3   , t+dt  ,v0,mask,mask_dx);
    !     u1 = u0+dt/6 * (k1 + 2*k2 +2*k3+k4) ;
    !
    ! end

    ! compute k_1, k_2, .... (coefficients for final stage)
    do j = 2, size(rk_coeffs, 1)-1


        call set_RK_input(dt, params, rk_coeffs(j,:), j, hvy_block, hvy_work, hvy_active, hvy_n)


        ! synchronize ghost nodes for new input
        ! further ghost nodes synchronization, fixed grid
        call synchronize_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .false., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer)

        t = time + dt*rk_coeffs(j,1)
        call RHS_wrapper(t, params, hvy_block(:,:,:,1:neq,:), hvy_work(:,:,:,j*neq+1:(j+1)*neq,:), &
        lgt_block, hvy_active, hvy_n)
    end do

    ! final stage
    call final_stage_RK(params, dt, hvy_work, hvy_block, hvy_active, hvy_n, rk_coeffs)


    !
    ! ! increase time variable after all RHS substeps
    time = time + dt
    !
    ! ! timings
    ! sub_t1   = MPI_Wtime()
    ! t_sum = t_sum + (sub_t1 - t0)
    ! call toc( params, "time_step (w/o ghost synch.)", t_sum)

    deallocate(rk_coeffs )

end subroutine time_stepper
