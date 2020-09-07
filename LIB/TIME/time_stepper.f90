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

subroutine time_stepper(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, hvy_mask, hvy_tmp, &
    hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist)
    implicit none

    !> time varible
    real(kind=rk), intent(inout)        :: time, dt
    integer(kind=ik), intent(inout)     :: iteration
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)

    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
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

    ! loop variables
    integer(kind=ik)                    :: k, j, Neqn, g, z1, z2, hvy_id
    integer(kind=ik), dimension(3) :: Bs
    ! time step, dx
    real(kind=rk)                       :: t
    ! array containing Runge-Kutta coefficients
    real(kind=rk), allocatable, save    :: rk_coeffs(:,:)

!---------------------------------------------------------------------------------------------
! variables initialization

    Neqn = params%n_eqn
    Bs    = params%Bs
    g     = params%n_ghosts

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

!---------------------------------------------------------------------------------------------
! main body

    if ( .not. All(params%periodic_BC) ) then
        !!! if we have NON-PERIODIC boundary conditions it is important to reset hvy_work.
        !!! this is important because hvy_work saves the RHS also in the ghost node layer of the
        !!! boundary blocks which is not synchronized. if RHS would be not 0 in the ghost node layer
        !!! then the integrator would change the values in the ghost node layer.
        do k = 1, hvy_n(tree_ID_flow)
            hvy_work(:, :, :, :, hvy_active(k,tree_ID_flow), :) = 0.0_rk
        enddo
    endif


    select case (params%time_step_method)
    case ("Krylov")
        !-----------------------------------------------------------------------
        ! krylov scheme
        !-----------------------------------------------------------------------
        ! use krylov time stepping
        call krylov_time_stepper(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, &
            hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist)


    case("RungeKuttaChebychev")
        !-----------------------------------------------------------------------
        ! runge-kutta chebychev scheme
        !-----------------------------------------------------------------------
        call RungeKuttaChebychev(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, &
            hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist)


    case("RungeKuttaGeneric")
        !-----------------------------------------------------------------------
        ! runge-kutta scheme
        !-----------------------------------------------------------------------

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

    case default
        call abort(19101816, "time_step_method is unkown: "//trim(adjustl(params%time_step_method)))

    end select

    ! increase time variable after all RHS substeps
    time = time + dt
    iteration = iteration + 1

end subroutine time_stepper
