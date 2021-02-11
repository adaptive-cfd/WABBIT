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

    integer(kind=ik) :: k

    ! currently not working (Thomas, 02-2021)
    ! call update_neighbors(params, lgt_block, hvy_neighbor, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), &
    ! lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), skip_diagonal_neighbors=.true.)

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
        call RungeKuttaGeneric(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, &
            hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist)

    case default
        call abort(19101816, "time_step_method is unkown: "//trim(adjustl(params%time_step_method)))

    end select


! currently not working (Thomas, 02-2021)
! call update_neighbors(params, lgt_block, hvy_neighbor, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), &
! lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), skip_diagonal_neighbors=.false.)

    ! increase time variable after all RHS substeps
    time = time + dt
    iteration = iteration + 1

end subroutine time_stepper
