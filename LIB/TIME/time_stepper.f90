!! Runge-Kutta: \n
!! data_field(t) = data_field(t) + sum(b_j*k_j)
!! with k_j = RHS(t+dt*c_j, datafield(t) + dt*sum(a_ji*k_i))         \n
!!
!! \image html time_step.svg "Time-stepper" width=300
!!
!! input:    - time variable
!!           - params
!!           - light and heavy data
!!           - neighbor list
!! output:   - time variable
!!           - heavy data array
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
! ********************************************************************************************

subroutine time_stepper(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, hvy_mask, hvy_tmp, &
    hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist, tree_id)
    implicit none

    real(kind=rk), intent(inout)        :: time, dt
    integer(kind=ik), intent(inout)     :: iteration
    type (type_params), intent(in)      :: params                       !> user defined parameter structure
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)              !> light data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)   !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in), optional :: tree_id                   !> is only needed for testing purposes
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)            !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_active(:,:)              !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)              !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_n(:)                     !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)                     !> number of active blocks (light data)
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:,:)     !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik)                    :: k,  tree_id_evolv

    ! currently not working (Thomas, 02-2021)
    ! call update_neighbors(params, lgt_block, hvy_neighbor, lgt_active(:,tree_ID_evolv), lgt_n(tree_ID_evolv), &
    ! lgt_sortednumlist(:,:,tree_ID_evolv), hvy_active(:,tree_ID_evolv), hvy_n(tree_ID_evolv), skip_diagonal_neighbors=.true.)


    if (present(tree_id)) then
          tree_id_evolv = tree_id
    else
          tree_id_evolv = tree_ID_flow
    end if

    if ( .not. All(params%periodic_BC) ) then
        !!! if we have NON-PERIODIC boundary conditions it is important to reset hvy_work.
        !!! this is important because hvy_work saves the RHS also in the ghost node layer of the
        !!! boundary blocks which is not synchronized. if RHS would be not 0 in the ghost node layer
        !!! then the integrator would change the values in the ghost node layer.
        do k = 1, hvy_n(tree_ID_evolv)
            hvy_work(:, :, :, :, hvy_active(k,tree_ID_evolv), :) = 0.0_rk
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
            hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist, tree_id_evolv)

    case("RungeKuttaGeneric-FSI")
        ! FSI versions of RK schemes advance a solid model simultaneously with the fluid. They
        ! are applicable only for ACM module currently
        call RungeKuttaGeneric_FSI(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, &
            hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist)

    case("RungeKuttaChebychev-FSI")
        ! FSI versions of RK schemes advance a solid model simultaneously with the fluid. They
        ! are applicable only for ACM module currently
        call RungeKuttaChebychev_FSI(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, &
            hvy_mask, hvy_tmp, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, lgt_sortednumlist)
    case default
        call abort(19101816, "time_step_method is unkown: "//trim(adjustl(params%time_step_method)))

    end select


! currently not working (Thomas, 02-2021)
! call update_neighbors(params, lgt_block, hvy_neighbor, lgt_active(:,tree_ID_evolv), lgt_n(tree_ID_evolv), &
! lgt_sortednumlist(:,:,tree_ID_evolv), hvy_active(:,tree_ID_evolv), hvy_n(tree_ID_evolv), skip_diagonal_neighbors=.false.)

    ! increase time variable after all RHS substeps
    time = time + dt
    iteration = iteration + 1

end subroutine time_stepper
