
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \file
!> \callgraph
!> \brief wrapper for RHS call in time step function, computes RHS in work array
!! (inplace)
!> \version 0.5
!> \author sm
!! \date 23/05/17 - create
!!
!
!>\details
!! calls RHS depending on physics
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
!
!**********************************************************************************************

subroutine RHS_wrapper(time, params, hvy_state, hvy_rhs, lgt_block, hvy_active, hvy_n)

!----------------------------------------------------------------------------------------------
! modules

!----------------------------------------------------------------------------------------------
! variables

   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure, hvy_active
    type (type_params), intent(in)      :: params
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_rhs(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(in)           :: hvy_state(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> global integral
    real(kind=rk), dimension(3)         :: volume_int

    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, neqn, lgt_id
    ! grid parameter, error variable
    integer(kind=ik)                    :: Bs, g


!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes
!---------------------------------------------------------------------------------------------
! main body
    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! performs initializations in the RHS module, such as resetting integrals
    select case(params%physics_type)
    case ("ACM-new")
      ! this call is not done for all blocks, but only once, globally.
      call RHS_ACM( time, hvy_state(:,:,:,:,hvy_active(1)), g, x0, dx, hvy_rhs(:,:,:,:,hvy_active(1)), "init_stage" )

    case ("ConvDiff-new")
      ! this call is not done for all blocks, but only once, globally.
      call RHS_convdiff( time, hvy_state(:,:,:,:, hvy_active(1)), g, x0, dx, hvy_rhs(:,:,:,:,hvy_active(1)), "init_stage" )
    case ("navier_stokes")
      ! this call is not done for all blocks, but only once, globally.
      call RHS_NStokes( time, hvy_state(:,:,:,:, hvy_active(1)), g, x0, dx, hvy_rhs(:,:,:,:,hvy_active(1)), "init_stage" )


    case default
      call abort(2152000, "[RHS_wrapper.f90]: physics_type is unknown"//params%physics_type)

    end select

    !-------------------------------------------------------------------------
    ! 2nd stage: integral_stage. (called for all blocks)
    !-------------------------------------------------------------------------
    ! For some RHS, the eqn depend not only on local, block based qtys, such as
    ! the state vector, but also on the entire grid, for example to compute a
    ! global forcing term (e.g. in FSI the forces on bodies). As the physics
    ! modules cannot see the grid, (they only see blocks), in order to encapsulate
    ! them nicer, two RHS stages have to be defined: integral / local stage.
    do k = 1, hvy_n
      ! convert given hvy_id to lgt_id for block spacing routine
      call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
      ! get block spacing for RHS
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

      !---------- different physics modules ----------
      select case(params%physics_type)
      case ("ACM-new")
        ! input state vector: hvy_block, output RHS vector: hvy_work
        call RHS_ACM( time, hvy_state(:,:,:,:,hvy_active(k)), g, &
             x0, dx, hvy_rhs(:,:,:,:, hvy_active(k)), "integral_stage" )

      case ("ConvDiff-new")
        ! input state vector: hvy_block, output RHS vector: hvy_work
        call RHS_convdiff( time, hvy_state(:,:,:,:,hvy_active(k)), g, &
             x0, dx, hvy_rhs(:,:,:,:,hvy_active(k)), "integral_stage" )

      case ("navier_stokes")
        ! input state vector: hvy_block, output RHS vector: hvy_work
        call RHS_NStokes( time, hvy_state(:,:,:,:,hvy_active(k)), g, &
             x0, dx, hvy_rhs(:,:,:,:,hvy_active(k)), "integral_stage" )

      case default
        call abort(2152000, "[RHS_wrapper.f90]: physics_type is unknown"//params%physics_type)

      end select
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used ror example for MPI_REDUCES
    select case(params%physics_type)
    case ("ACM-new")
      ! this call is not done for all blocks, but only once, globally.
      call RHS_ACM( time, hvy_state(:,:,:,:,hvy_active(1)), g, &
      x0, dx, hvy_rhs(:,:,:,:,hvy_active(1)), "post_stage" )

    case ("ConvDiff-new")
      ! this call is not done for all blocks, but only once, globally.
      call RHS_ConvDiff( time, hvy_state(:,:,:,:,hvy_active(1)), g, &
      x0, dx, hvy_rhs(:,:,:,:,hvy_active(1)), "post_stage" )

    case ("navier_stokes")
      ! this call is not done for all blocks, but only once, globally.
      call RHS_NStokes( time, hvy_state(:,:,:,:,hvy_active(1)), g, &
      x0, dx, hvy_rhs(:,:,:,:,hvy_active(1)), "post_stage" )

    case default
      call abort(2152000, "[RHS_wrapper.f90]: physics_type is unknown"//params%physics_type)
    end select


    !-------------------------------------------------------------------------
    ! 4th stage: local evaluation of RHS on all blocks (called for all blocks)
    !-------------------------------------------------------------------------
    ! the second stage then is what you would usually do: evaluate local differential
    ! operators etc.
    do k = 1, hvy_n
      ! convert given hvy_id to lgt_id for block spacing routine
      call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
      ! get block spacing for RHS
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

      !---------- different physics modules ----------
      select case(params%physics_type)
      case ("ACM-new")
        ! input state vector: hvy_block, output RHS vector: hvy_work
        call RHS_ACM( time, hvy_state(:,:,:,:, hvy_active(k)), g, &
             x0, dx, hvy_rhs(:,:,:,:, hvy_active(k)), "local_stage" )

      case ("ConvDiff-new")
        ! input state vector: hvy_block, output RHS vector: hvy_work
        call RHS_ConvDiff( time, hvy_state(:,:,:,:, hvy_active(k)), g, &
             x0, dx, hvy_rhs(:,:,:,:, hvy_active(k)), "local_stage" )

      case ("navier_stokes")
        ! input state vector: hvy_block, output RHS vector: hvy_work
        call RHS_NStokes( time, hvy_state(:,:,:,:, hvy_active(k)), g, &
             x0, dx, hvy_rhs(:,:,:,:, hvy_active(k)), "local_stage" )

      case default
        call abort(2152000, "Error [RHS_wrapper.f90]: physics_type is unknown"//params%physics_type)
      end select
    enddo


end subroutine RHS_wrapper
