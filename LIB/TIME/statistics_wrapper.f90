
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

subroutine statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, lgt_block, hvy_active, hvy_n)

!----------------------------------------------------------------------------------------------
! modules

!----------------------------------------------------------------------------------------------
! variables

   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time, dt
    !> user defined parameter structure, hvy_active
    type (type_params), intent(in)      :: params
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n


    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k,  lgt_id
    ! grid parameter
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3) :: Bs


!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%Bs
    g     = params%n_ghosts
!---------------------------------------------------------------------------------------------
! main body

    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! performs initializations in the RHS module, such as resetting integrals
    call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, hvy_active(1)), g, x0, dx,&
        hvy_tmp(:,:,:,:,hvy_active(1)), "init_stage")

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

      call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, hvy_active(k)), g, x0, dx,&
          hvy_tmp(:,:,:,:,hvy_active(k)), "integral_stage")
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used ror example for MPI_REDUCES
    call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, hvy_active(1)), g, x0, dx,&
        hvy_tmp(:,:,:,:,hvy_active(1)), "post_stage")


end subroutine statistics_wrapper
