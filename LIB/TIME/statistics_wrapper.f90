!> \brief wrapper for RHS call in time step function, computes RHS in work array
!! (inplace)
!!
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
!**********************************************************************************************

subroutine statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, &
    lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor)

   implicit none

    real(kind=rk), intent(in)           :: time, dt
    type (type_params), intent(in)      :: params                     !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)     !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)            !> light data array
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)    !> hvy_mask are qty that depend on the grid and not explicitly on time
    integer(kind=ik), intent(inout)     :: hvy_active(:,:)            !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n(:)                   !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)            !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)                   !> number of active blocks (light data)
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:,:)   !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)          !> heavy data array - neighbor data
    real(kind=rk), dimension(3)         :: dx, x0                     !> spacing and origin of a block
    integer(kind=ik)                    :: k,  lgt_id                 ! loop variables
    integer(kind=ik)                    :: g                          ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs

    ! grid parameter
    Bs    = params%Bs
    g     = params%n_ghosts

    call create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
    hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist)

    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! performs initializations in the RHS module, such as resetting integrals
    call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, hvy_active(1,tree_ID_flow)), g, x0, dx,&
        hvy_tmp(:,:,:,:,hvy_active(1,tree_ID_flow)), "init_stage", hvy_mask(:,:,:,:, hvy_active(1,tree_ID_flow)))

    !-------------------------------------------------------------------------
    ! 2nd stage: integral_stage. (called for all blocks)
    !-------------------------------------------------------------------------
    ! For some RHS, the eqn depend not only on local, block based qtys, such as
    ! the state vector, but also on the entire grid, for example to compute a
    ! global forcing term (e.g. in FSI the forces on bodies). As the physics
    ! modules cannot see the grid, (they only see blocks), in order to encapsulate
    ! them nicer, two RHS stages have to be defined: integral / local stage.
    do k = 1, hvy_n(tree_ID_flow)
      ! convert given hvy_id to lgt_id for block spacing routine
      call hvy2lgt( lgt_id, hvy_active(k,tree_ID_flow), params%rank, params%number_blocks )
      ! get block spacing for RHS
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

      call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, hvy_active(k,tree_ID_flow)), g, x0, dx,&
          hvy_tmp(:,:,:,:,hvy_active(k,tree_ID_flow)), "integral_stage", hvy_mask(:,:,:,:, hvy_active(k,tree_ID_flow)))
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used for example for MPI_REDUCES
    call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, hvy_active(1,tree_ID_flow)), g, x0, dx,&
        hvy_tmp(:,:,:,:,hvy_active(1,tree_ID_flow)), "post_stage", hvy_mask(:,:,:,:, hvy_active(1,tree_ID_flow)))


end subroutine statistics_wrapper
