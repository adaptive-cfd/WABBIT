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

subroutine statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, tree_ID)

    implicit none

    real(kind=rk), intent(in)           :: time, dt
    type (type_params), intent(inout)   :: params                     !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)     !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)    !> hvy_mask are qty that depend on the grid and not explicitly on time
    integer(kind=ik), intent(in)        :: tree_ID
    real(kind=rk), dimension(3)         :: dx, x0                     !> spacing and origin of a block
    integer(kind=ik)                    :: k,  lgt_id, hvy_id, hvy_id_mask         ! loop variables
    integer(kind=ik)                    :: g                          ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs

    ! grid parameter
    Bs     = params%Bs
    g      = params%g
    dx     = 0.0_rk
    x0     = 0.0_rk
    hvy_id = 1
    

    call createMask_tree(params, time, hvy_mask, hvy_tmp)

    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! performs initializations in the RHS module, such as resetting integrals
    ! NOTE: this is block-independent and therefore we pass dx=x0=0
    call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, 1), g, x0, dx,&
    hvy_tmp(:,:,:,:,1), "init_stage", hvy_mask(:,:,:,:, 1))

    !-------------------------------------------------------------------------
    ! 2nd stage: integral_stage. (called for all blocks)
    !-------------------------------------------------------------------------
    ! For some RHS, the eqn depend not only on local, block based qtys, such as
    ! the state vector, but also on the entire grid, for example to compute a
    ! global forcing term (e.g. in FSI the forces on bodies). As the physics
    ! modules cannot see the grid, (they only see blocks), in order to encapsulate
    ! them nicer, two RHS stages have to be defined: integral / local stage.
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
        ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
        ! Note: hvy_mask is not used if in this case by the RHS routines...
        if (size(hvy_mask,5) == 1) then
            hvy_id_mask = 1
        else 
            hvy_id_mask = hvy_id
        endif

        call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, hvy_id), g, x0, dx,&
        hvy_tmp(:,:,:,:,hvy_id), "integral_stage", hvy_mask(:,:,:,:, hvy_id_mask))
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used for example for MPI_REDUCES
    ! NOTE: this is block-independent and therefore we pass dx=x0=0
    call STATISTICS_meta(params%physics_type, time, dt, hvy_block(:,:,:,:, 1), g, x0, dx,&
    hvy_tmp(:,:,:,:,1), "post_stage", hvy_mask(:,:,:,:, 1))


end subroutine statistics_wrapper
