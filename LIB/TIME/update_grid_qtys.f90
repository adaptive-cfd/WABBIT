!-------------------------------------------------------------------------------
! While the state vector and many work variables (such as the mask function for penalization)
! are explicitly time dependent, some other quantities are not. They are rather grid-dependent
! but need not to be updated in every RK or krylov substep. Hence, those quantities are updated
! after the mesh is changed (i.e. after refine_mesh) and then kept constant during the evolution
! time step.
! An example for such a quantity would be geometry factors on non-cartesian grids, but also the
! body of an insect in tethered (=fixed) flight. In the latter example, only the wings need to be
! generated at every time t. This example generalizes to any combination of stationary and moving
! obstacle, i.e. insect behind fractal tree.
!-------------------------------------------------------------------------------
subroutine update_grid_qyts( time, params, lgt_block, hvy_gridQ, hvy_active, hvy_n )
    !> even though it is a bit odd, since those qtys shall be TIME INDEPENDENT, we pass time for debugging
    real(kind=rk), intent(in)              :: time
    !> user defined parameter structure
    type (type_params), intent(in)         :: params
    !> light data array
    integer(kind=ik), intent(inout)        :: lgt_block(:, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)           :: hvy_gridQ(:, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_n

    ! grid parameter
    integer(kind=ik) :: Bs, g, k, lgt_id
    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0


    ! grid parameter
    Bs = params%Bs
    g  = params%n_ghosts

    ! stage 1: initialization. called once and not for every block
    call UPDATE_GRID_QTYS_meta( time, params%physics_type, hvy_gridQ(:, :, :, :, hvy_active(1)), &
    g, x0, dx, "init_stage" )

    ! stage 2: local stage, called for each block.
    do k = 1, hvy_n
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        call UPDATE_GRID_QTYS_meta( time, params%physics_type, hvy_gridQ(:, :, :, :, hvy_active(k)), &
        g, x0, dx, "main_stage" )
    enddo

end subroutine
