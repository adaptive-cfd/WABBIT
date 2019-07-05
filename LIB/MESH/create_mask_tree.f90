subroutine create_mask_tree(params, time, lgt_block, hvy_mask, &
    hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist)
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    real(kind=rk), intent(in)           :: time
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
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
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)  :: lgt_sortednumlist(:,:)

    integer(kind=ik) :: k, lgt_id, Bs(1:3), g
    real(kind=rk) :: x0(1:3), dx(1:3)

    Bs    = params%Bs
    g     = params%n_ghosts

    ! create "time-dependent-part" here, add the existing "time-independent-part"
    ! if it is available, return the complete mask incl. all parts

    do k = 1, hvy_n
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, hvy_mask(:,:,:,:,hvy_active(k)), &
        "all-parts" )
    enddo


end subroutine


! subroutine create_mask_block(, requested_mask)
!
!     requested_mask == "time-dependent-part"
!     requested_mask == "time-independent-part"
!     requested_mask == "all-parts"
!
! end subroutine



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
! Updating those grid-depend quantities is a task for the physics modules: they should provide
! interfaces, if they require such qantities. In many cases, the grid_qtys are probably not used.
! Please note that in the current implementation, hvy_tmp also plays the role of a work array
!-------------------------------------------------------------------------------
