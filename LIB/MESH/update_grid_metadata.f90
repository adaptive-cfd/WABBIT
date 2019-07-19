! If the grid has changed (i.e. by adding and removing blocks, or loadbalancing)
! the meta-lists have to be updated on all mpiranks
!
!   This routine updates:
!   - the lists of currently active blocks
!   - the list used for block finding
!   - the neighboring information
!
! NOTE: we require the light data lgt_block to be synchronized BEFORE calling this routine.
!       No synchronization step is required afterwards.
!
subroutine update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, tree_ID)
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(out)       :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(out)       :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(out)    :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(out)       :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(out)       :: hvy_n
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk) :: t0

    t0 = MPI_wtime()

    call create_active_and_sorted_lists_tree( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID)

    call update_neighbors(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n)

    call toc( "update_grid_metadata (lists+neighbors)", MPI_wtime()-t0 )
end subroutine
