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
subroutine updateMetadata_tree(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)          !> light data array
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)        !> heavy data array - neighbor data
    integer(kind=ik), intent(out)       :: lgt_active(:,:)          !> list of active blocks (light data)
    integer(kind=ik), intent(out)       :: lgt_n(:)                 !> number of active blocks (light data)
    integer(kind=tsize), intent(out)    :: lgt_sortednumlist(:,:,:) !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(out)       :: hvy_active(:,:)          !> list of active blocks (heavy data)
    integer(kind=ik), intent(out)       :: hvy_n(:)                 !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk)                       :: t0
    t0 = MPI_wtime()

    call createActiveSortedLists_tree( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID)

    call updateNeighbors_tree(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID)

    call toc( "updateMetadata_tree (lists+neighbors)", MPI_wtime()-t0 )
end subroutine
