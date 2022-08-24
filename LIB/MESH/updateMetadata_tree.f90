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
subroutine updateMetadata_tree(params, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk)                       :: t0
    t0 = MPI_wtime()

    call createActiveSortedLists_tree(params, tree_ID)

    call updateNeighbors_tree(params, tree_ID)

    call toc( "updateMetadata_tree (lists+neighbors)", MPI_wtime()-t0 )
end subroutine
