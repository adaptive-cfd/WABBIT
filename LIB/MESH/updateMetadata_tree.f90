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
subroutine updateMetadata_tree(params, tree_ID, update_neighbors, update_family)
    implicit none

    type (type_params), intent(in) :: params            !< good ol' params
    integer(kind=ik), intent(in)   :: tree_ID           !< Tree to look at
    logical, optional, intent(in)  :: update_neighbors  !< flag if neighbours should be updated, defaults to .true.
    logical, optional, intent(in)  :: update_family     !< flag if family should be updated, defaults to .true.

    real(kind=rk) :: t0
    logical u_n, u_f
    t0 = MPI_wtime()

    u_n = .true.
    u_f = .true.
    if (present(update_neighbors)) u_n = update_neighbors
    if (present(update_family)) u_f = update_family

    call createActiveSortedLists_tree(params, tree_ID)
    ! call createActiveSortedLists_tree_old(params, tree_ID)

    if (u_n) call updateNeighbors_tree(params, tree_ID)
    if (u_f) call updateFamily_tree(params, tree_ID)

    call toc( "updateMetadata_tree (lists+neighbors)", MPI_wtime()-t0 )
end subroutine
