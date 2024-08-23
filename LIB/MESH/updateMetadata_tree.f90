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
subroutine updateMetadata_tree(params, tree_ID, update_lists, update_neighbors, update_family, search_overlapping, verbose_check)
    implicit none

    type (type_params), intent(in) :: params              !< good ol' params
    integer(kind=ik), intent(in)   :: tree_ID             !< Tree to look at
    logical, optional, intent(in)  :: update_lists        !< flag if active and sorted lists should be updated, defaults to .true.
    logical, optional, intent(in)  :: update_neighbors    !< flag if neighbours should be updated, defaults to .true.
    logical, optional, intent(in)  :: update_family       !< flag if family should be updated, defaults to .true.
    logical, optional, intent(in)  :: search_overlapping  !< for CVS multiple neighbors can coexist, so we search all of them
    logical, optional, intent(in)  :: verbose_check       !< Output verbose flag

    integer(kind=ik) :: k_b, lgt_ID
    real(kind=rk) :: t_full, t_part
    logical u_l, u_n, u_f, searchOverlapping
    t_full = MPI_wtime()

    u_l = .true.
    u_n = .true.
    u_f = .true.
    searchOverlapping = .false.
    if (present(update_lists)) u_l = update_lists
    if (present(update_neighbors)) u_n = update_neighbors
    if (present(update_family)) u_f = update_family
    if (present(search_overlapping)) searchOverlapping = search_overlapping

    if (u_l) call createActiveSortedLists_tree(params, tree_ID)

    ! if (present(verbose_check) .and. params%rank == 0) then
    !     if (params%rank == 0) then
    !         do k_b = 1, lgt_n(tree_ID)
    !             lgt_ID = lgt_active(k_b, tree_ID)
    !             write(*, '("2 - R0 - Exists BL-", i0, " L-", i0, " Ref-", i0, " TC-", i0, " - ", b32.32)') lgt_ID, &
    !                 lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2), lgt_block(lgt_id, IDX_TC_2)
    !         enddo
    !         do k_b = 1, lgt_n(tree_ID)
    !             lgt_ID = lgt_sortednumlist(1, k_b, tree_ID)
    !             write(*, '("3 - R0 - Exists BL-", i0, " L-", i0, " Ref-", i0, " TC-List-", i0, " - ", b32.32, " TC-", i0, " - ", b32.32)') lgt_ID, &
    !                 lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), &
    !                 lgt_sortednumlist(3, k_b, tree_ID), lgt_sortednumlist(3, k_b, tree_ID), &
    !                 lgt_block(lgt_id, IDX_TC_2), lgt_block(lgt_id, IDX_TC_2)
    !         enddo
    !     endif
    ! endif

    if (u_n) then
        t_part = MPI_wtime()
        call updateNeighbors_tree(params, tree_ID, search_overlapping=searchOverlapping, verbose_check=verbose_check)
        call toc( "updateMetadata_tree (neighbors)", 57, MPI_wtime()-t_part )
    endif

    if (u_f) then
        t_part = MPI_wtime()
        call updateFamily_tree(params, tree_ID)
        call toc( "updateMetadata_tree (family)", 58, MPI_wtime()-t_part )
    endif

    call toc( "updateMetadata_tree (TOTAL)", 59, MPI_wtime()-t_full )
end subroutine
