subroutine addSecurityZone_tree( time, params, level_this, tree_ID, hvy_block, hvy_tmp )

    use module_indicators

    implicit none
    real(kind=rk), intent(in)      :: time
    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    real(kind=rk), intent(inout)   :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)   :: level_this, tree_ID

    integer(kind=ik), parameter :: TMP_STATUS = 17
    integer(kind=ik) :: level_me, refinement_status, neighborhood,refinement_status_neighbor, &
    hvyID, lgtID, k, level_neighbor, lgtID_neighbor
    logical :: refinement

    !--------------------------------------------------------------------------------------------------------------------
    ! This routine performs two tasks:
    !   1) if a block has a significant neighbor on the SAME LEVEL, its (possible) coarsening flag is revoked to 0 (STAY)
    !   2) if a block has a significant neighbor on the FINER LEVEL, it must be refined (this includes gradedness, etc)
    !
    ! NOTES
    !   2) We usually refine everywhere, but here we do not. it's a good thing refine_tree is prepared for that.
    !--------------------------------------------------------------------------------------------------------------------

    ! note this is a parallel loop (because hvy_neighbor is not lgt_neighbor)
    ! note coarseningIndicator_tree synchronizes the refinement status
    ! note it does not matter that this is a hvy_loop: still all relevant blocks get the TMP_STATUS
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        level_me          = lgt_block( lgtID, IDX_MESH_LVL )
        refinement_status = lgt_block( lgtID, IDX_REFINE_STS )

        ! first task
        ! is the block on the level we look at?
        ! is this block assigned -1 (it wants to coarsen)?
        if ((level_me == level_this).and.(refinement_status == -1)) then
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor             = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor             = lgt_block( lgtID_neighbor, IDX_MESH_LVL )
                    refinement_status_neighbor = lgt_block( lgtID_neighbor, IDX_REFINE_STS )

                    ! this block is significant
                    if ((refinement_status_neighbor == 0).and.(level_neighbor == level_this)) then
                        ! revoke coarsening status. note we must use a temp status or otherwise all blocks
                        ! will revoke their coarsening
                        lgt_block( lgtID, IDX_REFINE_STS ) = TMP_STATUS
                    endif
                endif
            enddo
        endif

        ! second task: possible refinement - concerned blocks are one level coarser (J-1).
        if (level_me == level_this-1) then
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor             = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor             = lgt_block( lgtID_neighbor, IDX_MESH_LVL )
                    refinement_status_neighbor = lgt_block( lgtID_neighbor, IDX_REFINE_STS )

                    ! the finer neighboring block is significant:
                    if ((refinement_status_neighbor == 0).and.(level_neighbor == level_this)) then
                        ! this block shall be refined. it may entail other blocks to be refined as well
                        lgt_block( lgtID, IDX_REFINE_STS ) = +1
                        ! write(*,*) "It happened: a new block emerges"
                    endif
                endif
            enddo
        endif
    enddo

    ! finally lop over all blocks and convert 17 to 0
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        if (lgt_block( lgtID, IDX_REFINE_STS ) == TMP_STATUS) then
            lgt_block( lgtID, IDX_REFINE_STS ) = 0
        endif
    enddo

    ! after modifying all refinement flags, we need to synchronize light data
    call synchronize_lgt_data( params,  refinement_status_only=.true. )

    ! the following code is only required if new blocks have to be created (which is rare or even never occurs, it seems...)
    refinement = .false.
    do k = 1, lgt_n(tree_ID)
        if (lgt_block( lgt_active(k, tree_ID), IDX_REFINE_STS ) == +1) then
            refinement = .true.
        endif
    enddo

    ! create new blocks, if necessary.
    ! Note: we bypass "refine_tree" because the refinementIndicator_tree removes coarsening flags
    ! This is incredibly tricky, I need to think about some potential issues with this
    if (refinement) then
#ifdef DEV
        write(*, '("It happened! A new block emerges on rank ", i0)') params%rank
#endif

        call abort(197, "I don't believe in wonders so this case is abolished")

        ! ! creating new blocks is not always possible without creating even more blocks to ensure gradedness
        ! call ensureGradedness_tree( params, tree_ID )

        ! ! actual refinement of new blocks (Note: afterwards, new blocks have refinement_status=0)
        ! call refinementExecute_lvl2MallatWD( params, hvy_block, tree_ID, level_this )

        ! ! grid has changed...
        ! call updateMetadata_tree(params, tree_ID)

        ! ! When blocks on lower levels are created these outdate the neighbours there. If by chance this routine here
        ! ! will be called again and new blocks are created on lower levels, we need synched ghost patches
        ! call sync_ghosts_all( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )        
    endif
end subroutine
