!> This adds a full buffer zone around significant blocks
subroutine addSecurityZone_level( time, params, level_this, tree_ID, hvy_block, hvy_tmp )

    use module_indicators

    implicit none
    real(kind=rk), intent(in)      :: time
    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)    !> heavy data array - block data in wavelet decomposed form
    real(kind=rk), intent(inout)   :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)   :: level_this, tree_ID

    integer(kind=ik), parameter :: TMP_STATUS = 17
    integer(kind=ik) :: level_me, refinement_status, neighborhood,refinement_status_neighbor, &
    hvyID, lgtID, k, level_neighbor, lgtID_neighbor, idx(2,3)
    logical :: refinement

    ! We are not sure if security zone refinement is valid with coarse extension, as cells can be refined which did not experience coarse extension
    ! Therefore this can be en- and disabled here with a simple switch for now
    logical :: enable_refinement
    enable_refinement = .false.

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
        if (enable_refinement) then
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
    if (enable_refinement) then
        refinement = .false.
        do k = 1, lgt_n(tree_ID)
            if (lgt_block( lgt_active(k, tree_ID), IDX_REFINE_STS ) == +1) then
                refinement = .true.
            endif
        enddo

        ! create new blocks, if necessary.
        ! Note: we bypass "refine_tree" because the refinementIndicator_tree removes coarsening flags
        ! This should be wrong, as it bypasses the CE on the un-refined blocks and directly refines them
        ! As it is really rare, we could either accept the error from this refinement or from WC deletion
        ! Currently, we accept the WC deletion as this only takes out energy instead of skipping CE
        if (refinement) then
            write(*, '("It happened! A new block wanted to emerge on rank ", i0)') params%rank
#ifdef DEV
            ! call abort(197, "I don't believe in wonders, you should check up why this block wanted to be created")
#endif

            ! When blocks on lower levels are refined these need the boundary information from neighbouring patches
            ! However, these are not available as we only sync level-wise
            ! Since this case is super rare we just take the costs and do a full sync here
            call sync_ghosts_all( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )    

            ! creating new blocks is not always possible without creating even more blocks to ensure gradedness
            call ensureGradedness_tree( params, tree_ID )

            ! actual refinement of new blocks (Note: afterwards, new blocks have refinement_status=0)
            call refinementExecute_lvl2SpaghettiWD( params, hvy_block, tree_ID, level_this )

            ! grid has changed...
            call updateMetadata_tree(params, tree_ID)    
        endif
    endif
end subroutine



!> This adds a buffer zone around blocks if the CE would remove valuable WC
!  Problem for now: how do we tell neighbouring blocks on different processors that they revoke their status?
subroutine addSecurityZone_CE_level( time, params, level_this, tree_ID, hvy_block, hvy_tmp, indicator, norm, input_is_WD)

    use module_indicators

    implicit none
    real(kind=rk), intent(in)      :: time
    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)    !< heavy data array - block data might be wavelet decomposed form
    real(kind=rk), intent(inout)   :: hvy_tmp(:, :, :, :, :)      !< not used but we need to pass it to coarseningIdicator_block
    character(len=*), intent(in)   :: indicator                   !< how to choose blocks for refinement
    integer(kind=ik), intent(in)   :: level_this, tree_ID
    real(kind=rk), intent(inout)   :: norm(:)                     !< the computed norm for each component of the vector
    logical, intent(in)            :: input_is_WD                 !< flag if hvy_block is already wavelet decomposed

    integer(kind=ik), parameter :: TMP_STATUS = 17
    integer(kind=ik) :: level_me, ref_status, neighborhood, ref_status_neighbor, ref_check, &
    hvyID, lgtID, k, level_neighbor, lgtID_neighbor, idx(2,3), nx, ny, nz, g, Nwcl, Nwcr
    logical :: refinement, inputIsWD

    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    g = params%g
    Nwcl    = params%Nwcl
    Nwcr    = params%Nwcr
    idx(:, :) = 1

    !--------------------------------------------------------------------------------------------------------------------
    ! This routine performs checks if a block has a significant neighbor on the SAME LEVEL,
    ! then the patch where WC will be deleted will be checked again with coarseningIndicator_block,
    ! if it is signifant then the coarsening flag is revoked to 0 (= STAY)
    !--------------------------------------------------------------------------------------------------------------------

    ! note this is a parallel loop (because hvy_neighbor is not lgt_neighbor)
    ! note coarseningIndicator_tree synchronizes the refinement status
    ! note it does not matter that this is a hvy_loop: still all relevant blocks get the TMP_STATUS
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        level_me          = lgt_block( lgtID, IDX_MESH_LVL )
        ref_status = lgt_block( lgtID, IDX_REFINE_STS )

        ! is the block on the level we look at?
        ! is this block assigned 0 (it wants to stay and is significant)?
        if ((level_me == level_this).and.(ref_status == 0)) then
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor             = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor             = lgt_block( lgtID_neighbor, IDX_MESH_LVL )
                    ref_status_neighbor = lgt_block( lgtID_neighbor, IDX_REFINE_STS )

                    ! the neighbor wants to coarsen
                    if ((ref_status_neighbor == -1).and.(level_neighbor == level_this)) then
                        ! check for the patch where wc will be deleted, make sure to exclude ghost patches
                        call get_indices_of_modify_patch(params, neighborhood, idx, (/ nx, ny, nz/), (/Nwcl, Nwcl, Nwcl/), (/Nwcr, Nwcr, Nwcr/), &
                            X_s=(/ g, g, g/), X_e=(/ g, g, g/))
                        ref_check = -1

                        call coarseningIndicator_block( params, hvy_block(:,:,:,:,hvyID), &
                        hvy_tmp(:,:,:,:,hvyID), indicator, ref_check, norm, level_this, input_is_WD, indices=idx)

                        write(*, '("B", i0, " N", i0, " idx ", 6(i0, 1x), " significant patch: ", l1)') lgtID, neighborhood, idx, ref_check

                        ! revoke coarsening status. note we must use a temp status or otherwise all blocks
                        ! will revoke their coarsening
                        ! if (ref_check == 0) then
                        !     lgt_block( lgtID_neighbor, IDX_REFINE_STS ) = TMP_STATUS
                        ! endif
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
end subroutine
