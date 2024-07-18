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
            call sync_ghosts_tree( params, hvy_block, tree_ID )    

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
!> It is usually called after coarsening indicator where all blocks are WDed and have flag -1, 0 or REF_TMP_TREATED_COARSEN
subroutine addSecurityZone_CE_tree( time, params, tree_ID, hvy_block, hvy_tmp, indicator, norm, ignore_maxlevel, input_is_WD)

    use module_indicators

    implicit none
    real(kind=rk), intent(in)      :: time
    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)    !< heavy data array - block data might be wavelet decomposed form
    real(kind=rk), intent(inout)   :: hvy_tmp(:, :, :, :, :)      !< not used but we need to pass it to coarseningIdicator_block
    character(len=*), intent(in)   :: indicator                   !< how to choose blocks for refinement
    integer(kind=ik), intent(in)   :: tree_ID                     !< tree to look at
    real(kind=rk), intent(inout)   :: norm(:)                     !< the computed norm for each component of the vector
    logical, intent(in)            :: input_is_WD                 !< flag if hvy_block is already wavelet decomposed, for coarseningIndicatorBlock
    !> for the mask generation (time-independent mask) we require the mask on the highest
    !! level so the "force_maxlevel_dealiasing" option needs to be overwritten. Life is difficult, at times.
    logical, intent(in)                 :: ignore_maxlevel

    integer(kind=ik), parameter :: TMP_STATUS = 17
    integer(kind=ik) :: level_me, ref_status, neighborhood, i_neighborhood, ref_status_neighbor, ref_check, &
    hvyID, lgtID, k, level_neighbor, lgtID_neighbor, idx(2,3), nx, ny, nz, g, N_buffer_l, N_buffer_r
    logical :: refinement, inputIsWD

    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    g = params%g
    ! how large is the patch to be checked?
    ! this parameter could be set conservative (Nrec) or aggressive (Nwc+2)
    ! since for unlifted wavelets we also always coarsen aggressively over block boundaries, we do the same here and do not deploy any extra buffer
    N_buffer_l    = params%Nwcl+2
    N_buffer_r    = params%Nwcr+2
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

        level_me   = lgt_block( lgtID, IDX_MESH_LVL )
        ref_status = lgt_block( lgtID, IDX_REFINE_STS )

        ! ignore blocks on Jmax when they are forced to coarsen, it is quite unnecessary but better be safe then sorry
        if (level_me == params%JMax .and. params%force_maxlevel_dealiasing .and. .not. ignore_maxlevel) cycle

        ! is this block assigned 0 (it wants to stay and is significant)?
        if (ref_status == 0) then
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor             = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor             = lgt_block( lgtID_neighbor, IDX_MESH_LVL )
                    ref_status_neighbor = lgt_block( lgtID_neighbor, IDX_REFINE_STS )

                    ! if (lgtID == 32160) then
                    !     write(*, '("BH-", i0, " BL-", i0, " N-", i0, " NBL-", i0, " NR-", i0)') hvyID, lgtID, neighborhood, lgtID_neighbor, ref_status_neighbor
                    ! endif

                    ! the neighbor wants to coarsen
                    if ((ref_status_neighbor == -1 .or. ref_status_neighbor == REF_TMP_TREATED_COARSEN).and.(level_neighbor == level_me)) then
                        ! check for the patch where wc will be deleted, make sure to exclude ghost patches
                        call get_indices_of_modify_patch(params, neighborhood, idx, (/ nx, ny, nz/), (/N_buffer_l, N_buffer_l, N_buffer_l/), (/N_buffer_r, N_buffer_r, N_buffer_r/), &
                            X_s=(/ g, g, g/), X_e=(/ g, g, g/))
                        ref_check = -1

                        call coarseningIndicator_block( params, hvy_block(:,:,:,:,hvyID), &
                        hvy_tmp(:,:,:,:,hvyID), indicator, ref_check, norm, level_me, input_is_WD, indices=idx)

                        ! encode into refinement status if the neighboring block can coarsen or not (reasoning explained below)
                        if (ref_check == 0) then
                            ! write(*, '("SZ1 I-", i0, " R-", i0, " BL-", i0, " BH-", i0, " L-", i0, " Ref-", i0, " TC-", i0, " BLN-", i0, " RefN-", i0)') &
                            !     9, params%rank, lgtID, hvyid, level_me, ref_status, lgt_block(lgtid, IDX_TC_2), lgtID, ref_status_neighbor

                            call lgt_encode_significant_patch(lgtID, neighborhood, params%dim)
                        endif
                    endif
                endif
            enddo
        endif
    enddo

    ! We look at one block which is signfiant and want to change the refinement status of neighbours, however these neighbors might reside on a different proc
    ! The next synchronize_lgt_data call will then overwrite those values. We therefore need to communicate with the owners of the neighboring block
    ! if the status has changed or not. This is however a bit tricky as we don't know how many cells have to be updated and information send
    ! We circumvent this by encoding the patch which is significant into the refinement status and then decoding this information, knowing which blocks
    ! need to be updated.
    
    ! After encoding the significancy, we need to synch those between all blocks
    call synchronize_lgt_data( params,  refinement_status_only=.true. )


    ! loop over all of my blocks and check if they have a significant neighbor with significant patch in the neighborhood
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        level_me   = lgt_block( lgtID, IDX_MESH_LVL )
        ref_status = lgt_block( lgtID, IDX_REFINE_STS )

        ! is the block on the level we look at?
        ! is this block assigned -1 (it wants to coarsen)?
        if (ref_status == -1 .or. ref_status == REF_TMP_TREATED_COARSEN) then
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! make sure to invert neighborhood to access the correct patchIDs, this is used only for the singificant patch bit
                ! as before we looked b_significant->b_wants2coarsen and now we look b_wants2coarsen->b_significant
                i_neighborhood = inverse_neighbor(neighborhood, params%dim)

                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor      = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor      = lgt_block( lgtID_neighbor, IDX_MESH_LVL )
                    ref_status_neighbor = lgt_block( lgtID_neighbor, IDX_REFINE_STS )

                    ! does the neighbor has significancy flag and significant patch in this neighborhood relation?
                    ! Is neighbor on the same level so we can actually check the patches?
                    if (ref_status_neighbor > 0 .and. level_neighbor == level_me) then  ! need to check elsewise all significancy encodings are wrong
                        if (lgt_decode_significant_flag(lgtID_neighbor) .and. lgt_decode_significant_patch(lgtID_neighbor, i_neighborhood, params%dim)) then
                            ! write(*, '("SZ2 I-", i0, " R-", i0, " BL-", i0, " BH-", i0, " L-", i0, " Ref-", i0, " TC-", i0, " BLN-", i0, " RefN-", i0)') &
                            !     9, params%rank, lgtID, hvyid, level_me, ref_status, lgt_block(lgtid, IDX_TC_2), lgtID, ref_status_neighbor

                            ! revoke coarsening status. note we must use a temp status or otherwise all blocks
                            ! will revoke their coarsening
                            lgt_block( lgtID, IDX_REFINE_STS ) = TMP_STATUS
                        endif
                    endif
                endif
            enddo
        endif
    enddo

    ! After decoding the significancy, we need to synch changed stati between all blocks
    ! JB: Missing second synchronize was a bug I searched for long time
    call synchronize_lgt_data( params,  refinement_status_only=.true. )

    ! synchronous lgt_loop in order to set all encoded significancies and TMP_STATUS to 0
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        if (lgt_block( lgtID, IDX_REFINE_STS ) > 0) then  ! need to check elsewise all significancy encodings are wrong
            if (lgt_block( lgtID, IDX_REFINE_STS ) == TMP_STATUS .or. lgt_decode_significant_flag(lgtID)) then
                lgt_block( lgtID, IDX_REFINE_STS ) = 0
            endif
        endif
    enddo


end subroutine



!> \brief Encode into ref status that a block has a significant patch and wants a neighbor to stay
!> This encodes the significance of any of the 8 (2D, 4 sides 4 corners) or 26 (3D, 6 sides, 12 edges, 8 corners) patches into the refinement status
!> We assume that the refinement status can only contain -1 or 1 at that point, so we set one chosen bit to one to indicate that
!> binary refinement status of block that wants to coarsen (-1): 1111.1111.1111.1111:1111.1111.1111.1111
!> binary refinement status of block that wants to stay (0):     0000.0000.0000.0000:0000.0000.0000.0000
!> binary range r for any of the 26 patches with flag f, r2l:    0frr.rrrr.rrrr.rrrr.rrrr.rrrr.rrrr.0000
subroutine lgt_encode_significant_patch(lgtID, neighborhood, dim)
    implicit none

    integer(kind=ik), intent(in)   :: lgtID         !< ID of block where information will be encoded
    integer(kind=ik), intent(in)   :: neighborhood  !< neighborhood of significant patch
    integer(kind=ik), intent(in)   :: dim           !< params%dim

    integer(kind=ik)    :: patchID, ref_encoded

    call neighborhood2patchID(neighborhood, patchID, dim)

    ! get current ref status, should be zero except if already some value was encoded, then we will encode it on top of it
    ref_encoded = lgt_block(lgtID, IDX_REFINE_STS)

    ! safety feature, lgtID should never want to coarsen -1 (or smaller than 0 for any reason), then something went wrong
    if (ref_encoded < 0) then
        call abort(197002, "Block wants to coarsen but at the same time has a significant patch? That makes no sense!")
    endif
    
    ! In order to only set one bit, we do binary-or with a number were only the corresponding bit is set
    ! addition does not work as we might previously already have selected this patch to be significant and the addition would invalidate the results
    ! +3 as patchID is 1-based and we start at 4th bit from the right
    ref_encoded = IOR(ref_encoded, ISHFT(1, patchID+3))

    ! set flag on second-highest bit, this is used for looping, for example to wipe the encoding at the end back to 0 easily
    ref_encoded = IOR(ref_encoded, ISHFT(1, bit_size(ref_encoded)-2))

    lgt_block(lgtID, IDX_REFINE_STS) = ref_encoded

end subroutine



!> \brief Decode from refinement status if a block at a specific neighborhood has a significant patch
!> More explanation can be found at function lgt_encode_significant_patch
function lgt_decode_significant_patch(lgtID, neighborhood, dim)
    implicit none

    integer(kind=ik), intent(in)   :: lgtID                         !< ID of block where information will be encoded
    integer(kind=ik), intent(in)   :: neighborhood                  !< neighborhood of significant patch
    integer(kind=ik), intent(in)   :: dim                           !< params%dim
    logical                        :: lgt_decode_significant_patch  !< Flag if this neighborhood corresponds to an encoded significant patch

    integer(kind=ik)    :: patchID, ref_encoded, bit_extracted

    call neighborhood2patchID(neighborhood, patchID, dim)

    ! get current ref status, should be zero except if already some value was encoded, then we will encode it on top of it
    ref_encoded = lgt_block(lgtID, IDX_REFINE_STS)

    ! safety feature, lgtID should never want to coarsen -1 (or smaller than 0 for any reason), then something went wrong
    if (ref_encoded < 0) then
        call abort(197003, "Block wants to coarsen but at the same time has a significant patch? That makes no sense!")
    endif
    
    ! extract correct bit, +3 as patchID is 1-based and we start at 4th bit from the right
    bit_extracted = IBITS(ref_encoded, patchID+3, 1)
    lgt_decode_significant_patch = (bit_extracted .eq. 1)

end function



!> \brief Decode from refinement status if a block at a specific neighborhood has any significant patch, so it has the flag set
!> More explanation can be found at function lgt_encode_significant_patch
function lgt_decode_significant_flag(lgtID)
    implicit none

    integer(kind=ik), intent(in)   :: lgtID                         !< ID of block where information will be encoded
    logical                        :: lgt_decode_significant_flag  !< Flag if this neighborhood corresponds to an encoded significant patch

    integer(kind=ik)    :: patchID, ref_encoded, bit_extracted

    ! get current ref status, should be zero except if already some value was encoded, then we will encode it on top of it
    ref_encoded = lgt_block(lgtID, IDX_REFINE_STS)

    ! safety feature, lgtID should never want to coarsen -1 (or smaller than 0 for any reason), then something went wrong
    if (ref_encoded < 0) then
        call abort(197004, "Block wants to coarsen but at the same time has a significant patch? That makes no sense!")
    endif
    
    ! extract flag from second-highest bit
    bit_extracted = IBITS(ref_encoded, bit_size(ref_encoded)-2, 1)
    lgt_decode_significant_flag = (bit_extracted .eq. 1)

end function