!> \brief Modify the SC and WC of a wavelet decomposed blocks at fine/coarse interfaces
!> This routine assumes that the input is already wavelet decomposed in spaghetti form
!> It will work on all blocks or which have the REF_TMP_UNTREATED flag in refinement status for leaf-wise operation
subroutine coarse_extension_modify_tree(params, hvy_data, hvy_tmp, tree_ID, sc_skip_ghosts)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)     !< heavy data array, WDed in Spaghetti form
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)      !< heavy work data array - block data.
    integer(kind=ik), intent(in)        :: tree_ID                     !< Tree to be investigated
    logical, optional, intent(in)       :: sc_skip_ghosts              !< for second CE modify we can skip ghosts points

    integer(kind=ik)                    :: iteration, k, neighborhood, lgtID, hvyID, tree_me
    integer(kind=ik)                    :: nx,ny,nz,nc, level_me, level_neighbor, lgtID_neighbor
    logical                             :: toBeManipulated, scSkipGhosts

    nx = size(hvy_data, 1)
    ny = size(hvy_data, 2)
    nz = size(hvy_data, 3)
    nc = size(hvy_data, 4)

    scSkipGhosts = .false.
    if (present(sc_skip_ghosts)) scSkipGhosts = sc_skip_ghosts

    do k = 1, hvy_n(tree_ID)
        ! Set toBeManipulated - This is one of those sneaky errors I searched 10hours for - JB
        toBeManipulated = .false.

        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgtID, IDX_MESH_LVL )
        tree_me        = lgt_block( lgtID, IDX_TREE_ID )

        ! check if this block is to be modified
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                ! we proceed level-wise
                if ((level_neighbor < level_me) .and. (tree_me==tree_ID)) then
                    toBeManipulated = .true.
                    ! nnn = nnn + 1
                    ! its enough if one neighborhood is true
                    exit
                endif
            endif
        enddo

        if (toBeManipulated) then
            ! loop over all relevant neighbors
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )


                    if (level_neighbor < level_me) then
                        ! manipulation of coeffs
                        ! call coarseExtensionManipulateWC_block(params, wc, neighborhood)
                        ! call coarseExtensionManipulateSC_block(params, wc(:, :, :, :, 1), hvy_tmp(:,:,:,:,hvyID), neighborhood, scSkipGhosts)
                        call coarseExtensionManipulateWC_block(params, hvy_data(:,:,:,:,hvyID), neighborhood)
                        call coarseExtensionManipulateSC_block(params, hvy_data(:,:,:,:,hvyID), hvy_tmp(:,:,:,:,hvyID), neighborhood, scSkipGhosts)
                    elseif (level_neighbor > level_me) then
                        ! It is actually possible for a block to have both finer and coarser neighbors ranging from level J-1 to J+1
                        ! In order to reconstruct values on level J, we need the WC and SC which are on this level.
                        ! Finer neighbor has its decomposition ready for level J+1 but we cannot sync it for level J
                        ! So it is impossible to reconstruct values if the filters range into J+1 ghost nodes
                        
                        ! Let's ensure this crashes just to show that filters should not range into finer level ghost nodes
                        ! call coarseExtensionManipulateWC_block(params, wc, neighborhood, params%g, params%g, set_garbage=.true.)
                        call coarseExtensionManipulateWC_block(params, hvy_data(:,:,:,:,hvyID), neighborhood, params%g, params%g, set_garbage=.true.)
                    endif
                endif
            enddo
        endif
    enddo
end subroutine



!> \brief Apply CE for all blocks in a tree. This copies back the old values and then overwrites those on affected patches
subroutine coarse_extension_reconstruct_tree(params, hvy_data, hvy_tmp, tree_ID, REF_TMP_CHECK)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)     !< heavy data array
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)      !< heavy work data array - block data.
    integer(kind=ik), intent(in)        :: tree_ID                     !< which tree to study
    integer(kind=ik), intent(in), optional :: REF_TMP_CHECK               !< if this is passed, only blocks with this ref status will be affected at all

    integer(kind=ik)                    :: iteration, k, neighborhood, lgtID, hvyID, Nreconl, Nreconr
    integer(kind=ik)                    :: nx,ny,nz,nc, level_me, level_neighbor, lgtID_neighbor, idx(2,3), ref_me
    logical                             :: toBeManipulated
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: tmp_reconst

    integer(kind=ik) :: iy

    nx = size(hvy_data, 1)
    ny = size(hvy_data, 2)
    nz = size(hvy_data, 3)
    nc = size(hvy_data, 4)

    Nreconl = params%Nreconl
    Nreconr = params%Nreconr 

    if (allocated(tmp_reconst)) then
        if (size(tmp_reconst, 4) < nc) deallocate(tmp_reconst)
    endif
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )

    do k = 1, hvy_n(tree_ID)
        ! Set toBeManipulated - This is one of those sneaky errors I searched 10hours for - JB
        toBeManipulated = .false.

        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgtID, IDX_MESH_LVL )

        ! in some cases not all blocks were wavelet decomposed
        ! we want to retransform / copy them only if they have the passed refinement status, used for CE_update_tree
        if (present(REF_TMP_CHECK)) then
            ref_me = lgt_block( lgtID, IDX_REFINE_STS )
            if (ref_me /= REF_TMP_CHECK) cycle
        endif

        ! check if this block is to be modified
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                ! check if this patch is a CE patch
                if ((level_neighbor < level_me)) then
                    toBeManipulated = .true.
                    exit
                endif
            endif
        enddo

        if (toBeManipulated) then
        
            ! reconstruct from the manipulated coefficients
            tmp_reconst = hvy_data(:,:,:,1:nc,hvyID)
            call waveletReconstruction_block(params, tmp_reconst)

            ! if (lgt_block(lgtID, IDX_TC_2) == 4032 .and. lgt_block(lgtID, IDX_MESH_LVL)==5) then
            !     write(*, '("113 lvl ", i0)') level
            !     do iy = 1, 19
            !         write(*, '(19(es8.1))') tmp_reconst(1:19, iy, 1, 1)
            !     enddo
            ! endif

            hvy_data(:,:,:,1:nc,hvyID) = hvy_tmp(:,:,:,1:nc,hvyID)

            ! reconstruction part. We manipulated the data and reconstructed them on the entire block with modified coeffs.
            ! Now, we copy those reconstructed data back to the original block - this is
            ! the actual coarseExtension.
            ! JB: In theory we synched the SC and WC of the same level so we could overwrite the whole domain
            !     however, as with blocks with both finer and coarser neighbours the finer neighbours do not synch WC and SC
            !     the WR reconstructs false values at the border there
            !     For CVS this problem should solve itself as every block can find same-level neighbours for all directions
            do neighborhood = 1, size(hvy_neighbor,2)
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                    if (level_neighbor < level_me) then
                        ! coarse extension case (neighbor is coarser)
                        idx(:, :) = 1
                        call get_indices_of_modify_patch(params, neighborhood, idx, (/ nx, ny, nz/), (/Nreconl, Nreconl, Nreconl/), (/Nreconr, Nreconr, Nreconr/))

                        hvy_data(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc, hvyID) = &
                            tmp_reconst(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc)
                    endif
                endif
            enddo

            ! if (lgt_block(lgtID, IDX_TC_2) == 4032 .and. lgt_block(lgtID, IDX_MESH_LVL)==5) then
            !     write(*, '("113 data lvl ", i0)') level
            !     do iy = 1, 19
            !         write(*, '(19(es8.1))') hvy_data(1:19, iy, 1, 1, hvyID)
            !     enddo
            ! endif

        else  ! block is not modified, rewrite old values
            hvy_data(:,:,:,1:nc,hvyID) = hvy_tmp(:,:,:,1:nc,hvyID)
        endif
    enddo
end subroutine

