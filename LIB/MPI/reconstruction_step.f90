!> \brief Modify the SC and WC of a wavelet decomposed blocks at fine/coarse interfaces
!> This routine assumes that the input is already wavelet decomposed in spaghetti form
subroutine coarse_extension_modify(params, hvy_data, hvy_tmp, tree_ID, CE_case, s_val, clear_wc_only)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)     !< heavy data array, WDed in Spaghetti form
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)      !< heavy work data array - block data.
    integer(kind=ik), intent(in)        :: tree_ID                     !< Tree to be 
    character(len=*)                    :: CE_case                     !< String representing which kind of CE we want to do, act on tree, level or ref

    integer(kind=ik), intent(in), optional  :: s_val                   !< if CE_modify acts on level or ref, it is restricted to this value
    logical, optional, intent(in)       :: clear_wc_only               !< for second CE modify we only want to delete WC (actually only in ghost patch)

    integer(kind=ik)                    :: iteration, k, i_n, lgt_ID, hvy_ID, tree_me, CE_case_id
    integer(kind=ik)                    :: nx,ny,nz,nc, level_me, ref_me, level_n, lgt_ID_n
    logical                             :: ClearWcOnly

    nx = size(hvy_data, 1)
    ny = size(hvy_data, 2)
    nz = size(hvy_data, 3)
    nc = size(hvy_data, 4)

    ClearWcOnly = .false.
    if (present(clear_wc_only)) ClearWcOnly = clear_wc_only

    select case(CE_case)
    case("tree")
        CE_case_id = 0
    case("level")
        CE_case_id = 1
    case("ref")
        CE_case_id = 2
    case default
        call abort(240809, "My language does not have so many cases, so I have no idea what you want from me.")
    end select
    

    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgt_ID, IDX_MESH_LVL )
        ref_me         = lgt_block( lgt_ID, IDX_REFINE_STS)
        tree_me        = lgt_block( lgt_ID, IDX_TREE_ID )

        ! for some calls we don't want to work on whole tree so skip some blocks
        if (CE_case_ID == 1) then
            if (level_me /= s_val) cycle
        elseif (CE_case_ID == 2) then
            if (ref_me /= s_val) cycle
        endif

        ! loop over all relevant neighbors
        do i_n = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvy_ID, i_n) /= -1 ) then
                ! neighbor light data id
                lgt_ID_n = hvy_neighbor( hvy_ID, i_n )
                level_n = lgt_block( lgt_ID_n, IDX_MESH_LVL )

                ! coarse extension for coarser neighbors
                ! select only blocks which are leaf-blocks, additionally check if there is a same-lvl neighbor for the same patch
                if (level_n < level_me .and. block_is_leaf(params, hvy_ID) &
                .and. .not. block_has_valid_neighbor(params, hvy_id, i_n, 0)) then
                    ! manipulation of coeffs
                    call coarseExtensionManipulateWC_block(params, hvy_data(:,:,:,:,hvy_ID), i_n)
                    if (.not. clearWcOnly) then
                        call coarseExtensionManipulateSC_block(params, hvy_data(:,:,:,:,hvy_ID), hvy_tmp(:,:,:,:,hvy_ID), i_n)
                    endif
                elseif (level_n > level_me) then
                    ! It is actually possible for a block to have both finer and coarser neighbors ranging from level J-1 to J+1
                    ! In order to reconstruct values on level J, we need the WC and SC which are on this level.
                    ! Finer neighbor has its decomposition ready for level J+1 but we cannot sync it for level J
                    ! So it is impossible to reconstruct values if the filters range into J+1 ghost nodes
                    
                    ! Disabled for CVS testing as it is there not the case
! #ifdef DEV
!                         ! Let's ensure this crashes just to show that filters should not range into finer level ghost nodes
!                         call coarseExtensionManipulateWC_block(params, hvy_data(:,:,:,:,hvyID), neighborhood, params%g, params%g, set_garbage=.true.)
! #endif
                endif
            endif
        enddo
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

    integer(kind=ik)                    :: iteration, k, i_n, lgt_ID, hvy_ID, Nreconl, Nreconr
    integer(kind=ik)                    :: nx,ny,nz,nc, level_me, level_n, lgt_ID_n, idx(2,3), ref_me
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

        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgt_ID, IDX_MESH_LVL )

        ! in some cases not all blocks were wavelet decomposed
        ! we want to retransform / copy them only if they have the passed refinement status, used for CE_update_tree
        if (present(REF_TMP_CHECK)) then
            ref_me = lgt_block( lgt_ID, IDX_REFINE_STS )
            if (ref_me /= REF_TMP_CHECK) cycle
        endif

        ! check if this block is to be modified
        do i_n = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvy_ID, i_n) /= -1 ) then
                ! neighbor light data id
                lgt_ID_n = hvy_neighbor( hvy_ID, i_n )
                level_n = lgt_block( lgt_ID_n, IDX_MESH_LVL )

                ! select only blocks which are leaf-blocks, additionally check if there is a same-lvl neighbor for the same patch
                if ((level_n < level_me) .and. block_is_leaf(params, hvy_ID) &
                    .and. .not. block_has_valid_neighbor(params, hvy_id, i_n, 0)) then
                    toBeManipulated = .true.
                    exit
                endif
            endif
        enddo

        if (toBeManipulated) then
        
            ! reconstruct from the manipulated coefficients
            tmp_reconst = hvy_data(:,:,:,1:nc,hvy_ID)
            call waveletReconstruction_block(params, tmp_reconst)

            ! if (lgt_block(lgtID, IDX_TC_2) == 4032 .and. lgt_block(lgtID, IDX_MESH_LVL)==5) then
            !     write(*, '("113 lvl ", i0)') level
            !     do iy = 1, 19
            !         write(*, '(19(es8.1))') tmp_reconst(1:19, iy, 1, 1)
            !     enddo
            ! endif

            hvy_data(:,:,:,1:nc,hvy_ID) = hvy_tmp(:,:,:,1:nc,hvy_ID)

            ! reconstruction part. We manipulated the data and reconstructed them on the entire block with modified coeffs.
            ! Now, we copy those reconstructed data back to the original block - this is
            ! the actual coarseExtension.
            ! JB: In theory we synched the SC and WC of the same level so we could overwrite the whole domain
            !     however, as with blocks with both finer and coarser neighbours the finer neighbours do not synch WC and SC
            !     the WR reconstructs false values at the border there
            !     For CVS this problem should solve itself as every block can find same-level neighbours for all directions
            do i_n = 1, size(hvy_neighbor,2)
                if ( hvy_neighbor(hvy_ID, i_n) /= -1 ) then
                    ! neighbor light data id
                    lgt_ID_n = hvy_neighbor( hvy_ID, i_n )
                    level_n = lgt_block( lgt_ID_n, IDX_MESH_LVL )

                    ! coarse extension for coarser neighbors
                    ! select only blocks which are leaf-blocks, additionally check if there is a same-lvl neighbor for the same patch
                    if (level_n < level_me .and. block_is_leaf(params, hvy_ID) &
                        .and. .not. block_has_valid_neighbor(params, hvy_id, i_n, 0)) then
                        ! coarse extension case (neighbor is coarser)
                        idx(:, :) = 1
                        call get_indices_of_modify_patch(params%g, params%dim, i_n, idx, (/ nx, ny, nz/), &
                        (/Nreconl, Nreconl, Nreconl/), (/Nreconr, Nreconr, Nreconr/), &
                        g_m=(/params%g, params%g, params%g/), g_p=(/params%g, params%g, params%g/))

                        hvy_data(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc, hvy_ID) = &
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
            hvy_data(:,:,:,1:nc,hvy_ID) = hvy_tmp(:,:,:,1:nc,hvy_ID)
        endif
    enddo
end subroutine

