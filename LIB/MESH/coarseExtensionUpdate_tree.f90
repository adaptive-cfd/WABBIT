
!> \brief do coarse extensions for all coarse-fine interfaces on a block
!> This is very similar to adapt_tree but only affects blocks on coarse-fine interfaces and does not need to loop
!  !!! ATTENTION !!! This routine is finished but has not been cross-checked yet
subroutine coarseExtensionUpdate_tree( params, hvy_block, hvy_tmp, tree_ID, input_is_synced)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    logical, intent(in), optional       :: input_is_synced   !< in case input data was already synced we can pass this here

    real(kind=rk)                       :: t0, t1
    integer(kind=ik)                    :: k_b, k_n, g_spaghetti, g_this, REF_TMP, lgtID_neighbor, level_neighbor, ref_neighbor, hvyID, lgtID, level_me, ref_me
    logical                             :: inputIsSynced

    t0          = MPI_Wtime()
    t1          = MPI_Wtime()
    REF_TMP = 197

    ! it turns out, when the coefficients are spaghetti-ordered,
    ! we can sync only even numbers of points and save one for odd numbered
    g_spaghetti = params%g/2*2

    inputIsSynced = .false.
    if (present(input_is_synced)) inputIsSynced = input_is_synced


    ! To avoid that the incoming hvy_neighbor array and active lists are outdated
    ! we synchronize them.
    t0 = MPI_Wtime()
    call updateMetadata_tree(params, tree_ID)
    call toc( "coarseExntesionsUpdate_tree (update metadata)", 151, MPI_Wtime()-t0 )

    ! mark all blocks where to apply coarse-extension
    t0 = MPI_Wtime()
    do k_b = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k_b, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgtID, IDX_MESH_LVL )

        ! check if this block is to be modified
        do k_n = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, k_n) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, k_n )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                ! we proceed level-wise
                if ((level_neighbor < level_me)) then
                    lgt_block(lgtID, IDX_REFINE_STS) = REF_TMP_UNTREATED
                    exit
                endif
            endif
        enddo
    end do

    ! we also need to wavelet decompose all neighbors of affected blocks so that we can sync the SC and WC
    ! in theory some are not needed (far off the interface) but this is not implemented currently
    call synchronize_lgt_data( params, refinement_status_only=.true. )
    do k_b = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k_b, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgtID, IDX_MESH_LVL )

        ! check if neighbor wants to coarse, set temporary flag
        do k_n = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            lgtID_neighbor = hvy_neighbor( hvyID, k_n )
            if ( lgtID_neighbor /= -1 ) then
                ! neighbor light data id
                ref_neighbor = lgt_block( lgtID_neighbor, IDX_REFINE_STS )

                ! we proceed level-wise
                if (ref_neighbor == REF_TMP_UNTREATED) then
                    lgt_block(lgtID, IDX_REFINE_STS) = REF_TMP
                    exit
                endif
            endif
        enddo
    end do
    ! change temporary flag to temporary untreated flag
    do k_b = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k_b, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        ref_me       = lgt_block( lgtID, IDX_REFINE_STS )

        if (ref_me == REF_TMP) lgt_block( lgtID, IDX_REFINE_STS ) = REF_TMP_UNTREATED
    end do
    ! one more sync so that all ref flags are everywhere, probably this one is not needed anymore but lets keep it clean
    call synchronize_lgt_data( params, refinement_status_only=.true. )

    call toc( "coarseExtensionUpdate_tree (temp flag)", 152, MPI_Wtime()-t0 )

    ! sync values to blocks that want to wavelet decompose
    if (.not. inputIsSynced) then
        t0 = MPI_Wtime()
        g_this = max(ubound(params%HD,1),ubound(params%GD,1))
        call sync_TMP_from_MF( params, hvy_block, tree_ID, REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
        call toc( "coarseExtensionUpdate_tree (sync lvl <- MF)", 153, MPI_Wtime()-t0 )
    endif

    ! Wavelet-transform all blocks on this level
    ! From now on until wavelet retransform hvy_block will hold the wavelet decomposed values in spaghetti form for affected blocks
    t0 = MPI_Wtime()
    do k_b = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k_b, tree_ID)

        ! We compute detail coefficients on the fly here, for all blocks
        ! on the level.
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        ref_me = lgt_block( lgtID, IDX_REFINE_STS )

        ! FWT required for a block that is on the level
        if (ref_me == REF_TMP_UNTREATED) then
            ! hvy_tmp now is a copy with sync'ed ghost points.
            hvy_tmp(:,:,:,1:size(hvy_block, 4),hvyID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvyID)
            level_me = lgt_block( lgtID, IDX_MESH_LVL )

            ! Compute wavelet decomposition
            call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
        endif
    end do
    call toc( "coarseExtensionUpdate_tree (FWT)", 154, MPI_Wtime()-t0 )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! coarseExtension: remove wavelet coefficients near a fine/coarse interface
    ! on the fine block. Does nothing in the case of CDF60, CDF40 or CDF20.
    ! This is the first coarse extension before removing blocks
    ! As every block is assumed to be WDed here we can do it on the whole tree
    t0 = MPI_Wtime()
    if (params%useCoarseExtension .and. params%isLiftedWavelet) then
        call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="tree")
    endif
    call toc( "coarseExtensionUpdate_tree (coarse_extension_modify)", 155, MPI_Wtime()-t0 )

    ! synch SC and WC from new coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
    ! Attention1: For finer neighbours this is not possible, so near fine interfaces we cannot reconstruct correct values.
    ! Attention2: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
    ! Attention3: Neighboring blocks have fulfilled their role already so what they sync is not important, we could skip it but do it anyways for now
    t0 = MPI_Wtime()
    call sync_SCWC_from_MC( params, hvy_block, tree_ID, hvy_tmp, g_minus=g_spaghetti, g_plus=g_spaghetti)
    call toc( "coarseExtensionUpdate_tree (sync all <- MC)", 156, MPI_Wtime()-t0 )

    ! After synching with coarser neighbors, the spaghetti form of the whole ghost patch is filled with values
    ! We now need to delete all values in spots of WC, leaving only the correct SC values from coarse neighbours
    ! This uses the fact that values refined on the coarse grid points with wc=0 are the copied SC values
    ! We do this in order tu bundle up the synchronizations in the step before as the modification is really cheap
    t0 = MPI_Wtime()
    if (params%useCoarseExtension .and. params%isLiftedWavelet) then
        call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="tree", copy_sc=.false.)
    endif
    call toc( "coarseExtensionUpdate_tree (coarse_extension_modify)", 155, MPI_Wtime()-t0 )

    ! Wavelet-reconstruct blocks
    ! Copy back old values if no coarse extension is applied and elsewise just overwrite the affected patches inside the domain
    t0 = MPI_Wtime()
    call coarse_extension_reconstruct_tree(params, hvy_block, hvy_tmp, tree_ID, REF_TMP_CHECK=REF_TMP_UNTREATED)
    call toc( "coarseExtensionUpdate_tree (RWT)", 157, MPI_Wtime()-t0 )

    ! synchronize ghost nodes - final synch to update all neighbours with the new values
    ! Just once here at the end after reconstruct so everything is in order
    t0 = MPI_Wtime()
    call sync_ghosts_tree( params, hvy_block, tree_ID )
    call toc( "coarseExtensionUpdate_tree (sync post)", 158, MPI_Wtime()-t0 )

    call toc( "coarseExtensionUpdate_tree (TOTAL)", 150, MPI_wtime()-t1)

end subroutine