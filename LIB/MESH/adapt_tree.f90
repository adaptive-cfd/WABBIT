!> \author JB
!> \brief This routine performs the adaption of the mesh, where possible, with the full wavelet transformation. \n
!! 1. Grid for full wavelet transformation is prepared and mother blocks are created
!! 2. All blocks are wavelet transformed from fine to coarse in leaf-wise fashion, CE is applied on leaf-layer only
!! 3. CVS or denoising iterative loops are applied to compute threshold
!! 4. Grid adaption is computed and Blocks will be deleted where possible dependend on indicator
!! 5. All blocks will be reconstructed from JMin to JMax level-wise, updating daughter SC along the way
!
!> \note It is well possible to start with a very fine mesh and end up with only one active
!! block after this routine. You do *NOT* have to call it several times.
subroutine adapt_tree( time, params, hvy_block, tree_ID, indicator, hvy_tmp, hvy_work, hvy_mask, ignore_coarsening, ignore_maxlevel, std_est)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params  !< good ol' params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy tmp data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout), optional        :: hvy_work(:, :, :, :, :, :)
    !> mask data. we can use different trees (forest module) to generate time-dependent/indenpedent
    !! mask functions separately. This makes the mask routines tree-level routines (and no longer
    !! block level) so the physics modules have to provide an interface to create the mask at a tree
    !! level. All parts of the mask shall be included: chi, boundary values, sponges.
    !! Optional: if the grid is not adapted to the mask, passing hvy_mask is not required.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    character(len=*), intent(in)        :: indicator
    !> Sometimes we don't want to coarsen
    logical, intent(in), optional       :: ignore_coarsening
    !> during mask generation it can be required to ignore the maxlevel coarsening....life can suck, at times.
    logical, intent(in), optional       :: ignore_maxlevel
    integer(kind=ik), intent(in)        :: tree_ID
    !> Sometimes we want to pass in and out the estimated std - this probably will be made more smart at some point
    real(kind=rk), optional, intent(inout) :: std_est(:)
    
    ! loop variables
    integer(kind=ik)                    :: iteration, k, lgt_id
    real(kind=rk)                       :: t_block, t_all, t_loop
    integer(kind=ik)                    :: Jmax_active, Jmin_active, level, ierr, k1, hvy_id
    logical                             :: ignore_coarsening_apply, ignore_maxlevel_apply, iterate, toBeManipulated
    integer(kind=ik)                    :: level_me, ref_stat, Jmin, lgt_n_old, g_this
    character(len=clong)                :: toc_statement
    integer(kind=tsize)                 :: treecode

    real(kind=rk) :: thresh(1:params%n_eqn), thresh_old(1:params%n_eqn), thresh_save(1:params%n_eqn+1)

    real(kind=rk) :: norm(1:params%n_eqn)

    ! testing
    integer(kind=ik)                    :: ix, iy
    real(kind=rk)                       :: n_points


    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas
   
    t_block     = MPI_Wtime()
    t_all       = MPI_Wtime()
    iteration   = 0
    iterate     = .true.
    Jmin        = params%Jmin
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)
    level       = Jmax_active ! algorithm starts on maximum *active* level
    norm(:)     = 1.0_rk


    ignore_coarsening_apply = .false.
    ignore_maxlevel_apply = .false.
    if (present(ignore_coarsening)) ignore_coarsening_apply = ignore_coarsening
    if (present(ignore_maxlevel)) ignore_maxlevel_apply = ignore_maxlevel


    ! To avoid that the incoming hvy_neighbor array and active lists are outdated
    ! we synchronize them.
    t_block = MPI_Wtime()
    call updateMetadata_tree(params, tree_ID, search_overlapping=.false.)
    call toc( "adapt_tree (updateMetadata_tree)", 101, MPI_Wtime()-t_block )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !      Wavelet decomposition
    ! This rountine first create the full tree grid and then starts from leaf-layer and continues downwards with wavelet decomposition
    ! Blocks always pass their SC to their mother (pendant to excuteCoarsening of old function), new blocks only synch from medium neighbors to avoid CE
    ! This is repeated until all blocks on the layer JMin are present with wavelet decomposed values
    t_block = MPI_Wtime()
    call wavelet_decompose_full_tree(params, hvy_block, tree_ID, hvy_tmp)
    call toc( "adapt_tree (decompose_tree)", 102, MPI_Wtime()-t_block )
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !      Iterative loop to estimate thresholding
    ! For CVS and image denoising we iteratively estimate the thresholding later need for deletion of WC and actual coarsening
    if (indicator == "threshold-cvs" .or. indicator == "threshold-image-denoise") then

        ! call cvs_thresholding_local(params, hvy_block, hvy_tmp, hvy_work, indicator, tree_ID)

        call cvs_thresholding_global(params, hvy_block, hvy_tmp, hvy_work, indicator, tree_ID, thresh(:))

        ! debug: print thresholding factors to file
        ! JB: I don't know how to do it cleverly for now so I save it in a t-file
        if (params%rank==0) then
            thresh_save(1) = time
            thresh_save(2:) = thresh(:)
            call append_t_file( "thresholding.t", thresh_save(:) )
            ! write(*, '(A, 10(es10.3, 1x))') "Thresholding factors: ", thresh(:)
        endif
    endif
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ! Testing for de-aliasing - compute energy in JMax blocks before and after decomposition to compute what is dissipated
    ! call componentWiseNorm_tree(params, hvy_tmp, tree_ID, "L2", norm, norm_case="level", n_val=params%Jmax)
    ! if (params%rank == 0) write(*, '(A, 10(1x, es12.4))') "Norm on JMax  before de-aliasing:", norm(:)
    ! call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm, norm_case="SC_level", n_val=params%Jmax)
    ! if (params%rank == 0) write(*, '(A, 10(1x, es12.4))') "Norm on JMax-1 after de-aliasing:", norm(:)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !      Coarsening / Grid adaption
    ! Grid is adapted according to the indicator
    if (.not. ignore_coarsening_apply) then
        t_loop = MPI_Wtime()
        if (indicator == "threshold-cvs" .or. indicator == "threshold-image-denoise") then
            ! write back wavelet decomposed and normal values that have been backuped
            do k = 1, hvy_n(tree_ID)
                hvy_ID = hvy_active(k, tree_ID)
                  hvy_tmp(:, :, : ,1:size(hvy_block, 4), hvy_id) = hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 1)  ! normal values
                hvy_block(:, :, : ,1:size(hvy_block, 4), hvy_id) = hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 2)  ! WDed values
            enddo

            if (params%threshold_wc) then
                ! delete all WC where they are too small - this is the actual thresholding
                call wc_threshold_tree(params, hvy_block, tree_ID, thresh, Jmax_active)
            endif

            ! set norm to be used as threshold indicator as threshold
            norm(:) = thresh(:)
        elseif (params%threshold_wc) then
            ! compute norm
            if (params%eps_normalized) then
                call componentWiseNorm_tree(params, hvy_tmp, tree_ID, params%eps_norm, norm)
            else
                norm(:) = 1.0_rk
            endif

            ! delete all WC where they are too small - this is the actual thresholding
            call wc_threshold_tree(params, hvy_block, tree_ID, params%eps * norm(:), Jmax_active)
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Coarsening indicator: flag blocks for coarsening
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Check the entire grid where to coarsen. Note this is a wrapper for coarseningIndicator_block, which
        ! acts on a single block only.We distinguish two cases: wavelet or no wavelet (Shakespeare!). The wavelet case is the 
        ! default, other indicators are used mostly for testing. In the wavelet case, we check the blocks for their largest 
        ! detail (=wavelet coeff). The routine assigns -1 (COARSEN) to a block, if it matches the criterion. 
        ! This status may however be revoked below.
        t_block = MPI_Wtime()
        ! Mask will only be passed on if it is present
        call coarseningIndicator_tree( time, params, hvy_block, hvy_tmp, tree_ID, indicator, &
            ignore_maxlevel=ignore_maxlevel_apply, input_is_WD=.true., hvy_mask=hvy_mask, norm_inout=norm)
        call toc( "adapt_tree (coarseningIndicator_tree)", 108, MPI_Wtime()-t_block )

        ! After coarseningIndicator_tree, the situation is:
        ! coarseningIndicator_tree works on all blocks (for wavelet cases).
        ! blocks that are significant now have status 0, others have -1

        if (params%useSecurityZone .and. indicator/="everywhere" .and. indicator/="random" .and. indicator/='nothing (external)') then
            ! if we want to add a security zone, we check for every significant block if a neighbor wants to coarsen
            ! if this is the case, we check if any significant WC would be deleted (basically checking the thresholding for this patch)
            ! in that case we set the neighbouring block to be important as well (with a temporary flag)
            t_block = MPI_Wtime()
            call addSecurityZone_CE_tree( time, params, tree_ID, hvy_block, hvy_tmp, indicator, norm, ignore_maxlevel=ignore_maxlevel_apply, input_is_WD=.true.)
            call toc( "adapt_tree (security_zone_check)", 109, MPI_Wtime()-t_block)
        endif

        ! check if block has reached minimal level, if so, remove refinement flags
        call respectJmaxJmin_tree( params, tree_ID )

        ! unmark blocks that cannot be coarsened due to gradedness and completeness, in one go this should converge to the final grid
        call ensureGradedness_tree( params, tree_ID, check_daughters=.true., stay_value=REF_UNSIGNIFICANT_STAY)

        ! we do not need to move any cell data anymore, any block with -1 can simply be deleted
        do k = 1, lgt_n(tree_ID)
            lgt_ID = lgt_active(k, tree_ID)
            if ( lgt_block(lgt_ID, IDX_REFINE_STS) == -1) then
                ! delete blocks
                lgt_block(lgt_ID, :) = -1
            endif
        enddo

        ! update grid lists: active list, neighbor relations, etc
        t_block = MPI_Wtime()
        call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)
        call toc( "adapt_tree (updateMetadata_tree)", 101, MPI_Wtime()-t_block )

        ! This is the actual important coarse extension after all cells have been deleted, which modifies now the lasting c-f interfaces
        ! As all decompositions have been correct, we do not need to copy the SC and only delete the WC in order to keep the interfaces clean
        if (params%useCoarseExtension) then
            call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="tree", copy_sc=.false.)

            ! After CE we need to update bound cond for symmetry
            call bound_cond_generic(params, hvy_block, tree_ID, "tree", .true., edges_only=.false.)
        endif

        call toc( "adapt_tree (coarsening)", 104, MPI_Wtime()-t_loop )
    endif
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !      Wavelet reconstruction of blocks
    ! Reconstruct all blocks, this is done level-wise and scaling coefficients of daughter blocks are updated along the way
    ! If no SC or WC are modified (usually for unlifted wavelets and normal thresholding) we can simply copy back all values
    t_block = MPI_Wtime()
    if (.not. params%useCoarseExtension .and. indicator /= "threshold-cvs" .and. indicator /= "threshold-image-denoise") then
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            hvy_block(:, :, : ,1:size(hvy_block, 4), hvy_id) = hvy_tmp(:, :, : ,1:size(hvy_block, 4), hvy_id)  ! restore original values
        enddo
    else
        ! In a special case reconstruction can be skipped as it has already been thone, therefore the if-condition
        if (.not. ((indicator == "threshold-cvs" .or. indicator == "threshold-image-denoise") .and. ignore_coarsening_apply)) then
            call wavelet_reconstruct_full_tree(params, hvy_block, hvy_tmp, tree_ID)
        endif
    endif
    call toc( "adapt_tree (reconstruct_tree)", 105, MPI_Wtime()-t_block )
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! delete all non-leaf blocks with daughters as we for now do not have any use for them
    call prune_fulltree2leafs(params, tree_ID)

    ! At this point the coarsening is done. All blocks that can be coarsened are coarsened
    ! they may have passed several level also and non-leafs have been deleted. Rebalance for optimal loadbalancing
    ! ToDo: we sync directly afterwards so we could only transfer inner points, but this needs the buffering from xfer_bock_data
    t_block = MPI_Wtime()
    call balanceLoad_tree( params, hvy_block, tree_ID)
    call toc( "adapt_tree (balanceLoad_tree)", 106, MPI_Wtime()-t_block )

    ! synchronize ghost nodes - final synch to update all neighbours with the new values
    ! Just once here at the end after reconstruct so everything is in order
    t_block = MPI_Wtime()
    call sync_ghosts_tree( params, hvy_block, tree_ID )
    call toc( "adapt_tree (sync post)", 107, MPI_Wtime()-t_block )

    call toc( "adapt_tree (TOTAL)", 100, MPI_wtime()-t_all)
end subroutine



!> \brief function for full wavelet decomposition algorithm
!> \details This function implements the full favelet decomposition.
!! First, the full tree grid is prepared, then the leaf-layer is decomposed in one go
!! in order to have load-balanced work there, then the mothers are lvl-wise updated and consequently
!! decomposed until all blocks have the decomposed values. Afterwards, the grid will be in full tree format.
subroutine wavelet_decompose_full_tree(params, hvy_block, tree_ID, hvy_tmp, init_full_tree_grid, compute_SC_only, scaling_filter, verbose_check)
    implicit none

    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    !> heavy tmp data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    logical, intent(in), optional       :: init_full_tree_grid  !< if false, grid is already in full tree format, defaults to true
    logical, intent(in), optional       :: compute_SC_only  !< only apply scaling function filter, defaults to false
    real(kind=rk), optional :: scaling_filter(:)  !< filter to be used for the SC, only used with compute_SC_only=T, defaults to params%HD
    logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

    integer(kind=ik)     :: k, lgt_ID, hvy_ID, lgt_n_old, g_this, level_me, ref_stat, iteration, level
    real(kind=rk)        :: t_block, t_loop
    character(len=clong) :: toc_statement
    logical              :: iterate, use_leaf_first, initFullTreeGrid, computeSCOnly
    real(kind=rk), dimension(:), allocatable :: scalingFilter
    integer(kind=ik), dimension(:), allocatable, save :: lgt_refinementStatus_backup

    initFullTreeGrid = .true.
    if (present(init_full_tree_grid)) initFullTreeGrid = init_full_tree_grid
    computeSCOnly = .false.
    if (present(compute_SC_only)) computeSCOnly = compute_SC_only
    if (present(scaling_filter)) then
        ! passed filters might not have correct dimensions, so we set it here
        k = (ubound(scaling_filter, dim=1)-lbound(scaling_filter, dim=1))/2
        allocate(scalingFilter(-k:k))
        scalingFilter(:) = scaling_filter(:)
    else
        allocate(scalingFilter(lbound(params%HD, dim=1):ubound(params%HD, dim=1)))
        scalingFilter(:) = params%HD(:)
    endif

    ! the refinement_status is used in this routine for things other than the refinement_status
    ! but that may be problematic if we have set the refinement_status externally and like to keep it.
    ! This happens for example in postprocessing (--grid1-to-grid2). In this case, the indicator is
    ! "nothing (external)" and we should not overwrite it. Hence, we make a backup, and copy the original
    ! status (when entering this routine) back to lgt_block
    if (.not. allocated(lgt_refinementStatus_backup)) then
        allocate(lgt_refinementStatus_backup(1:size(lgt_block,1)))
    endif

    do k = 1, lgt_n(tree_ID)
        lgt_refinementStatus_backup(lgt_active(k, tree_ID)) = lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS)
    enddo


    ! we can only use the more parallel optimized leaf-first option with a specific minimum blocksize, so we test that here
    use_leaf_first = all(params%Bs(1:params%dim) >= 3* max(abs(lbound(params%HD, dim=1)), abs(ubound(params%HD, dim=1))))

    ! In order to start investigation, we set correct ref flag, being -1
    !    leaf-first : all blocks
    !    level-wise : only blocks on highest level
    level       = maxActiveLevel_tree(tree_ID)
    if (initFullTreeGrid) then
        do k = 1, lgt_n(tree_ID)
            lgt_ID = lgt_active(k, tree_ID)
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
            if (level_me == level .or. use_leaf_first) then
                lgt_block(lgt_ID, IDX_REFINE_STS) = -1
            endif
        end do
        ! do backup here so we need less logic in synching neighbors
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
        enddo

        ! at first, initialize all mothers without any values yet
        t_block = MPI_Wtime()
        call init_full_tree(params, tree_ID)
        call toc( "decompose_tree (init_mothers_tree)", 110, MPI_Wtime()-t_block )
        ! now, all mothers have been created and share the ref status REF_TMP_EMPTY, leafs share status -1
    else
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
            if (level_me == level .or. (use_leaf_first .and. block_is_leaf(params, hvy_ID))) then
                lgt_block(lgt_ID, IDX_REFINE_STS) = -1
            else
                lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_EMPTY
            endif
            ! do backup here so we need less logic in synching neighbors
            if (block_is_leaf(params, hvy_ID)) then
                hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
            endif
        end do
        ! we need to synch lgt data as loop to check for leaf-blocks is 
        call synchronize_lgt_data( params, refinement_status_only=.true.)
    endif

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !      leaf-decompose level-update loop
    ! This might sound confusing, but goes the following:
    !     first step decomposes all leaf-blocks, this will be the most expensive action per definition
    !     every update goes level-wise and starts from finest level, new blocks will therefore be ready to decompose leaf-wise
    ! With this, we keep the heavy computations leaf-wise, so load balanced on the first step, afterwise the level-wise loops needs no refinement flags and lgt_data update
    level       = maxActiveLevel_tree(tree_ID)
    iterate     = .true.
    iteration   = 0
    do while (iterate)
        t_loop = MPI_Wtime()
        lgt_n_old = lgt_n(tree_ID)

        ! ! I often need to check block values when debugging so I leave it here
        ! do k = 1, hvy_n(tree_ID)
        !     hvy_ID = hvy_active(k, tree_ID)
        !     call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        !     ! if (params%rank == 0 .and. lgt_block(lgt_id, IDX_REFINE_STS) == -1 .and. iteration >= 2) then
        !     ! if (lgt_block(lgt_id, IDX_MESH_LVL) == 3 .and. lgt_block(lgt_id, IDX_TC_2) == 54525952) then
        !         write(*, '("1 I-", i0, " R-", i2, " B-", i6, " L-", i2, " Ref-", i3, " TC-", i0)') iteration, params%rank, lgt_ID, &
        !             lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2)
        !         ! call write_neighborhood_info(lgt_id, hvy_neighbor(hvy_id, :), params%dim)
        !         ! call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), "block_cvs_9.txt", params%Bs, params%g)
        !     ! endif
        ! enddo

        ! synchronize ghost nodes - required to apply wavelet filters
        ! block only needs information from medium and fine neighbors as CE will cut dependency to coarse neighbors
        t_block = MPI_Wtime()
        ! first point is SC and last WC, so we need left(HD) or right(GD) for sync
        g_this = max(abs(lbound(params%HD,1)),ubound(params%GD,1))

        ! ! for coarse extension we are not dependend on coarser neighbors so lets skip the syncing
        ! if (params%isLiftedWavelet) then
        !     call sync_TMP_from_MF( params, hvy_block, tree_ID, -1, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
        !     call toc( "decompose_tree (sync TMP <- MF)", 111, MPI_Wtime()-t_block )

        !     write(toc_statement, '(A, i0, A)') "decompose_tree (it ", iteration, " sync lvl <- MF)"
        !     call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        ! ! unlifted wavelets need coarser neighbor values for their WC so we need to sync them too
        ! else
            call sync_TMP_from_all( params, hvy_block, tree_ID, -1, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
            call toc( "decompose_tree (sync TMP <- all)", 111, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "decompose_tree (it ", iteration, " sync lvl <- all)"
            call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        ! endif

        ! Wavelet-transform all blocks which are untreated
        ! From now on until wavelet retransform hvy_block will hold the wavelet decomposed values in spaghetti form for these blocks
        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            ref_stat = lgt_block( lgt_ID, IDX_REFINE_STS )
    
            ! FWT required for a block that is on the level
            if (ref_stat == -1) then
                ! hvy_tmp now is a copy with sync'ed ghost points.
                ! We actually copy here a second time but this is to ensure we have the correct ghost points saved, and copying is not too expensive anyways
                hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
                level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

                ! Compute wavelet decomposition
                ! For Jmax if dealiasing, just compute H filter
                ! Data SC/WC now in Spaghetti order
                if (computeSCOnly) then
                    ! if we only want to compute the SC, we need to apply the filter to the data
                    call blockFilterXYZ_vct( params, hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID), hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID), scalingFilter, &
                        lbound(scalingFilter, dim=1), ubound(scalingFilter, dim=1), do_restriction=.true.)
                elseif (level_me == params%Jmax .and. params%force_maxlevel_dealiasing) then
                    call blockFilterXYZ_vct( params, hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID), hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID), params%HD, &
                        lbound(params%HD, dim=1), ubound(params%HD, dim=1), do_restriction=.true.)
                else
                    call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvy_ID))
                endif
            endif
        end do
        call toc( "decompose_tree (FWT)", 112, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "decompose_tree (it ", iteration, " FWT)"
        call toc( toc_statement, 1150+iteration, MPI_Wtime()-t_block )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! coarseExtension: remove wavelet coefficients near a fine/coarse interface
        ! on the fine block. Does nothing in the case of CDF60, CDF40 or CDF20.
        ! This is the first coarse extension before removing blocks
        ! As no blocks are deleted, modifying newly decomposed blocks is sufficient, they are non-leafs anyways and will be skipped
        if (params%useCoarseExtension) then
            call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="ref", s_val=-1)

            ! After CE we need to update bound cond for symmetry
            call bound_cond_generic(params, hvy_block, tree_ID, "ref", .false., s_val=-1, edges_only=.false.)
        endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Upate refinement flags before synching:
        !    1. Blocks that have beend decomposed with CE are finished -> They get 0
        !    2. Mothers on level-1 with status REF_TMP_EMPTY get -1 now to enable sync to them
        !       they will be decomposed next loop
        do k = 1, lgt_n(tree_ID)
            lgt_ID = lgt_active(k, tree_ID)
            ! update that all blocks with -1 have been decomposed
            if (lgt_block(lgt_ID, IDX_REFINE_STS) == -1) then
                lgt_block(lgt_ID, IDX_REFINE_STS) = 0
            endif
            ! update all updated mother blocks which are empty (leaf-blocks will never be updated for leaf-first optimization) and set their refinement flag to -1
            if ((lgt_block(lgt_ID, IDX_REFINE_STS) == REF_TMP_EMPTY .or. .not. use_leaf_first) .and. lgt_block( lgt_ID, IDX_MESH_LVL ) == level-1) then
                lgt_block(lgt_ID, IDX_REFINE_STS) = -1
            endif
        end do

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !     Update mother blocks
        ! We update mother blocks level-wise, with this we always know that all blocks on this level are ready and need no lgt_data update
        ! This uses the already wavelet decomposed blocks and copies their SC into the new mother block
        t_block = MPI_Wtime()
        call sync_D2M(params, hvy_block, tree_ID, sync_case="level", s_val=level)
        call toc( "decompose_tree (sync level 2 mothers)", 113, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "decompose_tree (it ", iteration, " sync level 2 mothers)"
        call toc( toc_statement, 1200+iteration, MPI_Wtime()-t_block )

        ! iteration counter
        iteration = iteration + 1
        level = level - 1
        ! loop continues until we are on the lowest level.
        iterate = (level >= params%Jmin)

        write(toc_statement, '(A, i0, A)') "decompose_tree (it ", iteration, " TOTAL)"
        call toc( toc_statement, 1250+iteration, MPI_Wtime()-t_loop )

        if (present(verbose_check)) then
            write(*, '(A, i0, A, i0, A)') "Loop ", iteration, " with ", lgt_n(tree_ID), " blocks"
        endif

        ! ! for debugging purposes - savin does a sync if any ghost values are NaN
        ! write( toc_statement,'(A, i6.6, A)') 'TestWD_', iteration, "000000.h5"
        ! call saveHDF5_tree(toc_statement, dble(iteration), iteration, 1, params, hvy_block, tree_ID)
        ! write( toc_statement,'(A, i6.6, A)') 'TestN_', iteration, "000000.h5"
        ! call saveHDF5_tree(toc_statement, dble(iteration), iteration, 1, params, hvy_tmp, tree_ID)
    end do

    ! copy the original refinemtn_status (when entering this routine) back to lgt_block
    do k = 1, lgt_n(tree_ID)
        lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS) = lgt_refinementStatus_backup(lgt_active(k, tree_ID))
    enddo

    ! ! for debugging purposes - savin does a sync if any ghost values are NaN
    ! write( toc_statement,'(A, i6.6, A)') 'TestWD_', 100, "000000.h5"
    ! call saveHDF5_tree(toc_statement, dble(100), 100, 1, params, hvy_block, tree_ID)
    ! write( toc_statement,'(A, i6.6, A)') 'TestN_', 100, "000000.h5"
    ! call saveHDF5_tree(toc_statement, dble(100), 100, 1, params, hvy_tmp, tree_ID)
end subroutine



!> \brief function for full wavelet reconstruction algorithm
!> \details This function implements the full favelet reconstruction.
!! The input has to be in full grid format in order for the mother-update cycle to work.
!! Blocks are reconstructed level-wise starting and JMin and mothers are updated until all blocks are reconstructed.
!! Afterwards the grid is still in full tree format, but can be pruned by deleting all non-leafs
subroutine wavelet_reconstruct_full_tree(params, hvy_block, hvy_tmp, tree_ID)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy tmp data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    ! loop variables
    integer(kind=ik)                    :: iteration, k, Jmax_active, Jmin_active, level, hvy_ID, lgt_ID, level_me, Jmin
    real(kind=rk)                       :: t_block, t_all, t_loop
    character(len=clong)                :: toc_statement


    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    Jmin        = params%Jmin
    Jmax_active = maxActiveLevel_tree(tree_ID)
    iteration = 0
    do level = Jmin, Jmax_active
        t_loop = MPI_Wtime()

        ! ! Normal adapt_tree loop does another CE_modify here to copy the SC, but for full WT this is not necessary
        ! ! as interior blocks only get values from their medium neighbors and no SC need to be copied
        ! call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="level", s_val=level)

        ! synch SC and WC from coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
        ! Attention: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
        t_block = MPI_Wtime()
        call sync_SCWC_from_MC( params, hvy_block, tree_ID, hvy_tmp, g_minus=params%g, g_plus=params%g, level=level)
        call toc( "reconstruct_tree (sync lvl <- MC)", 115, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "reconstruct_tree (it ", iteration, " sync level <- MC)"
        call toc( toc_statement, 1300+iteration, MPI_Wtime()-t_block )

        ! In theory, for full WT this CE_modify is only here for two reasons:
        !    1. We need to clear the WC from coarser neighbors, as the values in there are interpolated SC and we need to set those to the correct WC being 0
        !    2. interior WC are deleted so that the SC copy for syncing later is exact with the original values, we assume those WC are not significant
        if (params%useCoarseExtension) then
            call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="level", s_val=level, copy_sc=.false.)

            ! After CE we need to update bound cond for symmetry
            call bound_cond_generic(params, hvy_block, tree_ID, "level", .true., s_val=level, edges_only=.false.)
        endif

        ! Wavelet-reconstruct blocks on level
        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

            ! RWT required for a block that is on the level
            if (level_me == level) then
                ! Compute wavelet reconstruction
                call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvy_ID))

                ! make copy of data in hvy_tmp as we sync from coarser neighbors using that and I do not want to write another logic currently
                ! One might be concerned if we sync from hvy_tmp and do prediction/interpolation, as it does not have it's ghost layers synced
                ! however, we discard all interpolated values with the CE so we cleverly circumvent problems there
                hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
            endif
        end do
        call toc( "reconstruct_tree (RWT)", 116, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "reconstruct_tree (it ", iteration, " RWT)"
        call toc( toc_statement, 1350+iteration, MPI_Wtime()-t_block )

        ! now we need to update the SC of all daughters
        t_block = MPI_Wtime()
        call sync_M2D(params, hvy_block, tree_ID, sync_case="level", s_val=level)
        call toc( "reconstruct_tree (sync level 2 daughters)", 117, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "reconstruct_tree (it ", iteration, " sync level 2 daughters)"
        call toc( toc_statement, 1400+iteration, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "reconstruct_tree (it ", iteration, " TOTAL)"
        call toc( toc_statement, 1450+iteration, MPI_Wtime()-t_loop )

        iteration = iteration + 1
    enddo

    ! ! for debugging purposes - savin does a sync if any ghost values are NaN
    ! write( toc_statement,'(A, i6.6, A)') 'TestWD_', 100, "000000.h5"
    ! call saveHDF5_tree(toc_statement, dble(100), 100, 1, params, hvy_block, tree_ID)
    ! write( toc_statement,'(A, i6.6, A)') 'TestN_', 100, "000000.h5"
    ! call saveHDF5_tree(toc_statement, dble(100), 100, 1, params, hvy_tmp, tree_ID)
end subroutine



!> \brief Function to init all mother blocks from leaf-only grid down to JMin
!> \details This functions inits all mother blocks from a leaf-only grid down to JMin.
!! This will be used for wavelet_decompose_full_tree in order to prepare the grid, which it then fills with values.
!! The benefit is, that the grid topology needs to be updated only once instead of once every loop cycle.
subroutine init_full_tree(params, tree_ID, set_ref, verbose_check)
    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik), intent(in), optional :: set_ref     !< if not set, ref will be set to REF_TMP_EMPTY, else to this number
    logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

    integer(kind=ik)     :: i_level, j, k, N, lgt_ID, hvy_ID, level_b, iteration, level, data_rank, lgt_merge_id, hvy_merge_id, hvy_n_old, digit, refset
    real(kind=rk)        :: t_block, t_loop
    character(len=clong) :: toc_statement
    logical              :: iterate
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: lgt_daughters(1:8), rank_daughters(1:8)
    integer(kind=tsize)                 :: treecode

    refset = REF_TMP_EMPTY
    if (present(set_ref)) refset = set_ref

    ! number of blocks to merge, 4 or 8
    N = 2**params%dim
    ! loop downwards until all blocks are on Jmin
    level       = maxActiveLevel_tree(tree_ID) ! leaf-wise blocks can have this as maximum level
    iteration   = 0

    ! J: I have decided for the following way of creating mother blocks
    !   - main loop goes refine-evolve-coarsen, so first iteration will work on all available block
    !     it is the one with the most data to send and it's nice that we can choose any of the ranks here (as they are all the same)
    !   - after first iteration, not all sister-information are available, I just create the following mothers after a deterministic way
    !     and choose the rank of the block with the highest digit, this ensures some form of spread of ranks, after all I do not really care anymore about the amount of data

    ! just loop until hvy_n does not change anymore, that means all possible mothers have been created
    do while( hvy_n_old /= hvy_n(tree_ID) )
        !---------------------------------------------------------------------------
        ! create new empty blocks on rank with block with digit = 2**dim-1
        !---------------------------------------------------------------------------
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
            level_b = lgt_block( lgt_ID, IDX_MESH_LVL )

            ! check if block does not have a mother yet and blocks are not on minimum level
            if ( hvy_family(hvy_ID, 1) == -1 .and. level_b /= params%Jmin) then
                ! Get digit
                treecode = get_tc(lgt_block( lgt_ID, IDX_TC_1:IDX_TC_2 ))
                digit = tc_get_digit_at_level_b(treecode, params%dim, level_b, params%Jmax)
                ! The merging will be done on the rank of the block with digit 2**dim-1
                if (digit == 2**params%dim -1) then
                    ! construct new mother block if on my rank and create light data entry for the new block
                    call get_free_local_light_id(params, params%rank, lgt_merge_id, message="init_mothers")
                    call lgt2hvy(hvy_merge_id, lgt_merge_id, params%rank, params%number_blocks)
                    ! set meta_data of mother block
                    lgt_block( lgt_merge_id, : ) = -1
                    call set_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2), tc_clear_until_level_b(treecode, &
                        dim=params%dim, level=level_b-1, max_level=params%Jmax))
                    lgt_block( lgt_merge_id, IDX_MESH_LVL ) = level_b-1
                    lgt_block( lgt_merge_id, IDX_REFINE_STS ) = refset
                    lgt_block( lgt_merge_id, IDX_TREE_ID ) = tree_ID
                    ! make sure family array for mother block is cleared
                    hvy_family(hvy_merge_id, :) = -1
                    
                    ! update myself that I have created my own mother
                    ! sisters will not be updated but as only one daughter will surely create the mother there will be no doubles
                    hvy_family(hvy_ID, 1) = lgt_merge_id
                endif
            endif
        enddo

        ! we need to update hvy_active locally to be able to loop correctly
        hvy_n_old = hvy_n(tree_ID)
        call createHvyActive_tree(params, tree_ID)

        iteration = iteration + 1

        if (present(verbose_check)) then
            write(*, '(A, i0, A, i0, A, i0)') "R", params%rank, " Deterministic mother addition it ", iteration, ", hvy_n= ", hvy_n(tree_ID)
        endif
    enddo

    ! the active lists are outdated, so lets resynch
    call synchronize_lgt_data( params, refinement_status_only=.false.)
    ! At last, we update all full neighbor relations only once
    call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)
end subroutine



!> \brief Threshold wavelet components
!> \details This function hard thresholds all wavelet components smaller than a specific value.
!! It is used for CVS and image denoising in order to remove noise from the field.
subroutine wc_threshold_tree(params, hvy_block, tree_ID, thresh, JMax_active)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)  :: params
    !> heavy data array, WDed values in spaghetti form
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)    :: tree_ID
    real(kind=rk), intent(in)       :: thresh(:)   !< Threshold for CVS to be applied to all WC
    integer(kind=ik), intent(in)    :: JMax_active !< if not provided it will be computed

    integer(kind=ik) :: hvy_ID, lgt_ID, ix, iy, iz, even_odd, i_b, i_p, level_me, l, g_dim(1:params%dim)
    real(kind=rk)    :: level_fac, wc_fac
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: u_wc
    integer, dimension(:), allocatable  :: mask_i, mask


    ! for threshold state vector components we do some slicing magic
    ! fortran does not allow logical mask slicing, so we build a mask of indices to slice with
    allocate(mask_i(1:params%n_eqn))
    do l = 1, params%n_eqn
        mask_i(l) = l
    end do

    ! is the first point a SC or WC?
    even_odd = mod(params%g + 1, 2)
    g_dim(:) = 0
    g_dim(1:params%dim) = params%g

    if (allocated(u_wc)) then
        if (size(u_wc, 4) < params%n_eqn) deallocate(u_wc)
    endif
    if (.not. allocated(u_wc)) allocate(u_wc(1:size(hvy_block, 1), 1:size(hvy_block, 2), 1:size(hvy_block, 3), 1:params%n_eqn ) )

    do i_b = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(i_b, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

        call wavelet_renorm_block(params, hvy_block(:, :, :, 1:params%n_eqn, hvy_ID), u_wc, level_me, level_ref=JMax_active)

        ! for threshlding components which are on their own
        do i_p = 1, params%n_eqn
            if (params%threshold_state_vector_component(i_p) == 1) cycle

            ! make sure to check if only the interior should be thresholded (syncing needed afterwards) or the whole block
            ! do iz = 1, params%Bs(3)+2*g_dim(3)
            !     do iy = 1, params%Bs(2)+2*g_dim(2)
            !         do ix = 1, params%Bs(1)+2*g_dim(1)
            do iz = g_dim(3)+1, params%Bs(3)+g_dim(3); do iy = g_dim(2)+1, params%Bs(2)+g_dim(2); do ix = g_dim(1)+1, params%Bs(1)+g_dim(1)
                ! skip SC, for 3D skip the check in z-direction
                if (mod(ix, 2) == even_odd .and. mod(iy, 2) == even_odd .and. mod(iz, 2) == merge(1, even_odd, params%dim==2)) cycle

                ! apply thresholding for thresholding components which are on their own
                if (abs(u_wc(ix, iy, iz, i_p)) < thresh(i_p)) then
                    hvy_block(ix, iy, iz, i_p, hvy_ID) = 0.0
                endif
            enddo; enddo; enddo
        enddo

        ! some thresholding components want to be treated together, as for example the velocity, so we compute a norm-equivalent if wanted
        ! every value >1 in thresholding_component are treated as belonging together
        ! so we could set 2 2 3 3 and actually have those treated as independent vector fields, where we compute the norm
        ! maybe a bit overkill because usually we only have the velocity, but why not have the capacity that we can build up on
        do l = 2, maxval(params%threshold_state_vector_component(:))
            ! convert logical mask to index mask
            mask = pack(mask_i, params%threshold_state_vector_component(:)==l)

            ! make sure to check if only the interior should be thresholded (syncing needed afterwards) or the whole block
            ! do iz = 1, params%Bs(3)+2*g_dim(3)
            !     do iy = 1, params%Bs(2)+2*g_dim(2)
            !         do ix = 1, params%Bs(1)+2*g_dim(1)
            do iz = g_dim(3)+1, params%Bs(3)+g_dim(3); do iy = g_dim(2)+1, params%Bs(2)+g_dim(2); do ix = g_dim(1)+1, params%Bs(1)+g_dim(1)
                ! skip SC, for 3D skip the check in z-direction
                if (mod(ix, 2) == even_odd .and. mod(iy, 2) == even_odd .and. mod(iz, 2) == merge(1, even_odd, params%dim==2)) cycle

                ! apply thresholding for thresholding components which are on their own
                ! all threshs should be similar, so maxval just takes one of them
                if (sqrt(sum(u_wc(ix, iy, iz, mask)**2)) < maxval(thresh(mask))) then
                    hvy_block(ix, iy, iz, mask, hvy_ID) = 0.0
                endif
            enddo; enddo; enddo
        enddo
    enddo

end subroutine


!> \brief Threshold wavelet components
!> \details This function hard thresholds all wavelet components smaller than a specific value.
!! It is used for CVS and image denoising in order to remove noise from the field.
subroutine wc_threshold_block(params, hvy_block, lgt_id, thresh)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)  :: params
    !> heavy data array, WDed values in spaghetti form
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :)
    integer(kind=ik), intent(in)    :: lgt_ID   !< lgt_id of this block
    real(kind=rk), intent(in)       :: thresh   !< Threshold for CVS to be applied to all WC

    integer(kind=ik) :: ix, iy, iz, even_odd, nc, level_me
    real(kind=rk)    :: level_fac, wc_fac
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: u_wc

    nc     = params%n_eqn
    level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

    ! is the first point a SC or WC?
    even_odd = mod(params%g + 1, 2)

    if (allocated(u_wc)) then
        if (size(u_wc, 4) < nc) deallocate(u_wc)
    endif
    if (.not. allocated(u_wc)) allocate(u_wc(1:size(hvy_block, 1), 1:size(hvy_block, 2), 1:size(hvy_block, 3), 1:nc ) )
    ! u_wc = 0.0_rk  ! set to zero, as maybe the first layer is not used and normed but will be thresholded nonetheless

    call wavelet_renorm_block(params, hvy_block(:, :, :, 1:nc), u_wc, level_me, level_ref=params%JMax, indices=(/1,params%Bs(1)+2*params%g,1,params%Bs(2)+2*params%g,1,merge(1,params%Bs(3)+2*params%g,params%dim==2)/))

    ! first point can always be ignored for spaghetti ordered values - if it is a WC, then it will not be used so we don't have to threshold it
    ! this seems unnecessary, but wavelet_renorm_block cannot affect the first point by setting indices and setting u_wc=0 is unnecessary resources spent
    do iz = 2, merge(1, params%Bs(3)+2*params%g, params%dim==2)
        do iy = 2, params%Bs(2)+2*params%g
            do ix = 2, params%Bs(1)+2*params%g
                ! skip SC, for 3D skip the check in z-direction
                if (mod(ix, 2) == even_odd .and. mod(iy, 2) == even_odd .and. mod(iz, 2) == merge(1, even_odd, params%dim==2)) cycle

                ! apply thresholding for norm of velocity only
                if (sqrt(sum(u_wc(ix, iy, iz, 1:params%dim)**2)) < thresh) then
                    hvy_block(ix, iy, iz, :) = 0.0
                endif
            enddo
        enddo
    enddo

end subroutine


!> \brief Prune tree by deleting all blocks which are not leafs
!> \details This function takes a whole tree and prunes it so that only leaf-blocks are left over.
subroutine prune_fulltree2leafs(params, tree_ID)
    implicit none

    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(in)    :: tree_ID

    integer(kind=ik) :: hvy_ID, lgt_ID, k
    real(kind=rk)    :: t_block

    ! delete all non-leaf blocks with daughters as we for now do not have any use for them
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )

        if ( .not. block_is_leaf(params, hvy_id) ) then
            ! we can savely delete blocks as long as we do not update family relations
            lgt_block( lgt_ID, : ) = -1
        endif
    end do

    ! deletion was only made locally (as family relations are hvy), so we need to sync lgt data
    call synchronize_lgt_data( params, refinement_status_only=.false. )

    ! update grid lists: active list, neighbor relations, etc
    call updateMetadata_tree(params, tree_ID, search_overlapping=.false.)

end subroutine


subroutine cvs_thresholding_local(params, hvy_block, hvy_tmp, hvy_work, indicator, tree_ID)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy tmp data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)
    character(len=*), intent(in)        :: indicator
    integer(kind=ik), intent(in)        :: tree_ID

    ! loop variables
    integer(kind=ik)                    :: iteration, k, k1, Jmax_active, Jmin_active, level, hvy_ID, lgt_ID, level_me, Jmin, g_spaghetti
    real(kind=rk)                       :: t_block, t_all, t_loop
    character(len=clong)                :: toc_statement

    ! testing
    integer(kind=ik)                    :: ix, iy
    real(kind=rk)                       :: n_points

    ! array for the norms and thresholds
    real(kind=rk), allocatable, save :: norm(:), thresh(:), thresh_old(:)

    if (.not. allocated(norm)) allocate(norm(params%number_blocks))
    if (.not. allocated(thresh)) allocate(thresh(params%number_blocks))
    if (.not. allocated(thresh_old)) allocate(thresh_old(params%number_blocks))

    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    ! it turns out, when the coefficients are spaghetti-ordered,
    ! we can sync only even numbers of points and save one for odd numbered
    g_spaghetti = params%g/2*2

    ! synch SC and WC from coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
    ! Attention: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
    t_block = MPI_Wtime()
    call sync_SCWC_from_MC( params, hvy_block, tree_ID, hvy_tmp, g_minus=g_spaghetti, g_plus=g_spaghetti)
    call toc( "reconstruct_tree (sync lvl <- MC)", 115, MPI_Wtime()-t_block )

    ! In theory, for full WT this CE_modify is only here for two reasons:
    !    1. We need to clear the WC from coarser neighbors, as the values in there are interpolated SC and we need to set those to the correct WC being 0
    !    2. interior WC are deleted so that the SC copy for syncing later is exact with the original values, we assume those WC are not significant
    if (params%useCoarseExtension) then
        call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="tree", copy_sc=.false.)
    endif

    ! now we are doing CVS loop, start by backing up important data
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)

        hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 1) =   hvy_tmp(:, :, : ,1:size(hvy_block, 4), hvy_id)  ! normal values
        hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 2) = hvy_block(:, :, : ,1:size(hvy_block, 4), hvy_id)  ! WDed values

        ! Compute wavelet reconstruction only from WC and save this one
        call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvy_ID), SC_reconstruct=.false., WC_reconstruct=.true.)
        hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 3) = hvy_block(:, :, : ,1:size(hvy_block, 4), hvy_id) ! WC reconstructed
    enddo

    ! ! debug what the wavelet reconstruction does
    ! call saveHDF5_tree("test-WCR.h5", 0.0_rk, 0, 1, params, hvy_block, tree_ID )
    ! call saveHDF5_tree("test-normal.h5", 0.0_rk, 0, 1, params, hvy_work(:, :, :, 1:size(hvy_block, 4), :, 1), tree_ID, no_sync=.true.)
    ! call saveHDF5_tree("test-WD.h5", 0.0_rk, 0, 1, params, hvy_work(:, :, :, 1:size(hvy_block, 4), :, 2), tree_ID, no_sync=.true.)
    ! ! return

    n_points = dble(product(params%bs(1:params%dim)))  ! block-wise, so all points for one block
    thresh(:) = 0.0_rk

    ! Wavelet-reconstruct only WC of all blocks
    do iteration = 0, min(100, params%azzalini_iterations)
        thresh_old(:) = thresh(:)

        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

            ! if (thresh(k) > 1e-2) then
            !     write(*, '(A, i0, 4(A, es12.3))') "Block ", k," thresh - ",thresh(k)," norm - ",norm(k)," : sum 1 - ", sum(abs(hvy_block(params%g+1:params%Bs(1)+params%g, params%g+1:params%Bs(2)+params%g, params%g+1:params%Bs(3)+params%g, 1:size(hvy_block, 4), hvy_id))), &
            !         " sum 2 - ", sum(abs(hvy_work(params%g+1:params%Bs(1)+params%g, params%g+1:params%Bs(2)+params%g, params%g+1:params%Bs(3)+params%g, 1:size(hvy_block, 4), hvy_id, 3)))
            ! endif

            ! subtract to get incoherent field
            if (iteration > 0) then
                hvy_tmp(:, :, :, 1:size(hvy_block, 4), hvy_id) = hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 3) - hvy_block(:, :, :, 1:size(hvy_block, 4), hvy_id)
            else
                hvy_tmp(:, :, :, 1:size(hvy_block, 4), hvy_id) = hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 3)
            endif

            call componentWiseNorm_block(params, hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID), lgt_ID, 'Energy_Prototype', norm(k))

        !     ! write the float value into 
        !     lgt_block(lgt_id, IDX_REFINE_STS) = transfer(real(norm(k), kind=kind(1.0)), 1_ik)

        ! enddo

        ! ! Sync the norm which is saved (for now) in the refinement status
        ! call synchronize_lgt_data( params, refinement_status_only=.true. )

        ! do k = 1, hvy_n(tree_ID)
        !     hvy_ID = hvy_active(k, tree_ID)
        !     call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
        !     level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

            ! compute thresholding according to Eqn (6) in Kadoch et al "On the role of vortical structures for turbulent mixing using direct
            ! numerical simulation and wavelet-based coherent vorticity extraction"
            ! J. Turb.(12)20: 1-17 (2011)
            ! thresh = sqrt( 4/3 Z/N log(N)) where N is the resolution - scales with level?
            ! thresh(k) = 0.0_rk
            ! do k1 = 1, 2**params%dim
            !     thresh(k) = thresh(k) + real(transfer(lgt_block(hvy_family(hvy_ID, 1+k1), IDX_REFINE_STS), 1.0), kind=rk)
            ! enddo
            ! thresh(k) = sqrt( (4.0_rk/3.0_rk) * thresh(k)/2.0_rk**params%dim * log(n_points*2**params%dim) )  ! block with all family
            thresh(k) = sqrt( (4.0_rk/3.0_rk) * norm(k) * log(n_points) )  ! norm already normalized by area/volume

            ! ! compute thresholding according to donoho1994 universal threshold T=\sqrt{2log(n)*sigma^2}
            ! thresh(k) = sqrt( 2.0_rk * norm(k) * log(n_points))  ! norm already normalized by area/volume

            ! recopy WDed values from backup
            hvy_block(:, :, : ,1:size(hvy_block, 4), hvy_ID) = hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_ID, 2) ! WDed values

            if (iteration < min(100, params%azzalini_iterations)) then
                ! delete all WC where they are too small
                call wc_threshold_block(params, hvy_block(:,:,:,:,hvy_ID), lgt_ID, thresh(k))

                ! Compute wavelet reconstruction only from WC
                call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvy_ID), SC_reconstruct=.false., WC_reconstruct=.true.)
            endif

            ! ! make copy of data in hvy_tmp as we sync from coarser neighbors using that and I do not want to write another logic currently
            ! hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
            
        end do

        if (params%rank == 0) then
            write(*, '(A, i2, 2(A, es12.4), A, 50(es12.4))') "CVS convergence it ", iteration, " thresh_mean= ", sum(thresh)/dble(hvy_n(tree_ID)), " thresh_max= ", maxval(thresh), " thresh = ", thresh(1:hvy_n(tree_ID))
        endif

        call toc( "CVS_loop", 10116, MPI_Wtime()-t_block )

        ! JB ToDo - set norm or threshk in refinement indicator? but there are only integer values, maybe set as float32?
        !           this way we can sync those betwee the blocks so that they have access to all norms, for prototyping that may be simplest
    enddo

    ! ! HACK: overwrite the thresholding if it is larger than wanted, possibly deleting too many details of the flow
    ! call componentWiseNorm_tree(params, hvy_work(:,:,:,:,:,1), tree_ID, indicator, norm)
    ! do k = 1, hvy_n(tree_ID)
    !     if (thresh(k) > params%eps * norm(k)) thresh(k) = params%eps * norm(k)
    ! enddo
    ! if (params%rank == 0) then
    !     write(*, '(A, i2, 2(A, es12.4))') "CVS convergence it ", -1, " thresh_mean= ", sum(thresh)/dble(hvy_n(tree_ID)), " thresh_max= ", maxval(thresh)
    ! endif

    ! now do the actual thresholding and reconstruction
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )

        ! write the float value into 
        lgt_block(lgt_id, IDX_REFINE_STS) = transfer(real(thresh(k), kind=kind(1.0)), 1_ik)

        call wc_threshold_block(params, hvy_block(:,:,:,:,hvy_ID), lgt_ID, thresh(k))
    enddo

    call wavelet_reconstruct_full_tree(params, hvy_block, hvy_tmp, tree_ID)

    ! ! debug: print thresholding factors to file
    ! ! JB: I don't know how to do it cleverly for now so I save it in a t-file
    ! if (params%rank==0) then
    !     thresh_save(1) = time
    !     thresh_save(2:) = thresh(:)
    !     call append_t_file( "thresholding.t", thresh_save(:) )
    !     ! write(*, '(A, 10(es10.3, 1x))') "Thresholding factors: ", thresh(:)
    ! endif
end subroutine


subroutine cvs_thresholding_global(params, hvy_block, hvy_tmp, hvy_work, indicator, tree_ID, thresh)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy tmp data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)
    character(len=*), intent(in)        :: indicator
    integer(kind=ik), intent(in)        :: tree_ID
    real(kind=rk), intent(inout)        :: thresh(1:params%n_eqn)

    ! loop variables
    integer(kind=ik)                    :: iteration, k, k1, Jmax_active, Jmin_active, level, hvy_ID, lgt_ID, level_me, Jmin, g_spaghetti
    real(kind=rk)                       :: t_block, t_all, t_loop
    character(len=clong)                :: toc_statement

    real(kind=rk) :: thresh_old(1:params%n_eqn), thresh_save(1:params%n_eqn+1), norm(1:params%n_eqn)

    ! testing
    integer(kind=ik)                    :: ix, iy
    real(kind=rk)                       :: n_points

    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    t_loop = MPI_Wtime()
    ! now we are doing CVS loop, start by backing up important data
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 1) =   hvy_tmp(:, :, : ,1:size(hvy_block, 4), hvy_id)  ! normal values
        hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 2) = hvy_block(:, :, : ,1:size(hvy_block, 4), hvy_id)  ! WDed values
    enddo

    ! compute first threshold, treat whole field as incoherent
    ! image-denoise: compute Q1^2 / 2 (done on all variables but only first one used?)
    ! CVS: assume we do ACM and compute enstrophy from computing vorticity on the fly from velocity components
    call componentWiseNorm_tree(params, hvy_tmp, tree_ID, indicator, norm)
    
    ! n_points = 2.0_rk**dble(params%dim*params%Jmax)*dble(product(params%bs(1:params%dim)))
    ! n_points = 2.0_rk**dble(params%dim*Jmax_active)*dble(product(params%bs(1:params%dim)))
    n_points = dble(product(params%bs(1:params%dim)))
    ! write(*, '(A, 2(i0, 1x))') "Points: ", params%Jmax, Jmax_active
    if (indicator == "threshold-cvs") then
        ! compute thresholding according to Eqn (6) in Kadoch et al "On the role of vortical structures for turbulent mixing using direct
        ! numerical simulation and wavelet-based coherent vorticity extraction"
        ! J. Turb.(12)20: 1-17 (2011)
        ! thresh = sqrt( 4/3 Z/N log(N)) where N is the resolution - scales with level?
        do k1 = 1, params%n_eqn
            thresh(k1) = sqrt( (4.0_rk/3.0_rk) * norm(k1) * log(n_points) )  ! norm already normalized by area/volume
        enddo
    elseif (indicator == "threshold-image-denoise") then
        ! compute thresholding according to donoho1994 universal threshold T=\sqrt{2log(n)*sigma^2}
        do k1 = 1, params%n_eqn
            thresh(k1) = sqrt( 2.0_rk * norm(k1) * log(n_points))  ! norm already normalized by area/volume
        enddo
    else
        ! no iterative loops for others, we probably compute the filter somewhere anyways
        thresh(:) = -1.0_rk
    endif
    thresh_old(:) = -1.0_rk

    if (params%rank == 0) then
        write(*, '(2(A, es12.4))') "CVS before loop       thresh= ", thresh(1), " est std= ", sqrt(norm(1))
    endif

    ! convergence loop
    do iteration = 1, max(100, params%azzalini_iterations)
        thresh_old(:) = thresh(:)
        ! recopy WDed values from backup
        if (iteration > 1) then
            do k = 1, hvy_n(tree_ID)
                hvy_ID = hvy_active(k, tree_ID)
                hvy_block(:, :, : ,1:size(hvy_block, 4), hvy_id) = hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 2) ! WDed values
            enddo
        endif

        ! delete all WC where they are too small
        call wc_threshold_tree(params, hvy_block, tree_ID, thresh, Jmax_active)

        ! now do the CVS loops
        call wavelet_reconstruct_full_tree(params, hvy_block, hvy_tmp, tree_ID)

        ! subtract to get incoherent field
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            hvy_tmp(:, :, :, 1:size(hvy_block, 4), hvy_id) = hvy_work(:, :, :, 1:size(hvy_block, 4), hvy_id, 1) - hvy_block(:, :, :, 1:size(hvy_block, 4), hvy_id)
        enddo

        ! compute new threshold, treat whole field as incoherent
        ! for now, we use first variable and compute Q1^2 / 2, luckily this is just the norm squared
        call componentWiseNorm_tree(params, hvy_tmp, tree_ID, indicator, norm)

        if (indicator == "threshold-cvs") then
            ! compute thresholding according to Eqn (6) in Kadoch et al "On the role of vortical structures for turbulent mixing using direct
            ! numerical simulation and wavelet-based coherent vorticity extraction"
            ! J. Turb.(12)20: 1-17 (2011)
            ! thresh = sqrt( 4/3 Z/N log(N)) where N is the resolution - scales with level?
            do k1 = 1, params%n_eqn
                thresh(k1) = sqrt( (4.0_rk/3.0_rk) * norm(k1) * log(n_points) )  ! norm already normalized by area/volume
            enddo
        elseif (indicator == "threshold-image-denoise") then
            ! compute thresholding according to donoho1994 universal threshold T=\sqrt{2log(n)*sigma^2}
            do k1 = 1, params%n_eqn
                thresh(k1) = sqrt( 2.0_rk * norm(k1) * log(n_points))  ! norm already normalized by area/volume
            enddo
        endif

        if (params%rank == 0) then
            write(*, '(A, i2, 3(A, es12.4))') "CVS convergence it ", iteration, " thresh= ", thresh(1), " est std= ", sqrt(norm(1)), " thresh_diff= ", abs(thresh(1) - thresh_old(1))
        endif

        ! exit early
        if (indicator == "threshold-image-denoise" .and. maxval(abs(thresh(:) - thresh_old(:))) < 1e-12) exit
        if (indicator == "threshold-cvs" .and. (iteration == params%azzalini_iterations .or. maxval(abs(thresh(:) - thresh_old(:))) == 0.0_rk)) exit  ! CVS computes one iteration and then uses new thresholding, totalling to 2 iterations
    enddo

    ! HACK: overwrite the thresholding if it is larger than wanted, possibly deleting too many details of the flow
    call componentWiseNorm_tree(params, hvy_work(:,:,:,:,:,1), tree_ID, indicator, norm)
    do k1 = 1, params%n_eqn
        if (thresh(k1) > params%eps * norm(k1)) thresh(k1) = params%eps * norm(k1)
    enddo
    if (params%rank == 0) then
        write(*, '(A, i2, 1(A, es12.4))') "CVS convergence it ", -1, " thresh= ", thresh(1)
    endif

    call toc( "adapt_tree (threshold iterative loop)", 103, MPI_Wtime()-t_loop )

end subroutine