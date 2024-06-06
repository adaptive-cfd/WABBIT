!! \author  engels, JB
!
!> \brief This routine performs the coarsing of the mesh, where possible. For the given mesh
!! we compute the details-coefficients on all blocks. If four sister blocks have maximum
!! details below the specified tolerance, (so they are insignificant), they are merged to
!! one coarser block one level below. This process is repeated until the grid does not change
!! anymore.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!
!> \note It is well possible to start with a very fine mesh and end up with only one active
!! block after this routine. You do *NOT* have to call it several times.
subroutine adapt_tree( time, params, hvy_block, tree_ID, indicator, hvy_tmp, hvy_mask, ignore_maxlevel)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    ! mask data. we can use different trees (forest module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! Optional: if the grid is not adapted to the mask, passing hvy_mask is not required.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    character(len=*), intent(in)        :: indicator
    ! during mask generation it can be required to ignore the maxlevel coarsening....life can suck, at times.
    logical, intent(in), optional       :: ignore_maxlevel

    integer(kind=ik), intent(in)        :: tree_ID
    ! loop variables
    integer(kind=ik)                    :: iteration, k, lgt_id
    real(kind=rk)                       :: t0, t1
    integer(kind=ik)                    :: Jmax_active, Jmin_active, level, ierr, k1, hvy_id
    logical                             :: ignore_maxlevel2, iterate, toBeManipulated
    integer(kind=ik)                    :: level_me, ref_stat, Jmin, lgt_n_old, g_this, g_spaghetti
    real(kind=rk), allocatable, dimension(:,:,:,:,:), save :: tmp_wd  ! used for data juggling


    ! testing
    integer(kind=ik)                    :: ix, iy


    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    t0          = MPI_Wtime()
    t1          = MPI_Wtime()
    iteration   = 0
    iterate     = .true.
    Jmin        = params%Jmin
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)
    level       = Jmax_active ! algorithm starts on maximum *active* level

    if (Jmin<1) call abort(2202243, "Currently, setting Jmin<1 is not possible")

    ignore_maxlevel2 = .false.
    if (present(ignore_maxlevel)) ignore_maxlevel2 = ignore_maxlevel

    if (allocated(tmp_wd)) then
        if (size(tmp_wd, 4) < size(hvy_block, 4)) deallocate(tmp_wd)
    endif
    if (.not. allocated(tmp_wd)) allocate(tmp_wd(1:size(hvy_block, 1), 1:size(hvy_block, 2), 1:size(hvy_block, 3), 1:size(hvy_block, 4), 1) )

    ! it turns out, when the coefficients are spaghetti-ordered,
    ! we can sync only even numbers of points and save one for odd numbered
    g_spaghetti = params%g/2*2


    ! To avoid that the incoming hvy_neighbor array and active lists are outdated
    ! we synchronize them.
    t0 = MPI_Wtime()
    call updateMetadata_tree(params, tree_ID)
    call toc( "adapt_tree (update neighbors)", MPI_Wtime()-t0 )

    ! Wavelet decomposition can be done for each block individually
    ! We iterate leaf-wise, meaning that each block is investigated and coarsened and resulting new blocks are then investigated
    ! until the number of blocks is constant (no new blocks are created in the process anyways) or iteration=JMax-JMin

    ! For fine and medium neighbors of a block first sync correctly sets boundary values and all operations (CE, coarsening, even CVS) on
    ! the neighbour do not change the wavelet decomposition on this block
    ! For coarse neighbors we apply coarse extension, this effectively decouples the wavelet decomposition from this blocks

    ! In order to not investigate blocks again, we will give everyone a temporary flag that they are not wavelet decomposed yet
    t0 = MPI_Wtime()
    do k = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k, tree_ID)
        lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_UNTREATED            
    end do
    call toc( "adapt_tree (temp flag)", MPI_Wtime()-t0 )

    do while (iterate)
        lgt_n_old = lgt_n(tree_ID)


        ! ! I often need to check block values when debugging so I leave it here
        ! do k = 1, hvy_n(tree_ID)
        !     hvy_ID = hvy_active(k, tree_ID)
        !     call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        !     if (lgt_block(lgt_id, IDX_MESH_LVL) == 3 .and. lgt_block(lgt_id, IDX_TC_2) == 58720256) then
        !         write(*, '("1 I-", i0, " R-", i0, " B-", i0, " L-", i0, " Ref-", i0, " TC-", i0)') iteration, params%rank, lgt_ID, &
        !             lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2)
        !         do iy = 1,34
        !             write(*, '(34(es8.1))') hvy_block(1:34, iy, 1, 1, hvy_id)
        !         enddo
        !     endif
        ! enddo


        ! synchronize ghost nodes - required to apply wavelet filters
        ! block only needs information from medium and fine neighbors as CE will cut dependency to coarse neighbors
        t0 = MPI_Wtime()
        g_this = max(ubound(params%HD,1),ubound(params%GD,1))
        call sync_TMP_from_MF( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_ID), hvy_n(tree_ID), &
            REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
        call toc( "adapt_tree (sync lvl <- MF)", MPI_Wtime()-t0 )


        ! Wavelet-transform all blocks on this level
        ! From now on until wavelet retransform hvy_block will hold the wavelet decomposed values in spaghetti form
        t0 = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
    
            ! We compute detail coefficients on the fly here, for all blocks
            ! on the level.
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            ref_stat = lgt_block( lgt_ID, IDX_REFINE_STS )
    
            ! FWT required for a block that is on the level
            if (ref_stat == REF_TMP_UNTREATED) then
                ! hvy_tmp now is a copy with sync'ed ghost points.
                hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
                level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

                ! Compute wavelet decomposition
                ! For Jmax if dealiasing, just compute H filter
                ! Data SC/WC now in Spaghetti order
                if (level_me == params%Jmax .and. params%force_maxlevel_dealiasing) then
                    call blockFilterXYZ_vct( params, hvy_block(:,:,:,:,hvy_ID), tmp_wd(:,:,:,:,1), params%HD, lbound(params%HD, dim=1), ubound(params%HD, dim=1))
                    call inflatedMallat2spaghetti_block(params, tmp_wd, hvy_block(:,:,:,:,hvy_ID), sc_only=.true.)
                else
                    call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvy_ID))
                endif
            endif
        end do
        call toc( "adapt_tree (FWT)", MPI_Wtime()-t0 )


        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! coarseExtension: remove wavelet coefficients near a fine/coarse interface
        ! on the fine block. Does nothing in the case of CDF60, CDF40 or CDF20.
        ! This is the first coarse extension before removing blocks
        ! As every block is assumed to be WDed here we can do it on the whole tree
        t0 = MPI_Wtime()
        if (params%useCoarseExtension) then
            call coarse_extension_modify_tree(params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
            hvy_n(tree_ID), lgt_n(tree_ID), tree_ID)
        endif
        call toc( "adapt_tree (coarse_extension_modify)", MPI_Wtime()-t0 )


        !> coarseningIndicator_tree resets ALL refinement_status to 0 (all blocks, not only level)
        ! Check the entire grid where to coarsen. Note this is a wrapper for coarseningIndicator_block, which
        ! acts on a single block only.
        ! We distinguish two cases: wavelet or no wavelet (Shakespeare!). The wavelet case is the default, other indicators
        ! are used mostly for testing.
        ! The routine first sets all blocks (regardless of level) to 0 (STAY).
        ! In the wavelet case, we check the blocks for their largest detail (=wavelet coeff).
        ! The routine assigns -1 (COARSEN) to a block, if it matches the criterion. This status may however be revoked below.
        t0 = MPI_Wtime()
        if (present(hvy_mask)) then
            ! if present, the mask can also be used for thresholding (and not only the state vector). However,
            ! as the grid changes within this routine, the mask will have to be constructed in coarseningIndicator_tree
            call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, &
                ignore_maxlevel=ignore_maxlevel2, input_is_WD=.true., check_ref_TMP=.true., hvy_mask=hvy_mask)
        else
            call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, &
                ignore_maxlevel=ignore_maxlevel2, input_is_WD=.true., check_ref_TMP=.true.)
        endif
        call toc( "adapt_tree (coarseningIndicator_tree)", MPI_Wtime()-t0 )

        ! After coarseningIndicator_tree, the situation is:
        ! coarseningIndicator_tree works on LEVEL (for wavelet cases).
        ! blocks that are significant on that level now have status 0, others (on this level) have -1
        ! Any blocks on other levels have status 0.

        ! ! Coarse Extension sets WC to zero inside a block regardless of if it is kept or not
        ! ! This can mean that we delete WC which are important, filtering out parts of the flow
        ! ! In order to prevent this, neighbouring blocks are kept if it would delete necessary WC
        ! if (params%useSecurityZone) then
        !     if ((indicator=="threshold-state-vector") .or. (indicator=="primary-variables")) then
        !         ! Note: we can add the security zone also for non-lifted wavelets (although this 
        !         ! does not make much sense, but for development...)
        !         t0 = MPI_Wtime()
        !         call addSecurityZone_level( time, params, level, tree_ID, hvy_block, hvy_tmp )
        !         call toc( "adapt_tree (addSecurityZone_level)", MPI_Wtime()-t0 )
        !     endif
        ! endif
        ! ! In addSecurityZone_tree, some blocks on level J have revoked their -1 status to 0, some
        ! ! new blocks may have been created and they have the status 0 as well.
        ! ! Note: as the algorithm proceeds level-wise, a block on level J is not checked again - it
        ! ! is not possible to 'accidentally' delete the newly created blocks later on.


        ! check if block has reached maximal level, if so, remove refinement flags
        t0 = MPI_Wtime()
        if (ignore_maxlevel2 .eqv. .false.) then
            call respectJmaxJmin_tree( params, tree_ID )
        endif
        call toc( "adapt_tree (respectJmaxJmin_tree)", MPI_Wtime()-t0 )

        ! unmark blocks that cannot be coarsened due to gradedness and completeness
        t0 = MPI_Wtime()
        call ensureGradedness_tree( params, tree_ID, mark_TMP_flag=.true. )
        call toc( "adapt_tree (ensureGradedness_tree)", MPI_Wtime()-t0 )

        ! Adapt the mesh, i.e. actually merge blocks
        ! This uses the already wavelet decomposed blocks and coarsens them by effectively copying their SC into the new mother block
        t0 = MPI_Wtime()
        call executeCoarsening_WD_tree( params, hvy_block, tree_ID, mark_TMP_flag=.true.)
        call toc( "adapt_tree (executeCoarsening_level)", MPI_Wtime()-t0 )

        ! update grid lists: active list, neighbor relations, etc
        ! JB: Why is this not in executeCoarsening? This might make more sense
        t0 = MPI_Wtime()
        call updateMetadata_tree(params, tree_ID)
        call toc( "adapt_tree (update metadata)", MPI_Wtime()-t0 )

        ! iteration counter
        iteration = iteration + 1
        level = level - 1
        ! loop continues until we are on the lowest level.
        iterate = (level >= Jmin)
        ! if at Jmin_active nothing happens anymore, then we can escape the loop now.
        if ((level <= Jmin_active).and.(lgt_n(tree_ID)==lgt_n_old)) iterate = .false.

    end do

    ! After the last coarsening step, some blocks were (possibly) coarsened.
    ! If that happened (and it is the usual case), suddenly new blocks have coarser neighbors, and
    ! on those, the coarseExt needs to be done again.
    t0 = MPI_Wtime()
    if (params%useCoarseExtension) then
        call coarse_extension_modify_tree(params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
        hvy_n(tree_ID), lgt_n(tree_ID), tree_ID)
    endif
    call toc( "adapt_tree (coarse_extension_modify)", MPI_Wtime()-t0 )


    ! synch SC and WC from new coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
    ! Attention1: For finer neighbours this is not possible, so near fine interfaces we cannot reconstruct correct values.
    ! Attention2: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
    t0 = MPI_Wtime()
    call sync_SCWC_from_MC( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active(:, tree_ID), hvy_n(tree_ID), hvy_tmp, g_minus=g_spaghetti, g_plus=g_spaghetti)

    call toc( "adapt_tree (sync all <- MC)", MPI_Wtime()-t0 )


    ! After synching with coarser neighbors, the spaghetti form of the whole ghost patch is filled with values
    ! We now need to delete all values in spots of WC, leaving only the correct SC values from coarse neighbours
    ! This uses the fact that values refined on the coarse grid points with wc=0 are the copied SC values
    ! We do this in order tu bundle up the synchronizations in the step before as the modification is really cheap
    t0 = MPI_Wtime()
    if (params%useCoarseExtension) then
        call coarse_extension_modify_tree(params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
        hvy_n(tree_ID), lgt_n(tree_ID), tree_ID, sc_skip_ghosts=.true.)
    endif
    call toc( "adapt_tree (coarse_extension_modify)", MPI_Wtime()-t0 )


    ! Wavelet-reconstruct all blocks in one go
    ! Copy back old values if no coarse extension is applied and elsewise just overwrite the affected patches inside the domain
    t0 = MPI_Wtime()
    if (params%useCoarseExtension) then
        call coarse_extension_reconstruct_tree(params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
        hvy_n(tree_ID), lgt_n(tree_ID))
    else  ! just set back all values to the original ones
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)    
            hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID)
        end do
    endif
    call toc( "adapt_tree (RWT)", MPI_Wtime()-t0 )

    ! synchronize ghost nodes - final synch to update all neighbours with the new values
    ! Just once here at the end after reconstruct so everything is in order
    t0 = MPI_Wtime()
    call sync_ghosts_all( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )
    call toc( "adapt_tree (sync post)", MPI_Wtime()-t0 )


    ! At this point the coarsening is done. All blocks that can be coarsened are coarsened
    ! they may have passed several level also. Now, the distribution of blocks may no longer
    ! be balanced, so we have to balance load now
    t0 = MPI_Wtime()
    call balanceLoad_tree( params, hvy_block, tree_ID )
    call toc( "adapt_tree (balanceLoad_tree)", MPI_Wtime()-t0 )

    call toc( "adapt_tree (TOTAL)", MPI_wtime()-t1)
end subroutine



!> \brief This routine performs the coarsing of the mesh, where possible. \n
!! 1. All blocks are wavelet transformed from fine to coarse and mother blocks are created if not present. \n
!! 2. CVS Filtering and coarse extension will be applied were necessary
!! 3. Blocks will be deleted where possible if all WC are zero
!! 4. All blocks will be reconstructed from bottom to top, updating daughter SC along the way
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!
!> \note It is well possible to start with a very fine mesh and end up with only one active
!! block after this routine. You do *NOT* have to call it several times.
subroutine adapt_tree_cvs( time, params, hvy_block, tree_ID, indicator, hvy_tmp, hvy_mask, ignore_maxlevel)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    ! mask data. we can use different trees (forest module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! Optional: if the grid is not adapted to the mask, passing hvy_mask is not required.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    character(len=*), intent(in)        :: indicator
    ! during mask generation it can be required to ignore the maxlevel coarsening....life can suck, at times.
    logical, intent(in), optional       :: ignore_maxlevel

    integer(kind=ik), intent(in)        :: tree_ID
    ! loop variables
    integer(kind=ik)                    :: iteration, k_block
    real(kind=rk)                       :: t0, t1
    integer(kind=ik)                    :: ierr, k1, hvyId, lgtid, nc
    logical                             :: ignore_maxlevel2, iterate
    ! level iterator loops from Jmax_active to Jmin_active for the levelwise coarsening
    integer(kind=ik)                    :: Jmax_active, Jmin_active, i_level, level_me, Jmin, lgt_n_old, g_this

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    t0          = MPI_Wtime()
    t1          = t0
    iteration   = 0
    iterate     = .true.
    nc          = size(hvy_block, 4)
    Jmin        = params%Jmin
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    if (Jmin<1) call abort(2202243, "Currently, setting Jmin<1 is not possible")

    if (present(ignore_maxlevel)) then
        ignore_maxlevel2 = ignore_maxlevel
    else
        ignore_maxlevel2 = .false.
    endif

    ! To avoid that the incoming hvy_neighbor array and active lists are outdated
    ! we synchronize them.
    t0 = MPI_Wtime()
    call updateMetadata_tree(params, tree_ID)
    call toc( "adapt_tree (update neighbors)", MPI_Wtime()-t0 )

    !***********************************************************************
    ! Wavelet decomposition - iterate from finest level (Jmax) to coarsest (JMin)
    !    1. Synch lvl J2J and J2J-1 (ToDo: restrictive synch)
    !    2. Wavelet decompose lvl J
    !    3. Synch SC from lvl J to mother lvl J-1
    !***********************************************************************
    do i_level= Jmax_active, Jmin_active, -1
        ! this call turned out useful for biorthogonal wavelets - the coarse extension deteriorates
        ! the loadbalancing substantially, thus the best balancing possible is advised
        call balanceLoad_tree( params, hvy_block, tree_ID )

        ! synchronize ghost nodes - required to apply wavelet filters
        ! note it can NOT be merged into coarseExtensionUpdate_tree because this is empty for CDFX0 wavelets
        t0 = MPI_Wtime()
        g_this = max(ubound(params%HD,1),ubound(params%GD,1))
        call sync_ghosts_all( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_ID), hvy_n(tree_ID), g_minus=g_this, g_plus=g_this)
        call toc( "adapt_tree (sync_WD_pre)", MPI_Wtime()-t0 )

        ! Wavelet decomposition
        t0 = MPI_Wtime()
        do k_block = 1, hvy_n(tree_ID)
            hvyID = hvy_active(k_block, tree_ID)
    
            ! We compute detail coefficients on the fly here, for all blocks
            ! on the level.
            call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
            level_me = lgt_block( lgtID, IDX_MESH_LVL )
    
            ! FWT required for a block that is on the level
            if (level_me == i_level) then
                ! hvy_tmp now is a copy with sync'ed ghost points.
                hvy_tmp(:,:,:,1:nc,hvyID) = hvy_block(:,:,:,1:nc,hvyID)
    
                ! data WC/SC now in Spaghetti order
                call waveletDecomposition_block(params, hvy_block(:,:,:,1:nc,hvyID))
            endif
        end do
        call toc( "adata_tree (FWT)", MPI_Wtime()-t0 )

        ! Synch SC to mother blocks, create this block if not present
        t0 = MPI_Wtime()
        g_this = max(ubound(params%HD,1),ubound(params%GD,1))
        call sync_ghosts_all( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_ID), hvy_n(tree_ID), g_minus=g_this, g_plus=g_this)
        call toc( "adapt_tree (sync_SC2Mother)", MPI_Wtime()-t0 )
    end do


    !     !> coarseningIndicator_tree resets ALL refinement_status to 0 (all blocks, not only level)
    !     ! Check the entire grid where to coarsen. Note this is a wrapper for coarseningIndicator_block, which
    !     ! acts on a single block only.
    !     ! We distinguish two cases: wavelet or no wavelet (Shakespeare!). The wavelet case is the default, other indicators
    !     ! are used mostly for testing.
    !     ! The routine first sets all blocks (regardless of level) to 0 (STAY).
    !     ! In the wavelet case, we check the blocks for their largest detail (=wavelet coeff). If the coarseExtension is used,
    !     ! the FWT of all blocks (on this level) is performed there (even for blocks that are not modified) - this reduces load
    !     ! imbalancing. If coarseExtension is not used, the actual FWT of blocks (on this level) is performed in coarseningIndicator_tree.
    !     ! The routine assigns -1 (COARSEN) to a block, if it matches the criterion. This status may however be revoked below.
    !     t0 = MPI_Wtime()
    !     if (present(hvy_mask)) then
    !         ! if present, the mask can also be used for thresholding (and not only the state vector). However,
    !         ! as the grid changes within this routine, the mask will have to be constructed in coarseningIndicator_tree
    !         call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, ignore_maxlevel2, hvy_mask)
    !     else
    !         call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, ignore_maxlevel2)
    !     endif
    !     call toc( "adapt_tree (coarseningIndicator_tree)", MPI_Wtime()-t0 )

    !     ! After coarseningIndicator_tree, the situation is:
    !     ! coarseningIndicator_tree works on LEVEL.
    !     ! blocks that are significant on that level now have status 0, others (on this level) have -1
    !     ! Any blocks on other levels have status 0.

    !     ! here, we should refine the coarse frontier blocks

    !     ! afterwards, the newly created blocks should have status 0 (check that!!)
    !     ! can we call refinement if the refinement_status incudes "-1" ? ===> to be checked.

    !     if (params%useSecurityZone) then
    !         if ((indicator=="threshold-state-vector") .or. (indicator=="primary-variables")) then
    !             ! Note: we can add the security zone also for non-lifted wavelets (although this 
    !             ! does not make much sense, but for development...)
    !             call addSecurityZone_tree( time, params, level, tree_ID, hvy_block, hvy_tmp )
    !         endif
    !     endif

    !     ! In addSecurityZone_tree, some blocks on level J have revoked their -1 status to 0, some
    !     ! new blocks may have been created and they have the status 0 as well.
    !     ! Note: as the algorithm proceeds level-wise, a block on level J is not checked again - it
    !     ! is not possible to 'accidentally' delete the newly created blocks later on.


    !     ! check if block has reached maximal level, if so, remove refinement flags
    !     t0 = MPI_Wtime()
    !     if (ignore_maxlevel2 .eqv. .false.) then
    !         call respectJmaxJmin_tree( params, tree_ID )
    !     endif
    !     call toc( "adapt_tree (respectJmaxJmin_tree)", MPI_Wtime()-t0 )

    !     ! unmark blocks that cannot be coarsened due to gradedness and completeness
    !     t0 = MPI_Wtime()
    !     call ensureGradedness_tree( params, tree_ID )
    !     call toc( "adapt_tree (ensureGradedness_tree)", MPI_Wtime()-t0 )

    !     ! adapt the mesh, i.e. actually merge blocks
    !     ! this applies the wavelet low-pass filter as well (h*h) before decimation
    !     t0 = MPI_Wtime()
    !     call executeCoarsening_tree( params, hvy_block, tree_ID )
    !     call toc( "adapt_tree (executeCoarsening_tree)", MPI_Wtime()-t0 )

    !     ! update grid lists: active list, neighbor relations, etc
    !     t0 = MPI_Wtime()
    !     call updateMetadata_tree(params, tree_ID)
    !     call toc( "adapt_tree (update neighbors)", MPI_Wtime()-t0 )

    !     ! After the coarsening step on this level, the some blocks were (possibly) coarsened.
    !     ! If that happened (and it is the usual case), suddenly new blocks have coarser neighbors, and
    !     ! on those, the coarseExt needs to be done. Note performing the coarseExt twice on a block does not
    !     ! alter the data, but is of course not for free. The usual workflow in adapt_tree is that many blocks
    !     ! can be coarsened, and thus the 2nd call to coarseExtensionUpdate_level is much cheaper.
    !     t0 = MPI_Wtime()
    !     if (params%useCoarseExtension) then
    !         call coarseExtensionUpdate_level( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
    !         hvy_n(tree_ID),lgt_n(tree_ID), inputDataSynced=.false., level=level, hvy_details=hvy_details )
    !     endif
    !     call toc( "adapt_tree (coarse_extension2)", MPI_Wtime()-t0 )

    !     ! iteration counter (used for random coarsening criterion)
    !     iteration = iteration + 1
    !     level = level - 1
    !     ! loop continues until we are on the lowest level.
    !     iterate = (level >= Jmin)
    !     ! if at Jmin_active nothing happens anymore, then we can escape the loop now.
    !     if ((level <= Jmin_active).and.(lgt_n(tree_ID)==lgt_n_old)) iterate = .false.
    ! end do

    ! At this point the coarsening is done. All blocks that can be coarsened are coarsened
    ! they may have passed several level also. Now, the distribution of blocks may no longer
    ! be balanced, so we have to balance load now
    t0 = MPI_Wtime()
    call balanceLoad_tree( params, hvy_block, tree_ID )
    call toc( "adapt_tree (balanceLoad_tree)", MPI_Wtime()-t0 )


    call toc( "adapt_tree (TOTAL)", MPI_wtime()-t1)
end subroutine