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
subroutine adapt_tree( time, params, hvy_block, tree_ID, indicator, hvy_tmp, hvy_mask, ignore_maxlevel, log_blocks, log_iterations)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> mask data. we can use different trees (forest module) to generate time-dependent/indenpedent
    !! mask functions separately. This makes the mask routines tree-level routines (and no longer
    !! block level) so the physics modules have to provide an interface to create the mask at a tree
    !! level. All parts of the mask shall be included: chi, boundary values, sponges.
    !! Optional: if the grid is not adapted to the mask, passing hvy_mask is not required.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    character(len=*), intent(in)        :: indicator
    !> during mask generation it can be required to ignore the maxlevel coarsening....life can suck, at times.
    logical, intent(in), optional       :: ignore_maxlevel
    !> some information that we can log so that we know what happened in the loops
    integer, intent(out), optional       :: log_blocks(:), log_iterations


    integer(kind=ik), intent(in)        :: tree_ID
    ! loop variables
    integer(kind=ik)                    :: iteration, k, lgt_id
    real(kind=rk)                       :: t_block, t_all, t_loop
    integer(kind=ik)                    :: Jmax_active, Jmin_active, level, ierr, k1, hvy_id
    logical                             :: ignore_maxlevel2, iterate, toBeManipulated
    integer(kind=ik)                    :: level_me, ref_stat, Jmin, lgt_n_old, g_this, g_spaghetti
    character(len=clong)                :: toc_statement

    real(kind=rk) :: norm(1:params%n_eqn)

    ! testing
    integer(kind=ik)                    :: ix, iy, lgt_id_sisters(2**params%dim)


    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    t_block          = MPI_Wtime()
    t_all          = MPI_Wtime()
    iteration   = 0
    iterate     = .true.
    Jmin        = params%Jmin
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)
    level       = Jmax_active ! algorithm starts on maximum *active* level
    norm(:)     = 1.0_rk

    if (Jmin<1) call abort(2202243, "Currently, setting Jmin<1 is not possible")

    ignore_maxlevel2 = .false.
    if (present(ignore_maxlevel)) ignore_maxlevel2 = ignore_maxlevel

    ! it turns out, when the coefficients are spaghetti-ordered,
    ! we can sync only even numbers of points and save one for odd numbered
    g_spaghetti = params%g/2*2


    ! To avoid that the incoming hvy_neighbor array and active lists are outdated
    ! we synchronize them.
    call updateMetadata_tree(params, tree_ID)

    ! Wavelet decomposition can be done for each block individually
    ! We iterate leaf-wise, meaning that each block is investigated and coarsened and resulting new blocks are then investigated
    ! until the number of blocks is constant (no new blocks are created in the process anyways) or iteration=JMax-JMin

    ! For fine and medium neighbors of a block first sync correctly sets boundary values and all operations (CE, coarsening, even CVS) on
    ! the neighbour do not change the wavelet decomposition on this block
    ! For coarse neighbors we apply coarse extension, this effectively decouples the wavelet decomposition from this blocks

    ! In order to not investigate blocks again, we will give everyone a temporary flag that they are not wavelet decomposed yet
    do k = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k, tree_ID)
        lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_UNTREATED            
    end do

    do while (iterate)
        t_loop = MPI_Wtime()
        lgt_n_old = lgt_n(tree_ID)


        ! ! I often need to check block values when debugging so I leave it here
        ! do k = 1, hvy_n(tree_ID)
        !     hvy_ID = hvy_active(k, tree_ID)
        !     call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        !     if (lgt_block(lgt_id, IDX_MESH_LVL) == 3 .and. lgt_block(lgt_id, IDX_TC_2) == 52428800) then
        !         write(*, '("1 I-", i0, " R-", i0, " B-", i0, " L-", i0, " Ref-", i0, " TC-", i0)') iteration, params%rank, lgt_ID, &
        !             lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2)
        !         do iy = 1,34
        !             write(*, '(34(es8.1))') hvy_block(1:34, iy, 1, 1, hvy_id)
        !         enddo
        !     endif
        ! enddo


        ! synchronize ghost nodes - required to apply wavelet filters
        ! block only needs information from medium and fine neighbors as CE will cut dependency to coarse neighbors
        t_block = MPI_Wtime()
        g_this = max(ubound(params%HD,1),ubound(params%GD,1))
        ! for coarse extension we are not dependend on coarser neighbors so lets skip the syncing
        if (params%useCoarseExtension .and. params%isLiftedWavelet) then
            call sync_TMP_from_MF( params, hvy_block, tree_ID, REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
            call toc( "adapt_tree (sync lvl <- MF)", 103, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " sync lvl <- MF)"
            call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        ! unlifted wavelets need coarser neighbor values for their WC so we need to sync them too
        else
            call sync_TMP_from_all( params, hvy_block, tree_ID, REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
            call toc( "adapt_tree (sync lvl <- all)", 103, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " sync lvl <- all)"
            call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        endif


        ! Wavelet-transform all remaining non-decomposed blocks
        ! From now on until wavelet retransform hvy_block will hold the wavelet decomposed values in spaghetti form
        t_block = MPI_Wtime()
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
                    call blockFilterXYZ_vct( params, hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID), hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID), params%HD, &
                        lbound(params%HD, dim=1), ubound(params%HD, dim=1), do_restriction=.true.)
                else
                    call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvy_ID))
                endif
            endif
        end do
        call toc( "adapt_tree (FWT)", 104, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " FWT)"
        call toc( toc_statement, 1200+iteration, MPI_Wtime()-t_block )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! coarseExtension: remove wavelet coefficients near a fine/coarse interface
        ! on the fine block. Does nothing in the case of CDF60, CDF40 or CDF20.
        ! This is the first coarse extension before removing blocks
        ! As every block is assumed to be WDed here we can do it on the whole tree
        t_block = MPI_Wtime()
        if (params%useCoarseExtension .and. params%isLiftedWavelet) then
            call coarse_extension_modify_tree(params, hvy_block, hvy_tmp, tree_ID)
        endif
        call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )     

        !> coarseningIndicator_tree resets ALL refinement_status to 0 (all blocks, not only level)
        ! Check the entire grid where to coarsen. Note this is a wrapper for coarseningIndicator_block, which
        ! acts on a single block only.
        ! We distinguish two cases: wavelet or no wavelet (Shakespeare!). The wavelet case is the default, other indicators
        ! are used mostly for testing.
        ! The routine first sets all blocks (regardless of level) to 0 (STAY).
        ! In the wavelet case, we check the blocks for their largest detail (=wavelet coeff).
        ! The routine assigns -1 (COARSEN) to a block, if it matches the criterion. This status may however be revoked below.
        t_block = MPI_Wtime()
        if (present(hvy_mask)) then
            ! if present, the mask can also be used for thresholding (and not only the state vector). However,
            ! as the grid changes within this routine, the mask will have to be constructed in coarseningIndicator_tree
            call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, &
                ignore_maxlevel=ignore_maxlevel2, input_is_WD=.true., leaf_loop_TMP=.true., hvy_mask=hvy_mask, norm_inout=norm)
        else
            call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, &
                ignore_maxlevel=ignore_maxlevel2, input_is_WD=.true., leaf_loop_TMP=.true., norm_inout=norm)
        endif
        call toc( "adapt_tree (coarseningIndicator_tree)", 106, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " coarseningIndicator_tree)"
        call toc( toc_statement, 1300+iteration, MPI_Wtime()-t_block )

        ! After coarseningIndicator_tree, the situation is:
        ! coarseningIndicator_tree works on all blocks (for wavelet cases).
        ! blocks that are significant now have status 0, others have -1

        if (params%useSecurityZone .and. indicator/="everywhere" .and. indicator/="random" .and. params%useCoarseExtension .and. params%isLiftedWavelet) then
            ! if we want to add a security zone, we check for every significant block if a neighbor wants to coarsen
            ! if this is the case, we check if any significant WC would be deleted (basically checking the thresholding for this patch)
            ! in that case we set the neighbouring block to be important as well (with a temporary flag)
            t_block = MPI_Wtime()
            call addSecurityZone_CE_tree( time, params, tree_ID, hvy_block, hvy_tmp, indicator, norm, .true.)
            call toc( "adapt_tree (security_zone_check)", 107, MPI_Wtime()-t_block)

            write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " security_zone_check)"
            call toc( toc_statement, 1400+iteration, MPI_Wtime()-t_block )
        endif

        ! check if block has reached maximal level, if so, remove refinement flags
        t_block = MPI_Wtime()
        if (ignore_maxlevel2 .eqv. .false.) then
            call respectJmaxJmin_tree( params, tree_ID )
        endif
        call toc( "adapt_tree (respectJmaxJmin_tree)", 108, MPI_Wtime()-t_block )

        ! unmark blocks that cannot be coarsened due to gradedness and completeness
        t_block = MPI_Wtime()
        call ensureGradedness_tree( params, tree_ID, mark_TMP_flag=.true. )
        call toc( "adapt_tree (ensureGradedness_tree)", 109, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " ensureGradedness_tree)"
        call toc( toc_statement, 1500+iteration, MPI_Wtime()-t_block )

        ! Adapt the mesh, i.e. actually merge blocks
        ! This uses the already wavelet decomposed blocks and coarsens them by effectively copying their SC into the new mother block
        t_block = MPI_Wtime()
        call executeCoarsening_WD_tree( params, hvy_block, tree_ID, mark_TMP_flag=.true.)
        call toc( "adapt_tree (executeCoarsening_tree)", 110, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " executeCoarsening_tree)"
        call toc( toc_statement, 1600+iteration, MPI_Wtime()-t_block )

        ! if (params%rank == 0) then
        !     do k = 1, lgt_n(tree_ID)
        !         lgt_ID = lgt_active(k, tree_ID)
        !         write(*, '("2 - R0 - Exists BL-", i0, " L-", i0, " Ref-", i0, "TC-", i0, " - ", b32.32)') lgt_ID, &
        !             lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2), lgt_block(lgt_id, IDX_TC_2)
        !     enddo
        ! endif

        ! update grid lists: active list, neighbor relations, etc
        ! JB: Why is this not in executeCoarsening? This might make more sense
        call updateMetadata_tree(params, tree_ID, verbose_check=.true.)

        ! iteration counter
        iteration = iteration + 1
        level = level - 1
        ! loop continues until we are on the lowest level.
        iterate = (level >= Jmin)
        ! if at Jmin_active nothing happens anymore, then we can escape the loop now.
        if ((level <= Jmin_active).and.(lgt_n(tree_ID)==lgt_n_old)) iterate = .false.

        ! log some statistics - amount of blocks, after iteration is increased so that it's 1-based
        if (present(log_blocks)) log_blocks(iteration) = lgt_n(tree_ID)
        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " TOTAL)"
        call toc( toc_statement, 1700+iteration, MPI_Wtime()-t_loop )
    end do

    ! log some statistics - amount of iterations
    if (present(log_iterations)) log_iterations = iteration

    if (params%useCoarseExtension .and. params%isLiftedWavelet) then
        ! If iteration hit Jmax-Jmin, some blocks might have been coarsened in last iteration.
        ! If that happened, suddenly new blocks have coarser neighbors, and on those, the coarseExt needs to be done again.
        t_block = MPI_Wtime()
        call coarse_extension_modify_tree(params, hvy_block, hvy_tmp, tree_ID)
        call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )

        ! synch SC and WC from new coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
        ! Attention1: For finer neighbours this is not possible, so near fine interfaces we cannot reconstruct correct values.
        ! Attention2: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
        t_block = MPI_Wtime()
        call sync_SCWC_from_MC( params, hvy_block, tree_ID, hvy_tmp, g_minus=g_spaghetti, g_plus=g_spaghetti)
        call toc( "adapt_tree (sync all <- MC)", 112, MPI_Wtime()-t_block )

        ! After synching with coarser neighbors, the spaghetti form of the whole ghost patch is filled with values
        ! We now need to delete all values in spots of WC, leaving only the correct SC values from coarse neighbours
        ! This uses the fact that values refined on the coarse grid points with wc=0 are the copied SC values
        ! We do this in order to bundle up the synchronizations in the step before as the modification is really cheap
        t_block = MPI_Wtime()
        call coarse_extension_modify_tree(params, hvy_block, hvy_tmp, tree_ID, sc_skip_ghosts=.true.)
        call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )


        ! Wavelet-reconstruct all blocks in one go
        ! Copy back old values if no coarse extension is applied and elsewise just overwrite the affected patches inside the domain
        t_block = MPI_Wtime()
        call coarse_extension_reconstruct_tree(params, hvy_block, hvy_tmp, tree_ID)
        call toc( "adapt_tree (reset or RWT)", 113, MPI_Wtime()-t_block )
    else 
        ! If there is no coarse extension to be done then we can simply set back all values directly
        ! just set back all values to the original ones
        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)    
            hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID)
        end do
        call toc( "adapt_tree (reset or RWT)", 113, MPI_Wtime()-t_block )
    endif


    ! At this point the coarsening is done. All blocks that can be coarsened are coarsened
    ! they may have passed several level also. Rebalance for optimal loadbalancing
    ! ToDo: we sync directly afterwards so we could only transfer inner points, but this needs the buffering from xfer_bock_data
    t_block = MPI_Wtime()
    call balanceLoad_tree( params, hvy_block, tree_ID)
    call toc( "adapt_tree (balanceLoad_tree)", 115, MPI_Wtime()-t_block )


    ! synchronize ghost nodes - final synch to update all neighbours with the new values
    ! Just once here at the end after reconstruct so everything is in order
    t_block = MPI_Wtime()
    call sync_ghosts_tree( params, hvy_block, tree_ID )
    call toc( "adapt_tree (sync post)", 114, MPI_Wtime()-t_block )


    call toc( "adapt_tree (TOTAL)", 100, MPI_wtime()-t_all)
end subroutine



!> \brief This routine performs the coarsing of the mesh, where possible. \n
!! 1. All blocks are wavelet transformed from fine to coarse and mother blocks are created if not present. \n
!! 2. coarse extension will be applied for SC were necessary
!! 3. Blocks will be deleted where possible if all WC are zero
!! 4. CVS filtering and coarse extension will be applied for WC
!! 5. All blocks will be reconstructed from JMin to JMax, updating daughter SC along the way
!
!> \note It is well possible to start with a very fine mesh and end up with only one active
!! block after this routine. You do *NOT* have to call it several times.
subroutine adapt_tree_cvs( time, params, hvy_block, tree_ID, indicator, hvy_tmp, hvy_mask, ignore_maxlevel, log_blocks, log_iterations)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> mask data. we can use different trees (forest module) to generate time-dependent/indenpedent
    !! mask functions separately. This makes the mask routines tree-level routines (and no longer
    !! block level) so the physics modules have to provide an interface to create the mask at a tree
    !! level. All parts of the mask shall be included: chi, boundary values, sponges.
    !! Optional: if the grid is not adapted to the mask, passing hvy_mask is not required.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    character(len=*), intent(in)        :: indicator
    !> during mask generation it can be required to ignore the maxlevel coarsening....life can suck, at times.
    logical, intent(in), optional       :: ignore_maxlevel
    !> some information that we can log so that we know what happened in the loops
    integer, intent(out), optional       :: log_blocks(:), log_iterations


    integer(kind=ik), intent(in)        :: tree_ID
    ! loop variables
    integer(kind=ik)                    :: iteration, k, lgt_id
    real(kind=rk)                       :: t_block, t_all, t_loop
    integer(kind=ik)                    :: Jmax_active, Jmin_active, level, ierr, k1, hvy_id
    logical                             :: ignore_maxlevel2, iterate, toBeManipulated
    integer(kind=ik)                    :: level_me, ref_stat, Jmin, lgt_n_old, g_this, g_spaghetti
    character(len=clong)                :: toc_statement

    real(kind=rk) :: norm(1:params%n_eqn)

    ! testing
    integer(kind=ik)                    :: ix, iy, lgt_id_sisters(2**params%dim)


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

    if (Jmin<1) call abort(2202243, "Currently, setting Jmin<1 is not possible")

    ignore_maxlevel2 = .false.
    if (present(ignore_maxlevel)) ignore_maxlevel2 = ignore_maxlevel

    ! it turns out, when the coefficients are spaghetti-ordered,
    ! we can sync only even numbers of points and save one for odd numbered
    g_spaghetti = params%g/2*2


    ! To avoid that the incoming hvy_neighbor array and active lists are outdated
    ! we synchronize them.
    call updateMetadata_tree(params, tree_ID)

    ! Wavelet decomposition can be done for each block individually
    ! We iterate leaf-wise, meaning that each block is investigated and coarsened and resulting new blocks are then investigated
    ! until the number of blocks is constant (no new blocks are created in the process anyways) or iteration=JMax-JMin

    ! For fine and medium neighbors of a block first sync correctly sets boundary values and all operations (CE, coarsening, even CVS) on
    ! the neighbour do not change the wavelet decomposition on this block
    ! For coarse neighbors we apply coarse extension, this effectively decouples the wavelet decomposition from this blocks

    ! In order to not investigate blocks again, we will give everyone a temporary flag that they are not wavelet decomposed yet
    do k = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k, tree_ID)
        lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_UNTREATED            
    end do

    ! iteration for wavelet decomposition - no blocks will be deleted during the loop!
    do while (iterate)
        t_loop = MPI_Wtime()
        lgt_n_old = lgt_n(tree_ID)


        ! ! I often need to check block values when debugging so I leave it here
        ! do k = 1, hvy_n(tree_ID)
        !     hvy_ID = hvy_active(k, tree_ID)
        !     call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        !     if (lgt_block(lgt_id, IDX_MESH_LVL) == 3 .and. lgt_block(lgt_id, IDX_TC_2) == 52428800) then
        !         write(*, '("1 I-", i0, " R-", i0, " B-", i0, " L-", i0, " Ref-", i0, " TC-", i0)') iteration, params%rank, lgt_ID, &
        !             lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2)
        !         do iy = 1,34
        !             write(*, '(34(es8.1))') hvy_block(1:34, iy, 1, 1, hvy_id)
        !         enddo
        !     endif
        ! enddo


        ! synchronize ghost nodes - required to apply wavelet filters
        ! block only needs information from medium and fine neighbors as CE will cut dependency to coarse neighbors
        t_block = MPI_Wtime()
        g_this = max(ubound(params%HD,1),ubound(params%GD,1))
        ! for coarse extension we are not dependend on coarser neighbors so lets skip the syncing
        if (params%isLiftedWavelet) then
            call sync_TMP_from_MF( params, hvy_block, tree_ID, REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
            call toc( "adapt_tree (sync lvl <- MF)", 103, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " sync lvl <- MF)"
            call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        ! unlifted wavelets need coarser neighbor values for their WC so we need to sync them too
        else
            call sync_TMP_from_all( params, hvy_block, tree_ID, REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
            call toc( "adapt_tree (sync lvl <- all)", 103, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " sync lvl <- all)"
            call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        endif


        ! Wavelet-transform all remaining non-decomposed blocks
        ! From now on until wavelet retransform hvy_block will hold the wavelet decomposed values in spaghetti form
        t_block = MPI_Wtime()
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
                    call blockFilterXYZ_vct( params, hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID), hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID), params%HD, &
                        lbound(params%HD, dim=1), ubound(params%HD, dim=1), do_restriction=.true.)
                else
                    call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvy_ID))
                endif
            endif
        end do
        call toc( "adapt_tree (FWT)", 104, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " FWT)"
        call toc( toc_statement, 1200+iteration, MPI_Wtime()-t_block )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! coarseExtension: remove wavelet coefficients near a fine/coarse interface
        ! on the fine block. Does nothing in the case of CDF60, CDF40 or CDF20.
        ! This is the first coarse extension before removing blocks
        ! As every block is assumed to be WDed here we can do it on the whole tree
        t_block = MPI_Wtime()
        call coarse_extension_modify_tree(params, hvy_block, hvy_tmp, tree_ID)
        call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Prepare refinement status for updating the mother blocks
        do k = 1, lgt_n(tree_ID)
            lgt_ID = lgt_active(k, tree_ID)
            if (any(lgt_block(lgt_ID, IDX_REFINE_STS) == (/ REF_TMP_UNTREATED , REF_TMP_TREATED_COARSEN /))) then
                lgt_block(lgt_ID, IDX_REFINE_STS) = -1
            endif
        end do

        ! respect lower level boundary
        t_block = MPI_Wtime()
        call respectJmaxJmin_tree( params, tree_ID )
        call toc( "adapt_tree (respectJmaxJmin_tree)", 108, MPI_Wtime()-t_block )

        ! watch for completeness
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call ensure_completeness( params, lgt_id, hvy_family(hvy_ID, 2:1+2**params%dim), mark_TMP_flag=.true. )
        enddo

        ! after locally modifying refinement statusses, we need to synchronize light data
        call synchronize_lgt_data( params, refinement_status_only=.true. )

        write(*, '("R", i0, " Loop ", i0, " before sending data")') params%rank, iteration

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Create or update mother blocks with SC of daughter blocks
        ! This uses the already wavelet decomposed blocks and copies their SC into the new mother block
        t_block = MPI_Wtime()
        call executeCoarsening_WD_tree( params, hvy_block, tree_ID, mark_TMP_flag=.true., no_deletion=.true.)
        call toc( "adapt_tree (update_mothers)", 111, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " update_mothers)"
        call toc( toc_statement, 1650+iteration, MPI_Wtime()-t_block )

        write(*, '("R", i0, " Loop ", i0, " after sending data")') params%rank, iteration

        ! update grid lists: active list, neighbor relations, etc
        ! Attention: Here we have an overfull grid, with the way updateNeighbors is build it still chooses neighbors as from a normal grid
        ! Meaning: First same-lvl J, then coarser J-1, then finer J+1
        ! For finished full WDed grids, no cell should have a finer neighbor but only the same-lvl mothers of those as neighbors
        call updateMetadata_tree(params, tree_ID, verbose_check=.true.)

        write(*, '("R", i0, " Loop ", i0, " after updateMetadata")') params%rank, iteration

        ! iteration counter
        iteration = iteration + 1
        level = level - 1
        ! loop continues until we are on the lowest level.
        iterate = (level >= Jmin)

        ! log some statistics - amount of blocks, after iteration it shoul be increased, the difference showing how many blocks we treated
        if (present(log_blocks)) log_blocks(iteration) = lgt_n(tree_ID)
        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " TOTAL)"
        call toc( toc_statement, 1700+iteration, MPI_Wtime()-t_loop )
    end do

    write(*, '("Finished first loop")')

    ! log some statistics - amount of iterations
    if (present(log_iterations)) log_iterations = iteration

    !> coarseningIndicator_tree resets ALL refinement_status to 0 (all blocks, not only level)
    ! Check the entire grid where to coarsen. Note this is a wrapper for coarseningIndicator_block, which
    ! acts on a single block only.
    ! We distinguish two cases: wavelet or no wavelet (Shakespeare!). The wavelet case is the default, other indicators
    ! are used mostly for testing.
    ! The routine first sets all blocks (regardless of level) to 0 (STAY).
    ! In the wavelet case, we check the blocks for their largest detail (=wavelet coeff).
    ! The routine assigns -1 (COARSEN) to a block, if it matches the criterion. This status may however be revoked below.
    t_block = MPI_Wtime()
    if (present(hvy_mask)) then
        ! if present, the mask can also be used for thresholding (and not only the state vector). However,
        ! as the grid changes within this routine, the mask will have to be constructed in coarseningIndicator_tree
        call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, &
            ignore_maxlevel=ignore_maxlevel2, input_is_WD=.true., leaf_loop_TMP=.true., hvy_mask=hvy_mask, norm_inout=norm)
    else
        call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, &
            ignore_maxlevel=ignore_maxlevel2, input_is_WD=.true., leaf_loop_TMP=.true., norm_inout=norm)
    endif
    call toc( "adapt_tree (coarseningIndicator_tree)", 106, MPI_Wtime()-t_block )

    ! After coarseningIndicator_tree, the situation is:
    ! coarseningIndicator_tree works on all blocks (for wavelet cases).
    ! blocks that are significant now have status 0, others have -1

    if (indicator/="everywhere" .and. indicator/="random") then
        ! if we want to add a security zone, we check for every significant block if a neighbor wants to coarsen
        ! if this is the case, we check if any significant WC would be deleted (basically checking the thresholding for this patch)
        ! in that case we set the neighbouring block to be important as well (with a temporary flag)
        t_block = MPI_Wtime()
        call addSecurityZone_CE_tree( time, params, tree_ID, hvy_block, hvy_tmp, indicator, norm, .true.)
        call toc( "adapt_tree (security_zone_check)", 107, MPI_Wtime()-t_block)
    endif

    ! check if block has reached maximal level, if so, remove refinement flags
    t_block = MPI_Wtime()
    if (ignore_maxlevel2 .eqv. .false.) then
        call respectJmaxJmin_tree( params, tree_ID )
    endif
    call toc( "adapt_tree (respectJmaxJmin_tree)", 108, MPI_Wtime()-t_block )

    ! unmark blocks that cannot be coarsened due to gradedness and completeness, in one go this should converge to the final grid
    t_block = MPI_Wtime()
    call ensureGradedness_tree( params, tree_ID, mark_TMP_flag=.true. )
    call toc( "adapt_tree (ensureGradedness_tree)", 109, MPI_Wtime()-t_block )

    ! we do not need to move any cell data anymore, any block with -1 can simply be deleted
    do k = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k, tree_ID)
        if ( lgt_block(lgt_ID, IDX_REFINE_STS) == -1) then
            ! delete blocks
            lgt_block(lgt_ID, :) = -1
        endif
        ! reset all refinement flags
        lgt_block(lgt_ID, IDX_REFINE_STS) = 0
    enddo

    ! update grid lists: active list, neighbor relations, etc
    call updateMetadata_tree(params, tree_ID, verbose_check=.true.)

    ! This is the actual important coarse extension after all cells have been deleted, which modifies now the lasting c-f interfaces
    t_block = MPI_Wtime()
    call coarse_extension_modify_tree(params, hvy_block, hvy_tmp, tree_ID)
    call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )

    ! Apply CVS filtering here if wanted

    Jmax_active = maxActiveLevel_tree(tree_ID)
    iteration = 0
    do level = Jmin, Jmax_active

        ! synch SC and WC from coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
        ! At this point no cells should have fine-level neighbors anyways
        ! Attention: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
        t_block = MPI_Wtime()
        call sync_SCWC_from_MC( params, hvy_block, tree_ID, hvy_tmp, g_minus=g_spaghetti, g_plus=g_spaghetti, level=level)
        call toc( "adapt_tree (sync level <- MC)", 112, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " sync level <- MC)"
        call toc( toc_statement, 1750+iteration, MPI_Wtime()-t_block )

        ! After synching with coarser neighbors, the spaghetti form of the whole ghost patch is filled with values
        ! We now need to delete all values in spots of WC, leaving only the correct SC values from coarse neighbours
        ! This uses the fact that values refined on the coarse grid points with wc=0 are the copied SC values
        ! We do this in order to bundle up the synchronizations in the step before as the modification is really cheap
        t_block = MPI_Wtime()
        call coarse_extension_modify_tree(params, hvy_block, hvy_tmp, tree_ID, sc_skip_ghosts=.true.)
        call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )

        ! Wavelet-reconstruct blocks on level
        ! Copy back old values if no coarse extension is applied and elsewise just overwrite the affected patches inside the domain
        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

            ! RWT required for a block that is on the level
            if (level_me == level) then
                ! Compute wavelet reconstruction
                call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvy_ID))

                ! make copy of data in hvy_tmp as we sync from coarser neighbors using that and I do not want to write another logic currently
                hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
            endif
        end do
        call toc( "adapt_tree (RWT)", 113, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " RWT)"
        call toc( toc_statement, 1850+iteration, MPI_Wtime()-t_block )

        ! now we need to update the SC of all daughters
        t_block = MPI_Wtime()
        call sync_2_daughters_level(params, hvy_block, tree_ID, level)
        call toc( "adapt_tree (sync level 2 daughters)", 114, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " sync level 2 daughters)"
        call toc( toc_statement, 1900+iteration, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (it ", iteration, " TOTAL)"
        call toc( toc_statement, 1950+iteration, MPI_Wtime()-t_loop )

        iteration = iteration + 1
    enddo

    ! delete all blocks with daughters
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
        level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

        ! RWT required for a block that is on the level
        if ( any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) then
            ! we can savely delete blocks as long as we do not update family relations
            lgt_block( lgt_ID, : ) = -1
        endif
    end do


    ! update grid lists: active list, neighbor relations, etc
    call updateMetadata_tree(params, tree_ID, verbose_check=.true.)


    ! At this point the coarsening is done. All blocks that can be coarsened are coarsened
    ! they may have passed several level also. Rebalance for optimal loadbalancing
    ! ToDo: we sync directly afterwards so we could only transfer inner points, but this needs the buffering from xfer_bock_data
    t_block = MPI_Wtime()
    call balanceLoad_tree( params, hvy_block, tree_ID)
    call toc( "adapt_tree (balanceLoad_tree)", 115, MPI_Wtime()-t_block )


    ! synchronize ghost nodes - final synch to update all neighbours with the new values
    ! Just once here at the end after reconstruct so everything is in order
    t_block = MPI_Wtime()
    call sync_ghosts_tree( params, hvy_block, tree_ID )
    call toc( "adapt_tree (sync post)", 114, MPI_Wtime()-t_block )


    call toc( "adapt_tree (TOTAL)", 100, MPI_wtime()-t_all)
end subroutine