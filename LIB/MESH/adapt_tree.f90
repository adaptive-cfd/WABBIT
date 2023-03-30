!! \author  engels
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
    integer(kind=ik)                    :: ierr, k1, hvy_id
    logical                             :: ignore_maxlevel2, iterate
    ! level iterator loops from Jmax_active to Jmin_active for the levelwise coarsening
    integer(kind=ik)                    :: Jmax_active, Jmin_active, level, Jmin

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    ! start time
    t0 = MPI_Wtime()
    t1 = t0
    iteration = 0
    iterate   = .true.
    Jmin      = params%Jmin
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)
    level       = Jmax_active ! algorithm starts on maximum *active* level

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

    !! we iterate from the highest current level to the lowest current level and then iterate further
    !! until the number of blocks is constant (note: as only coarsening
    !! is done here, no new blocks arise that could compromise the number of blocks -
    !! if it's constant, its because no more blocks are coarsened)
    do while (iterate)
        !> (a) check where coarsening is possible
        ! ------------------------------------------------------------------------------------
        ! first: synchronize ghost nodes - thresholding on block with ghost nodes
        ! synchronize ghostnodes, grid has changed, not in the first one, but in later loops
        t0 = MPI_Wtime()
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_ID), hvy_n(tree_ID) )
        call toc( "adapt_tree (sync_ghosts)", MPI_Wtime()-t0 )

        ! coarse extension: remove wavelet coefficients near a fine/coarse interface
        ! on the fine block. Does nothing in the case of CDF40 or CDF20.
        t0 = MPI_Wtime()
        call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
        hvy_n(tree_ID), inputDataSynced=.true. )
        call toc( "adapt_tree (coarse_extension)", MPI_Wtime()-t0 )

        !! calculate detail on the entire grid. Note this is a wrapper for coarseningIndicator_block, which
        !! acts on a single block only
        t0 = MPI_Wtime()
        if (present(hvy_mask)) then
            ! if present, the mask can also be used for thresholding (and not only the state vector). However,
            ! as the grid changes within this routine, the mask will have to be constructed in coarseningIndicator_tree
            call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, ignore_maxlevel2, hvy_mask)
        else
            call coarseningIndicator_tree( time, params, level, hvy_block, hvy_tmp, tree_ID, indicator, iteration, ignore_maxlevel2)
        endif
        call toc( "adapt_tree (coarseningIndicator_tree)", MPI_Wtime()-t0 )


        !> (b) check if block has reached maximal level, if so, remove refinement flags
        t0 = MPI_Wtime()
        if (ignore_maxlevel2 .eqv. .false.) then
            call respectJmaxJmin_tree( params, tree_ID )
        endif
        call toc( "adapt_tree (respectJmaxJmin_tree)", MPI_Wtime()-t0 )

        !> (c) unmark blocks that cannot be coarsened due to gradedness and completeness
        t0 = MPI_Wtime()
        call ensureGradedness_tree( params, tree_ID )
        call toc( "adapt_tree (ensureGradedness_tree)", MPI_Wtime()-t0 )

        !> (d) adapt the mesh, i.e. actually merge blocks
        ! this applies the wavelet low-pass filter as well (h*h) before decimation
        t0 = MPI_Wtime()
        call executeCoarsening_tree( params, hvy_block, tree_ID, .false. )
        call toc( "adapt_tree (executeCoarsening_tree)", MPI_Wtime()-t0 )

        ! update grid lists: active list, neighbor relations, etc
        t0 = MPI_Wtime()
        call updateMetadata_tree(params, tree_ID)
        call toc( "adapt_tree (update neighbors)", MPI_Wtime()-t0 )

        ! iteration counter (used for random coarsening criterion)
        iteration = iteration + 1
        level = level - 1
        ! loop continues until we are on the lowest level.
        iterate = (level >= Jmin)
    end do

    !> At this point the coarsening is done. All blocks that can be coarsened are coarsened
    !! they may have passed several level also. Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now
    t0 = MPI_Wtime()
    call balanceLoad_tree( params, hvy_block, tree_ID )
    call toc( "adapt_tree (balanceLoad_tree)", MPI_Wtime()-t0 )

    ! final coarse extension step
    t0 = MPI_Wtime()
    call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
    hvy_n(tree_ID), inputDataSynced=.false. )
    call toc( "adapt_tree (coarse_extension)", MPI_Wtime()-t0 )

    call toc( "adapt_tree (TOTAL)", MPI_wtime()-t1)
end subroutine
