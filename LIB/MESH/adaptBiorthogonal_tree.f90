subroutine adaptBiorthogonal_tree(time, params, hvy_block, tree_ID, hvy_tmp)
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)

    integer(kind=ik), intent(in)        :: tree_ID
    ! loop variables
    integer(kind=ik)                    :: k, lgt_id
    real(kind=rk)                       :: t0, t1
    integer(kind=ik)                    :: ierr, k1, hvy_id
    logical                             :: iterate
    ! level iterator loops from Jmax_active to Jmin_active for the levelwise coarsening
    integer(kind=ik)                    :: Jmax_active, Jmin_active, level, Jmin

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    if (params%penalization) call abort(2303051, "as threshold_mask is not implemented this currently works only without penalizatin")

    ! start time
    t0 = MPI_Wtime()
    t1 = t0
    iterate   = .true.
    Jmin      = params%Jmin
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)
    level       = Jmax_active ! algorithm starts on maximum *active* level


    ! 03/2023:  The new version of this routine
    ! * perform a wavelet decomposition for the blocks
    ! * apply the coarse extension to the WC/SC


    ! To avoid that the incoming hvy_neighbor array and active lists are outdated
    ! we synchronize them.
    call updateMetadata_tree(params, tree_ID)

    !! we iterate from the highest current level to the lowest current level and then iterate further
    !! until the number of blocks is constant (note: as only coarsening
    !! is done here, no new blocks arise that could compromise the number of blocks -
    !! if it's constant, its because no more blocks are coarsened)
    do while (iterate)
        ! call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_ID), hvy_n(tree_ID) )
        call substitution_step( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
        hvy_n(tree_ID) )

        call coarseningIndicatorBiorthogonal_tree( time, params, level, hvy_block, hvy_tmp, tree_ID)

        ! ?????now we know which blocks want to coarsen. now we must delete their WC

        !> (b) check if block has reached maximal level, if so, remove refinement flags
        call respectJmaxJmin_tree( params, tree_ID )

        !> (c) unmark blocks that cannot be coarsened due to gradedness and completeness
        call ensureGradedness_tree( params, tree_ID )

        !> (d) adapt the mesh, i.e. actually merge blocks
        ! also applies the HD filter
        call executeCoarsening_tree( params, hvy_block, tree_ID, ignorePrefilter=.false. )

        ! update grid lists: active list, neighbor relations, etc
        call updateMetadata_tree(params, tree_ID)

        level = level - 1
        iterate = (level >= Jmin)
    end do

    !> At this point the coarsening is done. All blocks that can be coarsened are coarsened
    !! they may have passed several level also. Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now
    call balanceLoad_tree( params, hvy_block, tree_ID )

    call substitution_step( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
    hvy_n(tree_ID) )

    call toc( "adaptBiorthogonal_tree (TOTAL)", MPI_wtime()-t1)
end subroutine
