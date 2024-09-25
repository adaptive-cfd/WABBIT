!> \brief Refine the mash (tag / interpolate / update lists)
!
!> \details This routine first sets the refinement flag for all blocks to +1
!! and then executes the refinement directly. Blocks that cannot be refined because they
!! are already on the finest allowed level are unaltered.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!!
!! Note we assume, on call, that active lists / neighbor are up-to-date
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
! ********************************************************************************************

subroutine refine_tree( params, hvy_block, hvy_tmp, indicator, tree_ID, hvy_mask  )
    use module_indicators

    implicit none

    type (type_params), intent(inout) :: params                   !> user defined parameter structure
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    real(kind=rk), intent(inout)   :: hvy_tmp(:, :, :, :, :)
    character(len=*), intent(in)   :: indicator                   !> how to choose blocks for refinement
    integer(kind=ik), intent(in)   :: tree_ID
    real(kind=rk), intent(inout), optional   :: hvy_mask(:, :, :, :, :)

    ! cpu time variables for running time calculation
    real(kind=rk)                  :: t0, t1, t2, t_misc
    integer(kind=ik)               :: k, hvy_n_afterRefinement, lgt_id, hvy_id
    real(kind=rk) :: norm(1:params%n_eqn)

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    ! start time
    t0 = MPI_Wtime()
    t_misc = 0.0_rk

    !> (a) loop over the blocks and set their refinement status.
    t1 = MPI_Wtime()

    ! JB: test if it works at all, this is not optimized and should be a proof of concept if we can skip refining non-significant blocks and it has an impact
    ! in order to work on all blocks and work with coarseningIndicator_tree logic, we set all ref status as unset to work as in a leaf-wise loop
    ! no need for syncing afterwards as coarseningIndicator works locally
    if (indicator == "significant") then
        ! significant needs access to functions which are not available for refinemendIndicator so for now I have it here

        ! JB: Time set to 0 for testing
        if (present(hvy_mask)) call createMask_tree(params, 0.0_rk, hvy_mask, hvy_tmp, all_parts=.false.)

        ! JB: Time set to 0 for testing
        if (params%eps_normalized .and. params%coarsening_indicator/="everywhere" .and. params%coarsening_indicator/="random" &
            .and. params%coarsening_indicator/="threshold-cvs" .and. params%coarsening_indicator/="threshold-image-denoise") then
            call componentWiseNorm_tree(params, hvy_block, tree_ID, params%eps_norm, norm)
            do k = 1, params%n_eqn
                if (norm(k) <= 1.0e-9_rk) norm(k) = 1.0_rk
            enddo
        else
            norm(:) = 1.0_rk
        endif

        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
            
            ! evaluate the criterion on this block.
            call coarseningIndicator_block( params, hvy_block(:,:,:,:,hvy_ID), hvy_tmp(:,:,:,:,hvy_ID), params%coarsening_indicator, &
            lgt_block(lgt_ID, IDX_REFINE_STS), norm(1:size(hvy_block, 4)), lgt_block(lgt_ID, IDX_MESH_LVL), input_is_WD=.false., block_mask=hvy_mask(:,:,:,:,hvy_ID))
        enddo

        ! very important: CPU1 cannot decide if blocks on CPU0 have to be refined.
        ! therefore we have to sync the lgt data
        call synchronize_lgt_data( params, refinement_status_only=.true. )

        ! for lifted wavelets we need the security zone being refined as well elsewise we get problems with CE
        if (params%useSecurityZone .and. params%coarsening_indicator/="everywhere" .and. params%coarsening_indicator/="random" .and. params%useCoarseExtension .and. params%isLiftedWavelet) then
            ! if we want to add a security zone, we check for every significant block if a neighbor wants to coarsen
            ! if this is the case, we check if any significant WC would be deleted (basically checking the thresholding for this patch)
            ! in that case we set the neighbouring block to be important as well (with a temporary flag)

            ! ATTENTION: Input_is_WD=.false. here is incredibly expensive, time set to 0 for testing
            call addSecurityZone_CE_tree( 0.0_rk, params, tree_ID, hvy_block, hvy_tmp, params%coarsening_indicator, norm, ignore_maxlevel=.false., input_is_WD=.false.)
        endif

        ! coarseningIndicator and Securityzone work for coarsening, so -1, we now have to translate this to refinement range
        do k = 1, lgt_n(tree_ID)
            lgt_ID = lgt_active(k, tree_ID)

            ! set 0 (significant) to +1 (refine), -1 (non-significant) to 0 (stay as it is)
            if (lgt_block( lgt_id, IDX_REFINE_STS ) == 0) lgt_block( lgt_id, IDX_REFINE_STS ) = +1
            if (lgt_block( lgt_id, IDX_REFINE_STS ) == -1) lgt_block( lgt_id, IDX_REFINE_STS ) = 0
        enddo

    else
        call refinementIndicator_tree( params, hvy_block, tree_ID, indicator )
    endif
    call toc( "refine_tree (refinementIndicator_tree)", 141, MPI_Wtime()-t1 )


    !> (b) remove refinement flag for blocks that are on the finest level and thus
    !! cannot be refined anymore.
    t1 = MPI_Wtime()
    call respectJmaxJmin_tree( params, tree_ID )
    call toc( "refine_tree (respectJmaxJmin_tree)", 142, MPI_Wtime()-t1 )


    !> (c) ensure gradedness of mesh. If the refinement is done everywhere, there is
    !! no way gradedness can be damaged, so we skip the call in this case. However,
    !! in all other indicators, this step is very important.
    t1 = MPI_Wtime()
    if ( indicator /= "everywhere" ) then
      call ensureGradedness_tree( params, tree_ID )
    endif
    call toc( "refine_tree (ensureGradedness_tree)", 143, MPI_Wtime()-t1 )


    !---------------------------------------------------------------------------
    ! check if the refinement step will succeed or if we will run out of memory
    !---------------------------------------------------------------------------
    ! NOTE: in the simulation part, a second check is performed BEFORE this routine
    ! is called. It then aborts the time loop and dumps a backup for late resuming.
    ! If the following check fails, then this will be in postprocessing, and it will kill the code.
    hvy_n_afterRefinement = 0
    do k = 1, hvy_n(tree_ID)
        ! loop over hvy is safer, because we can tell for each mpirank if it fails or not
        ! so we detect also problems arising from load-imbalancing.
        call hvy2lgt(lgt_id, hvy_active(k, tree_ID), params%rank, params%number_blocks)

        if ( lgt_block(lgt_id, IDX_REFINE_STS) == +1) then
            ! this block will be refined, so it creates 2**D new blocks but is deleted
            hvy_n_afterRefinement = hvy_n_afterRefinement + (2**params%dim - 1)
        else
            ! this block will not be refined, so it remains.
            hvy_n_afterRefinement = hvy_n_afterRefinement + 1
        endif
    enddo

    ! oh-oh case.
    if ( hvy_n_afterRefinement > params%number_blocks ) then
        write(*,'("On rank:",i5," hvy_n=",i5," will refine to ",i5," but limit is ",i5)') &
        params%rank, hvy_n(tree_ID), hvy_n_afterRefinement, params%number_blocks

        call abort (1909181827,"[refine_tree.f90]: The refinement step will fail (not enough memory).")
    endif
    !---------------------------------------------------------------------------


    !> (d) execute refinement, interpolate the new mesh. All blocks go one level up
    !! except if they are already on the highest level.
    t1 = MPI_Wtime()
    call refinement_execute_tree( params, hvy_block, tree_ID )
    call toc( "refine_tree (refinement_execute)", 144, MPI_Wtime()-t1 )


    !> (e) as the grid changed now with the refinement, we have to update the list of
    !! active blocks so other routines can loop just over these active blocks
    !! and do not have to ensure that the active list is up-to-date
    ! update list of sorted nunmerical treecodes, used for finding blocks
    t2 = MPI_wtime()
    call updateMetadata_tree(params, tree_ID)
    t_misc = MPI_wtime() - t2


    !> (f) At this point the refinement is done. Since not all blocks are refined, namely only those
    !! that were not on Jmax, Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now. EXCEPTION: if the flag force_maxlevel_dealiasing is set
    !! then we force blocks on Jmax to coarsen, even if their details are large. Hence, each refinement
    !! step is a true "everywhere". Then, there is no need for balancing, as all mpiranks automatically
    !! hold the same number of blocks, if they started with a balanced distribution (which is the
    !! output of adapt_tree.)
    !! This of course implies any other indicator than "everywhere" requires balancing here.
    if ((params%force_maxlevel_dealiasing .eqv. .false.) .or. (indicator/="everywhere")) then
        t1 = MPI_Wtime()
        call balanceLoad_tree( params, hvy_block, tree_ID )
        call toc( "refine_tree (balanceLoad_tree)", 145, MPI_Wtime()-t1 )
    endif

    ! call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, &
    ! hvy_active(:,tree_ID), hvy_n(tree_ID), inputDataSynced=.false. )

    call toc( "refine_tree (lists+neighbors)", 146, t_misc )
    call toc( "refine_tree (TOTAL)", 140, MPI_wtime()-t0 )

end subroutine refine_tree
