
!##############################################################
! This routine refines/coarsens to given target_level.
! If target_level is not passed, then max_treelevel is assumed
subroutine refineToEquidistant_tree(params, hvy_block, hvy_tmp, tree_ID, target_level, verbosity)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), intent(in)      :: tree_ID
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    integer(kind=ik), intent(in), optional :: target_level
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    integer(kind=ik):: level,k,treecode_size, lgt_id, rank, hvy_id, d, newBlocks, level_this
    logical :: verbose = .false.

    if (present(verbosity)) verbose=verbosity
    if (present(target_level)) then
        level = target_level
    else
        level = params%Jmax
    endif

    d = params%dim

    ! refine/coarse to attain desired level, respectively


    ! coarsen
    do while (maxActiveLevel_tree( tree_ID )>level)
        ! check where coarsening is actually needed and set refinement status to -1 (coarsen)
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)

            level_this = lgt_block(lgt_id, IDX_MESH_LVL)

            if (level_this > level) then
                lgt_block(lgt_id, IDX_REFINE_STS) = -1
            endif
        end do

        ! this might not be necessary since we start from an admissible grid
        call ensureGradedness_tree( params, tree_ID )

        call executeCoarsening_tree( params, hvy_block, tree_ID)

        call updateMetadata_tree(params, tree_ID)
    end do

    ! refine
    do while (minActiveLevel_tree( tree_ID )<level)

        call balanceLoad_tree( params, hvy_block, tree_ID )

        ! check where refinement is actually needed
        ! Even if the blocks are equally distributed (balance_load), some CPU now create more blocks, and
        ! those might no longer fit into the memory. This routine would then fail, even though globally, the 
        ! refined data would still fit into the memory. This is relevant in postprocessing, mostly.
        ! The trick is thus to tag just as many blocks for refinement as we can allow, and not all of them. 
        ! Then, the next iteration will balance the load again and so forth. Loop continues until we are indeed equidistant.
        ! Price to pay: we may do more iterations than absolutely necessary.
        ! Gain: optimal memory exploitation.
        newBlocks = 0
        ! heavy loop: done locally on all CPU
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            level_this = lgt_block(lgt_id, IDX_MESH_LVL)

            if ((level_this < level).and.(hvy_n(tree_ID)+newBlocks < params%number_blocks)) then
                lgt_block(lgt_id, IDX_REFINE_STS) = 1
                ! this is conservative, as one block gets deleted, so only (2**d - 1) new blocks
                newBlocks = newBlocks + 2**d
            endif
        end do

        ! each CPU tags their individual blocks, so sync'ing is required.
        call synchronize_lgt_data(params, refinement_status_only=.true.)

        call ensureGradedness_tree( params, tree_ID )

        if ( params%dim == 3 ) then
            call refinementExecute3D_tree( params, hvy_block, tree_ID )
        else
            call refinementExecute2D_tree( params, hvy_block(:,:,1,:,:), tree_ID )
        end if

        call updateMetadata_tree(params, tree_ID)

        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_ID), hvy_n(tree_ID) )
    end do

    ! final balancing: this routine is often used in postprocessing (and during simulations
    ! only for initialization), it's good to have the load balanced now
    call balanceLoad_tree( params, hvy_block, tree_ID )
end subroutine
!##############################################################
