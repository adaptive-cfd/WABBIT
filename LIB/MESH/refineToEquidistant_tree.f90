
!##############################################################
! This routine refines/coarsens to given target_level.
! If target_level is not passed, then max_treelevel is assumed
subroutine refineToEquidistant_tree(params, hvy_block, hvy_tmp, tree_ID, target_level, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), intent(in)      :: tree_ID
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    integer(kind=ik), intent(in), optional :: target_level
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    integer(kind=ik):: level,k,treecode_size, lgt_id
    logical :: verbose = .false.

    if (present(verbosity)) verbose=verbosity
    if (present(target_level)) then
        level = target_level
    else
        level = params%Jmax
    endif

    ! refine/coarse to attain desired level, respectively
    ! coarsen
    do while (maxActiveLevel_tree( tree_ID )>level)
        ! check where coarsening is actually needed and set refinement status to -1 (coarsen)
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)

            if (treecode_size(lgt_block(lgt_id,:), params%Jmax) > level) then
                lgt_block(lgt_id, params%Jmax + IDX_REFINE_STS) = -1
            endif
        end do

        ! this might not be necessary since we start from an admissible grid
        call ensureGradedness_tree( params, tree_ID )

        call executeCoarsening_tree( params, hvy_block, tree_ID)

        call updateMetadata_tree(params, tree_ID)
    end do

    ! refine
    do while (minActiveLevel_tree( tree_ID )<level)
        ! check where refinement is actually needed
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)

            if (treecode_size(lgt_block(lgt_id,:), params%Jmax) < level) then
                lgt_block(lgt_id, params%Jmax + IDX_REFINE_STS) = 1
            endif
        end do

        call ensureGradedness_tree( params, tree_ID )

        if ( params%dim == 3 ) then
            call refinementExecute3D_tree( params, hvy_block, tree_ID )
        else
            call refinementExecute2D_tree( params, hvy_block(:,:,1,:,:), tree_ID )
        end if

        call updateMetadata_tree(params, tree_ID)

        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:, tree_ID), hvy_n(tree_ID) )
    end do

    call balanceLoad_tree( params, hvy_block, tree_ID )
end subroutine
!##############################################################
