!> \brief module for all mesh subroutines
! ********************************************************************************************

module module_mesh

    use mpi
    use module_forestMetaData
    use module_params               ! global parameters
    use module_timing               ! debug module
    use module_interpolation        ! interpolation routines
    ! use MPI module, since threshold_block needs to synch ghosts
    use module_MPI
    use module_treelib              ! module with evrything related to treecodes (encoding, decoding, neighbors, etc)
    use module_operators, only: component_wise_tree_norm
    ! used in executeCoarsening_tree
    use module_helpers, only: most_common_element
    ! if the threshold_mask option is used, then the mesh module needs to create the mask function here
    ! hence we require the metamodule to be used.
    use module_physics_metamodule

    implicit none

    interface set_desired_num_blocks_per_rank
        module procedure set_desired_num_blocks_per_rank1, &
                         set_desired_num_blocks_per_rank2
    end interface


contains

    ! create all active (lgt/hvy) lists, create also sorted lgt data list
#include "create_active_and_sorted_lists.f90"
#include "createMask_tree.f90"

#include "block_xfer_nonblocking.f90"

#include "updateNeighbors_tree.f90"

#include "find_neighbors.f90"
#include "doesBlockExist_tree.f90"

    ! block refinement subroutine
#include "refine_tree.f90"

    ! check treelevel restrictions
#include "respectJmaxJmin_tree.f90"

    ! check treelevel restrictions
#include "refinementExecute2D_tree.f90"
#include "refinementExecute3D_tree.f90"

    ! adapt the mesh
#include "adapt_tree.f90"
#include "coarseningIndicator_tree.f90"

    ! gradedness check
#include "ensureGradedness_tree.f90"

    ! completeness check
#include "ensure_completeness.f90"

    ! coarse mesh
#include "executeCoarsening_tree.f90"
#include "merge_blocks.f90"

    ! balance the load
#include "balanceLoad_tree.f90"

    ! create list with number of blocks per rank
#include "set_desired_num_blocks_per_rank.f90"

    ! treecode to 2D z-sfc position
#include "treecode_to_sfc_id_2D.f90"

    ! treecode to 3D z-sfc position
#include "treecode_to_sfc_id_3D.f90"

    ! transfer treecode to 2D hilbert code
#include "treecode_to_hilbertcode_2D.f90"

    ! transfer treecode to 3D hilbert code
#include "treecode_to_hilbertcode_3D.f90"

    ! goes back from a treecode to xyz cartesian coordinates
#include "get_block_spacing_origin.f90"

    ! find sisters to a given block
#include "findSisters_tree.f90"

    ! find globally coarsest / finest levels
#include "ActiveLevel_tree.f90"
#include "get_free_local_light_id.f90"
#include "quicksort.f90"
#include "updateMetadata_tree.f90"

    ! routines for creation of an initial grid
#include "createEquidistantGrid_tree.f90"
#include "createRandomGrid_tree.f90"

    ! allocate and reset all memory required for the gird
#include "reset_tree.f90"
#include "allocate_forest.f90"

#include "write_block_distribution.f90"

    ! lgt_block synchronization
#include "check_lgt_block_synchronization.f90"
#include "remove_nonperiodic_neighbors.f90"
#include "forest.f90"


! check if we still have enough memory left: for very large simulations
! we cannot affort to have them fail without the possibility to resume them
logical function notEnoughMemoryToRefineEverywhere_tree(params, tree_ID)
    implicit none
    type (type_params), intent(inout) :: params
    integer(kind=ik), intent(in)      :: tree_ID

    integer :: lgt_n_max, k, lgt_n_afterRefinement

    notEnoughMemoryToRefineEverywhere_tree = .false.
    params%out_of_memory = .false.

    ! without adaptivity, this routine makes no sense, as the memory is constant
    ! the run either crashes right away or never
    if (params%adapt_tree .eqv. .false.) then
        notEnoughMemoryToRefineEverywhere_tree = .false.
        return
    endif

    ! this is the available maximum number of active blocks.
    lgt_n_max = params%number_blocks*params%number_procs

    ! remove blocks already used for mask etc
    lgt_n_max = lgt_n_max - sum(lgt_n(2:size(lgt_n)))

    ! safety margin. inhomogenoeus distribution
    ! of blocks can make trouble. Therefore, we raise the alarm earlier.
    lgt_n_max = lgt_n_max * 9 / 10


    ! however, this routine is called BEFORE refinement. figure out how many blocks
    ! that we be after refinement
    if (params%force_maxlevel_dealiasing) then
        ! with dealiasing, the relation is very simple.
        lgt_n_afterRefinement = lgt_n(tree_ID) * (2**params%dim)
    else
        lgt_n_afterRefinement = 0
        do k = 1, lgt_n(tree_ID)
            if (lgt_block( lgt_active(k, tree_ID), params%max_treelevel + IDX_MESH_LVL ) < params%max_treelevel) then
                ! this block can be refined (and will be) (it is refined to it has 8 children
                ! but disappears, so 7 new blocks)
                lgt_n_afterRefinement = lgt_n_afterRefinement + (2**params%dim)-1
            else
                ! this block is at maximum refinement and will not be refined, but not be deleted neither
                lgt_n_afterRefinement = lgt_n_afterRefinement + 1
            endif
        enddo

    endif


    if ( lgt_n_afterRefinement > lgt_n_max) then
        ! oh-oh.
        notEnoughMemoryToRefineEverywhere_tree = .true.

        if (params%rank==0) then
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("Nblocks_total=",i7," Nblocks_fluid=",i7)') sum(lgt_n), lgt_n(tree_ID)
            write(*,'("Nblocks_allocated=",i7," Nblocks_fluid_available=",i7)') params%number_blocks*params%number_procs, lgt_n_max
            write(*,'("Nblocks_after_refinement=",i7)') lgt_n_afterRefinement
            write(*,'("Nblocks_after_refinement=",i7)') lgt_n_afterRefinement
            write(*,'("-------------_> oh-oh.")')
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("-----------OUT OF MEMORY--------------")')
        endif

        params%out_of_memory = .true.
    endif


end function


!##############################################################
! This routine refines/coarsens to given target_level.
! If target_level is not passed, then max_treelevel is assumed
subroutine toEquidistant_tree(params, hvy_block, hvy_tmp, tree_ID, target_level, verbosity)

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
        level = params%max_treelevel
    endif

    ! refine/coarse to attain desired level, respectively
    ! coarsen
    do while (maxActiveLevel_tree( tree_ID )>level)
        ! check where coarsening is actually needed and set refinement status to -1 (coarsen)
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)

            if (treecode_size(lgt_block(lgt_id,:), params%max_treelevel) > level) then
                lgt_block(lgt_id, params%max_treelevel + IDX_REFINE_STS) = -1
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

            if (treecode_size(lgt_block(lgt_id,:), params%max_treelevel) < level) then
                lgt_block(lgt_id, params%max_treelevel + IDX_REFINE_STS) = 1
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


end module module_mesh
