!> \file
!
!> \brief module for all mesh subroutines
!
!> \details
!> \author msr
!! \date 24/11/16 - create
! ********************************************************************************************

module module_mesh

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_timing
    ! interpolation routines
    use module_interpolation
    ! use MPI module, since thrshold_block needs to synch ghosts
    use module_MPI
    ! module with evrything related to treecodes (encoding, decoding, neighbors, etc)
    use module_treelib
    !
    use module_boundary_conditions
    use module_operators, only: component_wise_tree_norm
    ! used in coarse_mesh
    use module_helpers, only: most_common_element
    ! if the threshold_mask option is used, then the mesh module needs to create the mask function here
    ! hence we require the metamodule to be used.
    use module_physics_metamodule

    implicit none

    ! interface generalising the mono-tree to multi-trees
    interface create_active_and_sorted_lists
        module procedure create_active_and_sorted_lists_tree_old, &
                         create_active_and_sorted_lists_forest
    end interface

    interface set_desired_num_blocks_per_rank
        module procedure set_desired_num_blocks_per_rank1, &
                         set_desired_num_blocks_per_rank2
    end interface

    interface allocate_grid
        module procedure allocate_forest, allocate_tree
    end interface

    interface deallocate_grid
        module procedure deallocate_forest, deallocate_tree
    end interface

contains

    ! create all active (lgt/hvy) lists, create also sorted lgt data list
#include "create_active_and_sorted_lists.f90"

#include "block_xfer_nonblocking.f90"

    ! update neighbors, 2D/3D
#include "update_neighbors.f90"
#include "update_neighbors_2D.f90"
#include "update_neighbors_3D.f90"
    ! neighbor search, edge
#include "find_neighbor_edge_2D.f90"
#include "find_neighbor_edge_3D.f90"

    ! neighbor search, corner
#include "find_neighbor_corner_2D.f90"
#include "find_neighbor_corner_3D.f90"
    ! neighbor search, face
#include "find_neighbor_face_3D.f90"

    ! search routine to find light data id corresponding to treecode
#include "does_block_exist.f90"

    ! block refinement subroutine
#include "refine_mesh.f90"

    ! check treelevel restrictions
#include "respect_min_max_treelevel.f90"

    ! check treelevel restrictions
#include "refinement_execute_2D.f90"
#include "refinement_execute_3D.f90"

    ! adapt the mesh
#include "adapt_mesh.f90"
#include "grid_coarsening_indicator.f90"

    ! gradedness check
#include "ensure_gradedness.f90"

    ! completeness check
#include "ensure_completeness.f90"

    ! coarse mesh
#include "coarse_mesh.f90"
#include "merge_blocks.f90"

    ! balance the load
#include "balance_load.f90"

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
#include "find_sisters.f90"

    ! find globally coarsest / finest levels
#include "max_active_level.f90"
#include "min_active_level.f90"
#include "get_free_local_light_id.f90"
#include "quicksort.f90"
#include "update_grid_metadata.f90"

    ! routines for creation of an initial grid
#include "create_equidistant_grid.f90"
#include "create_random_grid.f90"

    ! allocate and reset all memory required for the gird
#include "reset_grid.f90"
#include "allocate_grid.f90"

#include "write_block_distribution.f90"

    ! lgt_block synchronization
#include "check_lgt_block_synchronization.f90"


! check if we still have enough memory left: for very large simulations
! we cannot affort to have them fail without the possibility to resume them
logical function not_enough_memory(params, lgt_block, lgt_active, lgt_n)
    implicit none
    integer, intent(in) :: lgt_n(:)
    integer(kind=ik), intent(in) :: lgt_block(:, :)
    integer(kind=ik), intent(in) :: lgt_active(:,:)
    type (type_params), intent(inout) :: params

    integer :: lgt_n_max, k, lgt_n_afterRefinement

    not_enough_memory = .false.
    params%out_of_memory = .false.

    ! without adaptivity, this routine makes no sense, as the memory is constant
    ! the run either crashes right away or never
    if (params%adapt_mesh .eqv. .false.) then
        not_enough_memory = .false.
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
        lgt_n_afterRefinement = lgt_n(tree_ID_flow) * (2**params%dim)
    else
        lgt_n_afterRefinement = 0
        do k = 1, lgt_n(tree_ID_flow)
            if (lgt_block( lgt_active(k, tree_ID_flow), params%max_treelevel + IDX_MESH_LVL ) < params%max_treelevel) then
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
        not_enough_memory = .true.

        if (params%rank==0) then
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("Nblocks_total=",i7," Nblocks_fluid=",i7)') sum(lgt_n), lgt_n(tree_ID_flow)
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
subroutine to_dense_mesh(params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, target_level, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
    integer(kind=ik), intent(inout)   :: lgt_n !< number of light active blocks
    integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
    integer(kind=ik), intent(inout)   :: lgt_active(:), hvy_active(:) !< active lists
    integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    integer(kind=ik), intent(in), optional :: target_level
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    integer(kind=ik):: level, tree_id,k,treecode_size
    logical :: verbose=.false.

    if (present(verbosity)) verbose=verbosity
    if (present(target_level)) then
        level=target_level
    else
        level = params%max_treelevel
    endif

    ! get tree_id form lgt data
    tree_id = lgt_block(lgt_active(1),params%max_treelevel+IDX_TREE_ID)

    ! refine/coarse to attain desired level, respectively
    !coarsen
    do while (max_active_level( lgt_block, lgt_active, lgt_n )>level)
        ! check where coarsening is actually needed and set refinement status to -1 (coarsen)
        do k = 1, lgt_n
            if (treecode_size(lgt_block(lgt_active(k),:), params%max_treelevel) > level)&
                lgt_block(lgt_active(k), params%max_treelevel + IDX_REFINE_STS) = -1
        end do
        ! this might not be necessary since we start from an admissible grid
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n )

        call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_active, hvy_n, tree_ID=1)

        call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, tree_id)
    end do

    ! refine
    do while (min_active_level( lgt_block, lgt_active, lgt_n )<level)
        ! check where refinement is actually needed
        do k = 1, lgt_n
            if (treecode_size(lgt_block(lgt_active(k),:), params%max_treelevel) < level)&
                lgt_block(lgt_active(k), params%max_treelevel + IDX_REFINE_STS) = 1
        end do
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n )
        if ( params%dim == 3 ) then
            ! 3D:
            call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
        else
            ! 2D:
            call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:),&
                hvy_active, hvy_n )
        end if

        call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, tree_id)

        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    end do

    call balance_load( params, lgt_block, hvy_block, &
        hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_id )
end subroutine
!##############################################################










end module module_mesh
