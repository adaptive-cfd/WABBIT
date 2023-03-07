!> \brief module for all mesh subroutines
! ********************************************************************************************

module module_mesh

    use mpi
    use module_hdf5_wrapper
    use module_forestMetaData
    use module_params               ! global parameters
    use module_timing               ! debug module
    use module_interpolation        ! interpolation routines
    ! use MPI module, since threshold_block needs to synch ghosts
    use module_MPI
    use module_treelib              ! module with evrything related to treecodes (encoding, decoding, neighbors, etc)
    use module_operators
    ! used in executeCoarsening_tree
    use module_helpers
    ! if the threshold_mask option is used, then the mesh module needs to create the mask function here
    ! hence we require the metamodule to be used.
    use module_physics_metamodule

    implicit none

    interface set_desired_num_blocks_per_rank
        module procedure set_desired_num_blocks_per_rank1, &
                         set_desired_num_blocks_per_rank2
    end interface


contains


#include "unitTest_ghostSync.f90"
#include "unitTest_waveletDecomposition.f90"
#include "waveletDecomposition_tree.f90"
#include "refineToEquidistant_tree.f90"
#include "InputOutput_Flusi.f90"
#include "InputOutput.f90"
#include "create_active_and_sorted_lists.f90"
#include "createMask_tree.f90"
#include "block_xfer_nonblocking.f90"
#include "updateNeighbors_tree.f90"
#include "find_neighbors.f90"
#include "doesBlockExist_tree.f90"
#include "refine_tree.f90"
#include "respectJmaxJmin_tree.f90"
#include "refinementExecute2D_tree.f90"
#include "refinementExecute3D_tree.f90"
#include "adapt_tree.f90"
#include "coarseningIndicator_tree.f90"
#include "ensureGradedness_tree.f90"
#include "ensure_completeness.f90"
#include "executeCoarsening_tree.f90"
#include "merge_blocks.f90"
#include "balanceLoad_tree.f90"
#include "set_desired_num_blocks_per_rank.f90"
#include "treecode_to_sfc_id_2D.f90"
#include "treecode_to_sfc_id_3D.f90"
#include "treecode_to_hilbertcode_2D.f90"
#include "treecode_to_hilbertcode_3D.f90"
#include "get_block_spacing_origin.f90"
#include "findSisters_tree.f90"
#include "ActiveLevel_tree.f90"
#include "get_free_local_light_id.f90"
#include "quicksort.f90"
#include "updateMetadata_tree.f90"
#include "createEquidistantGrid_tree.f90"
#include "createRandomGrid_tree.f90"
#include "reset_tree.f90"
#include "allocate_forest.f90"
#include "write_block_distribution.f90"
#include "check_lgt_block_synchronization.f90"
#include "remove_nonperiodic_neighbors.f90"
#include "forest.f90"
#include "notEnoughMemoryToRefineEverywhere_tree.f90"

end module module_mesh
