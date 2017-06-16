!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_mesh.f90
!> \version 0.4
!> \author msr
!
!> \brief module for all mesh subroutines
!
!> \details
!! = log ======================================================================================
!! \n
!! 24/11/16 - create
! ********************************************************************************************

module module_mesh

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! interpolation routines
    use module_interpolation
    ! use MPI module, since thrshold_block needs to synch ghosts
    use module_MPI
    ! module with evrything related to treecodes (encoding, decoding, neighbors, etc)
    use module_treelib
!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! create all active (lgt/hvy) lists, create also sorted lgt data list
    include "create_active_and_sorted_lists.f90"

    ! update neighbors, 2D/3D
    include "update_neighbors.f90"
    include "update_neighbors_2D.f90"
    include "update_neighbors_3D.f90"

    ! neighbor search, edge
    include "find_neighbor_edge_2D.f90"
    include "find_neighbor_edge_3D.f90"

    ! neighbor search, corner
    include "find_neighbor_corner_2D.f90"
    include "find_neighbor_corner_3D.f90"

    ! neighbor search, face
    include "find_neighbor_face_3D.f90"

    ! search routine to find light data id corresponding to treecode
    include "does_block_exist.f90"

    ! block refinement subroutine
    include "refine_mesh.f90"

    ! check treelevel restrictions
    include "respect_min_max_treelevel.f90"

    ! check treelevel restrictions
    include "refinement_execute_2D.f90"
    include "refinement_execute_3D.f90"

    ! adapt the mesh
    include "adapt_mesh.f90"

    ! gradedness check
    include "ensure_gradedness.f90"

    ! completeness check
    include "ensure_completeness.f90"

    ! coarse mesh
    include "coarse_mesh.f90"
    include "merge_blocks.f90"

    ! balance the load
    include "balance_load_2D.f90"
    include "balance_load_3D.f90"

    ! create list with number of blocks per rank
    include "set_desired_num_blocks_per_rank.f90"

    ! create friends table
    include "compute_friends_table.f90"

    ! affinity list
    include "compute_affinity.f90"

    ! treecode to 2D z-sfc position
    include "treecode_to_sfc_id_2D.f90"

    ! treecode to 3D z-sfc position
    include "treecode_to_sfc_id_3D.f90"

    ! transfer treecode to 2D hilbert code
    include "treecode_to_hilbertcode_2D.f90"

    ! transfer treecode to 3D hilbert code
    include "treecode_to_hilbertcode_3D.f90"

    ! goes back from a treecode to xyz cartesian coordinates
    include "get_block_spacing_origin.f90"

    ! find sisters to a given block
    include "find_sisters.f90"

    ! find globally coarsest / finest levels
    include "max_active_level.f90"
    include "min_active_level.f90"

    !
    include "get_free_local_light_id.f90"
    include "gather_blocks_on_proc.f90"
    include "quicksort.f90"

end module module_mesh
