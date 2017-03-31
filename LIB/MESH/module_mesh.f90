! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_mesh.f90
! version: 0.4
! author: msr
!
! module for all mesh subroutines
!
! = log ======================================================================================
!
! 24/11/16 - create
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
    ! time step module
    use module_time_step

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! create list of active blocks (light data)
    include "create_lgt_active_list.f90"

    ! create list of active blocks (heavy data)
    include "create_hvy_active_list.f90"

    ! update neighbors, 2D/3D
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
    include "refine_everywhere.f90"

    ! check treelevel restrictions
    include "respect_min_max_treelevel.f90"

    ! check treelevel restrictions
    include "refine_mesh_2D.f90"
    include "refine_mesh_3D.f90"

    ! adapt the mesh
    include "adapt_mesh.f90"

    ! threshold the blocks
    include "threshold_block.f90"

    ! gradedness check
    include "ensure_gradedness.f90"

    ! completeness check
    include "ensure_completeness_2D.f90"
    include "ensure_completeness_3D.f90"

    ! coarse mesh
    include "coarse_mesh_2D.f90"
    include "coarse_mesh_3D.f90"

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
    include "treecode_to_sfc_id.f90"

    ! transfer treecode to 2D hilbert code
    include "treecode_to_hilbercode.f90"

    ! transfer treecode to 3D z-curve code
    include "treecode_to_3D_z_curve.f90"

    ! goes back from a treecode to xyz cartesian coordinates
    include "decoding.f90"
    include "get_block_spacing_origin.f90"

end module module_mesh
