! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_debug_functions.f90
! version: 0.4
! author: msr
!
! module for all debug subroutines
!
! = log ======================================================================================
!
! 29/11/16 - create
!
! TODO: union with debug data structure
! ********************************************************************************************

module module_debug_functions

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! lgt_block synchronization
    include "check_lgt_block_synchronization.f90"

    ! check ghost nodes
    include "check_ghost_nodes.f90"

    ! check future mesh gradedness
    include "write_future_mesh_lvl.f90"

end module module_debug_functions
