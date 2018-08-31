!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_initialization.f90
!> \version 0.4
!> \author msr
!
!> \brief module for all init subroutines
!
!!\details
!! \date 24/11/16 - create
!! \date 04/04/17 - rename to module_initialization (as specific iniconds are in the module module_initital_conditions.f90)
!!            here we provide interfaces for mesh creation and inicond setting
! ********************************************************************************************

module module_initialization

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! mesh module
    use module_mesh
    ! read routines
    use module_IO
    ! to set the initial condition depending on pysics, we have to include them here
    use module_physics_metamodule

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! init_data subroutine
    include "set_initial_grid.f90"
    include "set_inicond_blocks.f90"
    include "get_inicond_from_file.f90"

    ! subroutine to write new heavy block data
    include "new_block_heavy.f90"

end module module_initialization
