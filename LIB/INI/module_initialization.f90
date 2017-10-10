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
!>
!! = log ======================================================================================
!! \n
!! 24/11/16 - create \n
!! 04/04/17 - rename to module_initialization (as specific iniconds are in the module module_initital_conditions.f90)
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
    ! this module contains the routines to set the initial condition on blocks
    use module_initial_conditions
!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! init_data subroutine
    include "set_blocks_initial_condition.f90"
    include "set_inicond_all_blocks.f90"
    include "get_inicond_from_file.f90"

    ! allocate and reset all memotry requred for the gird
    include "allocate_grid.f90"
    include "reset_grid.f90"
    include "create_equidistant_base_mesh.f90"
    include "allocate_com_arrays.f90"


    ! initial block distribution - 2D case
    include "initial_block_distribution_2D.f90"

    ! subroutine to write new heavy block data
    include "new_block_heavy.f90"

    ! vorticity filaments testcase
    ! include "inicond_vorticity_filaments.f90"
    !
    ! ! initial zeros for all fields
    include "inicond_zeros.f90"

    ! initial block distribution - 3D case
    include "initial_block_distribution_3D.f90"

    ! ! start field, 3D sphere
    ! include "inicond_sphere.f90"
    !
    ! ! richtmyer meshkov instability setup
    ! include "inicond_richtmyer_meshkov.f90"
    !
    ! ! shear layer setup
    ! include "inicond_shear_layer.f90"
    !
    ! ! sinus initialization
    ! include "inicond_sinus_2D.f90"
    ! include "inicond_sinus_3D.f90"

end module module_initialization
