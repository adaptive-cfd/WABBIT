!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name: module_unit_test.f90
!> \version: 0.5
!> \author: msr
!
!> \brief module for all unit test subroutines
!
!> = log ======================================================================================
!! \n
!! 21/03/17 - create
!
! ********************************************************************************************

module module_unit_test

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! init module
    use module_initialization
    ! mesh module
    use module_mesh
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

    ! ghost nodes unit test
    include "unit_test_ghost_nodes_synchronization.f90"

    include "unit_test_wavelet_compression.f90"

    ! time stepper convergence unit test
    include "unit_test_time_stepper_convergence.f90"

    ! spatial convergence order test
    include "unit_test_spatial_convergence_order.f90"

    include "unit_test_treecode.f90"

end module module_unit_test
