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
    use mpi
    ! global parameters
    use module_params
    ! init module
    use module_initialization
    ! mesh module
    use module_mesh
    ! time step module
    use module_time_step

    implicit none


contains

#include "unit_test_ghost_nodes_synchronization.f90"
#include "unit_test_treecode.f90"

end module module_unit_test
