!> \brief module for all unit test subroutines
! ********************************************************************************************

module module_unit_test
    use mpi
    use module_params           ! global parameters
    use module_initialization   ! init module
    use module_mesh
    use module_time_step
    use module_treelib          ! for treecode test

    implicit none

contains

#include "unit_test_treecode.f90"
#include "unit_test_Sync.f90"
#include "unit_test_ghostSync.f90"
#include "unit_test_waveletDecomposition.f90"
#include "unit_test_waveletDecomposition_invertibility.f90"
#include "unit_test_refineCoarsen.f90"
#include "createTestGrids.f90"

end module module_unit_test
