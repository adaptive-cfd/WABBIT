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
#include "unitTest_Sync.f90"
#include "unitTest_ghostSync.f90"
#include "unitTest_waveletDecomposition.f90"
#include "unitTest_refineCoarsen.f90"

end module module_unit_test
