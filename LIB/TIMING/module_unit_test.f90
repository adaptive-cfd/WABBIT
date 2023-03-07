!> \brief module for all unit test subroutines
! ********************************************************************************************

module module_unit_test
    use mpi
    use module_params           ! global parameters
    use module_initialization   ! init module
    use module_mesh
    use module_time_step

    implicit none

contains

#include "unit_test_treecode.f90"

end module module_unit_test
