!> \brief This module contains routines for refinement and coarsening indicators, i.e. routines
!! to tag blocks that require refinement or can possibly be coarsened
!> \details Refinement/coarsening indicators are expected to grow to larger numbers in the foreseeable
!! future, which is why they are outsourced from the mesh module to a mmodule on their own.
! ********************************************************************************************

module module_indicators
    use mpi
    use module_params
    use module_timing
    use module_interpolation
    ! we now have an indicator which computes the vorticity, so include the operator module
    use module_operators
    use module_forestMetaData

    ! some operators may depend on the actual data (that is, heavy data), for example
    ! for shock or mask detection. These criteria are computed mpi-locally (because of course
    ! each block is associated to one CPU), afterwards, a lgt_data synchronization step is
    ! required, therefore use module_MPI.
    use module_MPI

    implicit none


contains

#include "refinementIndicator_tree.f90"
#include "coarseningIndicator_block.f90"
#include "threshold_block.f90"

end module module_indicators
