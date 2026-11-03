module module_saving

    ! this module is to disentangle module_IO from saving (ie preapre
    ! data for storing.) Otherwise, we had circular dependencies.
    use mpi
    use hdf5
    use module_params
    use module_timing
    use module_mpi
    ! use mesh module, since we want to compute origin/spacing of blocks
    use module_mesh
    ! use module operators for computation of the vorticity field
    use module_operators, only: compute_vorticity
    ! use physics modules to save the data
    use module_physics_metamodule
    ! actual storing and reading of HDF5 files:
    use module_helpers, only : check_file_exists, block_contains_NaN
    use module_forestMetaData

    ! HACK to save the beam equation backup when saving HDF5 output :(
    ! We currently have no clean interface for that.
    use module_ACM, only : additionalSaving_ACM

    implicit none

contains

#include "save_data.f90"

end module
