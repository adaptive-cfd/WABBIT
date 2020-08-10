module module_saving

    ! this module is to disentangle module_IO from saving (ie preapre
    ! data for storing.) Otherwise, we had circular dependencies.
    use mpi
    use hdf5
    ! global parameters
    use module_params
    ! timing module
    use module_timing
    ! own MPI module
    use module_mpi
    ! use mesh module, since we want to compute origin/spacing of blocks
    use module_mesh
    ! use module operators for computation of the vorticity field
    use module_operators, only: compute_vorticity
    ! use physics modules to save the data
    use module_physics_metamodule
    ! actual storing and reading of HDF5 files:
    use module_IO

    use module_helpers, only : check_file_exists, block_contains_NaN
    use module_mask


    implicit none


contains


#include "save_data.f90"

end module
