module module_IO
    use mpi
    use hdf5
    use module_params
    use module_timing
    use module_hdf5_wrapper
    use module_mpi
    ! use mesh module, since we want to compute origin/spacing of blocks
    use module_mesh
    ! use module operators for computation of the vorticity field
    use module_operators, only: compute_vorticity
    ! use physics modules to save the data
    use module_physics_metamodule
    ! for reading parameters and storing them in each hdf file
    use module_ini_files_parser_mpi
    use module_helpers, only : check_file_exists, block_contains_NaN
    use module_forest
    use module_treelib

    implicit none

contains

#include "write_field.f90"
#include "read_mesh.f90"
#include "read_field.f90"
#include "read_attributes.f90"
#include "read_file_flusi.f90"
#include "forest_IO.f90"

end module module_IO
