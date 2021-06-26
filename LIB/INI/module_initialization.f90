module module_initialization
    use mpi
    use module_params
    use module_timing
    use module_mesh
    ! read routines
    use module_IO
    ! to set the initial condition depending on pysics, we have to include them here
    use module_physics_metamodule

    use module_mask
    use module_treelib

    implicit none

contains

    ! init_data subroutine
#include "set_initial_grid.f90"
#include "set_inicond_blocks.f90"
#include "get_inicond_from_file.f90"
end module module_initialization
