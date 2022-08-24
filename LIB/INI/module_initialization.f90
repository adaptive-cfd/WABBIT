module module_initialization
    use mpi
    use module_params
    use module_timing
    use module_mesh
    ! read routines
    use module_IO
    ! to set the initial condition depending on pysics, we have to include them here
    use module_physics_metamodule

    use module_treelib
    use module_forestMetaData

    implicit none

contains

#include "setInitialCondition_tree.f90"
#include "setInicondBlocks_tree.f90"
end module module_initialization
