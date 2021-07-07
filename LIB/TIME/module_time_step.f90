module module_time_step
    use mpi
    use module_params     ! global parameters
    use module_timing
    use module_MPI
    ! use mesh module, since we need to compute dx and origin of blocks
    use module_mesh !, only : get_block_spacing_origin, max_active_level
    ! to call RHS routines:
    use module_physics_metamodule, only : RHS_meta, STATISTICS_meta, FILTER_meta
    use module_mask
    use module_treelib

    implicit none

    ! coefficient tables for RKC scheme
    real(kind=rk), dimension(1:40,1:40), private, save :: mu=0.0_rk, mu_tilde=0.0_rk, nu=0.0_rk, gamma_tilde=0.0_rk, c=0.0_rk

contains

#include "time_stepper.f90"
#include "RHS_wrapper.f90"
#include "filter_wrapper.f90"
#include "statistics_wrapper.f90"
#include "krylov.f90"
#include "calculate_time_step.f90"
#include "runge_kutta_chebychev.f90"
#include "runge_kutta_chebychev_FSI.f90"
#include "runge_kutta_generic.f90"
#include "runge_kutta_generic_FSI.f90"

end module module_time_step
