!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_time_step.f90
!> \version 0.4
!> \author msr
!
!> \brief time step module
!
!>
!! = log ======================================================================================
!! \n
!! 24/11/16 - create
! ********************************************************************************************

module module_time_step

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! MPI module
    use module_MPI
    ! use mesh module, since we need to compute dx and origin of blocks
    use module_mesh, only : get_block_spacing_origin, max_active_level

    use module_operators, only: volume_integral
    ! to call RHS routines:
    use module_physics_metamodule, only : RHS_meta, STATISTICS_meta, FILTER_meta

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! time step
    include "time_stepper.f90"
    include "set_RK_input.f90"
    include "RHS_wrapper.f90"
    include "filter_wrapper.f90"
    include "statistics_wrapper.f90"
    include "final_stage_RK.f90"

    ! dt calculation
    include "calculate_time_step.f90"

    ! velocity norm calculation
    include "get_block_max_velocity_norm.f90"

end module module_time_step
