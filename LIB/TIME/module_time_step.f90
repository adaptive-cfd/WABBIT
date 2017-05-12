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
    use module_mesh, only : get_block_spacing_origin

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! time step
    include "time_step_RK4.f90"

    ! filter
    include "filter_block.f90"
    include "filter_1D.f90"

    ! dt calculation
    include "calculate_time_step.f90"

end module module_time_step
