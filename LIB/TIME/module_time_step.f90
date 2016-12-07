! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_time_step.f90
! version: 0.4
! author: msr
!
! time step module
!
! = log ======================================================================================
!
! 24/11/16 - create
! ********************************************************************************************

module module_time_step

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! interpolation routines
    use module_interpolation
    ! navier stokes rhs test module
    use module_rhs_navier_stokes

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

    ! ghost nodes synchronization
    include "synchronize_ghosts.f90"

    ! coyp internal ghost nodes
    include "copy_ghost_nodes.f90"

    ! send/receive external ghost nodes
    include "send_receive_data.f90"

    ! test rk for navier stokes rhs
    include "time_step_RK4_2.f90"

end module module_time_step
