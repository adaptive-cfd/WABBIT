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

    ! write send buffer
    include "create_send_buffer.f90"

    ! write received buffer
    include "write_receive_buffer.f90"

end module module_time_step
