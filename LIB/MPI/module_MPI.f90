! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_MPI.f90
! version: 0.4
! author: msr
!
! MPI module
!
! = log ======================================================================================
!
! 12/01/17 - create
! ********************************************************************************************

module module_MPI

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

    ! ghost nodes synchronization
    include "synchronize_ghosts.f90"

    ! coyp internal ghost nodes
    include "copy_ghost_nodes.f90"

    ! write send buffer
    include "create_send_buffer.f90"

    ! write received buffer
    include "write_receive_buffer.f90"

end module module_MPI
