!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_MPI.f90
!> \version 0.4
!> \author msr
!> \brief MPI module
!> \details
!! = log ======================================================================================
!! \n
!! 12/01/17 - create
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

    ! ghost nodes synchronization
    include "synchronize_internal_nodes.f90"
    include "synchronize_external_nodes.f90"

    ! coyp internal ghost nodes
    include "copy_ghost_nodes_2D.f90"
    include "copy_ghost_nodes_3D.f90"

    ! write send buffer
    include "create_send_buffer_2D.f90"
    include "create_send_buffer_3D.f90"

    ! write received buffer
    include "write_receive_buffer_2D.f90"
    include "write_receive_buffer_3D.f90"

    ! maximal number of communications
    include "max_com_num.f90"

    ! fill send buffer
    include "fill_send_buffer.f90"

    ! fill receive buffer
    include "fill_receive_buffer.f90"

    ! soubroutine for get data with lock/unlock synchronization
    include "RMA_lock_unlock_get_data.f90"

    ! soubroutine to put data with lock/unlock synchronization
    include "RMA_lock_unlock_put_data.f90"

    ! non blocking data transfer
    include "isend_irecv_data.f90"

    include "blocks_per_mpirank.f90"

end module module_MPI
