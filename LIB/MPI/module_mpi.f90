!>\dir
!> MPI synchronization routines

!> \file 
!> \brief File of MPI Module 


!> \brief This module implements all MPI synchronization routines
!> \details
!>          * synchronization of ghost nodes
!>          * synchronization of light data
!>          * creation of rank block rank lists
!>          * copy, write send and recieve buffers
!> \version 0.4
!> \author msr
!! 12/01/17 - create
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
    include "create_external_com_list.f90"

    ! coyp internal ghost nodes
    include "copy_ghost_nodes_2D.f90"
    include "copy_ghost_nodes_3D.f90"

    ! redundant nodes subroutine - use to ensure correct order of ghost nodes synch
    include "copy_redundant_nodes_2D.f90"

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
    include "synchronize_lgt_data.f90"

    ! routine to reset ghost nodes to uniform number (for debuging)
    include "reset_ghost_nodes.f90"

end module module_MPI
