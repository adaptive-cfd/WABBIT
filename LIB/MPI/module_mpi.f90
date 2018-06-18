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
!> \version 0.6
!> \author engels, reiss, msr
!! 12/01/17 - create
!! 18/07/2018 - remove all old ghost node routines, build in JR's improved version of MSR's new ghost nodes
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

    include "synchronize_ghosts.f90"
    include "blocks_per_mpirank.f90"
    include "synchronize_lgt_data.f90"
    include "reset_ghost_nodes.f90"

end module module_MPI
