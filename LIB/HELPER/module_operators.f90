!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_operators.f90
!> \version 0.5
!> \author sm
!
!> \brief module for all operator routines
!
!>
!! = log ======================================================================================
!! \n
!! 28/7/17 - create
! *********************************************************************************************

module module_operators

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! own MPI module
    !use module_mpi
    ! use mesh module, since we want to compute origin/spacing of blocks
    use module_mesh, only : get_block_spacing_origin

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    include "volume_integral.f90"
    include "compute_vorticity.f90"
    include "compute_forcing.f90"

end module module_operators
