! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_mesh.f90
! version: 0.4
! author: msr
!
! module for all mesh subroutines
!
! = log ======================================================================================
!
! 24/11/16 - create
! ********************************************************************************************

module module_mesh

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

    ! create list of active blocks (light data)
    include "create_lgt_active_list.f90"

    ! create list of active blocks (heavy data)
    include "create_hvy_active_list.f90"

    ! update neighbors
    include "update_neighbors.f90"

end module module_mesh
