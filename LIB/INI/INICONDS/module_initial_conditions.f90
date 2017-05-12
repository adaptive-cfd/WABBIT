!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_initialization.f90
!> \version 0.4
!> \author engels
!
!> \brief module for all init subroutines
!
!>
!! = log ======================================================================================
!! \n
!! 03 Apr 2017 - create \n
!! 12/04/17 - add 3D sphere
!!
! ********************************************************************************************

module module_initial_conditions

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

  include "initial_condition_on_block_wrapper.f90"
  include "inicond_gauss_blob.f90"
  include "inicond_sinus_2D.f90"

  ! 3D sphere initialization
  include "inicond_sphere.f90"

end module module_initial_conditions
