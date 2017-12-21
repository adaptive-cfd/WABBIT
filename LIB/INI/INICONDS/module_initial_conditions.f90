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
    use module_ACM_new
    use module_convdiff_new
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
  include "inicond_constant_acm.f90"

  ! 3D sphere initialization
  include "inicond_sphere.f90"

  ! shear layer initialization
  include "inicond_shear_layer.f90"

end module module_initial_conditions
