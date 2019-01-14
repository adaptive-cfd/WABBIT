!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_indicators.f90
!> \version 0.5
!> \author engels
!
!> \brief This module contains routines for refinement and coarsening indicators, i.e. routines
!! to tag blocks that require refinement or can possibly be coarsened
!
!> \details Refinement/coarsening indicators are expected to grow to larger numbers in the foreseeable
!! future, which is why they are outsourced from the mesh module to a mmodule on their own.
!! = log ======================================================================================
!! \n
!! 23/05/2017 create
! ********************************************************************************************

module module_indicators

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! timing module
    use module_timing
    ! interpolation routines
    use module_interpolation
    ! we now have an indicator which computes the vorticity, so include the operator module
    use module_operators

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! interpolation fields
    real(kind=rk), SAVE, PRIVATE, allocatable :: u1(:,:,:), u2(:,:,:), u3(:,:,:)

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

  include "refinement_indicator.f90"
  include "block_coarsening_indicator.f90"
  ! threshold the blocks
  include "threshold_block.f90"
end module module_indicators
