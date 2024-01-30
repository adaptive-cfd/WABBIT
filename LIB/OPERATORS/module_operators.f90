!> \brief module for all operator routines
! *********************************************************************************************

module module_operators

use mpi
use module_params     ! global parameters
use module_timing
use module_treelib
use module_forestMetaData

implicit none

PRIVATE
!**********************************************************************************************
! These are the important routines that are visible to WABBIT:
!**********************************************************************************************
PUBLIC :: compute_vorticity, compute_vorticity_abs, divergence, gradient, compute_Qcriterion, componentWiseNorm_tree

contains

#include "compute_Qcriterion.f90"
#include "compute_vorticity.f90"
#include "divergence.f90"
#include "gradient.f90"
#include "componentWiseNorm_tree.f90"


end module module_operators
