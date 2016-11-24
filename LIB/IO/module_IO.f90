! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_IO.f90
! version: 0.4
! author: msr
!
! module for all IO routines
!
! = log ======================================================================================
!
! 24/11/16 - create
! ********************************************************************************************

module module_IO

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! hdf5 file wrapper
    use module_hdf5_wrapper

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! create list of active blocks (light data)
    include "save_data.f90"

    ! create list of active blocks (heavy data)
    include "write_field.f90"

end module module_IO
