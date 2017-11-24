!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_IO.f90
!> \version 0.4
!> \author msr
!
!> \brief module for all IO routines
!
!>
!! = log ======================================================================================
!! \n
!! 24/11/16 - create
! ********************************************************************************************

module module_IO

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    use hdf5
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! hdf5 file wrapper
    use module_hdf5_wrapper
    ! own MPI module
    use module_mpi
    ! use mesh module, since we want to compute origin/spacing of blocks
    use module_mesh
    ! use module operators for computation of the vorticity field
    use module_operators, only: compute_vorticity
    
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

    ! write fields to disk
    include "write_field.f90"

    ! compute vorticity to save it on disk
    include "write_vorticity.f90"

    ! write mask term to disk
    include "write_mask.f90"

    ! read mesh properties and time from input file
    include "read_mesh_and_attributes.f90"

    ! read field from input file
    include "read_field.f90"

    ! check if input file exists
    include "check_file_exists.f90"

end module module_IO
