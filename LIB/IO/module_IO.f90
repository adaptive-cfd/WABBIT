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
    ! timing module
    use module_timing
    ! hdf5 file wrapper
    use module_hdf5_wrapper
    ! own MPI module
    use module_mpi
    ! use mesh module, since we want to compute origin/spacing of blocks
    use module_mesh
    ! use module operators for computation of the vorticity field
    use module_operators, only: compute_vorticity
    ! use physics modules to save the data
    use module_physics_metamodule
    ! for reading parameters and storing them in each hdf file
    use module_ini_files_parser_mpi

    use module_helpers, only : check_file_exists, block_contains_NaN

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

    ! read mesh properties and time from input file
    include "read_mesh.f90"

    ! read field from input file
    include "read_field.f90"

    include "read_attributes.f90"

    include "read_file_flusi.f90"

end module module_IO
