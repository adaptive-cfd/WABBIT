! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_debug.f90
! version: 0.4
! author: msr
!
! module for all debug subroutines
!
! = log ======================================================================================
!
! 29/11/16 - create
!
! TODO: union with debug data structure
! ********************************************************************************************

module module_debug

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! global user defined debug structure
    type type_debug

        ! computing time measurement array
        ! row number: id corresponding to names list
        ! column 1: number of subroutine calls for one time loop
        ! column 2: sum (time) of all subroutine calls for one time loop
        ! column 3: number of subroutine calls over complete program
        ! column 4: sum (time) of all subroutine calls over complete program
        real(kind=rk), dimension(:,:), allocatable          :: comp_time

        ! names of time measurements
        ! row number: id
        ! column: name
        character(len=40), dimension(:), allocatable        :: name_comp_time

    end type type_debug

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (type_debug), save                                 :: debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! lgt_block synchronization
    include "check_lgt_block_synchronization.f90"

    ! check ghost nodes
    include "check_ghost_nodes.f90"

    ! check future mesh gradedness
    include "write_future_mesh_lvl.f90"

    ! write time measurements
    include "write_debug_times.f90"

    ! write block distribution
    include "write_block_distribution.f90"

    ! write communciations list
    include "write_com_list.f90"

    ! com matrix
    include "write_com_matrix.f90"

end module module_debug
