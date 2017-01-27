! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_init.f90
! version: 0.4
! author: msr
!
! module for all init subroutines
!
! = log ======================================================================================
!
! 24/11/16 - create
! ********************************************************************************************

module module_init

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! ini file parser module
    use module_ini_files_parser
    ! debug module
    use module_debug

!---------------------------------------------------------------------------------------------
! variables

    implicit none

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! init_data subroutine
    include "init_data.f90"

    ! block list allocation
    include "allocate_block_list.f90"

    ! block data allocation
    include "allocate_block_data.f90"

    ! block work data allocation
    include "allocate_work_data.f90"

    ! start field, gauss blob
    include "inicond_gauss_blob.f90"

    ! initial block distribution - 2D case
    include "initial_block_distribution_2D.f90"

    ! subroutine to write new heavy block data
    include "new_block_heavy.f90"

    ! vorticity filaments testcase
    include "inicond_vorticity_filaments.f90"

    ! read ini file
    include "ini_file_to_params.f90"

    ! initial zeros for all fields
    include "inicond_zeros.f90"

    ! initial block distribution - 3D case
    include "initial_block_distribution_3D.f90"

end module module_init
