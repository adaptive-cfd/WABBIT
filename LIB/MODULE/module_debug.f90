! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_debug.f90
! version: 0.4
! author: msr
!
! debug data structure
!
!
! = log ======================================================================================
!
! 24/11/16 - create
! ********************************************************************************************

module module_debug

!---------------------------------------------------------------------------------------------
! modules

    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! global user defined debug structure
    type type_debug

        ! computing time measurement array
        ! dim 1: proc rank
        ! dim 2: number of measured subroutines
        real(kind=rk), dimension(:,:), allocatable      :: comp_time

        ! names of time measurements
        character(len=40), dimension(:), allocatable    :: name_comp_time

    end type type_debug

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (type_debug), save                             :: debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

end module module_debug
