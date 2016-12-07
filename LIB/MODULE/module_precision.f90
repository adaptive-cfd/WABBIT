! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_precision.f90
! version: 0.4
! author: msr
!
! module for data precision
!
! = log ======================================================================================
!
! 06/12/16 - create
! ********************************************************************************************

module module_precision

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! define data precision parameters
    integer, parameter, public   :: sngl_prec=selected_real_kind(4)
    integer, parameter, public   :: dble_prec=selected_real_kind(8)

    integer, parameter, public   :: int_prec=selected_int_kind(8)

    integer, parameter, public   :: rk=dble_prec
    integer, parameter, public   :: ik=int_prec

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

end module module_precision
