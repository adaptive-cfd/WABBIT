!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name proc_to_lgt_data_start_id.f90
!> \version 0.4
!> \author msr
!
!> \brief return start index on light data corresponding to proc rank
!
!>
!! input:    - rank, number of blocks \n
!! output:   - start of light data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 23/11/16 - create
! ********************************************************************************************

subroutine proc_to_lgt_data_start_id( lgt_start, rank, N )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> light data start index
    integer(kind=ik), intent(out)       :: lgt_start

    !> rank of proc
    integer(kind=ik), intent(in)        :: rank

    !> number of blocks per proc
    integer(kind=ik), intent(in)        :: N

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    lgt_start = rank*N + 1

end subroutine proc_to_lgt_data_start_id
