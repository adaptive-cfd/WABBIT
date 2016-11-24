! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: lgt_id_to_proc_rank.f90
! version: 0.4
! author: msr
!
! return proc rank corresponding to given light id
!
! input:    - light id, number of blocks
! output:   - rank
!
! = log ======================================================================================
!
! 23/11/16 - create
! ********************************************************************************************

subroutine lgt_id_to_proc_rank( rank, lgt_id, N )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data start index
    integer(kind=ik), intent(out)       :: rank

    ! rank of proc
    integer(kind=ik), intent(in)        :: lgt_id

    ! number of blocks per proc
    integer(kind=ik), intent(in)        :: N

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    rank = (lgt_id-1) / N

end subroutine lgt_id_to_proc_rank
