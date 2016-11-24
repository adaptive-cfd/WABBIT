! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: lgt_id_to_hvy_id.f90
! version: 0.4
! author: msr
!
! convert light id into heavy id
!
! input:    - light id, proc rank, number of blocks per proc
! output:   - heavy data id
!
! = log ======================================================================================
!
! 23/11/16 - create
! ********************************************************************************************

subroutine lgt_id_to_hvy_id( hvy_id, lgt_id, rank, N )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! heavy id
    integer(kind=ik), intent(out)       :: hvy_id

    ! light id
    integer(kind=ik), intent(in)        :: lgt_id

    ! rank of proc
    integer(kind=ik), intent(in)        :: rank

    ! number of blocks per proc
    integer(kind=ik), intent(in)        :: N

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    hvy_id = lgt_id - rank*N

end subroutine lgt_id_to_hvy_id
