!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name lgt_id_to_hvy_id.f90
!> \version 0.4
!> \author msr
!
!> \brief convert light id into heavy id
!
!> \details
!! input:    - light id, proc rank, number of blocks per proc \n
!! output:   - heavy data id \n
!!
!!
!! = log ======================================================================================
!! \n
!! 23/11/16 - create
! ********************************************************************************************

subroutine lgt_id_to_hvy_id( hvy_id, lgt_id, rank, N )
    ! global parameters
    use module_params

    implicit none

    !> heavy id
    integer(kind=ik), intent(out)       :: hvy_id

    !> light id
    integer(kind=ik), intent(in)        :: lgt_id

    !> rank of proc
    integer(kind=ik), intent(in)        :: rank

    !> number of blocks per proc
    integer(kind=ik), intent(in)        :: N

    hvy_id = lgt_id - rank*N

end subroutine
