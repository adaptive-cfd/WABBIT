!> \brief return proc rank corresponding to given light id
!! input:    - light id, number of blocks \n
!! output:   - rank \n
! ********************************************************************************************

subroutine lgt2proc( rank, lgt_id, N )
    ! global parameters
    use module_params

    implicit none

    integer(kind=ik), intent(out)       :: rank       !> rank of proc
    integer(kind=ik), intent(in)        :: lgt_id     !> light data start index
    integer(kind=ik), intent(in)        :: N          !> number of blocks per proc

    rank = (lgt_id-1) / N

end subroutine
