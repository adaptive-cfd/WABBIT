!! input:    - light id, proc rank, number of blocks per proc \n
!! output:   - heavy data id \n
! ********************************************************************************************

subroutine lgt2hvy( hvy_id, lgt_id, rank, N )
    ! global parameters
    use module_params

    implicit none
    integer(kind=ik), intent(out)       :: hvy_id
    integer(kind=ik), intent(in)        :: lgt_id
    integer(kind=ik), intent(in)        :: rank       !> rank of proc
    integer(kind=ik), intent(in)        :: N          !> number of blocks per proc

    hvy_id = lgt_id - rank*N

end subroutine
