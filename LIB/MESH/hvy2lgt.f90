!> \brief convert heavy id into light id, requiring both to be on the same processor.
!! That means mpirank=1 has the heavy id hvy_id=17, which is the position of the
!! block in the heavy data array. The lgt_id is the corresponding position in the
!! global list of blocks, the light data.
!! input:     - heavy id
!!            - proc rank
!!            - number of blocks per proc
!! output:    - light data id
! ********************************************************************************************

subroutine hvy2lgt( lgt_id, hvy_id, rank, N )
    use module_params

    implicit none
    integer(kind=ik), intent(in)        :: hvy_id
    integer(kind=ik), intent(out)       :: lgt_id
    integer(kind=ik), intent(in)        :: rank       !> rank of proc
    integer(kind=ik), intent(in)        :: N          !> number of blocks per proc

    lgt_id = rank*N + hvy_id

end subroutine hvy2lgt
