!> \brief To a given light id "my_id", find the 3 (2D) or 7 (3D) sister block that have a common mother
!! block. They are returned in the sisters array. \n
!! input:    - light data array \n
!! output:   - light data array
! ********************************************************************************************

subroutine findSisters_tree(params, lgt_my_id, lgt_sisters_id, tree_ID)

    implicit none

    type (type_params), intent(in)      :: params                     !> user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_my_id                  !> this is the block whose sisters we look for
    !> here we will return the sisters. This array is allocated before calling
    !! this routine, and it can be either 4 or 8 or in length (2D / 3D), depending on whether
    !! you want to include the block whose sisters we look for or not.
    integer(kind=ik), intent(inout)     :: lgt_sisters_id(:)
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik)                    :: i                          ! loop variables
    integer(kind=tsize)                 :: tcb_sister, tcb_me         ! treecode variable
    integer(kind=ik)                    :: N_sisters
    integer(kind=ik)                    :: my_level
    logical                             :: exists

    ! check out how many sisters we look for. The number can be 4 or 8 in 2D or 3D. Note the
    ! block whose sisters we look for is returned as well
    N_sisters = size(lgt_sisters_id)

#ifdef DEV
    if ( N_sisters /= 4 .and. N_sisters /= 8 ) then
        call abort(123123, "findSisters_tree: you don't ask for a valid number of sisters")
    endif
#endif

    ! who am i?
    my_level = lgt_block( lgt_my_id, IDX_MESH_LVL )
    tcb_me = get_tc(lgt_block(lgt_my_id, IDX_TC_1 : IDX_TC_2))

    lgt_sisters_id = -1


    do i = 1, N_sisters
        if ( i-1 == tc_get_digit_at_level_b(tcb_me, dim=params%dim, level=my_level, max_level=params%Jmax) ) then
            ! this is the block itself, no need to look for it
            lgt_sisters_id(i) = lgt_my_id

        else
            ! sister, who are you?
            tcb_sister = tc_set_digit_at_level_b(tcb_me, i-1, dim=params%dim, level=my_level, max_level=params%Jmax)

            ! sister, are you there?
            call doesBlockExist_tree(tcb_sister, exists, lgt_sisters_id(i), dim=params%dim, level=my_level, tree_id=tree_ID, max_level=params%Jmax)
        end if
    enddo

end subroutine findSisters_tree
