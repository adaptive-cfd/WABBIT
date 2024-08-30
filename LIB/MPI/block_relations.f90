!> \brief Check if block is a leaf-block and has no daughters
function block_is_leaf(params, hvy_ID, check_empty)
    implicit none

    type (type_params), intent(in)      :: params      !< user defined parameter structure
    integer(kind=ik), intent(in)        :: hvy_ID      !< block to be investigated
    logical, intent(in), optional       :: check_empty !< if daughters should be checked for REF_TMP_EMPTY

    integer(kind=ik) :: i_d, lgt_ID_d
    logical :: block_is_leaf  !< output status

    block_is_leaf = all(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) == -1)

    if (.not. block_is_leaf .and. (present(check_empty))) then
        if (check_empty) then
            do i_d = 2+2**params%dim, 1+2**(params%dim+1)
                lgt_ID_d = hvy_family(hvy_ID, i_d)
                if (lgt_ID_d /= -1) then
                    if (lgt_block(lgt_ID_d, IDX_REFINE_STS) /= REF_TMP_EMPTY) return
                endif
            enddo
            block_is_leaf = .true.
        endif
    endif

end function block_is_leaf

!> \brief Check if block is a root-block and has no mother
function block_is_root(params, hvy_ID, check_empty)
    implicit none

    type (type_params), intent(in)      :: params      !< user defined parameter structure
    integer(kind=ik), intent(in)        :: hvy_ID      !< block to be investigated
    logical, intent(in), optional       :: check_empty !< if daughters should be checked for REF_TMP_EMPTY

    logical :: block_is_root  !< output status

    block_is_root = hvy_family(hvy_ID, 1) == -1

    if (.not. block_is_root .and. (present(check_empty))) then
        if (check_empty) block_is_root = lgt_block(hvy_family(hvy_ID, 1), IDX_REFINE_STS) == REF_TMP_EMPTY
    endif
end function block_is_root



!> \brief Check if block has a valid neighbor for same patch for the given level_difference
function block_has_valid_neighbor(params, hvy_ID, i_n, lvl_diff)
    implicit none

    type (type_params), intent(in)      :: params      !< user defined parameter structure
    integer(kind=ik), intent(in)        :: hvy_ID      !< block to be investigated
    integer(kind=ik), intent(in)        :: i_n         !< neighborhood relation
    integer(kind=ik), intent(in)        :: lvl_diff    !< lvl_difference, 0=medium, 1=coarser, -1=finer

    integer(kind=ik) :: lgt_id_n, ref_n, i_s
    logical :: block_has_valid_neighbor  !< output status

    block_has_valid_neighbor = .false.

    ! there are several different possibilities, once any is found we can exit
    do i_s = np_l(i_n, lvl_diff), np_u(i_n, lvl_diff)
        lgt_id_n = hvy_neighbor(hvy_ID, i_s)

        ! if it exists we have to make sure that its ref flag is not set to REF_TMP_EMPTY
        if (lgt_id_n /= -1) block_has_valid_neighbor = lgt_block(lgt_id_n, IDX_REFINE_STS) /= REF_TMP_EMPTY

        ! if we found a valid neighbor we can exit
        if (block_has_valid_neighbor) exit
    enddo
end function block_has_valid_neighbor