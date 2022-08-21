!> \brief To a given light id "my_id", find the 3 (2D) or 7 (3D) sister block that have a common mother
!! block. They are returned in the sisters array. \n
!! input:    - light data array \n
!! output:   - light data array
! ********************************************************************************************

subroutine findSisters_tree(params, lgt_my_id, lgt_sisters_id, lgt_block, lgt_n, lgt_sortednumlist, tree_ID)

    implicit none

    type (type_params), intent(in)      :: params                     !> user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_my_id                  !> this is the block whose sisters we look for
    !> here we will return the sisters. This array is allocated before calling
    !! this routine, and it can be either 4 or 8 or in length (2D / 3D), depending on whether
    !! you want to include the block whose sisters we look for or not.
    integer(kind=ik), intent(inout)     :: lgt_sisters_id(:)
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)            !> light data array
    integer(kind=ik), intent(in)        :: lgt_n(:)                   !> number of active blocks (light data)
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:,:)   !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik)                    :: i                          ! loop variables
    integer(kind=ik), allocatable, save :: all_treecodes(:,:)         ! treecode variable
    integer(kind=ik)                    :: N_sisters
    integer(kind=ik)                    :: mother_level, my_level
    logical                             :: exists

    ! check out how many sisters we look for. The number can be 4 or 8 in 2D or 3D. Note the
    ! block whose sisters we look for is returned as well
    N_sisters = size(lgt_sisters_id)

#ifdef DEV
    if ( N_sisters /= 4 .and. N_sisters /= 8 ) then
        call abort(123123, "findSisters_tree: you don't ask for a valid number of sisters")
    endif
#endif

    ! allocate an array for all treecodes (including all 4/8 sisters)
    if (.not.allocated(all_treecodes)) allocate( all_treecodes(1:N_sisters,1:params%max_treelevel) )
    ! initialize array as -1, since we do not use all of it, possibly (if we do not happen to
    ! be on the highest level)
    all_treecodes = -1

    my_level = lgt_block( lgt_my_id, params%max_treelevel + IDX_MESH_LVL )

    lgt_sisters_id = -1

    ! Find sisters. The sister blocks have the same mother, that means their treecode
    ! is idential up to the last entry
    mother_level = my_level - 1


    do i = 1, N_sisters
        if ( i-1 == lgt_block(lgt_my_id, my_level) ) then
            ! this is the block itself, no need to look for it
            lgt_sisters_id(i) = lgt_my_id

        else
            ! copy the identical mother level
            all_treecodes(i,1:mother_level) = lgt_block(lgt_my_id, 1:mother_level)
            ! the last treecode digit is (0..3) or (0..7)
            all_treecodes(i,mother_level+1) = i-1
!!!            tree_ID = lgt_block(lgt_my_id, params%max_treelevel + IDX_TREE_ID)
            ! look for the sisters in the list of blocks (light data), store their ID if found
            ! (-1 otherwise)
            call doesBlockExist_tree(all_treecodes(i,:), exists, lgt_sisters_id(i), lgt_sortednumlist, lgt_n, tree_ID)

        end if
    enddo

end subroutine findSisters_tree
