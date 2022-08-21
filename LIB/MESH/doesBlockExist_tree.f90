!> \brief given a treecode, look for the block, return .true. and its lgtID if found
! ********************************************************************************************

subroutine doesBlockExist_tree(tcBlock, exists, lgtID, lgt_sortednumlist, lgt_n, tree_ID)
    implicit none

    integer(kind=ik), intent(in)        :: tcBlock(:)                 !> block treecode we are looking for, array representation
    logical, intent(out)                :: exists                     !> true if block with treecode exists
    integer(kind=ik), intent(out)       :: lgtID                      !> light data id of block if found
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:,:)   !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(in)        :: lgt_n(:)                   !> it helps to know how many active light blocks we have in total
    integer(kind=ik), intent(in)        :: tree_ID                    !> index of the tree we are looking at
    integer(kind=ik)                    :: k, i1, i2, imid            ! loop variables
    integer(kind=tsize)                 :: num_treecode               ! numerical treecode

    exists   = .false.
    lgtID = -1

    !> 1st: given the array treecode, compute the numerical value of the treecode we're
    !! looking for. note these values are stored for easier finding in lgt_sortednumlist
    num_treecode = treecode2int(tcBlock, tree_ID)

    !> 2nd: binary search. start with the entire interval, then choose either right or left half
    i1 = 1
    i2 = lgt_n(tree_ID)

    ! if the interval is only two entries, check both and return the winner
    do while ( abs(i2-i1) >= 2 )
        ! cut interval in two parts
        imid = (i1+i2) / 2
        if (lgt_sortednumlist(imid,2,tree_ID) < num_treecode) then
            i1 = imid
        else
            i2 = imid
        end if
    end do

    if ( num_treecode == lgt_sortednumlist(i1,2,tree_ID) .and. lgt_sortednumlist(i1,1,tree_ID) > 0) then
        ! found the block we're looking for
        exists = .true.
        lgtID = int( lgt_sortednumlist(i1,1,tree_ID), kind=ik)
        return
    end if

    if ( num_treecode == lgt_sortednumlist(i2,2,tree_ID) .and. lgt_sortednumlist(i2,1,tree_ID) > 0) then
        ! found the block we're looking for
        exists = .true.
        lgtID = int( lgt_sortednumlist(i2,1,tree_ID), kind=ik)
        return
    end if
end subroutine doesBlockExist_tree
