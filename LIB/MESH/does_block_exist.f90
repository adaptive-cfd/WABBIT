!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name does_block_exist.f90
!> \version 0.5
!> \author engels
!
!> \brief given a treecode, look for the block, return .true. and its light_id if found
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine does_block_exist(treecode, exists, light_id, lgt_sortednumlist, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> block treecode we are looking for, array representation
    integer(kind=ik), intent(in)        :: treecode(:)
    !> logical, .true. if block with treecode exists
    logical, intent(out)                :: exists
    !> light data id of block if found
    integer(kind=ik), intent(out)       :: light_id
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)
    !> it helps to know how many active light blocks we have in total
    integer(kind=ik), intent(in)        :: lgt_n
    ! loop variables
    integer(kind=ik)                    :: k, i1, i2, imid
    ! numerical treecode
    integer(kind=tsize)                 :: num_treecode

    exists   = .false.
    light_id = -1

    !> 1st: given the array treecode, compute the numerical value of the treecode we're
    !! looking for. note these values are stored for easier finding in lgt_sortednumlist
    num_treecode = treecode2int(treecode)

    !> 2nd: binary search. start with the entire interval, then choose either right or left half
    i1 = 1
    i2 = lgt_n

    ! if the interval is only two entries, check both and return the winner
    do while ( abs(i2-i1) >= 2 )
        ! cut interval in two parts
        imid = (i1+i2) / 2
        ! if (num_treecode < lgt_sortednumlist(imid,2)) then
        if (lgt_sortednumlist(imid,2) < num_treecode) then
            i1 = imid
        else
            i2 = imid
        end if
    end do

    if ( num_treecode == lgt_sortednumlist(i1,2) .and. lgt_sortednumlist(i1,1) > 0) then
        ! found the block we're looking for
        exists = .true.
        light_id = int( lgt_sortednumlist(i1,1), kind=ik)
        return
    end if

    if ( num_treecode == lgt_sortednumlist(i2,2) .and. lgt_sortednumlist(i2,1) > 0) then
        ! found the block we're looking for
        exists = .true.
        light_id = int( lgt_sortednumlist(i2,1), kind=ik)
        return
    end if
end subroutine does_block_exist
