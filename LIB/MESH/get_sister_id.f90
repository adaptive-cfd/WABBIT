! ********************************
! 2D AMR prototype
! --------------------------------
!
! return id's of sister blocks
!
! name: get_sister_id.f90
! date: 18.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine get_sister_id(id1, id2, id3, treecode, level)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), dimension(10), intent(in) :: treecode
    integer(kind=ik), intent(in)                :: level
    integer(kind=ik), intent(out)               :: id1, id2, id3

    integer(kind=ik)                            :: k, i
    integer(kind=ik), dimension(10)             :: sister_treecode
    integer(kind=ik), dimension(3)              :: ids

    ! return id = -1, if a sister-block is not present
    ids             = -1
    sister_treecode = treecode
    i               = 0

    do k = 1, 4
        ! sister treecode differs only on last element
        if (treecode(level) /= k-1) then
            i                       = i + 1
            sister_treecode(level)  = k-1
            call find_block_id(sister_treecode, ids(i))
        end if
    end do
    id1 = ids(1)
    id2 = ids(2)
    id3 = ids(3)

end subroutine get_sister_id
