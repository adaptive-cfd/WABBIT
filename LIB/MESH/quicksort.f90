!> \brief Quicksort algorithm sorting after one or several numbers
!> \details This algorithm sorts the array a from position first
!! to position last. Sortsize defines the size (3 or 4 currently
!! supported). The actual comparison in tc_id_lower defines the
!! order of the elements. Currently two modi are implemented:
!!   -  sortsize=3, sort after a(:, 2) and a(:, 3) which are
!!      combined to one treecode
!!   -  sortsize=4, sort after treecode a(:, 2) and a(:, 3)
!!      as well as a(:, 4) which is level and tree_id together.
!!      Sorting will be done after tree_id, then level then TC
recursive subroutine quicksort(a, first, last, sortsize)
    use module_params

    implicit none

    integer(kind=ik), intent(inout)    :: a(:,:)
    integer(kind=ik), intent(in)       :: sortsize
    integer(kind=ik), dimension(sortsize)     :: x, t
    integer(kind=ik)                   :: first, last
    integer(kind=ik)                   :: i, j

    ! x = a( (first+last) / 2 , sortdim)
    x = a( (first+last) / 2 , :)
    i = first
    j = last

    ! if we've arrived at small lists, call interchange sort and return
    if ( j-i < 6) then
        call interchange_sort(a, first, last, sortsize)
        return
    endif

    ! otherwise do recursive quicksort
    do
        ! do while (a(i,sortdim) < x(sortdim))
        do while (tc_id_lower(a(i,2:sortsize), x(2:sortsize), sortsize>=4))
            i=i+1
        end do
        ! do while (x(sortdim) < a(j,sortdim))
        do while (tc_id_lower(x(2:sortsize), a(j,2:sortsize), sortsize>=4))
            j=j-1
        end do
        if (i >= j) exit
        t = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t
        i=i+1
        j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1, sortsize)
    if (j+1 < last)  call quicksort(a, j+1, last, sortsize)

end subroutine quicksort

!> \brief Interchange algorithm sorting after one or several numbers
!> \details This algorithm sorts the array a from position first
!! to position last. Sortsize defines the size (3 or 4 currently
!! supported). The actual comparison in tc_id_lower defines the
!! order of the elements. Currently two modi are implemented:
!!   -  sortsize=3, sort after a(:, 2) and a(:, 3) which are
!!      combined to one treecode
!!   -  sortsize=4, sort after treecode a(:, 2) and a(:, 3)
!!      as well as a(:, 4) which is level and tree_id together.
!!      Sorting will be done after tree_id, then level then TC
subroutine interchange_sort(a, left_end, right_end, sortsize)
    use module_params
    implicit none
    integer(kind=ik), intent(inout) ::  a(:,:)
    integer(kind=ik) :: left_end, right_end
    integer(kind=ik), intent(in) :: sortsize

    integer(kind=ik) :: i, j
    integer(kind=ik), dimension(sortsize) :: temp

    do i = left_end, right_end - 1
        do j = i+1, right_end
        ! if (a(j,sortdim) < a(i,sortdim)) then
        if (tc_id_lower(a(j,2:sortsize), a(i,2:sortsize), sortsize>=4)) then
            temp = a(i,:)
            a(i,:) = a(j,:)
            a(j,:) = temp
            end if
        end do
    end do

 end subroutine interchange_sort