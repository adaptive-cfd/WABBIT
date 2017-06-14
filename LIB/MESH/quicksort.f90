!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name quicksort.f90
!> \version 0.5
!> \author engels, msr
!
!> \brief quicksort subroutine
!
!>
!! input:    -
!! output:   -
!!
!!
!! = log ======================================================================================
!! \n
!! 14/06/17 - create
!!
! ********************************************************************************************
recursive subroutine quicksort(a, first, last, sortdim)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    integer(kind=tsize), intent(inout) ::  a(:,:)
    integer(kind=ik), intent(in) :: sortdim
    integer(kind=tsize), dimension(2) :: x, t
    integer(kind=ik) :: first, last
    integer(kind=ik) :: i, j

!---------------------------------------------------------------------------------------------
! variables initialization

    x = a( (first+last) / 2 , 2)
    i = first
    j = last

!---------------------------------------------------------------------------------------------
! main body

    ! if we've arrived at small lists, call interchange sort and return
    if ( j-i < 6) then
        call interchange_sort(a, first, last, sortdim)
        return
    endif

    ! otherwise do recursive quicksort
    do
        do while (a(i,sortdim) < x(sortdim))
            i=i+1
        end do
        do while (x(sortdim) < a(j,sortdim))
            j=j-1
        end do
        if (i >= j) exit
        t = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t
        i=i+1
        j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1, sortdim)
    if (j+1 < last)  call quicksort(a, j+1, last, sortdim)

end subroutine quicksort

subroutine interchange_sort(a, left_end, right_end, sortdim)
   use module_params
   implicit none
   integer(kind=tsize), intent(inout) ::  a(:,:)
   integer(kind=ik) :: left_end, right_end
   integer(kind=ik), intent(in) :: sortdim

   integer(kind=ik) :: i, j
   integer(kind=tsize), dimension(2) :: temp

   do i = left_end, right_end - 1
      do j = i+1, right_end
         if (a(i,sortdim) > a(j,sortdim)) then
            temp = a(i,:)
            a(i,:) = a(j,:)
            a(j,:) = temp
         end if
      end do
   end do

 end subroutine interchange_sort
