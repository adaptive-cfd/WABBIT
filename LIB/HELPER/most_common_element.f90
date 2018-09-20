!> \brief from a list returns the most common element
!
!> \details
!! input:    - list of integers \n
!! output:   - most common element from list \n
!! \author mtriebeck
!! \date 18/9/18 - create
! ********************************************************************************************

function most_common_element( list )

    implicit none

    !> input list of integers
    integer(kind=ik), intent(in)	:: list(:)
    !> array of acummulated number of each element in list
    integer(kind=ik)              :: list_quantity(size(list))
    integer(kind=ik)     		:: max_element_pos, i, j, N, most_common_element


    N = size(list)

    ! for each rank in the array finding out its quantity
    do i = 1, N
	     list_quantity(i) = 1
	      do j = 1, N
		        if ( i/=j .and. list(i)==list(j) ) then
			           list_quantity(i) = list_quantity(i) + 1
		        end if
	       end do
    end do

    ! finding out position of max number in list_quantity and element itself
    max_element_pos = 0
    do i = 1, N
	     if ( list_quantity(i) > max_element_pos ) then
			      max_element_pos = i
	     end if
    end do

  most_common_element = list(max_element_pos)

end function most_common_element
