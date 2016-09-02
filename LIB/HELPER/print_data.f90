! ********************************
! 2D AMR prototype
! --------------------------------
!
! print data field to screen
!
! name: print_data.f90
! date: 23.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine print_data(data, Bs)

    implicit none

    integer, intent(in)                       :: Bs
    real(8), dimension(Bs,Bs), intent(in)     :: data

    integer                                   :: i, j

    do i = 1, Bs
        do j = 1, Bs
            write(*,'(f5.2,1x)',advance='no') data(i,j)
        end do
        write(*,*) ""
    end do

    write(*,*)

end subroutine print_data
