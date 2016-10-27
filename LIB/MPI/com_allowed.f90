! ********************************
! WABBIT
! --------------------------------
!
! check if communication between two procs is possible
!
! name: com_allowed.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

logical function com_allowed(array, N)

    implicit none

    integer, intent(in)                   :: N
    logical, dimension(N), intent(in)     :: array

    integer                               :: i
    logical                               :: one_com_allowed

    com_allowed   = .false.
    one_com_allowed = .false.

    do i = 1, N

        if ( one_com_allowed == .true. ) then
            if ( array(i) == .true. ) then
                ! second allowed communication detected
                com_allowed = .true.
            end if
        end if

        if ( one_com_allowed == .false. ) then
            if ( array(i) == .true. ) then
                ! first allowed communication detected
                one_com_allowed = .true.
            end if
        end if

    enddo

end function com_allowed
