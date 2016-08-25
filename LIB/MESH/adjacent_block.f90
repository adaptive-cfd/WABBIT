! ********************************
! 2D AMR prototype
! --------------------------------
!
! give treecode for adjacent block
!
! name: adjacent_block.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

recursive subroutine adjacent_block(me, neighbor, direction)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), dimension(10), intent(in)        :: me
    integer(kind=ik), dimension(10), intent(out)       :: neighbor
    character(len=2), intent(in)                       :: direction

    integer(kind=ik)                                   :: i, n, treecode_size
    integer(kind=ik), dimension(10)                    :: neighbor2

    neighbor    = -1
    n           = 0

    ! treecode size
    n = treecode_size(me)

    select case(direction)

        case('NO')
        ! north
            neighbor(n) = modulo(me(n)+2, 4)
            i = n - 1

            do while (i /= 0)
                if ( (me(i+1)==2) .or. (me(i+1)==3) ) then
                    neighbor(1:i) = me(1:i)
                    i = 0
                else
                    neighbor(i) = modulo(me(i)+2, 4)
                    i = i - 1
                end if
            end do

        case('SO')
        ! south
            neighbor(n) = modulo(me(n)+2, 4)
            i = n - 1

            do while (i /= 0)
                if ( (me(i+1)==0) .or. (me(i+1)==1) ) then
                    neighbor(1:i) = me(1:i)
                    i = 0
                else
                    neighbor(i) = modulo(me(i)+2, 4)
                    i = i - 1
                end if
            end do

        case('EA')
        ! east
            if ( (me(n)==0) .or. (me(n)==2) ) then
                neighbor(n) = modulo(me(n)+1, 4)
            else
                neighbor(n) = modulo(me(n)-1, 4)
            end if

            i = n - 1

            do while (i /= 0)
                if ( (me(i+1)==0) .or. (me(i+1)==2) ) then
                    neighbor(1:i) = me(1:i)
                    i = 0
                else
                    if ( (me(i)==1) .or. (me(i)==3) ) then
                        neighbor(i) = modulo(me(i)-1, 4)
                    else
                        neighbor(i) = modulo(me(i)+1, 4)
                    end if
                    i = i - 1
                end if
            end do

        case('WE')
        ! west
            if ( (me(n)==0) .or. (me(n)==2) ) then
                neighbor(n) = modulo(me(n)+1, 4)
            else
                neighbor(n) = modulo(me(n)-1, 4)
            end if

            i = n - 1

            do while (i /= 0)
                if ( (me(i+1)==1) .or. (me(i+1)==3) ) then
                    neighbor(1:i) = me(1:i)
                    i = 0
                else
                    if ( (me(i)==1) .or. (me(i)==3) ) then
                        neighbor(i) = modulo(me(i)-1, 4)
                    else
                        neighbor(i) = modulo(me(i)+1, 4)
                    end if
                    i = i - 1
                end if
            end do

        case('NE')
        ! northeast
            call adjacent_block(me, neighbor2, 'NO')
            call adjacent_block(neighbor2, neighbor, 'EA')

        case('NW')
        ! northwest
            call adjacent_block(me, neighbor2, 'NO')
            call adjacent_block(neighbor2, neighbor, 'WE')

        case('SE')
        ! southeast
            call adjacent_block(me, neighbor2, 'SO')
            call adjacent_block(neighbor2, neighbor, 'EA')

        case('SW')
        ! southwest
            call adjacent_block(me, neighbor2, 'SO')
            call adjacent_block(neighbor2, neighbor, 'WE')

    end select

end subroutine adjacent_block
