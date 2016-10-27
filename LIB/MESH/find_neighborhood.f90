! ********************************
! WABBIT
! --------------------------------
!
! find neighborhood id
!
! name: find_neighborhood.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

integer function find_neighborhood(block_num, neighbor_dir)

    use module_blocks

    implicit none

    integer, intent(in)             :: block_num
    character(len=2), intent(in)    :: neighbor_dir

    integer                         :: k
    character(len=2)                :: my_dir

    select case(neighbor_dir)
        case('NO')
            my_dir = 'SO'
        case('SO')
            my_dir = 'NO'
        case('WE')
            my_dir = 'EA'
        case('EA')
            my_dir = 'WE'
        case('NE')
            my_dir = 'SW'
        case('NW')
            my_dir = 'SE'
        case('SE')
            my_dir = 'NW'
        case('SW')
            my_dir = 'NE'
    end select

    do k = 1, 8
        if ( blocks(block_num)%neighbor_dir(k) == my_dir ) then
            find_neighborhood = k
        end if
    enddo

end function find_neighborhood
