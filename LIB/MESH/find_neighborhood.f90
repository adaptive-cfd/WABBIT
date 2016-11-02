! ********************************
! WABBIT
! --------------------------------
!
! find neighborhood id for neighbor relation
! between block1 -> block2
!
! name: find_neighborhood.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

integer function find_neighborhood(dir)

    use module_blocks

    implicit none

    character(len=3), intent(in)    :: dir

    ! -------------------------------------------------------------------------------------------------------------------------
    ! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
    ! -------------------------------------------------------------------------------------------------------------------------

    ! list id for corresponding neighborhood
    select case(dir)
        case('__N')
            find_neighborhood = 3
        case('__S')
            find_neighborhood = 1
        case('__W')
            find_neighborhood = 2
        case('__E')
            find_neighborhood = 4
        case('_NE')
            find_neighborhood = 8
        case('_NW')
            find_neighborhood = 7
        case('_SE')
            find_neighborhood = 6
        case('_SW')
            find_neighborhood = 5
        case('NNE')
            find_neighborhood = 11
        case('NNW')
            find_neighborhood = 12
        case('SSE')
            find_neighborhood = 9
        case('SSW')
            find_neighborhood = 10
        case('ENE')
            find_neighborhood = 15
        case('ESE')
            find_neighborhood = 16
        case('WNW')
            find_neighborhood = 13
        case('WSW')
            find_neighborhood = 14
    end select

end function find_neighborhood
