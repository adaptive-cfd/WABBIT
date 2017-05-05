!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name    adjacent_block_2D.f90
!> \version  0.4
!> \author   msr
!
!> \brief give treecode for adjacent block
!
!> 
!! input:    
!!                    - treecode for block N
!!                    - direction for neighbor search
!!                    - max treelevel
!!
!! output:   
!!                    - neighbor treecode, for neighbor on same level
!!
!!
!! = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4
! ********************************************************************************************

recursive subroutine adjacent_block_2D(me, neighbor, direction, level, max_treelevel)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    !> mesh level
    integer(kind=ik), intent(in)        :: level
    !> block treecode
    integer(kind=ik), intent(in)        :: me(max_treelevel)
    !> direction for neighbor search
    character(len=3), intent(in)        :: direction

    !> neighbor treecode
    integer(kind=ik), intent(out)       :: neighbor(max_treelevel)

    ! treecode variable
    integer(kind=ik)                    :: neighbor2(max_treelevel)
    ! loop variable
    integer(kind=ik)                    :: i

!---------------------------------------------------------------------------------------------
! variables initialization

    neighbor = -1

!---------------------------------------------------------------------------------------------
! main body

    select case(direction)

        case('__N')
        ! north
            neighbor(level) = modulo(me(level)+2, 4)
            i = level - 1

            do while (i /= 0)
                if ( (me(i+1)==2) .or. (me(i+1)==3) ) then
                    neighbor(1:i) = me(1:i)
                    i = 0
                else
                    neighbor(i) = modulo(me(i)+2, 4)
                    i = i - 1
                end if
            end do

        case('__S')
        ! south
            neighbor(level) = modulo(me(level)+2, 4)
            i = level - 1

            do while (i /= 0)
                if ( (me(i+1)==0) .or. (me(i+1)==1) ) then
                    neighbor(1:i) = me(1:i)
                    i = 0
                else
                    neighbor(i) = modulo(me(i)+2, 4)
                    i = i - 1
                end if
            end do

        case('__E')
        ! east
            if ( (me(level)==0) .or. (me(level)==2) ) then
                neighbor(level) = modulo(me(level)+1, 4)
            else
                neighbor(level) = modulo(me(level)-1, 4)
            end if

            i = level - 1

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

        case('__W')
        ! west
            if ( (me(level)==0) .or. (me(level)==2) ) then
                neighbor(level) = modulo(me(level)+1, 4)
            else
                neighbor(level) = modulo(me(level)-1, 4)
            end if

            i = level - 1

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

        case('_NE')
        ! northeast
            call adjacent_block_2D(me, neighbor2, '__N', level, max_treelevel)
            call adjacent_block_2D(neighbor2, neighbor, '__E', level, max_treelevel)

        case('_NW')
        ! northwest
            call adjacent_block_2D(me, neighbor2, '__N', level, max_treelevel)
            call adjacent_block_2D(neighbor2, neighbor, '__W', level, max_treelevel)

        case('_SE')
        ! southeast
            call adjacent_block_2D(me, neighbor2, '__S', level, max_treelevel)
            call adjacent_block_2D(neighbor2, neighbor, '__E', level, max_treelevel)

        case('_SW')
        ! southwest
            call adjacent_block_2D(me, neighbor2, '__S', level, max_treelevel)
            call adjacent_block_2D(neighbor2, neighbor, '__W', level, max_treelevel)

    end select

end subroutine adjacent_block_2D
