!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name adjacent_block_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief give treecode for adjacent block in 3D \n
!
!> 
!! input:    
!!                    - treecode for block N
!!                    - direction for neighbor search
!!                    - max treelevel
!!
!!output: 
!!                    - neighbor treecode, for neighbor on same level 
!!
!! \n
!!
!  --------------------------------------------------------------------------------------------
!> neighbor codes: \n
!  --------------- 
!> for imagination:  
!!                   - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front 
!!                   - edge: boundary between two sides - use sides numbers for coding
!!                   - corner: between three sides - so use all three sides numbers
!!                   - block on higher/lower level: block shares face/edge and one unique corner,
!!                     so use this corner code in second part of neighbor code
!!
!! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
!! edges:  '_12/___', '_13/___', '_14/___', '_15/___' \n
!!         '_62/___', '_63/___', '_64/___', '_65/___' \n
!!         '_23/___', '_25/___', '_43/___', '_45/___' \n
!! corner: '123/___', '134/___', '145/___', '152/___' \n
!!         '623/___', '634/___', '645/___', '652/___' \n
!! \n
!! complete neighbor code array, 74 possible neighbor relations \n
!! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!!               '_43/134', '_43/634', '_45/145', '_45/645' /) \n
! --------------------------------------------------------------------------------------------
!
!> \details = log ====================================================================================== 
!! \n
!! 27/01/17 - start
!
! ********************************************************************************************

recursive subroutine adjacent_block_3D(me, neighbor, direction, level, max_treelevel)

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
    character(len=7), intent(in)        :: direction

    !> neighbor treecode
    integer(kind=ik), intent(out)       :: neighbor(max_treelevel)
    ! treecode variable
    integer(kind=ik)                    :: neighbor2(max_treelevel), neighbor3(max_treelevel)

    ! loop variable
    integer(kind=ik)                    :: i

!---------------------------------------------------------------------------------------------
! variables initialization

    neighbor = me

!---------------------------------------------------------------------------------------------
! main body

    select case(direction)

        case('__1/___')
        ! '1'-side
            neighbor(level) = modulo(me(level) + 4, 8)
            i = level - 1

            do while (i /= 0)
                select case( me(i+1) )
                    case (0,1,2,3)
                        ! nothing to do, leave loop
                        i = 0

                    case (4,5,6,7)
                        neighbor(i) = modulo(me(i) + 4, 8)
                        i = i - 1

                end select
            end do

        case('__6/___')
        ! '6'-side
            neighbor(level) = modulo(me(level) - 4, 8)
            i = level - 1

            do while (i /= 0)
                select case( me(i+1) )
                    case (4,5,6,7)
                        ! nothing to do, leave loop
                        i = 0

                    case (0,1,2,3)
                        neighbor(i) = modulo(me(i) - 4, 8)
                        i = i - 1

                end select
            end do

        case('__3/___')
        ! '3'-side
            select case(me(level))
                case(0,1,4,5)
                    neighbor(level) = modulo(me(level) + 2, 8)
                case(2,3,6,7)
                    neighbor(level) = modulo(me(level) + 6, 8)
            end select

            i = level - 1

            do while (i /= 0)
                select case( me(i+1) )
                    case (0,1,4,5)
                        ! nothing to do, leave loop
                        i = 0

                    case (2,3,6,7)
                        select case(me(i))
                            case (0,1,4,5)
                                neighbor(i) = modulo(me(i) + 2, 8)
                            case (2,3,6,7)
                                neighbor(i) = modulo(me(i) + 6, 8)
                        end select
                        i = i - 1

                end select
            end do

        case('__5/___')
        ! '5'-side
            select case(me(level))
                case(2,3,6,7)
                    neighbor(level) = modulo(me(level) - 2, 8)
                case(0,1,4,5)
                    neighbor(level) = modulo(me(level) - 6, 8)
            end select

            i = level - 1

            do while (i /= 0)
                select case( me(i+1) )
                    case (2,3,6,7)
                        ! nothing to do, leave loop
                        i = 0

                    case (0,1,4,5)
                        select case(me(i))
                            case (2,3,6,7)
                                neighbor(i) = modulo(me(i) - 2, 8)
                            case (0,1,4,5)
                                neighbor(i) = modulo(me(i) - 6, 8)
                        end select
                        i = i - 1

                end select
            end do

        case('__2/___')
        ! '2'-side
            select case(me(level))
                case(1,3,5,7)
                    neighbor(level) = modulo(me(level) - 1, 8)
                case(0,2,4,6)
                    neighbor(level) = modulo(me(level) - 7, 8)
            end select

            i = level - 1

            do while (i /= 0)
                select case( me(i+1) )
                    case(1,3,5,7)
                        ! nothing to do, leave loop
                        i = 0

                    case(0,2,4,6)
                        select case(me(i))
                            case(1,3,5,7)
                                neighbor(i) = modulo(me(i) - 1, 8)
                            case(0,2,4,6)
                                neighbor(i) = modulo(me(i) - 7, 8)
                        end select
                        i = i - 1

                end select
            end do

        case('__4/___')
        ! '4'-side
            select case(me(level))
                case(0,2,4,6)
                    neighbor(level) = modulo(me(level) + 1, 8)
                case(1,3,5,7)
                    neighbor(level) = modulo(me(level) + 7, 8)
            end select

            i = level - 1

            do while (i /= 0)
                select case( me(i+1) )
                    case(0,2,4,6)
                        ! nothing to do, leave loop
                        i = 0

                    case(1,3,5,7)
                        select case(me(i))
                            case(0,2,4,6)
                                neighbor(i) = modulo(me(i) + 1, 8)
                            case(1,3,5,7)
                                neighbor(i) = modulo(me(i) + 7, 8)
                        end select
                        i = i - 1

                end select
            end do

        case('_12/___')
        ! '12'-edge
            call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

        case('_13/___')
        ! '13'-edge
            call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

        case('_14/___')
        ! '14'-edge
            call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

        case('_15/___')
        ! '15'-edge
            call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

        case('_62/___')
        ! '62'-edge
            call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

        case('_63/___')
        ! '63'-edge
            call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

        case('_64/___')
        ! '64'-edge
            call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

        case('_65/___')
        ! '65'-edge
            call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

        case('_23/___')
        ! '23'-edge
            call adjacent_block_3D(me, neighbor2, '__2/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

        case('_25/___')
        ! '25'-edge
            call adjacent_block_3D(me, neighbor2, '__2/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

        case('_43/___')
        ! '43'-edge
            call adjacent_block_3D(me, neighbor2, '__4/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

        case('_45/___')
        ! '45'-edge
            call adjacent_block_3D(me, neighbor2, '__4/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

        case('123/___')
        ! '123'-corner
            call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__2/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

        case('134/___')
        ! '134'-corner
            call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__3/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

        case('145/___')
        ! '145'-corner
            call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__4/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

        case('152/___')
        ! '152'-corner
            call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__5/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

        case('623/___')
        ! '623'-corner
            call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__2/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

        case('634/___')
        ! '634'-corner
            call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__3/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

        case('645/___')
        ! '645'-corner
            call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__4/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

        case('652/___')
        ! '652'-corner
            call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
            call adjacent_block_3D(neighbor3, neighbor2, '__5/___', level, max_treelevel)
            call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

    end select

end subroutine adjacent_block_3D
