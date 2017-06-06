!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name treecode_to_hilbertcode_3D.f90
!> \version 0.5
!> \author sm
!> \brief convert given treecode to code in hilbert curve
!
!> 
!! algorithm from: Michael Bader. "Space-Filling Curves: An Introduction With Applications in Scientific Computing."
!! Springer Science & Business Media, 2012. p. 115
!! \n
!! input:   
!!           - treecode
!!           - size of treecode
!!
!! output:   
!!           - hilbert code
!!
!! = log ======================================================================================
!! \n
!! 02/05/17 - create 
! ********************************************************************************************
!
!> \image html 3dhilbert.svg "The 12 Basic Patterns for the Hilbert Curve in 3D" width=400
!> \image latex 3dhilbert.eps "The 12 Basic Patterns for the Hilbert Curve in 3D"

subroutine treecode_to_hilbertcode_3D(treecode, hilbertcode, n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> treecode size
    integer(kind=ik), intent(in)        :: n

    !> treecode
    integer(kind=ik), intent(in)        :: treecode(n)

    !> hilbert code
    integer(kind=ik), intent(out)       :: hilbertcode(n)

    ! loop variable
    integer(kind=ik)                    :: k

    ! treecode number
    integer(kind=ik)                    :: tree_i, prev_tree_i

    integer(kind=ik)                    :: hilbert_pattern, prev_hilbert_pattern

!---------------------------------------------------------------------------------------------
! variables initialization

    prev_tree_i          = 0
    prev_hilbert_pattern = 0
    hilbert_pattern      = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over treecode
    do k = 1, n

        ! read treecode number
        ! -1 is used as 0
        if ( treecode(k) == -1 ) then
            tree_i = 0
        else
            tree_i = treecode(k)
        end if

        ! set hilbert pattern, in first step: always start with first pattern
        if ( k == 1 ) then
            ! first step
            hilbert_pattern = 1
        else
            ! calculate pattern
            call prev_pattern_to_pattern_3D( prev_hilbert_pattern, hilbert_pattern, prev_tree_i )
        end if

        ! position in new pattern
        call pattern_pos_3D( hilbert_pattern, tree_i, hilbertcode(k) )

        ! save previous pattern and treecode
        prev_hilbert_pattern    = hilbert_pattern
        prev_tree_i             = tree_i

    end do

end subroutine treecode_to_hilbertcode_3D

!------------------------------------------------------------------------------------------------------

! subroutine to calculate new hilbert pattern
subroutine prev_pattern_to_pattern_3D( prev_hilbert_pattern, hilbert_pattern, pos )

    implicit none

    ! position number
    integer(kind=ik), intent(in)        :: pos

    ! hilbert pattern, note: use integer
    ! pattern : A  B  C  D  E  F  G  H  I  J  K  L
    ! integer : 1  2  3  4  5  6  7  8  9  10 11 12
    integer(kind=ik), intent(out)       :: hilbert_pattern
    integer(kind=ik), intent(in)        :: prev_hilbert_pattern

    select case(prev_hilbert_pattern)
        ! case A
        case(1)
            select case(pos)
                case(0)
                    hilbert_pattern = 3
                case(1)
                    hilbert_pattern = 5
                case(2)
                    hilbert_pattern = 8
                case(3)
                    hilbert_pattern = 5
                case(4)
                    hilbert_pattern = 4
                case(5)
                    hilbert_pattern = 11
                case(6)
                    hilbert_pattern = 8
                case(7)
                    hilbert_pattern = 11
            end select
        ! case B
        case(2)
            select case(pos)
                case(0)
                    hilbert_pattern = 12
                case(1)
                    hilbert_pattern = 3
                case(2)
                    hilbert_pattern = 12
                case(3)
                    hilbert_pattern = 7
                case(4)
                    hilbert_pattern = 6
                case(5)
                    hilbert_pattern = 4
                case(6)
                    hilbert_pattern = 6
                case(7)
                    hilbert_pattern = 7
            end select
        ! case C
        case(3)
            select case(pos)
                case(0)
                    hilbert_pattern = 5
                case(1)
                    hilbert_pattern = 12
                case(2)
                    hilbert_pattern = 1
                case(3)
                    hilbert_pattern = 2
                case(4)
                    hilbert_pattern = 9
                case(5)
                    hilbert_pattern = 9
                case(6)
                    hilbert_pattern = 1
                case(7)
                    hilbert_pattern = 2
            end select
        ! case D
        case(4)
            select case(pos)
                case(0)
                    hilbert_pattern = 10
                case(1)
                    hilbert_pattern = 10
                case(2)
                    hilbert_pattern = 1
                case(3)
                    hilbert_pattern = 2
                case(4)
                    hilbert_pattern = 11
                case(5)
                    hilbert_pattern = 6
                case(6)
                    hilbert_pattern = 1
                case(7)
                    hilbert_pattern = 2
            end select
        ! case E
        case(5)
            select case(pos)
                case(0)
                    hilbert_pattern = 1
                case(1)
                    hilbert_pattern = 6
                case(2)
                    hilbert_pattern = 7
                case(3)
                    hilbert_pattern = 6
                case(4)
                    hilbert_pattern = 3
                case(5)
                    hilbert_pattern = 3
                case(6)
                    hilbert_pattern = 10
                case(7)
                    hilbert_pattern = 10
            end select
        ! case F
        case(6)
            select case(pos)
                case(0)
                    hilbert_pattern = 4
                case(1)
                    hilbert_pattern = 4
                case(2)
                    hilbert_pattern = 9
                case(3)
                    hilbert_pattern = 9
                case(4)
                    hilbert_pattern = 5
                case(5)
                    hilbert_pattern = 2
                case(6)
                    hilbert_pattern = 5
                case(7)
                    hilbert_pattern = 8
            end select
        ! case G
        case(7)
            select case(pos)
                case(0)
                    hilbert_pattern = 2
                case(1)
                    hilbert_pattern = 5
                case(2)
                    hilbert_pattern = 10
                case(3)
                    hilbert_pattern = 5
                case(4)
                    hilbert_pattern = 2
                case(5)
                    hilbert_pattern = 11
                case(6)
                    hilbert_pattern = 9
                case(7)
                    hilbert_pattern = 11
            end select
        ! case H
        case(8)
            select case(pos)
                case(0)
                    hilbert_pattern = 12
                case(1)
                    hilbert_pattern = 1
                case(2)
                    hilbert_pattern = 12
                case(3)
                    hilbert_pattern = 10
                case(4)
                    hilbert_pattern = 6
                case(5)
                    hilbert_pattern = 1
                case(6)
                    hilbert_pattern = 6
                case(7)
                    hilbert_pattern = 9
            end select
        ! case I
        case(9)
            select case(pos)
                case(0)
                    hilbert_pattern = 7
                case(1)
                    hilbert_pattern = 8
                case(2)
                    hilbert_pattern = 3
                case(3)
                    hilbert_pattern = 3
                case(4)
                    hilbert_pattern = 7
                case(5)
                    hilbert_pattern = 8
                case(6)
                    hilbert_pattern = 11
                case(7)
                    hilbert_pattern = 6
            end select
        ! case J
        case(10)
            select case(pos)
                case(0)
                    hilbert_pattern = 7
                case(1)
                    hilbert_pattern = 8
                case(2)
                    hilbert_pattern = 5
                case(3)
                    hilbert_pattern = 12
                case(4)
                    hilbert_pattern = 7
                case(5)
                    hilbert_pattern = 8
                case(6)
                    hilbert_pattern = 4
                case(7)
                    hilbert_pattern = 4
            end select
        ! case K
        case(11)
            select case(pos)
                case(0)
                    hilbert_pattern = 4
                case(1)
                    hilbert_pattern = 4
                case(2)
                    hilbert_pattern = 9
                case(3)
                    hilbert_pattern = 9
                case(4)
                    hilbert_pattern = 1
                case(5)
                    hilbert_pattern = 12
                case(6)
                    hilbert_pattern = 7
                case(7)
                    hilbert_pattern = 12
            end select
        ! case L
        case(12)
            select case(pos)
                case(0)
                    hilbert_pattern = 11
                case(1)
                    hilbert_pattern = 2
                case(2)
                    hilbert_pattern = 11
                case(3)
                    hilbert_pattern = 8
                case(4)
                    hilbert_pattern = 3
                case(5)
                    hilbert_pattern = 3
                case(6)
                    hilbert_pattern = 10
                case(7)
                    hilbert_pattern = 10
            end select

    end select

end subroutine prev_pattern_to_pattern_3D

!------------------------------------------------------------------------------------------------------------

! position in given hilbert pattern
subroutine pattern_pos_3D( hilbert_pattern, pos, hilbert_pos )

    implicit none

    ! position number
    integer(kind=ik), intent(in)        :: pos

    ! hilbert pattern, note: use integer
    ! pattern : A  B  C  D  E  F  G  H  I  J  K  L
    ! integer : 1  2  3  4  5  6  7  8  9  10 11 12
    integer(kind=ik), intent(in)        :: hilbert_pattern

    ! hilbert position number
    ! note: position starts with 0 and ends with 7 (as treecode position)
    integer(kind=ik), intent(out)       :: hilbert_pos

    select case(hilbert_pattern)
        ! case A
        case(1)
            select case(pos)
                case(0)
                    hilbert_pos = 0
                case(1)
                    hilbert_pos = 1
                case(2)
                    hilbert_pos = 3
                case(3)
                    hilbert_pos = 2
                case(4)
                    hilbert_pos = 7
                case(5)
                    hilbert_pos = 6
                case(6)
                    hilbert_pos = 4
                case(7)
                    hilbert_pos = 5
            end select

        ! case B
        case(2)
           select case(pos)
                case(0)
                    hilbert_pos = 6
                case(1)
                    hilbert_pos = 7
                case(2)
                    hilbert_pos = 5
                case(3)
                    hilbert_pos = 4
                case(4)
                    hilbert_pos = 1
                case(5)
                    hilbert_pos = 0
                case(6)
                    hilbert_pos = 2
                case(7)
                    hilbert_pos = 3
            end select
        ! case C
        case(3)
           select case(pos)
                case(0)
                    hilbert_pos = 0
                case(1)
                    hilbert_pos = 7
                case(2)
                    hilbert_pos = 1
                case(3)
                    hilbert_pos = 6
                case(4)
                    hilbert_pos = 3
                case(5)
                    hilbert_pos = 4
                case(6)
                    hilbert_pos = 2
                case(7)
                    hilbert_pos = 5
            end select
        ! case D
        case(4)
            select case(pos)
                case(0)
                    hilbert_pos = 4
                case(1)
                    hilbert_pos = 3
                case(2)
                    hilbert_pos = 5
                case(3)
                    hilbert_pos = 2
                case(4)
                    hilbert_pos = 7
                case(5)
                    hilbert_pos = 0
                case(6)
                    hilbert_pos = 6
                case(7)
                    hilbert_pos = 1
            end select
        ! case E
        case(5)
            select case(pos)
                case(0)
                    hilbert_pos = 0
                case(1)
                    hilbert_pos = 3
                case(2)
                    hilbert_pos = 7
                case(3)
                    hilbert_pos = 4
                case(4)
                    hilbert_pos = 1
                case(5)
                    hilbert_pos = 2
                case(6)
                    hilbert_pos = 6
                case(7)
                    hilbert_pos = 5
            end select
        ! case F
        case(6)
            select case(pos)
                case(0)
                    hilbert_pos = 2
                case(1)
                    hilbert_pos = 1
                case(2)
                    hilbert_pos = 5
                case(3)
                    hilbert_pos = 6
                case(4)
                    hilbert_pos = 3
                case(5)
                    hilbert_pos = 0
                case(6)
                    hilbert_pos = 4
                case(7)
                    hilbert_pos = 7
            end select
        ! case G
        case(7)
            select case(pos)
                case(0)
                    hilbert_pos = 4
                case(1)
                    hilbert_pos = 5
                case(2)
                    hilbert_pos = 7
                case(3)
                    hilbert_pos = 6
                case(4)
                    hilbert_pos = 3
                case(5)
                    hilbert_pos = 2
                case(6)
                    hilbert_pos = 0
                case(7)
                    hilbert_pos = 1
            end select
        ! case H
        case(8)
            select case(pos)
                case(0)
                    hilbert_pos = 2
                case(1)
                    hilbert_pos = 3
                case(2)
                    hilbert_pos = 1
                case(3)
                    hilbert_pos = 0
                case(4)
                    hilbert_pos = 5
                case(5)
                    hilbert_pos = 4
                case(6)
                    hilbert_pos = 6
                case(7)
                    hilbert_pos = 7
            end select
        ! case I
        case(9)
            select case(pos)
                case(0)
                    hilbert_pos = 2
                case(1)
                    hilbert_pos = 5
                case(2)
                    hilbert_pos = 3
                case(3)
                    hilbert_pos = 4
                case(4)
                    hilbert_pos = 1
                case(5)
                    hilbert_pos = 6
                case(6)
                    hilbert_pos = 0
                case(7)
                    hilbert_pos = 7
            end select
        ! case J
        case(10)
            select case(pos)
                case(0)
                    hilbert_pos = 6
                case(1)
                    hilbert_pos = 1
                case(2)
                    hilbert_pos = 7
                case(3)
                    hilbert_pos = 0
                case(4)
                    hilbert_pos = 5
                case(5)
                    hilbert_pos = 2
                case(6)
                    hilbert_pos = 4
                case(7)
                    hilbert_pos = 3
            end select
        ! case K
        case(11)
            select case(pos)
                case(0)
                    hilbert_pos = 6
                case(1)
                    hilbert_pos = 5
                case(2)
                    hilbert_pos = 1
                case(3)
                    hilbert_pos = 2
                case(4)
                    hilbert_pos = 7
                case(5)
                    hilbert_pos = 4
                case(6)
                    hilbert_pos = 0
                case(7)
                    hilbert_pos = 3
            end select
        ! case L
        case(12)
            select case(pos)
                case(0)
                    hilbert_pos = 4
                case(1)
                    hilbert_pos = 7
                case(2)
                    hilbert_pos = 3
                case(3)
                    hilbert_pos = 0
                case(4)
                    hilbert_pos = 5
                case(5)
                    hilbert_pos = 6
                case(6)
                    hilbert_pos = 2
                case(7)
                    hilbert_pos = 1
            end select
    end select

end subroutine pattern_pos_3D
