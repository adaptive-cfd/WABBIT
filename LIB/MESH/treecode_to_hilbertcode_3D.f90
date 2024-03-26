!> \brief convert given treecode to code in hilbert curve
!! algorithm from: Michael Bader. "Space-Filling Curves: An Introduction With Applications in Scientific Computing."
!! Springer Science & Business Media, 2012. p. 115
!! input:    - treecode
!!           - size of treecode
!! output:   - hilbert code
! ********************************************************************************************
!> \image html 3dhilbert.svg "The 12 Basic Patterns for the Hilbert Curve in 3D" width=400
!> \image latex 3dhilbert.eps "The 12 Basic Patterns for the Hilbert Curve in 3D"

subroutine treecode_to_hilbertcode_3D(treecode, hilbertcode, dim, level, max_level)

    implicit none
    !> dimension (2 or 3), defaults to 3
    integer(kind=ik), optional, intent(in)    :: dim  
    !> Level at which to encode, can be negative to set from max_level, defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level
    !> Numerical treecode in
    integer(kind=tsize), intent(in)     :: treecode
    !> Position on 2D hilbert curve
    integer(kind=tsize), intent(out)    :: hilbertcode
    integer(kind=ik)                    :: k                                    ! loop variable
    integer(kind=ik)                    :: tree_i, prev_tree_i                  ! treecode number
    ! hilbert pattern, note: use integer
    ! pattern : A B C D
    ! integer : 1 2 3 4
    integer(kind=ik)                    :: hilbert_pattern, prev_hilbert_pattern
    integer(kind=ik)                    :: hilbert_i
    integer(kind=ik)                    :: n_dim, n_level, max_tclevel          ! for optional setting

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! init
    prev_tree_i          = 0
    prev_hilbert_pattern = 0
    hilbert_pattern      = 1
    hilbert_i            = 0
    hilbertcode          = 0_tsize

    ! loop over treecode
    do k = 1, n_level

        tree_i = tc_get_level_b(treecode, dim=n_dim, level=k, max_level=max_tclevel)

        ! set hilbert pattern, in first step: always start with first pattern
        if ( k > 1 ) then
            call prev_pattern_to_pattern_3D( prev_hilbert_pattern, hilbert_pattern, prev_tree_i )
        end if

        ! position in new pattern
        call pattern_pos_3D( hilbert_pattern, tree_i, hilbert_i )

        hilbertcode = tc_set_level_b(hilbertcode, hilbert_i, dim=n_dim, level=k, max_level=max_tclevel)

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
