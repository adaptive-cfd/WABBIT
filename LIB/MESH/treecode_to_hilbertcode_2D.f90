!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name treecode_to_hilbertcode_2D.f90
!> \version 0.4
!> \author msr
!
!> \brief convert given treecode to code in hilbert curve
!
!> 
!! input:   
!!           - treecode
!!           - size of treecode
!!
!! output:   
!!           - hilbert code
!!
!! hilbert pattern
! ---------------
!>
!!   |A: |B: |C: |D: |
!!   |---|---|---|---|
!!   |0 3|0 1|2 1|2 3|
!!   |1 2|3 2|3 0|1 0|
!!
!! one level up:
!!    |A: |B: |C: |D: |
!!    |---|---|---|---|
!!    |B D|A B|C C|D A|
!!    |A A|C B|B D|D C|
!!
!! = log ======================================================================================
!! \n
!! 25/01/17 - create
! ********************************************************************************************
!
!> \image html hilbert.svg "A Path taken by the Hilbert Curve in 2D" width=300
!> \image latex hilbert.eps "A Path taken by the Hilbert Curve in 2D"

subroutine treecode_to_hilbertcode_2D(treecode, hilbertcode, n)

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

    ! hilbert pattern, note: use integer
    ! pattern : A B C D
    ! integer : 1 2 3 4
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

        ! set hilbert pattern, in first step: always start with A pattern
        if ( k == 1 ) then
            ! first step
            hilbert_pattern = 1
        else
            ! calculate pattern
            call prev_pattern_to_pattern_2D( prev_hilbert_pattern, hilbert_pattern, prev_tree_i )
        end if

        ! position in new pattern
        call pattern_pos_2D( hilbert_pattern, tree_i, hilbertcode(k) )

        ! save previous pattern and treecode
        prev_hilbert_pattern    = hilbert_pattern
        prev_tree_i             = tree_i

    end do


end subroutine treecode_to_hilbertcode_2D

!---------------------------------------------------------------------------------------------

! subroutine to calculate new hilbert pattern
subroutine prev_pattern_to_pattern_2D( prev_hilbert_pattern, hilbert_pattern, pos )

    implicit none

    ! position number
    integer(kind=ik), intent(in)        :: pos

    ! hilbert pattern, note: use integer
    ! pattern : A B C D
    ! integer : 1 2 3 4
    integer(kind=ik), intent(out)       :: hilbert_pattern
    integer(kind=ik), intent(in)        :: prev_hilbert_pattern

    select case(prev_hilbert_pattern)
        ! case A
        case(1)
            select case(pos)
                case(0)
                    hilbert_pattern = 2
                case(1)
                    hilbert_pattern = 4
                case(2)
                    hilbert_pattern = 1
                case(3)
                    hilbert_pattern = 1
            end select

        ! case B
        case(2)
            select case(pos)
                case(0)
                    hilbert_pattern = 1
                case(1)
                    hilbert_pattern = 2
                case(2)
                    hilbert_pattern = 3
                case(3)
                    hilbert_pattern = 2
            end select

        ! case C
        case(3)
            select case(pos)
                case(0)
                    hilbert_pattern = 3
                case(1)
                    hilbert_pattern = 3
                case(2)
                    hilbert_pattern = 2
                case(3)
                    hilbert_pattern = 4
            end select

        ! case D
        case(4)
            select case(pos)
                case(0)
                    hilbert_pattern = 4
                case(1)
                    hilbert_pattern = 1
                case(2)
                    hilbert_pattern = 4
                case(3)
                    hilbert_pattern = 3
            end select

    end select

end subroutine prev_pattern_to_pattern_2D

! position in given hilbert pattern
subroutine pattern_pos_2D( hilbert_pattern, pos, hilbert_pos )

    implicit none

    ! position number
    integer(kind=ik), intent(in)        :: pos

    ! hilbert pattern, note: use integer
    ! pattern : A B C D
    ! integer : 1 2 3 4
    integer(kind=ik), intent(in)        :: hilbert_pattern

    ! hilbert position number
    ! note: position starts with 0 and ends with 3 (as treecode position)
    integer(kind=ik), intent(out)       :: hilbert_pos

    select case(hilbert_pattern)
        ! case A
        case(1)
            select case(pos)
                case(0)
                    hilbert_pos = 0
                case(1)
                    hilbert_pos = 3
                case(2)
                    hilbert_pos = 1
                case(3)
                    hilbert_pos = 2
            end select

        ! case B
        case(2)
            select case(pos)
                case(0)
                    hilbert_pos = 0
                case(1)
                    hilbert_pos = 1
                case(2)
                    hilbert_pos = 3
                case(3)
                    hilbert_pos = 2
            end select

        ! case C
        case(3)
            select case(pos)
                case(0)
                    hilbert_pos = 2
                case(1)
                    hilbert_pos = 1
                case(2)
                    hilbert_pos = 3
                case(3)
                    hilbert_pos = 0
            end select

        ! case D
        case(4)
            select case(pos)
                case(0)
                    hilbert_pos = 2
                case(1)
                    hilbert_pos = 3
                case(2)
                    hilbert_pos = 1
                case(3)
                    hilbert_pos = 0
            end select

    end select

end subroutine pattern_pos_2D
