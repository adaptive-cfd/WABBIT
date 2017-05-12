!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name max_com_num.f90
!> \version 0.4
!> \author msr
!
!> \brief find max communication number
!
!>
!! input:    
!!           - com matrix line
!!           - proc rank
!!
!! output:   
!!           - number of maximal communication to one other proc
!!
!! = log ======================================================================================
!! \n
!! 13/01/17 - create for v0.4
! ********************************************************************************************

subroutine max_com_num( N, P, com_matrix_line, rank )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> number of coms, number of procs
    integer(kind=ik), intent(out)                   :: N, P

    !> com matrix line
    integer(kind=ik), intent(in)                    :: com_matrix_line(:)

    !> proc rank
    integer(kind=ik), intent(in)                    :: rank

    ! loop variable
    integer(kind=ik)                                :: k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset N
    N = 0
    ! reste P
    P = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all line elements
    do k = 1, size(com_matrix_line,1)

        ! communication to other proc, do not count internal communications
        if ( (com_matrix_line(k) /= 0) .and. (k /= rank+1) ) then

            ! max number
            N = max(N, com_matrix_line(k))
            ! add one proc
            P = P + 1

        end if

    end do

end subroutine max_com_num
