! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: ensure_completeness.f90
! version: 0.4
! author: msr
!
! sets refinement status to -2 for all sister blocks, if coarsening is possible
!
! input:    - light data array
! output:   - light data array
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine ensure_completeness( block_list, max_treelevel )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel

    ! loop variables
    integer(kind=ik)                    :: k, N, l, i, j

    ! sister ids
    integer(kind=ik)                    :: id(3)

    ! treecode variable
    integer(kind=ik)                    :: treecode(max_treelevel)
    ! block level
    integer(kind=ik)                    :: level

    ! function to compare treecodes
    logical                             :: array_compare

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N  = size(block_list, 1)
    i  = 0
    id = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all blocks
    do k = 1, N

        ! block is active
        if ( block_list( k, 1 ) /= -1 ) then

            ! if the block want to coarsen, check refinement status of sister blocks
            if ( block_list( k, max_treelevel+2 ) == -1) then

                ! get sister id's
                treecode = block_list( k, 1:max_treelevel )
                level    = block_list( k, max_treelevel+1 )

                i = 0
                do l = 1, 4
                    ! sister treecode differs only on last element
                    if ( block_list( k, level ) /= l-1) then

                        i               = i + 1
                        treecode(level) = l-1
                        ! find block id
                        do j = 1, N
                            if ( block_list( j, 1 ) /= -1 ) then
                                if (array_compare( block_list( j, 1:max_treelevel ), treecode, max_treelevel)) then
                                    id(i) = j
                                end if
                            end if
                        end do

                    end if
                end do

                ! if all sisters exists
                if ( (id(1)>0) .and. (id(2)>0) .and. (id(3)>0) ) then

                    ! if all sister blocks want to coarsen, then coarsening is allowed
                    if ( ( block_list( id(1), max_treelevel+2 ) == -1) .and. ( block_list( id(2), max_treelevel+2 ) == -1) .and. ( block_list( id(3), max_treelevel+2 ) == -1) ) then
                        block_list( k, max_treelevel+2 )      = -2
                        block_list( id(1), max_treelevel+2 )  = -2
                        block_list( id(2), max_treelevel+2 )  = -2
                        block_list( id(3), max_treelevel+2 )  = -2
                    end if

                end if

            end if

        end if

    end do

end subroutine ensure_completeness
