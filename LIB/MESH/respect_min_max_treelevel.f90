! ********************************
! WABBIT
! --------------------------------
!
! unset refinement status in respect
! of min/max treelevel
!
! name: respect_min_max_treelevel.f90
! date: 27.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine respect_min_max_treelevel()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                            :: k, N

    ! loop over all heavy data
    N           = blocks_params%number_max_blocks

    ! loop over all blocks
    do k = 1, N

        if ( blocks(k)%active ) then

            if ( (blocks(k)%refinement == 1) .and. (blocks(k)%level >= params%max_treelevel) ) then
                ! can not refine
                blocks(k)%refinement = 0
            end if

            if ( (blocks(k)%refinement == -1) .and. (blocks(k)%level <= params%min_treelevel) ) then
                ! can not coarsen
                blocks(k)%refinement = 0
            end if

        end if

    end do

end subroutine respect_min_max_treelevel
