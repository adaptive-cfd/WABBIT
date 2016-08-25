! ********************************
! 2D AMR prototype
! --------------------------------
!
! unset refinement status in respect
! of min/max treelevel
!
! name: respect_min_max_treelevel.f90
! date: 18.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine respect_min_max_treelevel()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                            :: k, N, block_num

    N  = size(blocks_params%active_list, dim=1)

    ! loop over all active blocks
    do k = 1, N

        block_num = blocks_params%active_list(k)

        if ( (blocks(block_num)%refinement == 1) .and. (blocks(block_num)%level >= params%max_treelevel) ) then
            ! can not refine
            blocks(block_num)%refinement = 0
        end if

        if ( (blocks(block_num)%refinement == -1) .and. (blocks(block_num)%level <= params%min_treelevel) ) then
            ! can not coarsen
            blocks(block_num)%refinement = 0
        end if

    end do

end subroutine respect_min_max_treelevel
