! ********************************
! WABBIT
! --------------------------------
!
! sets refinement status to -2 for all
! sister blocks, if coarsening is
! possible
!
! name: ensure_completeness.f90
! date: 28.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine ensure_completeness()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik) :: k, N, id1, id2, id3

    N  = blocks_params%number_max_blocks

    ! loop over all blocks
    do k = 1, N

        if (blocks(k)%active) then

            ! if the block want to coarsen, check refinement status of sister blocks
            if (blocks(k)%refinement == -1) then

                ! get sister id's
                call get_sister_id(id1, id2, id3, blocks(k)%treecode, blocks(k)%level)

                ! if all sisters exists
                if ( (id1>0) .and. (id2>0) .and. (id3>0) ) then

                    ! if all sister blocks want to coarsen, then coarsening is allowed
                    if ( (blocks(id1)%refinement == -1) .and. (blocks(id2)%refinement == -1) .and. (blocks(id3)%refinement == -1) ) then
                        blocks(k)%refinement    = -2
                        blocks(id1)%refinement  = -2
                        blocks(id2)%refinement  = -2
                        blocks(id3)%refinement  = -2
                    end if

                end if

            end if

        end if

    end do

end subroutine ensure_completeness
