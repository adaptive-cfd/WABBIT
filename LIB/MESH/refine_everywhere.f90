! ********************************
! WABBIT
! --------------------------------
!
! refine every block to create the
! wavelet safety zone
!
! name: refine_everywhere.f90
! date: 30.09.2016
! author: engels, msr
! version: 0.2
!
! ********************************

subroutine refine_everywhere()

    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    integer(kind=ik)    :: k, N

    N           = blocks_params%number_max_blocks

    ! set status "refine" for all active blocks
    do k = 1, N
        if (blocks(k)%active) then
            blocks(k)%refinement = 1
        end if
    end do

    ! check if block has reached maximal level
    call respect_min_max_treelevel()

    ! interpolate the new mesh
    call interpolate_mesh()

end subroutine refine_everywhere
