! ********************************
! 2D AMR prototype
! --------------------------------
!
! return free block id
!
! name: get_free_block.f90
! date: 22.08.2016
! author: msr, engels
! version: 0.1
!
! ********************************

subroutine get_free_block(id)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(out)   :: id
    integer(kind=ik)                :: k

    id = -1

    ! loop over all blocks
    do k = 1, blocks_params%number_max_blocks
      ! check if the block is active. if it is not, then we found a free block
      ! to return
      if (blocks(k)%active .eqv. .false.) then
        id = k
        exit
      end if
    end do

    ! error catching: is there no more free blocks on the list?
    if (id == -1) then
      write(*,*) "We try to fetch a free block ID from the list but all blocks are used."
      write(*,*) "Stop."
      stop
    end if

end subroutine get_free_block
