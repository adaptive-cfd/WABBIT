! ********************************
! 2D AMR prototype
! --------------------------------
!
! check and count block number
!
! name: block_check.f90
! date: 25.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine block_check()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                    :: k, N, block_num
    integer(kind=ik), dimension(10)     :: block_count

    N               = size(blocks_params%active_list, dim=1)
    block_count     = 0

    do k = 1, N
        block_num                               = blocks_params%active_list(k)
        block_count( blocks(block_num)%level )  = block_count( blocks(block_num)%level ) + 1
    end do

    write(*, '(a)', advance='no') "block check: "

    do k = 1, 10
        if (block_count(k) > 0) then
            write(*, '(i4,a,i4,2x)', advance='no') block_count(k), " blocks on level ", k
        end if
    end do

    do k = 10, 1, -1
        if (block_count(k) > 0) then
            if ( modulo(block_count(k), 4) == 0 ) then
                if (k > 1) then
                    block_count(k-1) = block_count(k-1) + block_count(k)/4
                end if
            else
                print*, "error: block inconsistency"
                stop
            end if
        end if
    end do

    write(*,*)
    write(*,'(80("-"))')

end subroutine block_check
