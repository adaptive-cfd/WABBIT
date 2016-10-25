! ********************************
! WABBIT
! --------------------------------
!
! encoding treecode
!
! name: encoding.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine encoding(treecode, i, j, nx, ny)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                        :: i, j, nx, ny
    integer(kind=ik), dimension(10), intent(out)        :: treecode

    real(kind=rk)                                       :: Jx, Jy
    integer(kind=ik)                                    :: N, k, allocate_error
    integer(kind=ik), dimension(:), allocatable         :: c, d

    Jx = log(dble(nx)) / log(2.0_rk)
    Jy = log(dble(ny)) / log(2.0_rk)

    N = nint(Jx)

    treecode = -1

    allocate( c(N), stat=allocate_error)
    allocate( d(N), stat=allocate_error)

    call int_to_binary(i-1, N, c)
    call int_to_binary(j-1, N, d)

    do k = 1, N
        if (c(N-k+1)==0 .and. d(N-k+1)==0) then
            treecode(k) = 0
        end if
        if (c(N-k+1)==0 .and. d(N-k+1)==1) then
            treecode(k) = 1
        end if
        if (c(N-k+1)==1 .and. d(N-k+1)==0) then
            treecode(k) = 2
        end if
        if (c(N-k+1)==1 .and. d(N-k+1)==1) then
            treecode(k) = 3
        end if
    end do

    deallocate( c, stat=allocate_error)
    deallocate( d, stat=allocate_error)

end subroutine encoding
