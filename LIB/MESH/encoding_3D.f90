! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: encoding_3D.f90
! version: 0.5
! author: msr
!
! encoding 3D treecode
!
! input:    - block position varaibles i, j, k
!           - number of blocks
!           - length of treecode in light data
! output:   - treecode
!
! = log ======================================================================================
!
! 26/01/17 - create
!
! ********************************************************************************************

subroutine encoding_3D(treecode, i, j, k, block_num, treeN)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! block position coordinates
    integer(kind=ik), intent(in)    :: i, j, k
    ! number of blocks
    integer(kind=ik), intent(in)    :: block_num

    ! treecode size
    integer(kind=ik), intent(in)    :: treeN
    ! treecode
    integer(kind=ik), intent(out)   :: treecode(treeN)

    ! variables for calculate real treecode length N
    real(kind=rk)                   :: Jn
    integer(kind=ik)                :: N

    ! allocation error variable
    integer(kind=ik)                :: allocate_error

    ! auxiliary vectors
    integer(kind=ik), allocatable   :: c(:), d(:), e(:)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! real treecode length
    !Jn = log(dble(block_num)) / log(2.0_rk)
    Jn = log(dble(block_num)) / log(8.0_rk)
    N  = nint(Jn)

    ! reset output
    treecode = -1

!---------------------------------------------------------------------------------------------
! main body

    ! allocate auxiliary vectors
    allocate( c(N), stat=allocate_error)
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( d(N), stat=allocate_error)
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( e(N), stat=allocate_error)
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if



    ! convert block coordinates into binary numbers
    call int_to_binary(i-1, N, c)
    call int_to_binary(j-1, N, d)
    call int_to_binary(k-1, N, e)

    ! calculate treecode
    treecode(1:N) = e*4 + d*2 + c

    ! clean up
    deallocate( c, stat=allocate_error)
    deallocate( d, stat=allocate_error)
    deallocate( e, stat=allocate_error)

end subroutine encoding_3D
