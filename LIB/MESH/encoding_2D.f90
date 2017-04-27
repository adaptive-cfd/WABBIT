!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name encoding_2D.f90
!> \version 0.4
!> \author msr
!
!> \brief encoding 2D treecode
!
!> \details input:    
!!                   - block position varaibles i, j
!!                   - number of blocks
!!                   - length of treecode in light data
!!
!!          output: 
!!                   - treecode
!!
!> = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4 \n
!! 26/01/17 - rename to encoding 2D
!
! ********************************************************************************************

subroutine encoding_2D(treecode, i, j, nx, ny, treeN)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> block position coordinates
    integer(kind=ik), intent(in)    :: i, j
    !> number of blocks
    integer(kind=ik), intent(in)    :: nx, ny

    !> treecode size
    integer(kind=ik), intent(in)    :: treeN
    !> treecode
    integer(kind=ik), intent(out)   :: treecode(treeN)

    ! variables for calculate real treecode length N
    real(kind=rk)                   :: Jx, Jy
    integer(kind=ik)                :: N

    ! loop variables
    integer(kind=ik)                :: k

    ! allocation error variable
    integer(kind=ik)                :: allocate_error

    ! auxiliary vectors
    integer(kind=ik), allocatable   :: c(:), d(:)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! real treecode length
    Jx = log(dble(nx)) / log(2.0_rk)
    Jy = log(dble(ny)) / log(2.0_rk)
    N = nint(Jx)

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

    ! convert block coordinates into binary numbers
    call int_to_binary(i-1, N, c)
    call int_to_binary(j-1, N, d)

    ! loop over binary vectors to calculate treecode
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

    ! clean up
    deallocate( c, stat=allocate_error)
    deallocate( d, stat=allocate_error)

end subroutine encoding_2D
