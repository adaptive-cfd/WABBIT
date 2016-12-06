! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: write_com_list.f90
! version: 0.4
! author: msr
!
! write communication list to file
! note: existing file will be overridden
!
! input:    - current com list
! output:   -
!
! = log ======================================================================================
!
! 05/12/16 - create
! ********************************************************************************************

subroutine write_com_list( com_list )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! iteration
    integer(kind=ik), intent(in)        :: com_list(:, :)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! file existence variable
    logical                             :: file_exists
    ! file IO error variable
    integer(kind=ik)                    :: io_error

    ! loop variable
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

!---------------------------------------------------------------------------------------------
! main body

    ! check file existence, if not create file
    inquire(file="com_list.dat", exist=file_exists)

    if ( rank == 0 ) then
        if (file_exists) then
            ! open and overwrite
            open(unit=99,file="com_list.dat",status='old', action='write', iostat=io_error)
        else
            ! first opening
            open(unit=99,file="com_list.dat",status='new',action='write', iostat=io_error)
        end if
    end if

    ! write data
    k = 1
    if (rank == 0) then

        ! write file header
        write(99, '("rank", 1x, "neighbor-rank", 1x, "block", 1x, "neighbor-block", 1x, "direction")', advance='no')
        write(99,*)

        do while ( com_list(k,1) /= -1 )

            write(99, '(i3, 7x, i3, 8x, i3, 7x, i3, 7x, i3)', advance='no') com_list(k,2), com_list(k,3), com_list(k,4), com_list(k,5), com_list(k,6)
            k = k + 1
            ! next line
            write(99,*)

        end do

        ! close file
        close(unit=99)

    end if

end subroutine write_com_list
