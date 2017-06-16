!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name write_debug_times.f90
!> \version 0.4
!> \author msr
!
!> \brief write time measurements
!
!>
!! input:    - current iteration \n
!! output:   -                   \n
!!
!!
!! = log ======================================================================================
!! \n
!! 02/12/16 - create
! ********************************************************************************************

subroutine write_debug_times( iteration )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> iteration
    integer(kind=ik), intent(in)        :: iteration

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

    ! file name
    character(len=80)                   :: fname

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

!---------------------------------------------------------------------------------------------
! main body

    ! file name
    write( fname,'(i5.5, "times.dat")') rank

    ! check file existence, if not create file
    inquire(file=fname, exist=file_exists)

    if (file_exists) then
        ! open for append
        open(unit=99,file=fname,status='old', position="append", action='write', iostat=io_error)
    else
        ! first opening
        open(unit=99,file=fname,status='new',action='write', iostat=io_error)
    end if

    ! write data

    ! write file header
    write(99,'(80("_"))')
    write(99, '(42x, "calls", 2x, "sum", 5x, "time", 6x, "sum")', advance='no')
    write(99,*)

    ! write times
    k = 1
    do while ( debug%name_comp_time(k) /= "---" )

        ! write name
        write(99, '(a)', advance='no') debug%name_comp_time(k)
        ! write number of calls
        write(99, '(2x,i5)', advance='no') int(debug%comp_time(k,1))
        ! write global number of calls
        write(99, '(2x,i7)', advance='no') int(debug%comp_time(k,3))
        ! write time
        write(99, '(2x,f12.6)', advance='no') debug%comp_time(k,2)
        ! write global time
        write(99, '(2x,f12.6)', advance='no') debug%comp_time(k,4)
        ! next line
        write(99,*)
        ! loop variable
        k = k + 1

    end do

    write(99,'(80("-"))')
    write(99, '("iteration: ", i7)', advance='no') iteration
    write(99,*)

    ! close file
    close(unit=99)


end subroutine write_debug_times
