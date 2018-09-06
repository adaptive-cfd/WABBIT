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

subroutine write_debug_times( iteration, params )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> iteration
    integer(kind=ik), intent(in)        :: iteration

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    ! process rank
    integer(kind=ik)                    :: rank
    ! loop variable
    integer(kind=ik)                    :: k
    ! file name
    character(len=80)                   :: fname

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    rank = params%rank

!---------------------------------------------------------------------------------------------
! main body

    ! file name
    write( fname,'(i5.5, "times.dat")') rank

    ! we always destroy existing data and re-create the file - those files otherwise
    ! get really big (TB range)
    open(unit=99,file=fname, status='replace')


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
