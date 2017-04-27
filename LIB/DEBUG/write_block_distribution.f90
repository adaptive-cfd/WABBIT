!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name write_block_distribution.f90
!> \version 0.4
!> \author msr
!
!> \brief write current block distribution to file
!
!> \details
!! input:    - current distribution list \n
!! output:   -
!! \n
!! = log ======================================================================================
!! \n
!! 05/12/16 - create
! ********************************************************************************************

subroutine write_block_distribution( dist_list )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> iteration
    integer(kind=ik), intent(in)        :: dist_list(:)

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
    inquire(file="block_dist.dat", exist=file_exists)

    if ( rank == 0 ) then
        if (file_exists) then
            ! open for append
            open(unit=99,file="block_dist.dat",status='old', position="append", action='write', iostat=io_error)
        else
            ! first opening
            open(unit=99,file="block_dist.dat",status='new',action='write', iostat=io_error)
        end if
    end if

    ! write data
    if (rank == 0) then

        do k = 1, size(dist_list)

            write(99, '(i3,1x)', advance='no') dist_list(k)

        end do

        ! next line
        write(99,*)

        ! close file
        close(unit=99)

    end if

end subroutine write_block_distribution
