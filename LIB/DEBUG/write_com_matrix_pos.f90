!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name write_com_matrix_pos.f90
!> \version 0.4
!> \author msr
!
!> \brief write position communication matrix to file
!
!> \details
!! input:    - current position com matrix \n
!! output:   -
!! \n
!! = log ======================================================================================
!! \n
!! 16/12/16 - create
! ********************************************************************************************

subroutine write_com_matrix_pos( com_matrix )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> iteration
    integer(kind=ik), intent(in)        :: com_matrix(:, :)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! file existence variable
    logical                             :: file_exists
    ! file IO error variable
    integer(kind=ik)                    :: io_error

    ! loop variable
    integer(kind=ik)                    :: k, l

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

!---------------------------------------------------------------------------------------------
! main body

    ! check file existence, if not create file
    inquire(file="com_matrix_pos.dat", exist=file_exists)

    if ( rank == 0 ) then
        if (file_exists) then
            ! open and overwrite
            open(unit=99,file="com_matrix_pos.dat",status='old', position="append", action='write', iostat=io_error)
        else
            ! first opening
            open(unit=99,file="com_matrix_pos.dat",status='new',action='write', iostat=io_error)
        end if
    end if

    ! write data
    if (rank == 0) then

        do k = 1, size(com_matrix,1)

            do l = 1, size(com_matrix,1)

                write(99, '(i5, 2x)', advance='no') com_matrix(k,l)

            end do

            ! next line
            write(99,*)

        end do

        write(99,'(80("-"))')
        write(99,*)

        ! close file
        close(unit=99)

    end if

end subroutine write_com_matrix_pos
