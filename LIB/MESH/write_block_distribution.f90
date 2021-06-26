!> \brief write current block distribution to file
!========================================================================================

    subroutine write_block_distribution( params, dist_list, filename )
        implicit none

        integer(kind=ik), intent(in)        :: dist_list(:)   !< iteration
        type (type_params), intent(in)      :: params         !< user defined parameter structure
        character(len=*), intent(in)        :: filename
        logical                             :: file_exists    ! file existence variable
        integer(kind=ik)                    :: io_error, k    ! file IO error variable, counter

        if ( params%rank == 0 ) then
            ! check file existence, if not create file
            inquire(file=filename, exist=file_exists)

            ! if the file is not there, now we create it.
            if (file_exists) then
                ! open for append
                open(unit=99,file=filename,status='old', position="append", action='write', iostat=io_error)
            else
                ! first opening
                open(unit=99,file=filename,status='new',action='write', iostat=io_error)
            end if

            do k = 1, size(dist_list)
                write(99, '(i4,1x)', advance='no') dist_list(k)
            end do
            ! next line
            write(99,*)
            ! close file
            close(unit=99)
        end if

    end subroutine write_block_distribution
!========================================================================================
