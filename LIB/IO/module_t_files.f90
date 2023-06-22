module module_t_files
    use mpi
    use module_globals
    implicit none

    ! precision statement
    integer, parameter :: flush_frequency = 5
    integer, parameter :: max_parallel_files = 50
    integer, parameter :: max_columns = 160
    integer, save :: mpirank = 7
    ! variables
    real(kind=rk), save, allocatable :: data_buffer(:,:,:)
    character(len=cshort), save, allocatable :: filenames(:)
    integer, save, allocatable :: iteration(:), n_columns(:)
    logical :: disable_all_output = .false.
    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE


    PUBLIC :: init_t_file, append_t_file, close_t_file, close_all_t_files, disable_all_t_files_output

contains

    ! It may sometimes be useful not to have any output generated, for example during
    ! testing and development. calling this function effectively disables the entire module.
    subroutine disable_all_t_files_output()
        implicit none
        disable_all_output = .true.
    end subroutine

    subroutine init_t_file(fname, overwrite, header)
        implicit none
        character(len=*), intent(in) :: fname
        character(len=15), dimension(:), optional, intent(in) :: header
        logical, intent(in) :: overwrite

        integer :: fileID, mpicode
        logical :: exists
        character(len=cshort) :: format

        if (disable_all_output) return

        ! the first called routine will have to set this...
        call MPI_COMM_RANK (MPI_COMM_WORLD, mpirank, mpicode)

        if (mpirank/=0) return

        if (.not. allocated(data_buffer) ) then
            allocate(data_buffer(1:flush_frequency, 1:max_columns, 1:max_parallel_files))
            allocate(iteration(1:max_parallel_files))
            allocate(n_columns(1:max_parallel_files))
            allocate(filenames(1:max_parallel_files))

            data_buffer = 0.0_rk
            iteration = 1
            n_columns = 0
            filenames = "---"
        end if



        ! find a free or the corresponding slot in the array:
        do fileID = 1, max_parallel_files
            if ( filenames(fileID) == "---" .or. filenames(fileID) == fname ) then
                filenames(fileID) = fname
                exit
            endif
        enddo

        if (filenames(fileID) /= fname) then
            write(*,*) "Unfortunately, the buffered *.t files writer cannot handle more than ", max_parallel_files
            write(*,*) "files. You need to increase the number, then compile again."
            call MPI_ABORT(MPI_COMM_WORLD, 2409192, mpicode)
        endif


        ! does the file exist already?
        inquire ( file=fname, exist=exists )

        ! create empty file, write header
        if (.not. exists .or. overwrite) then
            open (55, file=fname, status='replace')
            if (present(header)) then
                ! make format string
                write(format, '(A,i3.3,A)') "('%',", size(header,1) ,"(A15,1x))"
                ! write header
                write(55, format) header
            endif
            close(55)
        endif

        ! now, the file exists, and we have assigned a buffer slot to it
    end subroutine



    subroutine append_t_file(fname, data_line)
        implicit none
        character(len=*), intent(in) :: fname
        real(kind=rk), dimension(:), intent(in) :: data_line

        integer :: fileID, mpicode

        if (disable_all_output) return

        ! the first called routine will have to set this...
        call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpicode)

        ! dump data to buffer, flush if it is time to do so
        if (mpirank/=0) return

        ! if the first call to append_t_file is before the first call to init_t_file
        ! then the arrays are not allocated => code can crash. hence: init_t_file (without overwriting)
        if (.not. allocated(filenames)) call init_t_file(fname, .false.)


        do fileID = 1, max_parallel_files-1
            ! entry for current file exists
            if ( filenames(fileID) == fname ) exit
        end do


        ! not found?
        if (filenames(fileID) /= fname) then
            ! the file is not open for writing: the user forgot to call init_t_file
            ! so do it now.
            call init_t_file(fname, .false.)

            ! .. and find the ID again
            do fileID = 1, max_parallel_files
                ! entry for current file exists
                if ( filenames(fileID) == fname ) exit
            end do
        endif

        if (size(data_line) > max_columns) then
            write(*,*) "Buffered *.t file writer: too many columns in file, increase parameter, compile again."
            call MPI_ABORT(MPI_COMM_WORLD, 2409198, mpicode)
        endif

        data_buffer( iteration(fileID), 1:size(data_line,1), fileID  ) = data_line
        iteration(fileID) = iteration(fileID) + 1
        n_columns(fileID) = size(data_line,1)

        if (iteration(fileID) == flush_frequency+1) then
            call flush_t_file( fileID )
        endif

    end subroutine



    subroutine flush_t_file( fileID )
        implicit none
        integer, intent(in) :: fileID
        character(len=cshort) :: format
        integer :: i

        if (disable_all_output) return

        if (mpirank/=0) return

        if (iteration(fileID)-1 == 0 ) return

        ! write(*,*) "flushing", fileID, trim(adjustl(filenames(fileID))), iteration(fileID)-1
        open(14, file=filenames(fileID), status='unknown', position='append')

        do i = 1, iteration(fileID)-1
            ! make format string
            write(format, '(A,i3.3,A)') "(", n_columns(fileID), "(es15.8,1x))"
            ! write header
            write(14,format) data_buffer( i, 1:n_columns(fileID), fileID  )
        enddo

        close(14)

        iteration(fileID) = 1
    end subroutine



    subroutine close_t_file(fname)
        implicit none
        character(len=*), intent(in) :: fname
        integer :: fileID

        if (disable_all_output) return
        if (mpirank/=0) return

        ! dump data to buffer, flush if it is time to do so

        fileID = 1
        do while (  filenames(fileID) /= "---" )
            ! entry for current subroutine exists
            if (  filenames(fileID) == fname ) exit
            fileID = fileID + 1
        end do

        call flush_t_file( fileID )

    end subroutine



    subroutine close_all_t_files()
        implicit none
        integer :: fileID

        if (disable_all_output) return
        if (mpirank/=0) return

        ! dump data to buffer, flush if it is time to do so

        fileID = 1
        do while (  filenames(fileID) /= "---" )
            call flush_t_file( fileID )
            fileID = fileID + 1
        end do

    end subroutine



end module
