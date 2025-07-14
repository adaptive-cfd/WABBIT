!> \brief Module to handle writing of so-called t-files where array entries can be written to (log-files).
!> For every t-file a buffer area will be created and data values are stored and only flushed with specific frequency
!! in order to avoid excessive file opening and closing
!> Why do we have this module ?
!>
!> On the supercomputers, we are allowed to do very little IO, and that in particular
!! extends to even very small operations, like dumping a log file to disk. This is only a few
!! bytes, but opening and closing the file stresses the machine. Therefore, this module is a buffer
!! for the very simple fortran write statement: We just collect some time steps (flush_frequency)
!! before dumping this data to disk.
module module_t_files
    use mpi
    use module_globals
    implicit none

    ! precision statement
    integer(kind=ik), save, public :: flush_frequency = 5 ! default value may be overwritten in ini_file_to_params
    integer(kind=ik), parameter :: max_parallel_files = 150
    integer(kind=ik), parameter :: max_columns = 160
    integer(kind=ik), save :: mpirank = 7

    ! variables
    real(kind=rk), save, allocatable :: data_buffer(:,:,:)
    character(len=cshort), save, allocatable :: filenames(:)
    integer(kind=ik), save, allocatable :: iteration(:), n_columns(:)
    logical :: disable_all_output = .false.

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE


    PUBLIC :: init_t_file, append_t_file, close_all_t_files, disable_all_t_files_output

contains

    ! It may sometimes be useful not to have any output generated, for example during
    ! testing and development. calling this function effectively disables the entire module.
    subroutine disable_all_t_files_output()
        implicit none
        disable_all_output = .true.
    end subroutine


    !> \brief Initialize a file written during simulation
    !> This initializes the buffer where values for the t_files will be temporarily stored before flushing
    subroutine init_t_file(fname, overwrite, header)
        implicit none
        character(len=*), intent(in) :: fname  !< Name of file, mostly should end in ".t"
        character(len=*), dimension(:), optional, intent(in) :: header  !< Optional array of headers to be written in first line
        logical, intent(in) :: overwrite  !< If file should be cleared if it already exists

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



    !> \brief Appends an array of real values to the given t-file
    !> This file has to be initialized with init_t_file before and afterwards closed appropriately
    !> Values are first stored in an array before being flushed to the file by flush_t_file depending on flush_frequency
    subroutine append_t_file(fname, data_line)
        implicit none
        character(len=*), intent(in) :: fname                   !< name of the file, mostly should end in ".t"
        real(kind=rk), dimension(:), intent(in) :: data_line    !< array of real values to append to the file

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



    !> \brief Flushes temporarily stored values to file
    subroutine flush_t_file( fileID )
        implicit none
        integer, intent(in) :: fileID   !< name of the file, mostly should end in ".t"
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
            ! write data
            write(14,format) data_buffer( i, 1:n_columns(fileID), fileID  )
        enddo

        close(14)

        iteration(fileID) = 1
    end subroutine



    !> \brief Write latest buffered values to files for all files.
    !> This does not actually close a file but rather ensures that the buffers are appropriately flushed
    subroutine close_all_t_files()
        implicit none
        integer :: fileID

        if (disable_all_output) return
        if (mpirank/=0) return

        ! flush all open t files to disk
        do fileID = 1, max_parallel_files
            if (  filenames(fileID) /= "---" ) then
                call flush_t_file( fileID )
            endif       
        enddo
    end subroutine



end module
