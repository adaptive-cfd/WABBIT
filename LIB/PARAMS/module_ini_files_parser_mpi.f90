! Ini files parser for MPI parallel codes
! This module is just a layer that includes the serial ini_files_parser
! if you read a parameter, only root rank extracts it from the file (which is also
! read by root only) and then broadcasts it to all other processes
module module_ini_files_parser_mpi
! use the serial ini files parser. the present module is just a wrapper for mpi codes
! since the original ini files parser is serial
use module_ini_files_parser
use mpi


! the generic call "read_param" redirects to these routines, depending on the data
! type and the dimensionality. vectors can be read without setting a default.
interface read_param_mpi
    module procedure param_dbl_mpi, param_int_mpi, param_vct_mpi, param_str_mpi, &
        param_bool_mpi, param_vct_str_mpi, param_boolvct_mpi, param_matrix_mpi
    end interface

    !!!!!!!!
contains
    !!!!!!!!



    !-----------------------------------------------------------------------------
    ! read an array from an ascii file (SERIAL version, to be executed only on root)
    ! note: array is assumed-shape and its size defines what we try to read
    !-----------------------------------------------------------------------------
    subroutine read_array_from_ascii_file_mpi(file, array, n_header)
        implicit none
        character(len=*), intent(in) :: file
        integer, intent(in) :: n_header
        real(kind=rk), intent(inout) :: array (1:,1:)
        integer :: nlines, ncols, mpicode, mpirank

        ! check if communicator is set
        if (WABBIT_COMM==-1) then
            call abort(3567632,"Error[module_ini_files_parser_mpi.f90]: Communicator not set")
        endif
        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        nlines = size(array,1)
        ncols = size(array,2)

        ! only root reads from file...
        if (mpirank==0) call read_array_from_ascii_file(file, array, n_header)
        ! ... then broadcast
        call MPI_BCAST(array,nlines*ncols,MPI_DOUBLE_PRECISION,0,WABBIT_COMM,mpicode)

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.



    end subroutine read_array_from_ascii_file_mpi

    subroutine read_intarray_from_ascii_file_mpi(file, array, n_header)
        implicit none
        character(len=*), intent(in) :: file
        integer, intent(in) :: n_header
        integer(kind=ik), intent(inout) :: array (1:,1:)
        integer :: nlines, ncols, mpicode, mpirank

        ! check if communicator is set
        if (WABBIT_COMM==-1) then
            call abort(3567632,"Error[module_ini_files_parser_mpi.f90]: Communicator not set")
        endif
        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        nlines = size(array,1)
        ncols = size(array,2)

        ! only root reads from file...
        if (mpirank==0) call read_intarray_from_ascii_file(file, array, n_header)
        ! ... then broadcast
        call MPI_BCAST(array,nlines*ncols,MPI_INTEGER4,0,WABBIT_COMM,mpicode)

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine read_intarray_from_ascii_file_mpi


    !-----------------------------------------------------------------------------
    ! count the number of lines in an ascii file, skip n_header lines
    !-----------------------------------------------------------------------------
    subroutine count_lines_in_ascii_file_mpi(file, num_lines, n_header)
        implicit none
        character(len=*), intent(in) :: file
        integer, intent(out) :: num_lines
        integer, intent(in) :: n_header
        integer :: mpicode, mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        ! only root reads from file...
        if (mpirank==0) call count_lines_in_ascii_file(file, num_lines, n_header)
        ! ... then broadcast
        call MPI_BCAST(num_lines,1,MPI_INTEGER,0,WABBIT_COMM,mpicode)

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine count_lines_in_ascii_file_mpi



    !-----------------------------------------------------------------------------
    ! count the number of columns in an ascii file, skip n_header lines
    !-----------------------------------------------------------------------------
    subroutine count_cols_in_ascii_file_mpi(file, num_cols, n_header)
        implicit none
        character(len=*), intent(in) :: file
        integer, intent(out) :: num_cols
        integer, intent(in) :: n_header
        integer :: mpicode, mpirank

        ! fetch my process id
        call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpicode)

        ! only root reads from file...
        if (mpirank==0) call count_cols_in_ascii_file(file, num_cols, n_header)
        ! ... then broadcast
        call MPI_BCAST(num_cols,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode)

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine count_cols_in_ascii_file_mpi

    !-------------------------------------------------------------------------------
    ! clean a previously read ini file, deallocate its string array, and reset
    ! verbosity to .true. (as a matter of precaution)
    !-------------------------------------------------------------------------------
    subroutine clean_ini_file_mpi(PARAMS)
        implicit none
        type(inifile), intent(inout) :: PARAMS
        integer :: mpirank, mpicode

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        if (mpirank == 0) then
            call clean_ini_file(PARAMS)
        endif
    end subroutine clean_ini_file_mpi


    !-------------------------------------------------------------------------------
    ! Read the file paramsfile, count the lines and put the
    ! text in PARAMS.
    !-------------------------------------------------------------------------------
    subroutine read_ini_file_mpi(PARAMS, file, verbose, remove_comments)
        implicit none

        type(inifile), intent(inout) :: PARAMS
        character(len=*) :: file ! this is the file we read the PARAMS from
        logical, intent(in) :: verbose
        logical, optional, intent(in) ::  remove_comments
        integer :: mpirank, mpicode
        logical :: exists

        ! check if communicator is set
        if (WABBIT_COMM==-1) then
            call abort(3567632,"Error[module_ini_files_parser_mpi.f90]: Communicator not set")
        endif

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        if (mpirank==0) then
            inquire ( file=file, exist=exists )

            ! check if the specified file exists
            if (.not. exists) then
              write (*,'("ERROR! file: ",A," not found")') trim(adjustl(file))
              call abort(30302020,"INI file not found")
            endif

            if (present(remove_comments)) then
                call read_ini_file( PARAMS, file, verbose, remove_comments )
            else
                call read_ini_file( PARAMS, file, verbose )
            endif
        endif
    end subroutine read_ini_file_mpi




    !-------------------------------------------------------------------------------
    ! Fetches a REAL VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find the parameter, we return this and warn
    ! Output:
    !       params_real: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_dbl_mpi (PARAMS, section, keyword, params_real, defaultvalue)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        real (kind=rk) :: params_real, defaultvalue
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call read_param (PARAMS, section, keyword, params_real, defaultvalue)
        endif

        ! And then broadcast
        call MPI_BCAST( params_real, 1, MPI_DOUBLE_PRECISION, 0, WABBIT_COMM, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_dbl_mpi




    !-------------------------------------------------------------------------------
    ! Fetches a STRING VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find the parameter, we return this and warn
    ! Output:
    !       params_string: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_str_mpi (PARAMS, section, keyword, params_string, defaultvalue)
        implicit none

        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! what section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        character(len=*), intent (inout) :: params_string
        character(len=*), intent (in) :: defaultvalue
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call read_param (PARAMS, section, keyword, params_string, defaultvalue)
        endif

        ! And then broadcast
        call MPI_BCAST( params_string, len(params_string), MPI_CHARACTER, 0, WABBIT_COMM, mpicode)

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_str_mpi



    !-------------------------------------------------------------------------------
    ! Fetches a VECTOR VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find a vector, we return this and warn
    !       n: length of vector
    ! Output:
    !       params_vector: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_vct_mpi (PARAMS, section, keyword, params_vector, defaultvalue)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        real(kind=rk), intent(inout) :: params_vector(1:)
        real(kind=rk), optional, intent(in) :: defaultvalue(1:)

        integer :: n
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        n = size(params_vector,1)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            if (present(defaultvalue)) then
                call read_param (PARAMS, section, keyword, params_vector, defaultvalue)
            else
                call read_param (PARAMS, section, keyword, params_vector)
            endif
        endif

        ! And then broadcast
        call MPI_BCAST( params_vector, n, MPI_DOUBLE_PRECISION, 0, WABBIT_COMM, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_vct_mpi


    !-------------------------------------------------------------------------------
    ! Fetches a STRING VALUED vector parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find the parameter, we return this and warn
    ! Output:
    !       params_string: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_vct_str_mpi (PARAMS, section, keyword, params_vector, defaultvalue)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout)    :: PARAMS
        character(len=*), intent(in)    :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in)    :: keyword ! what keyword do you look for? for example nx=128
        character(len=*), intent(inout) :: params_vector(1:)
        character(len=*), intent(in)    :: defaultvalue(1:)

        integer :: n
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        n = size(params_vector,1)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call read_param (PARAMS, section, keyword, params_vector, defaultvalue)
        endif

        call MPI_BCAST( params_vector, len(params_vector(1))*n, MPI_CHARACTER, 0, WABBIT_COMM, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_vct_str_mpi



    !-------------------------------------------------------------------------------
    ! Fetches a MATRIX VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find a matrix, we return this and warn
    !       n: length of matrix
    ! Output:
    !       params_matrix: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_matrix_mpi (PARAMS, section, keyword, matrix, defaultvalue)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        real(kind=rk), allocatable, intent(out) :: matrix(:,:)
        real(kind=rk), intent(in)               :: defaultvalue(:,:)

        integer :: n,m
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call read_param (PARAMS, section, keyword, matrix, defaultvalue)
            n = size(matrix,1)
            m = size(matrix,2)
        endif

        ! And then broadcast
        call MPI_BCAST(  n, 1, MPI_INTEGER, 0, WABBIT_COMM, mpicode )
        call MPI_BCAST(  m, 1, MPI_INTEGER, 0, WABBIT_COMM, mpicode )

        if ( .not. allocated(matrix) ) then
            allocate(matrix(1:n,1:m))
        endif

        call MPI_BCAST(  matrix, n*m, MPI_DOUBLE_PRECISION, 0, WABBIT_COMM, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_matrix_mpi

    !-------------------------------------------------------------------------------
    ! Fetches a INTEGER VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find the parameter, we return this and warn
    ! Output:
    !       params_int: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_int_mpi(PARAMS, section, keyword, params_int, defaultvalue)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        integer :: params_int, defaultvalue
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call read_param(PARAMS, section, keyword, params_int, defaultvalue)
        endif

        ! And then broadcast
        call MPI_BCAST( params_int, 1, MPI_INTEGER, 0, WABBIT_COMM, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_int_mpi


    !-------------------------------------------------------------------------------
    ! Fetches a BOOLEAN VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find the parameter, we return this and warn
    ! Output:
    !       params_int: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_bool_mpi(PARAMS, section, keyword, params_bool, defaultvalue)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        logical :: params_bool, defaultvalue
        integer :: mpirank, mpicode

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call read_param(PARAMS, section, keyword, params_bool, defaultvalue)
        endif

        ! And then broadcast
        call MPI_BCAST( params_bool, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_bool_mpi


    !-------------------------------------------------------------------------------
    ! Fetches a VECTOR VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    !       defaultvalue: if the we can't find a vector, we return this and warn
    !       n: length of vector
    ! Output:
    !       params_vector: this is the parameter you were looking for
    !-------------------------------------------------------------------------------
    subroutine param_boolvct_mpi (PARAMS, section, keyword, params_vector, defaultvalue)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        logical, intent(inout)       :: params_vector(1:)
        logical, intent(in):: defaultvalue(1:)

        integer :: n
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        n = size(params_vector,1)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call read_param (PARAMS, section, keyword, params_vector, defaultvalue)
        endif

        ! And then broadcast
        call MPI_BCAST( params_vector, n, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_boolvct_mpi

    !-----------------------------------------------------------------------------
    ! sometimes, it turned out to be useful to provide some values as multiples of
    ! dx, when performing convergence tests (eg thickness=5*dx;) the unit of dx is
    ! set here:
    !-----------------------------------------------------------------------------
    subroutine set_lattice_spacing_mpi(dx_unit)
        implicit none
        real(kind=rk), intent(in) :: dx_unit

        call set_lattice_spacing(dx_unit)
    end subroutine


    !-------------------------------------------------------------------------------
    ! Fetches a MATRIX VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    ! Output:
    !       matrixlines, matrixcols are the dimensions of the matrix
    !
    ! NOTE: Annoyingly, the fujitsu SXF90 compiler cannot handle allocatable arrays
    ! as arguments. so we have to split the routine in one part that returns the size
    ! of the array, then let the caller allocate, then read the matrix. very tedious.
    !
    ! EXAMPLE:
    !   call param_matrix_size_mpi(PARAMS,"Stuff","matrix",a,b)
    !   allocate(matrix(1:a,1:b))
    !   call param_matrix_read_mpi(PARAMS,"Stuff","matrix",matrix)
    !-------------------------------------------------------------------------------
    subroutine param_matrix_size_mpi (PARAMS, section, keyword, matrixlines, matrixcols)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        integer, intent(out) :: matrixcols, matrixlines

        integer :: n,m
        integer :: mpicode
        integer :: mpirank

        ! fetch my process id
        call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpicode)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call param_matrix_size (PARAMS, section, keyword, matrixlines, matrixcols)
        endif

        ! And then broadcast
        call MPI_BCAST(  matrixlines, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpicode )
        call MPI_BCAST(  matrixcols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_matrix_size_mpi

    !-------------------------------------------------------------------------------
    ! Fetches a MATRIX VALUED parameter from the PARAMS.ini file.
    ! Displays what it does on stdout (so you can see whats going on)
    ! Input:
    !       PARAMS: the complete *.ini file
    !       section: the section we're looking for
    !       keyword: the keyword we're looking for
    ! Output:
    !       matrixlines, matrixcols are the dimensions of the matrix
    !
    ! NOTE: Annoyingly, the fujitsu SXF90 compiler cannot handle allocatable arrays
    ! as arguments. so we have to split the routine in one part that returns the size
    ! of the array, then let the caller allocate, then read the matrix. very tedious.
    !
    ! EXAMPLE:
    !   call param_matrix_size_mpi(PARAMS,"Stuff","matrix",a,b)
    !   allocate(matrix(1:a,1:b))
    !   call param_matrix_read_mpi(PARAMS,"Stuff","matrix",matrix)
    !-------------------------------------------------------------------------------
    subroutine param_matrix_read_mpi (PARAMS, section, keyword, matrix)
        implicit none
        ! Contains the ascii-params file
        type(inifile), intent(inout) :: PARAMS
        real(kind=rk), intent(inout) :: matrix(1:,1:)
        character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
        integer :: matrixcols, matrixlines

        integer :: mpicode
        integer :: mpirank

        matrixlines = size(matrix,1)
        matrixcols = size(matrix,2)

        ! fetch my process id
        call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpicode)

        ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
        if (mpirank==0) then
            call param_matrix_read (PARAMS, section, keyword, matrix)
        endif

        ! And then broadcast
        call MPI_BCAST(  matrix, matrixlines*matrixcols, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(WABBIT_COMM, mpicode) ! note this is irrelevant for performance here.

    end subroutine param_matrix_read_mpi

end module
