!-------------------------------------------------------------------------------
!> \brief FORTRAN ini files parser module
!-------------------------------------------------------------------------------
!> \details
!! This module reads values from standard ini-files.\n
!! The first step is to read in the entire file into a derived datatype. this datatype
!! is provided here. After reading the file, you can call the generic routine
!! "read_param" to extract a value, e.g. parameter=10.0;\n
!! note we expect lines to end on a colon (;)
!! commented lines beginning with # ; ! are ignored.\n
!! The ini files are organized in sections and values:\n
!!              
!!        [Section]
!!        parameter=10;
!-------------------------------------------------------------------------------
!! \details This module is not MPI-aware. use the mpi layer in ini_files_parser_mpi for this
!-------------------------------------------------------------------------------
module module_ini_files_parser
  use module_precision

  ! it sometimes is useful, for codes with equidistant resolution, to specify
  ! values as multiples of grid spacing, mostly for convergence tests and the like
  ! the corresponding spacing would be this:
  real(kind=rk), private, save :: dx
  logical, private, save :: lattice_spacing_set = .false.

  ! maximum width of parameter file. note we have very long lines if we read long
  ! arrays, as it happens for example when we read fourier coefficients for insect
  ! kinematics
  integer, parameter :: maxcolumns=4096

  ! is set to true, we'll produce some output on the screen (for documentation of runs)
  ! the flag is set with the read_ini_file routine
  logical, private, save :: verbosity = .true.

  ! format of read output
  character(len=100), private,save ::FORMAT1='("read [", A,"]",T30,A,T60,"=",TR2,A)'
  ! type definition of an inifile type. it contains an allocatable string array
  ! of maxcolumns length. we will allocate the array, then fill it with what we
  ! read from the file.
  type inifile
    ! string array that contains the text file:
    character(len=maxcolumns),allocatable,dimension(:) :: PARAMS
    ! number of lines in file:
    integer :: nlines
  end type

  ! the generic call "read_param" redirects to these routines, depending on the data
  ! type and the dimensionality. vectors can be read without setting a default.
  interface read_param
    module procedure param_sgl, param_dbl, param_int, param_vct, param_str, param_bool, param_matrix, param_vct_str
  end interface


!!!!!!!!
contains
!!!!!!!!


  !-----------------------------------------------------------------------------
  ! read an array from an ascii file (SERIAL version, to be executed only on root)
  ! note: array is assumed-shape and its size defines what we try to read
  ! --> MPI wrapper in the MPI parser module
  !-----------------------------------------------------------------------------
  subroutine read_array_from_ascii_file(file, array, n_header)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(in) :: n_header
    real(kind=rk), intent(inout) :: array (1:,1:)
    integer :: nlines, ncols, i, io_error
    character(len=maxcolumns) :: dummy
    character(len=16) :: fmt
    character(len=3) :: ncols_str


    nlines = size(array,1)
    ncols = size(array,2)

    write(*,'(80("-"))')
    write(*,'("INFO: reading ",i5," lines with ",i3," colums from ",A)') nlines, ncols, file

    ! set up format string
    write(ncols_str,'(i3.3)') ncols
    fmt = '('//ncols_str//'(es12.4,1x))'

    io_error = 0
    i = 0

    open(unit=14,file=trim(adjustl(file)),action='read',status='old')
    do while (io_error==0)
      ! read a line from file
      read (14,'(A)',iostat=io_error) dummy
      i = i + 1
      ! if we're past the header AND the read worked (i.e. not end of file)
      if (i > n_header .and. io_error==0) then
        read(dummy,*) array(i-n_header,:)
!        write(*,fmt) array(i-n_header,:)
      endif
    enddo
    close (14)

    write(*,'("Done reading.")')
    write(*,'(80("-"))')
  end subroutine read_array_from_ascii_file


  !-----------------------------------------------------------------------------
  ! count the number of lines in an ascii file, skip n_header lines
  ! --> MPI wrapper in the MPI parser module
  !-----------------------------------------------------------------------------
  subroutine count_lines_in_ascii_file(file, num_lines, n_header)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(out) :: num_lines
    integer, intent(in) :: n_header
    integer :: io_error, i
    character(len=maxcolumns) :: dummy


    ! count the lines
    io_error = 0
    i = 0
    open(unit=14,file=trim(adjustl(file)),action='read',status='old')
    do while (io_error==0)
      read (14,'(A)',iostat=io_error) dummy
      if (io_error==0) i = i+1
    enddo
    close (14)
    num_lines = i - n_header
  end subroutine count_lines_in_ascii_file



  !-------------------------------------------------------------------------------
  ! clean a previously read ini file, deallocate its string array, and reset
  ! verbosity to .true. (as a matter of precaution)
  !-------------------------------------------------------------------------------
  subroutine clean_ini_file(PARAMS)
    implicit none
    type(inifile), intent(inout) :: PARAMS

    if (allocated(PARAMS%PARAMS)) deallocate(PARAMS%PARAMS)
    verbosity = .true.
  end subroutine clean_ini_file


  !-----------------------------------------------------------------------------
  ! sometimes, it turned out to be useful to provide some values as multiples of
  ! dx, when performing convergence tests (eg thickness=5*dx;) the unit of dx is
  ! set here:
  !-----------------------------------------------------------------------------
  subroutine set_lattice_spacing(dx_unit)
    implicit none
    real(kind=rk), intent(in) :: dx_unit

    dx = dx_unit
    lattice_spacing_set = .true.
  end subroutine


  !-------------------------------------------------------------------------------
  ! Read the file paramsfile, count the lines and put the
  ! text in PARAMS.
  !-------------------------------------------------------------------------------
  subroutine read_ini_file(PARAMS, file, verbose)
    implicit none

    type(inifile), intent(inout) :: PARAMS
    character(len=*) :: file ! this is the file we read the PARAMS from
    character(len=maxcolumns) :: line
    logical, optional, intent(in) :: verbose

    integer :: io_error, i
    logical :: exists

    ! check if the specified file exists
    inquire ( file=file, exist=exists )
    if ( exists .eqv. .false.) then
      write (*,'("ERROR! file: ",A," not found")') trim(adjustl(file))
      stop
    endif

    ! we set the module-global variable verbosity. if set to false, all routines
    ! will perform their task quietly.
    if (present(verbose)) then
      verbosity = verbose
    else
      verbosity = .true.
    endif

    if (verbosity) then
      write (*,'(A,i3)') " *** info: reading params from "//trim(file)
    endif

    !-----------------------------------------------------------------------------
    ! determine number of lines in file
    !-----------------------------------------------------------------------------
    call count_lines_in_ascii_file( file, PARAMS%nlines, 0)
    allocate( PARAMS%PARAMS(1:PARAMS%nlines) )

    !-----------------------------------------------------------------------------
    ! read actual file
    !-----------------------------------------------------------------------------
    io_error = 0
    i = 1
    open(unit=14,file=trim(adjustl(file)),action='read',status='old')
    do while (io_error==0)
      read (14,'(A)',iostat=io_error) line
      ! if we're not yet at EoF
      if (io_error == 0) then

        !-----------------------------------------------------------------------
        ! Preprocessing of ini file
        !-----------------------------------------------------------------------
        ! remove leading spaces
        line = adjustl(line)

        ! remove everthing after ";", if it occurs ( comments )
        if (index(line,";") /= 0) then
          line( index(line,";")+1:len_trim(line) ) = " "
        endif

        ! remove commented lines completely
        if (line(1:1) == ";" .or. line(1:1) == "!" .or. line(1:1) == "#" .or. line(1:1) == "%") then
          line = " "
        endif

        PARAMS%PARAMS(i) = line
        i = i+1
      endif
    enddo
    close (14)
  end subroutine read_ini_file






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
  subroutine param_dbl (PARAMS, section, keyword, params_real, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns) :: value    ! returns the value
    double precision :: params_real, defaultvalue

    call GetValue(PARAMS, section, keyword, value)

    if (value /= '') then
      if ( index(value,'*dx') == 0 ) then
        !-- value to be read is in absolute form (i.e. we just read the value)
        read (value, *) params_real
        write (value,'(g10.3)') params_real
      else
        !-- the value is given in gridpoints (e.g. thickness=5*dx)
        read (value(1:index(value,'*dx')-1),*) params_real
        params_real = params_real*dx
        write (value,'(g10.3,"(=",g10.3,"*dx)")') params_real, params_real/dx
        if ( lattice_spacing_set .eqv. .false.) then
          write(*,*) "INI FILES PARSER ERROR: you try to read relative values without setting dx first."
          stop
        endif
      endif
    else
      ! no value red, use default value
      write (value,'(g10.3," (THIS IS THE DEFAULT VALUE!)")') defaultvalue
      params_real = defaultvalue
    endif

    ! in verbose mode, inform about what we did
    if (verbosity) then
      write (*,FORMAT1) trim(section),trim(keyword),adjustl(trim(value))
      !write (*,FORMAT1) trim(section),trim(keyword),adjustl(trim(value))
    endif
  end subroutine param_dbl


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
  subroutine param_sgl (PARAMS, section, keyword, params_real, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns) :: value    ! returns the value
    real :: params_real, defaultvalue

    call GetValue(PARAMS, section, keyword, value)

    if (value /= '') then
      if ( index(value,'*dx') == 0 ) then
        !-- value to be read is in absolute form (i.e. we just read the value)
        read (value, *) params_real
        write (value,'(g10.3)') params_real
      else
        !-- the value is given in gridpoints (e.g. thickness=5*dx)
        read (value(1:index(value,'*dx')-1),*) params_real
        params_real = params_real*real(dx,kind=4)
        write (value,'(g10.3,"(=",g10.3,"*dx)")') params_real, params_real/dx

        if ( lattice_spacing_set .eqv. .false.) then
          write(*,*) "INI FILES PARSER ERROR: you try to read relative values without setting dx first."
          stop
        endif

      endif
    else
      ! no value red, use default value
      write (value,'(g10.3," (THIS IS THE DEFAULT VALUE!)")') defaultvalue
      params_real = defaultvalue
    endif

    ! in verbose mode, inform about what we did
    if (verbosity) then
      write (*,FORMAT1) trim(section),trim(keyword),adjustl(trim(value))
    endif
  end subroutine param_sgl

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
  subroutine param_str (PARAMS, section, keyword, params_string, defaultvalue)
    implicit none

    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! what section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns)  value    ! returns the value
    character(len=*), intent (inout) :: params_string
    character(len=*), intent (in) :: defaultvalue

    call GetValue(PARAMS, section, keyword, value)

    if (value .ne. '') then
      params_string = value
    else
      value = trim(adjustl(defaultvalue))//" (THIS IS THE DEFAULT VALUE!)"
      params_string = defaultvalue
    endif

    ! in verbose mode, inform about what we did
    if (verbosity) then
      write (*,FORMAT1)trim(section),trim(keyword),adjustl(trim(value))
    endif
  end subroutine param_str



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
  subroutine param_vct (PARAMS, section, keyword, params_vector, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    real(kind=rk) :: params_vector(1:)
    real(kind=rk), optional, intent(in) :: defaultvalue(1:)

    integer :: n,m
    character(len=maxcolumns) :: value
    character(len=14)::formatstring

    n = size(params_vector,1)
    ! empty vector??
    if (n==0) return

    if ( present(defaultvalue) ) then
      m = size(defaultvalue,1)
      if (n/=m) then
        write(*,*) "error: vector and default value are not of the same length"
      endif
    endif

    write(formatstring,'("(",i2.2,"(g10.3,1x))")') n

    call GetValue(PARAMS, section, keyword, value)

    if (value .ne. '') then
      ! read the n values from the vector string
      read (value, *) params_vector
      write (value,formatstring) params_vector
    else
      if (present(defaultvalue)) then
        ! return default
        write (value,formatstring) defaultvalue
        value = trim(adjustl(value))//" (THIS IS THE DEFAULT VALUE!)"
        params_vector = defaultvalue
      else
        ! return zeros
        params_vector = 0.d0
        write (value,formatstring) params_vector
        value = trim(adjustl(value))//" (RETURNING ZEROS - NO DEFAULT SET!)"
      endif
    endif

    ! in verbose mode, inform about what we did
    if (verbosity) then
      write (*,FORMAT1) trim(section),trim(keyword),adjustl(trim(value))
    endif
  end subroutine param_vct

  !-------------------------------------------------------------------------------
  !! \brief
  !! Fetches a STRING VALUED vector parameter from the PARAMS.ini file.
  !! Displays what it does on stdout (so you can see whats going on)
  !! \details
  !! Input:
  !!       PARAMS: the complete *.ini file
  !!       section: the section we're looking for
  !!       keyword: the keyword we're looking for
  !!       defaultvalue: if the we can't find the parameter, we return this and warn
  !! Output:
  !!       params_string: this is the parameter you were looking for
  !-------------------------------------------------------------------------------
  subroutine param_vct_str (PARAMS, section, keyword, params_vector, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout)    :: PARAMS
    character(len=*), intent(in)    :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in)    :: keyword ! what keyword do you look for? for example nx=128
    character(len=*), intent (inout) :: params_vector(1:)
    character(len=*), intent (in) :: defaultvalue(1:)

    integer :: n,m
    character(len=maxcolumns) :: value
    character(len=14)::formatstring

    n = size(params_vector,1)
    m = size(defaultvalue,1)
    if (n==0) return

    if (n/=m) then
      write(*,*) "error: vector and default value are not of the same length"
    endif

    write(formatstring,'("(",i2.2,"(g8.3,1x))")') n

    call GetValue(PARAMS, section, keyword, value)

    if (value .ne. '') then
      ! read the three values from the vector string
      read (value, *) params_vector
      write (value,formatstring) params_vector
      !params_vector = value
    else
      write (value,formatstring) defaultvalue
      value = trim(adjustl(value))//" (THIS IS THE DEFAULT VALUE!)"
      params_vector = defaultvalue
    endif

    ! in verbose mode, inform about what we did
    if (verbosity) then
      write (*,FORMAT1) trim(section),trim(keyword),adjustl(trim(value))
    endif
  end subroutine param_vct_str



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
  subroutine param_int(PARAMS, section, keyword, params_int, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns) ::  value    ! returns the value
    integer :: params_int, defaultvalue

    call GetValue(PARAMS, section, keyword, value)

    if (value .ne. '') then
      read (value, *) params_int
      write (value,'(i9)') params_int
    else
      write (value,'(i9," (THIS IS THE DEFAULT VALUE!)")') defaultvalue
      params_int = defaultvalue
    endif

    ! in verbose mode, inform about what we did
    if (verbosity) then
      write (*,FORMAT1) trim(section),trim(keyword),adjustl(trim(value))
    endif
  end subroutine param_int


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
  subroutine param_bool(PARAMS, section, keyword, params_bool, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns) ::  value    ! returns the value
    logical :: params_bool
    logical :: defaultvalue

    call GetValue(PARAMS, section, keyword, value)

    if (value .ne. '') then
      select case (value)
      case ("yes","1","true",".true.")
        params_bool = .true.
      case ("no","0","false",".false.")
        params_bool = .false.
      end select

      write (value,'(L)') params_bool
    else
      write (value,'(L," (THIS IS THE DEFAULT VALUE!)")') defaultvalue
      params_bool = defaultvalue
    endif

    ! in verbose mode, inform about what we did
    if (verbosity) then
      write (*,FORMAT1) trim(section),trim(keyword),adjustl(trim(value))
    endif
  end subroutine param_bool


  !-------------------------------------------------------------------------------
  ! Fetches a MATRIX VALUED parameter from the PARAMS.ini file.
  ! Displays what it does on stdout (so you can see whats going on)
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  !       defaultvalue: if the we can't find the parameter, we return this and warn
  ! Output:
  !       params_int: this is the parameter you were looking for
  !-------------------------------------------------------------------------------
  subroutine param_matrix(PARAMS, section, keyword, matrix, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    real(kind=rk), intent(in)    :: defaultvalue(:,:)
    real(kind=rk), allocatable, intent(out) :: matrix(:,:)
    character(len=maxcolumns) ::  value    ! returns the value
    character(len=29)         :: defaultmessage
    integer :: i, j, matrixcols, matrixlines
    integer :: index1, index2
    logical :: foundsection

    foundsection = .false.
    value = ''

    if (allocated(matrix)) then
      write(*,*) "INIFILES: ERROR: matrix already allocated"
      stop
    end if

    !-- loop over the lines of PARAMS.ini file
    do i=1, PARAMS%nlines
      !-- ignore commented lines completely, if first non-blank character is one of #,!,;,%
      if ((PARAMS%PARAMS(i)(1:1).ne.'#').and.(PARAMS%PARAMS(i)(1:1).ne.';').and.&
      (PARAMS%PARAMS(i)(1:1).ne.'!').and.(PARAMS%PARAMS(i)(1:1).ne.'%')) then

      !-- does this line contain the "[section]" statement?
      if (index(PARAMS%PARAMS(i),'['//section//']')==1) then
        ! yes, it does
        foundsection = .true.
      elseif (PARAMS%PARAMS(i)(1:1) == '[') then
        ! we're already at the next section mark, so we left the section we
        ! were looking for again
        foundsection = .false.
      endif

      !-- we're inside the section we want
      if (foundsection) then
        ! for a matrix, prototype is
        ! [SECTION]
        ! keyword=(/1 2 3 4
        ! 4 5 6 7
        ! 9 9 9 9/);

        ! does this line contain the beginning of the keyword we're looking for ?
        if (index(PARAMS%PARAMS(i),keyword//'=(/')==1) then
          ! yes, it does.
          index1 = index(PARAMS%PARAMS(i),'=(/')+3
          index2 = len_trim( PARAMS%PARAMS(i) )
          ! remove spaces between (/     and values
          value = adjustl(PARAMS%PARAMS(i)(index1:index2))

          ! first we check how many columns we have, by looping over the first line
          ! containing the =(/ substring
          matrixcols = 1
          do j = 1, len_trim(value)
            ! count elements in the line by counting the separating spaces
            if ( value(j:j) == " " ) then
              matrixcols = matrixcols + 1
            end if
          enddo


          ! now count lines in array, loop until you find /) substring (or = substring)
          matrixlines = 0
          do j = i, PARAMS%nlines
            matrixlines = matrixlines +1
            if (index(PARAMS%PARAMS(j),"/)") /= 0) then
              ! we found the terminal line of the matrix statement
              exit
            elseif (index(PARAMS%PARAMS(j),"=") /= 0 .and. j>i) then
              ! a = would mean we skipped past the matrix definition to the next variable..
              write(*,*) "INIFILES: ERROR: invalid ini matrix (code 767626201)"
              stop
            end if
          end do

          allocate( matrix(1:matrixlines,1:matrixcols) )


          do j = i, i+matrixlines-1
            if ( j == i ) then
              ! first line
              index1 = index(PARAMS%PARAMS(j),"(/")+2
              index2 = len_trim(PARAMS%PARAMS(j))
            elseif (j == i+matrixlines-1) then
              ! last line
              index1 = 1
              index2 = index(PARAMS%PARAMS(j),"/)")-1
            else
              ! interior lines
              index1 = 1
              index2 = len_trim(PARAMS%PARAMS(j))
            endif
            ! remove leading spaces, then read
            value = adjustl(PARAMS%PARAMS(j)(index1:index2))
            read( value, * ) matrix(j-i+1,:)
          enddo

          exit ! loop over lines
        endif ! found first line
      endif ! found section
    end if
    end do ! loop over lines

  ! section not found, use default value
  if (.not. allocated(matrix)) then
      allocate ( matrix(size(defaultvalue,1),size(defaultvalue,2)))
      !value = trim(adjustl(defaultvalue))//" (THIS IS THE DEFAULT VALUE!)"
      defaultmessage = " (THIS IS THE DEFAULT VALUE!)"
      matrix = defaultvalue
  else
      defaultmessage = ''
  end if

  ! in verbose mode, inform about what we did
  if (verbosity) then
         write(*,'(" read [",A,"]",T30,A,T60,"as Matrix of size ",i4," x ",i4, a)') trim(section), trim(keyword), size(matrix,1) , size(matrix,2), trim(defaultmessage)
  end if

end subroutine param_matrix


  !-------------------------------------------------------------------------------
  ! Extracts a value from the PARAMS.ini file, which is in "section" and
  ! which is named "keyword"
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  ! Output:
  !       value: is a string contaiing everything between '=' and ';'
  !              to be processed further, depending on the expected type
  !              of variable (e.g. you read an integer from this string)
  !-------------------------------------------------------------------------------
  subroutine GetValue (PARAMS, section, keyword, value)
    implicit none

    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword   ! what keyword do you look for? for example nx=128
    character(len=*), intent(inout) :: value   ! returns the value
    integer :: i
    integer :: index1,index2
    logical :: foundsection

    foundsection = .false.
    value = ''

    !-- loop over the lines of PARAMS.ini file
    do i=1, PARAMS%nlines
      !-- ignore commented lines completely, if first non-blank character is one of #,!,;,%
      ! shouldn't happen anymore since inifile-preprocessing while reading should
      ! directly skip these lines..
      if ((PARAMS%PARAMS(i)(1:1).ne.'#').and.(PARAMS%PARAMS(i)(1:1).ne.';').and.&
      (PARAMS%PARAMS(i)(1:1).ne.'!').and.(PARAMS%PARAMS(i)(1:1).ne.'%')) then

      !-- does this line contain the "[section]" statement?
      if (index(PARAMS%PARAMS(i),'['//section//']')==1) then
        ! yes, it does
        foundsection = .true.
      elseif (PARAMS%PARAMS(i)(1:1) == '[') then
        ! we're already at the next section mark, so we leave the section we
        ! were looking for again
        foundsection = .false.
      endif

      !-- we're inside the section we want
      if (foundsection) then
        if (index(PARAMS%PARAMS(i),keyword//'=')==1) then
          if (index(PARAMS%PARAMS(i),';')/=0) then
            ! found "keyword=" in this line, as well as the delimiter ";"
            index1 = index(PARAMS%PARAMS(i),'=')+1
            index2 = index(PARAMS%PARAMS(i),';')-1
            value = adjustl(PARAMS%PARAMS(i)(index1:index2))
            exit ! leave do loop
          else
            ! found "keyword=" in this line, but delimiter ";" is missing.
            ! proceed, but this may fail (e.g. value=999commentary or
            ! value=999\newline fails)
            if (verbosity) then
              write (*,'(A)') "Though found the keyword, I'm unable to find &
                            &value for variable --> "//trim(keyword)//" <-- &
                            &missing delimiter (;)"
            endif
            index1 = index(PARAMS%PARAMS(i),'=')+1
            index2 = index(PARAMS%PARAMS(i),' ')-1
            value = PARAMS%PARAMS(i)(index1:index2)
            if (verbosity) then
              write (*,'(A)') "As you forgot to set ; at the end of your value,&
                          & I try to extract a value nonetheless. Proceed, with&
                          & fingers crossed (there is no guarantee this will work)"
              write (*,'(A)') "Extracted --->"//trim(value)//"<-----"
            endif
            exit
          endif
        endif
      endif

    endif
  enddo
end subroutine getvalue




end module
