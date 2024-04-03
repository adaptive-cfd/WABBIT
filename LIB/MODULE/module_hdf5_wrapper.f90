!> \brief Read a field from a file
!> \note a single file can contain many arrays, if required
!!
!! |INPUT      |                                                                  |
!! |-----------|------------------------------------------------------------------|
!! |filename   |     file to write to/read from                                   |
!! |dsetname   |     datasetname, i.e. the name of the array in the file          |
!! |field      |     actual data (ALLOCATED!!!)                                   |
!! |rared      |     lower bounds of memory portion hold by the CPU               |
!! |rbred      |     upper bounds of memory portion hold by the CPU               |
!!                   \note rared and rbred are 1:3 arrays. if running on one proc
!!                   they are rared=(/0,0,0/) rbred=(/nx-1,ny-1,nz-1/). If the data
!!                   is distributed among procs, each proc has to indicate which
!!                   portion of the array it holds \n
!! |OUTPUT     |                                                                  |
!! |-----------|------------------------------------------------------------------|
!! |field      |     data read from file/written to file                          |
! ********************************************************************************************

module module_hdf5_wrapper
    use hdf5
    use mpi
    use module_params       ! global parameters
    use module_globals

    implicit none

    ! Precision of doubles
    character(len=cshort) :: field_precision = "double" !"single"
    integer(kind=hsize_t), parameter :: max_chunk = 128

    ! interface for writing attributes. an attribute is an object which is attached
    ! to a dataset, in our case a array saved in the file. we put useful information
    ! in attributes, for example the time, resolution and domain size
    ! both routines take vector values input
    interface write_attribute
        module procedure write_attrib_dble, write_attrib_int, write_attrib_str
    end interface

    ! we can also read attributes from existing files and datasets. thus, for example
    ! when reading a field from file, we check the attribute nxyz for the size,
    ! so we know how much memory to allocate
    interface read_attribute
        module procedure read_attrib_dble, read_attrib_dble2, read_attrib_int, read_attrib_int2
    end interface read_attribute

    interface read_dset_mpi_hdf5
      module procedure read_int_dset_mpi_hdf5_1D, read_dble_dset_mpi_hdf5_1D, &
        read_int_dset_mpi_hdf5_2D, read_dble_dset_mpi_hdf5_2D, &
        read_int_dset_mpi_hdf5_3D, read_dble_dset_mpi_hdf5_3D, &
        read_int_dset_mpi_hdf5_4D, read_dble_dset_mpi_hdf5_4D, &
        read_long_dset_mpi_hdf5_1D
    end interface read_dset_mpi_hdf5

    interface write_dset_mpi_hdf5
      module procedure write_int_dset_mpi_hdf5_1D, write_dble_dset_mpi_hdf5_1D, &
        write_int_dset_mpi_hdf5_2D, write_dble_dset_mpi_hdf5_2D, &
        write_int_dset_mpi_hdf5_3D, write_dble_dset_mpi_hdf5_3D, &
        write_int_dset_mpi_hdf5_4D, write_dble_dset_mpi_hdf5_4D, &
        write_long_dset_mpi_hdf5_1D
    end interface write_dset_mpi_hdf5

contains

  !-----------------------------------------------------------------------------
  ! open an hdf5 file, return handle file_id
  !-----------------------------------------------------------------------------
  subroutine open_file_hdf5(filename, file_id, ovrwrte)
    implicit none
    integer(hid_t), intent(out) ::  file_id
    character(len=*), intent (in) :: filename
    logical, intent(in) :: ovrwrte

    integer(hid_t) :: plist_id
    integer  :: error
    integer :: mpirank, mpicode
    logical :: exist
    
    ! Initialize HDF5 library and Fortran interfaces.
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    ! this sets up a property list (plist_id) with standard values for
    ! FILE_ACCESS
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    ! Modify the property list and store the MPI IO comminucator
    ! information in the file access property list
    call h5pset_fapl_mpio_f(plist_id, WABBIT_COMM, MPI_INFO_NULL, error)

    !---------------------------------------------------------------------------
    ! open the file
    !---------------------------------------------------------------------------
    ! check if the file already exists
    call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)
    if (mpirank == 0) then
      inquire ( file=filename, exist=exist )
    endif
    call MPI_BCAST( exist, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )


    if ((exist .eqv. .false. ) .or. (ovrwrte .eqv. .true.) ) then
        ! file does not exist, create the file collectively
        call h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    else
        ! file does exist, open it. if the dataset we want to write exists already
        ! it will be overwritten. however, if other datasets are present, they will not
        ! be erased
        call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F , file_id, error, access_prp = plist_id)
    endif

    ! this closes the property list plist_id (we'll re-use it)
    call h5pclose_f(plist_id, error)
  end subroutine open_file_hdf5

  !-----------------------------------------------------------------------------
  ! close a file handle.
  !-----------------------------------------------------------------------------
  subroutine close_file_hdf5(file_id)
    implicit none
    integer(hid_t),intent(inout) ::  file_id

    integer :: error
    !!! Close dataspaces:
    call h5fclose_f(file_id, error) ! Close the file.
    call h5close_f(error) ! Close Fortran interfaces and HDF5 library.
  end subroutine close_file_hdf5

  !-----------------------------------------------------------------------------
  ! open an hdf5 file, return handle file_id
  !-----------------------------------------------------------------------------
  subroutine open_file_hdf5_serial(filename, file_id, ovrwrte)
      implicit none
      integer(hid_t),intent(out) ::  file_id
      character(len=*), intent (in) :: filename
      logical, intent(in) :: ovrwrte

      integer(hid_t) :: plist_id
      integer :: error
      logical :: exist

      ! Initialize HDF5 library and Fortran interfaces.
      call h5open_f(error)

      call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F , file_id, error)

  end subroutine open_file_hdf5_serial


  !-----------------------------------------------------------------------------
  ! close a file handle.
  !-----------------------------------------------------------------------------
  subroutine close_file_hdf5_serial(file_id)
      implicit none
      integer(hid_t),intent(inout) ::  file_id

      integer :: error
      !!! Close dataspaces:
      call h5fclose_f(file_id, error) ! Close the file.
      call h5close_f(error) ! Close Fortran interfaces and HDF5 library.
  end subroutine close_file_hdf5_serial


  subroutine get_size_datafield(datarank, file_id, dsetname, dims_file)
    implicit none
    integer, intent(in)                       :: datarank ! Data dimensionality

    integer(hid_t), intent(in)                :: file_id
    character(len=*),intent(in)               :: dsetname
    integer(hsize_t), dimension(datarank)     :: dims_file, dims_dummy
    integer(hid_t)                            :: dset_id       ! dataset identifier
    integer(hid_t)                            :: filespace     ! dataspace identifier in file
    integer                                   :: error  ! error flags

    ! Open an existing dataset.
    call h5dopen_f(file_id, dsetname, dset_id, error)
    ! get its dataspace
    call h5dget_space_f(dset_id, filespace, error)
    ! get the dimensions of the field in the file
    call h5sget_simple_extent_dims_f(filespace, dims_file, dims_dummy, error)
    ! Close dataset
    call h5dclose_f(dset_id, error)
  end subroutine get_size_datafield


  subroutine get_rank_datafield(file_id, dsetname, datarank)
    implicit none
    ! returns dimensionality (datarank) of the array in the file (e.g. 4D array)
    integer, intent(out)             :: datarank ! Data dimensionality
    integer(hid_t), intent(in)       :: file_id
    character(len=*),intent(in)      :: dsetname
    integer(hid_t)                   :: dset_id       ! dataset identifier
    integer(hid_t)                   :: filespace     ! dataspace identifier in file
    integer                          :: error  ! error flags

    ! Open an existing dataset.
    call h5dopen_f(file_id, dsetname, dset_id, error)
    ! get its dataspace
    call h5dget_space_f(dset_id, filespace, error)
    ! get the dimensions of the field in the file
    call h5sget_simple_extent_ndims_f(filespace, datarank, error)
    ! Close dataset
    call h5dclose_f(dset_id, error)
  end subroutine get_rank_datafield


  ! This function call is massive, but since every read or write of dataset
  ! does the same thing this ensures it works for all
  !> \brief Common preparations for dataset read or writes with HDF5
  subroutine prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, memspace, plist_id, &
    file_precision, dims_global, dims_local, file_mode)

    implicit none
    integer, intent(in)                     :: datarank         !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU


    integer(hid_t), intent(out) :: dset_id   !< dataset identifier
    integer(hid_t), intent(out) :: filespace !< dataspace identifier in file
    integer(hid_t), intent(out) :: memspace  !< dataspace identifier in memory
    integer(hid_t), intent(out) :: plist_id  !< property list identifier

    integer(hid_t), intent(in)  :: file_precision  !< file precision for saving
    !> Prepare for reading or writing: "r" - read, "w" - write, "f" - force-write
    character, intent(in)       :: file_mode

    integer(hsize_t), dimension(datarank)              :: dims_file, dims_dummy
    integer(hsize_t), dimension(datarank), intent(out) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank), intent(out) :: dims_local  !< hyperslab dimensions
    integer(hsize_t), dimension(datarank)              :: chunk_dims  !< chunk dimensions
    integer(hsize_t), dimension(datarank)              :: h5_count  !< how many blocks to select from dataspace
    integer(hssize_t), dimension(datarank)             :: h5_offset  !< offset, set after lbounds
    !> stride is spacing between elements, this is 1 here. striding is done in the
    !! caller; here, we just write the entire (possibly downsampled) field to disk.
    integer(hsize_t), dimension(datarank)              :: h5_stride
    
    ! temporary variables
    integer :: error  ! error flag
    integer :: i_dim, mindim, maxdim
    logical :: exist  ! for file checking

    ! debugging: write what is written or saved
    ! write(*, '(a, " : ", a, " in dimension ", i0)') file_mode, dsetname, datarank

    ! ----------------------------------------------------------------------------
    ! Initializes the size of the data as well as hdf5 write settings
    ! ----------------------------------------------------------------------------
    do i_dim = 1, datarank
      ! ----------------------------------------------------------------------------
      ! Compute the dimension of the complete field (i.e. the union of all CPU's)
      ! which we will write to file.
      ! ----------------------------------------------------------------------------
      call MPI_ALLREDUCE ( lbounds(i_dim),mindim,1,MPI_INTEGER,MPI_MIN,WABBIT_COMM,error)
      call MPI_ALLREDUCE ( ubounds(i_dim),maxdim,1,MPI_INTEGER,MPI_MAX,WABBIT_COMM,error)
      ! size of the global array
      dims_global(i_dim) = int( maxdim-mindim+1, kind=hsize_t )
      ! size of array on this cpu
      dims_local(i_dim) = ubounds(i_dim)-lbounds(i_dim) + 1

      !-----------------------------------------------------------------------------
      ! chunking.
      ! HDF writes "chunks" of data at once, and their size can be used for tuning
      ! i/o performance. For example, when chunks fit the cache of the HDD, the performance
      ! may be better.
      !-----------------------------------------------------------------------------
      ! Each process knows how much data it has and where to store it.
      ! now, define the dataset chunking. Chunking is largest dimension in
      ! each direction, but no more than 128 points (so biggest possible chunk is 128^3
      ! which is about 16MB)
      call MPI_ALLREDUCE ( dims_local(i_dim),chunk_dims(i_dim),1,MPI_INTEGER8,MPI_MAX,WABBIT_COMM,error)
      chunk_dims(i_dim) = min(chunk_dims(i_dim), max_chunk )

      ! Tell HDF5 how our data is organized:
      h5_offset(i_dim) = lbounds(i_dim)
      h5_count(i_dim) = 1
      h5_stride(i_dim) = 1
    enddo

    !-----------------------------------------------------------------------------
    ! create dataspace "filespace" to write to
    !-----------------------------------------------------------------------------
    ! Create the data space for the  dataset.
    ! Dataspace in the file: contains all data from all procs
    call h5screate_simple_f(datarank, dims_global, filespace, error)
    !-----------------------------------------------------------------------------
    ! create dataspace "memspace" to be written
    !-----------------------------------------------------------------------------
    ! dataspace in memory: contains only local data
    call h5screate_simple_f(datarank, dims_local, memspace, error)

    ! Create chunked dataset.
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call h5pset_chunk_f(plist_id, datarank, chunk_dims, error)

    !-----------------------------------------------------------------------------
    ! Check about dataset existence
    !-----------------------------------------------------------------------------
    call h5lexists_f(file_id, dsetname, exist, error)
    if (exist .and. file_mode == "w") then
      write(*,'("Write_hdf5 error: You are trying to write to an existing dataset...this is not supported.")')
      call MPI_ABORT(WABBIT_COMM,4441,error)
    elseif (.not. exist .and. file_mode == "r") then
      write(*,'("Read_hdf5 error: You are trying to read a dataset that does not exist: ", a)') dsetname
      call MPI_ABORT(WABBIT_COMM,4441,error)
    endif

    !-----------------------------------------------------------------------------
    ! Open or create dataset
    !-----------------------------------------------------------------------------
    if (file_mode == "r") then
      ! Open an existing dataset.
      call h5dopen_f(file_id, dsetname, dset_id, error)
      call h5dget_space_f(dset_id, filespace, error)
    
      ! get the dimensions of the field in the file and check if they match
      call h5sget_simple_extent_dims_f(filespace, dims_file, dims_dummy, error)
      if ( any(dims_global(1:datarank)/=dims_file(1:datarank)) ) then
        write(*,'("Read_hdf5 error: File dimensions do not match.")')
        write(*,'("Anticipated size by simulation: ", 5(i0, 1x))') dims_global
        write(*,'("Actual size in file: ", 5(i0, 1x))') dims_file
        call MPI_ABORT(WABBIT_COMM, 100043, error)
      endif
    elseif (file_mode == "w" .or. file_mode == "f") then
      ! create the dataset
      call h5dcreate_f(file_id, dsetname, file_precision, filespace, dset_id, error, plist_id)
    
      call h5sclose_f(filespace, error)
      call h5dget_space_f(dset_id, filespace, error)
    endif

    ! Select hyperslab in the file.
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, h5_offset, h5_count, &
    error, h5_stride, dims_local)

    ! Create property list for collective dataset read or write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  end subroutine


  !-------------------------------------------------------------------------------
  !> \brief Read 1D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_int_dset_mpi_hdf5_1D(file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 1     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_INTEGER, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_int_dset_mpi_hdf5_1D


  !-------------------------------------------------------------------------------
  !> \brief Read 2D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_int_dset_mpi_hdf5_2D(file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 2     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_INTEGER, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_int_dset_mpi_hdf5_2D


  !-------------------------------------------------------------------------------
  !> \brief Read 3D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_int_dset_mpi_hdf5_3D(file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 3     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_INTEGER, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_int_dset_mpi_hdf5_3D


  !-------------------------------------------------------------------------------
  !> \brief Read 4D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_int_dset_mpi_hdf5_4D(file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 4     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3),lbounds(4):ubounds(4))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_INTEGER, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_int_dset_mpi_hdf5_4D


  !-------------------------------------------------------------------------------
  !> \brief Read 1D long integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_long_dset_mpi_hdf5_1D(file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 1     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=tsize), intent(inout)      :: field(lbounds(1):ubounds(1))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_STD_I64LE, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_STD_I64LE, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_long_dset_mpi_hdf5_1D


  !-------------------------------------------------------------------------------
  !> \brief Read 1D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_dble_dset_mpi_hdf5_1D( file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 1     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_dble_dset_mpi_hdf5_1D


  !-------------------------------------------------------------------------------
  !> \brief Read 2D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_dble_dset_mpi_hdf5_2D( file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 2     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_dble_dset_mpi_hdf5_2D


  !-------------------------------------------------------------------------------
  !> \brief Read 3D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_dble_dset_mpi_hdf5_3D( file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 3     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_dble_dset_mpi_hdf5_3D


  !-------------------------------------------------------------------------------
  !> \brief Read 4D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !! \attention Ensure field is already allocated
  !-------------------------------------------------------------------------------
  subroutine read_dble_dset_mpi_hdf5_4D( file_id, dsetname, lbounds, ubounds, field )
    implicit none
    integer, parameter                      :: datarank = 4     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3),lbounds(4):ubounds(4))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "r")

    ! Everything is prepared - we can read!
    ! Q: Why dims_local here?
    call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dims_local, error, &
    mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine read_dble_dset_mpi_hdf5_4D


  !-------------------------------------------------------------------------------
  !> \brief Write 1D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_int_dset_mpi_hdf5_1D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 1     !< data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_int_dset_mpi_hdf5_1D


  !-------------------------------------------------------------------------------
  !> \brief Write 2D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_int_dset_mpi_hdf5_2D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 2     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_int_dset_mpi_hdf5_2D


  !-------------------------------------------------------------------------------
  !> \brief Write 3D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_int_dset_mpi_hdf5_3D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 3     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_int_dset_mpi_hdf5_3D


  !-------------------------------------------------------------------------------
  !> \brief Write 4D integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_int_dset_mpi_hdf5_4D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 4     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=ik), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3),lbounds(4):ubounds(4))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_INTEGER, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_int_dset_mpi_hdf5_4D


  !-------------------------------------------------------------------------------
  !> \brief Write 1D long integer array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_long_dset_mpi_hdf5_1D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 1     !< data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    integer(kind=tsize), intent(inout)            :: field(lbounds(1):ubounds(1))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_STD_I64LE, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_STD_I64LE, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_long_dset_mpi_hdf5_1D


  !-------------------------------------------------------------------------------
  !> \brief Write 1D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_dble_dset_mpi_hdf5_1D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 1     !< data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_dble_dset_mpi_hdf5_1D


  !-------------------------------------------------------------------------------
  !> \brief Write 2D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_dble_dset_mpi_hdf5_2D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 2     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_dble_dset_mpi_hdf5_2D


  !-------------------------------------------------------------------------------
  !> \brief Write 3D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_dble_dset_mpi_hdf5_3D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 3     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_dble_dset_mpi_hdf5_3D


  !-------------------------------------------------------------------------------
  !> \brief Write 4D double array to dataset of a HDF5 file
  !> \note a single file can contain many arrays, if required
  !!
  !> `lbounds` and `ubounds` are arays of size 3/4. If running on n proc they are
  !! `lbounds=(0,0,0,0)`, `ubounds=(nx-1,ny-1,nz-1,n)`. If the data is distributed
  !!  among procs, each proc has to indicate which portion of the array it holds.
  !-------------------------------------------------------------------------------
  subroutine write_dble_dset_mpi_hdf5_4D(file_id, dsetname, lbounds, ubounds, field)
    implicit none
    integer, parameter                      :: datarank = 4     !< Data dimensionality
    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer,dimension(datarank), intent(in) :: lbounds, ubounds !< Bounds of memory portion hold by the CPU
    !> Field data to be read
    real(kind=rk), intent(inout)            :: field(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3),lbounds(4):ubounds(4))

    integer(hid_t) :: dset_id       !< dataset identifier
    integer(hid_t) :: filespace     !< dataspace identifier in file
    integer(hid_t) :: memspace      !< dataspace identifier in memory
    integer(hid_t) :: plist_id      !< property list identifier

    integer(hsize_t), dimension(datarank) :: dims_global  !< dataset dimensions in the file.
    integer(hsize_t), dimension(datarank) :: dims_local   !< hyperslab dimensions
    integer :: error  ! error flags

    !----------------------------------------------------------------------------
    ! init size of the data and hdf5 write settings
    ! init hdf5 dataset
    !----------------------------------------------------------------------------
    call prepare_dset_hdf5(file_id, dsetname, lbounds, ubounds, datarank, dset_id, filespace, &
      memspace, plist_id, H5T_NATIVE_DOUBLE, dims_global, dims_local, "w")

    !-----------------------------------------------------------------------------
    ! actual writing of heavy data
    !-----------------------------------------------------------------------------
    ! Write the dataset collectively, double precision in memory
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dims_global, &
    error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

    ! Make sure to close everythig properly (every create or open)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error) ! Close properties
    call h5dclose_f(dset_id, error)  ! Close dataset
  end subroutine write_dble_dset_mpi_hdf5_4D


  ! write an array of strings (of length slen)  to a separate dataset in an HDF5
  ! file. used to backup parameter files.
  subroutine write_string_dset_hdf5(file_id, dsetname, data, slen)
      ! use h5lt
      implicit none

      integer(hid_t), intent(in) :: file_id
      integer, intent(in) :: slen
      character(len=*), intent(in) :: dsetname
      character(len=slen), DIMENSION(:), intent(in) :: data

      INTEGER(HID_T)  :: c_type, dset_id, dspace_id
      integer :: hdferr
      integer(HSIZE_T), dimension(:), allocatable  :: dims
      integer(SIZE_T):: s_dim


      allocate(dims(1))
      dims  = size(data,1)
      s_dim = slen

      call h5tcopy_f(H5T_FORTRAN_S1, c_type, hdferr)

      call h5tset_size_f(c_type, s_dim, hdferr)

      call h5screate_simple_f(1, dims, dspace_id, hdferr)

      call h5dcreate_f(file_id, dsetname, c_type, dspace_id, dset_id, hdferr)

      call h5dwrite_f(dset_id, c_type, data, dims, hdferr)

      call h5tclose_f(c_type, hdferr)
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(dspace_id, hdferr)
  end subroutine


  !> \brief Check if a dataset is present in the data
  function dset_exists(file_id, dsetname)
    implicit none

    integer(hid_t), intent(in)              :: file_id          !< File identifier
    character(len=*),intent(in)             :: dsetname         !< Name of set
    integer :: hdferr

    logical   :: dset_exists  !< Output logical
    !-----------------------------------------------------------------------------
    ! Check about dataset existence
    !-----------------------------------------------------------------------------
    call h5lexists_f(file_id, dsetname, dset_exists, hdferr)
  end function


  !-------------------------------------------------------------------------------
  ! write an attribute
  ! INPUT:
  !   filename  what file to write to (e.g. hallo.h5)
  !   dsetname  what dataset to write to (e.g. stuff to append to hallo.h5:stuff)
  !   aname     the name of the attribute to write
  !   attribute the vector that will hold the attribute. note: this routine uses
  !             assumed shaped arrays: it will try to read as many values as the size of the vectors
  ! OUTPUT:
  !   none
  !-------------------------------------------------------------------------------
  subroutine write_attrib_dble(file_id,dsetname,aname,attribute)
    implicit none

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent (in) :: dsetname, aname
    real(kind=rk), DIMENSION(:), intent (in) :: attribute

    integer, parameter :: arank = 1
    integer :: dim
    integer :: error  ! error flags
    integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
    integer(hid_t) :: attr_id   ! Attribute identifier
    integer(hid_t) :: dset_id  ! dataset identifier
    integer(hsize_t) :: adims(1)  ! Attribute dimension
    logical :: exists

    ! convert input data for the attribute to the precision required by the HDF library
    dim = size(attribute)
    adims = int(dim, kind=hsize_t)

    ! open the dataset
    call h5dopen_f(file_id, dsetname, dset_id, error)

    ! check if attribute exists already
    call h5aexists_f(dset_id, aname, exists, error)

    if (exists) then
      ! open attribute (it exists already)
      call h5aopen_f(dset_id, aname, attr_id, error)
      ! Get dataspace
      call h5aget_space_f(attr_id, aspace_id, error)
      ! Write the attribute data attribute to the attribute identifier attr_id.
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attribute, adims, error)
    else
      ! Determine the dataspace identifier aspace_id
      call h5screate_simple_f(arank,adims,aspace_id,error)
      ! set attr_id, ie create an attribute attached to the object dset_id
      call h5acreate_f(dset_id,aname,H5T_NATIVE_DOUBLE,aspace_id,attr_id,error)
      ! Write the attribute data attribute to the attribute identifier attr_id.
      call h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,attribute,adims,error)
    endif

    call h5aclose_f(attr_id,error) ! Close the attribute.
    call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
    call h5dclose_f(dset_id,error)

  end subroutine write_attrib_dble

  !-------------------------------------------------------------------------------
  ! write an attribute
  ! INPUT:
  !   filename  what file to write to (e.g. hallo.h5)
  !   dsetname  what dataset to write to (e.g. stuff to append to hallo.h5:stuff)
  !   aname     the name of the attribute to write
  !   attribute the vector that will hold the attribute. note: this routine uses
  !             assumed shaped arrays: it will try to read as many values as the size of the vectors
  ! OUTPUT:
  !   none
  !-------------------------------------------------------------------------------
  subroutine write_attrib_int(file_id,dsetname,aname,attribute)

    implicit none

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent (in) :: dsetname, aname
    integer(kind=ik), DIMENSION(:), intent (in) :: attribute

    integer(hsize_t) :: adims(1)  ! Attribute dimension
    integer, parameter :: arank = 1
    integer :: dim
    integer :: error  ! error flags
    integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
    integer(hid_t) :: attr_id   ! Attribute identifier
    integer(hid_t) :: dset_id  ! dataset identifier
    logical :: exists


    ! convert input data for the attribute to the precision required by the HDF library
    dim = size(attribute)
    adims = int(dim, kind=hsize_t)

    ! open the dataset
    call h5dopen_f(file_id, dsetname, dset_id, error)

    ! check if attribute exists already
    call h5aexists_f(dset_id, aname, exists, error)

    if (exists) then
      ! open attribute (it exists already)
      call h5aopen_f(dset_id, aname, attr_id, error)
      ! Get dataspace
      call h5aget_space_f(attr_id, aspace_id, error)
      ! Write the attribute data attribute to the attribute identifier attr_id.
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attribute, adims, error)
    else
      ! Determine the dataspace identifier aspace_id
      call h5screate_simple_f(arank,adims,aspace_id,error)
      ! set attr_id, ie create an attribute attached to the object dset_id
      call h5acreate_f(dset_id,aname,H5T_NATIVE_INTEGER,aspace_id,attr_id,error)
      ! Write the attribute data attribute to the attribute identifier attr_id.
      call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attribute,adims,error)
    endif

    call h5aclose_f(attr_id,error) ! Close the attribute.
    call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
    call h5dclose_f(dset_id,error)

  end subroutine write_attrib_int


  subroutine write_attrib_str(file_id, dsetname, aname, attribute)

    implicit none

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname, aname
    character(len=*), intent(in) :: attribute

    integer(hsize_t) :: adims(1)  ! Attribute dimension
    integer, parameter :: arank = 1
    integer :: dim
    integer :: error  ! error flags
    integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
    integer(hid_t) :: attr_id   ! Attribute identifier
    integer(hid_t) :: dset_id  ! dataset identifier
    logical :: exists


    adims = 1

    ! open the dataset
    call h5dopen_f(file_id, dsetname, dset_id, error)

    ! check if attribute exists already
    call h5aexists_f(dset_id, aname, exists, error)

    if (exists) then
      ! open attribute (it exists already)
      call h5aopen_f(dset_id, aname, attr_id, error)
      ! Get dataspace
      call h5aget_space_f(attr_id, aspace_id, error)
      ! Write the attribute data attribute to the attribute identifier attr_id.
      call h5awrite_f(attr_id, H5T_FORTRAN_S1, attribute, adims, error)
    else
      ! Determine the dataspace identifier aspace_id
      call h5screate_simple_f(arank,adims,aspace_id,error)
      ! set attr_id, ie create an attribute attached to the object dset_id
      call h5acreate_f(dset_id,aname,H5T_FORTRAN_S1,aspace_id,attr_id,error)
      ! Write the attribute data attribute to the attribute identifier attr_id.
      call h5awrite_f(attr_id,H5T_FORTRAN_S1,attribute,adims,error)
    endif

    call h5aclose_f(attr_id,error) ! Close the attribute.
    call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
    call h5dclose_f(dset_id,error)

  end subroutine write_attrib_str



  !-------------------------------------------------------------------------------
  ! Read an attribute
  ! INPUT:
  !   filename  what file to read from (e.g. hallo.h5)
  !   dsetname  what dataset to read from within the file (e.g. stuff to read hallo.h5:stuff)
  !   aname     the name of the attribute to read
  !   attribute the vector that will hold the attribute. note: this routine uses
  !             assumed shaped arrays: it will try to read as many values as the size of the vectors
  ! OUTPUT:
  !   attribute the values read from file
  !-------------------------------------------------------------------------------
  subroutine read_attrib_dble(file_id,dsetname,aname,attribute)
    implicit none

    integer(hid_t), intent(in)                  :: file_id
    character(len=*), intent (in)               :: dsetname, aname
    real(kind=rk), DIMENSION(:), intent (inout) :: attribute

    integer            :: dim
    integer            :: error  ! error flags
    integer(hid_t)     :: aspace_id ! Attribute Dataspace identifier
    integer(hid_t)     :: attr_id   ! Attribute identifier
    integer(hid_t)     :: dset_id  ! dataset identifier
    integer(hsize_t)   :: adims(1)  ! Attribute dimension
    logical            :: exists

    ! convert input data for the attribute to the precision required by the HDF library
    dim = size(attribute)
    adims = int(dim, kind=hsize_t)

    ! open the dataset
    call h5dopen_f(file_id, dsetname, dset_id, error)

    ! check if attribute exists
    call h5aexists_f(dset_id, aname, exists, error)

    if (exists) then
      ! open attribute
      call h5aopen_f(dset_id, aname, attr_id, error)
      ! Get dataspace for attribute
      call h5aget_space_f(attr_id, aspace_id, error)
      ! read attribute data
      call h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attribute, adims, error)
      ! close attribute
      call h5aclose_f(attr_id,error) ! Close the attribute.
      call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
    else
      attribute = 0.0_rk
    endif

    call h5dclose_f(dset_id,error)

  end subroutine read_attrib_dble

  subroutine read_attrib_dble2(file_id, dsetname, aname, attribute)
    implicit none

    integer(hid_t), intent(in)                  :: file_id
    character(len=*), intent (in)               :: dsetname, aname
    real(kind=rk), intent (inout) :: attribute

    integer            :: dim
    integer            :: error  ! error flags
    integer(hid_t)     :: aspace_id ! Attribute Dataspace identifier
    integer(hid_t)     :: attr_id   ! Attribute identifier
    integer(hid_t)     :: dset_id  ! dataset identifier
    integer(hsize_t)   :: adims(1)  ! Attribute dimension
    logical            :: exists

    ! convert input data for the attribute to the precision required by the HDF library
    dim = 1
    adims = int(dim, kind=hsize_t)

    ! open the dataset
    call h5dopen_f(file_id, dsetname, dset_id, error)

    ! check if attribute exists
    call h5aexists_f(dset_id, aname, exists, error)

    if (exists) then
      ! open attribute
      call h5aopen_f(dset_id, aname, attr_id, error)
      ! Get dataspace for attribute
      call h5aget_space_f(attr_id, aspace_id, error)
      ! read attribute data
      call h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attribute, adims, error)
      ! close attribute
      call h5aclose_f(attr_id,error) ! Close the attribute.
      call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
    else
      attribute = 0.0_rk
    endif

    call h5dclose_f(dset_id,error)

  end subroutine read_attrib_dble2

  !-------------------------------------------------------------------------------
  ! Read an attribute
  ! INPUT:
  !   filename  what file to read from (e.g. hallo.h5)
  !   dsetname  what dataset to read from within the file (e.g. stuff to read hallo.h5:stuff)
  !   aname     the name of the attribute to read
  !   attribute the vector that will hold the attribute. note: this routine uses
  !             assumed shaped arrays: it will try to read as many values as the size of the vectors
  ! OUTPUT:
  !   attribute the values read from file
  !-------------------------------------------------------------------------------
  subroutine read_attrib_int(file_id, dsetname, aname, attribute, default_value)
    implicit none

    integer(hid_t), intent(in)                     :: file_id
    character(len=*), intent (in)                  :: dsetname, aname
    integer(kind=ik), DIMENSION(:), intent (inout) :: attribute
    integer(kind=ik), DIMENSION(1:size(attribute)), intent (in), optional :: default_value

    integer             :: dim
    integer             :: error  ! error flags
    integer(hid_t)      :: aspace_id ! Attribute Dataspace identifier
    integer(hid_t)      :: attr_id   ! Attribute identifier
    integer(hid_t)      :: dset_id  ! dataset identifier
    integer(hsize_t)    :: adims(1)  ! Attribute dimension
    logical             :: exists

    ! convert input data for the attribute to the precision required by the HDF library
    dim = size(attribute)
    adims = int(dim, kind=hsize_t)

    ! open the dataset
    call h5dopen_f(file_id, dsetname, dset_id, error)

    ! check if attribute exists
    call h5aexists_f(dset_id, aname, exists, error)

    if (exists) then
        ! open attribute
        call h5aopen_f(dset_id, aname, attr_id, error)
        ! Get dataspace for attribute
        call h5aget_space_f(attr_id, aspace_id, error)
        ! read attribute data
        call h5aread_f( attr_id, H5T_NATIVE_INTEGER, attribute, adims, error)
        ! close attribute
        call h5aclose_f(attr_id,error) ! Close the attribute.
        call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
    else
        if (present(default_value)) then
            attribute = default_value
        else
            attribute = 0
        endif
    endif

    call h5dclose_f(dset_id,error)

  end subroutine read_attrib_int

  !-------------------------------------------------------------------------------
  ! Read an attribute
  ! INPUT:
  !   filename  what file to read from (e.g. hallo.h5)
  !   dsetname  what dataset to read from within the file (e.g. stuff to read hallo.h5:stuff)
  !   aname     the name of the attribute to read
  !   attribute the vector that will hold the attribute. note: this routine uses
  !             assumed shaped arrays: it will try to read as many values as the size of the vectors
  ! OUTPUT:
  !   attribute the values read from file
  !-------------------------------------------------------------------------------
  subroutine read_attrib_int2(file_id, dsetname, aname, attribute, default_value)
    implicit none

    integer(hid_t), intent(in)                     :: file_id
    character(len=*), intent (in)                  :: dsetname, aname
    integer(kind=ik), intent (inout)               :: attribute
    integer(kind=ik), intent (in), optional        :: default_value

    integer             :: dim
    integer             :: error  ! error flags
    integer(kind=ik)    :: attribute_vct(1)
    integer(hid_t)      :: aspace_id ! Attribute Dataspace identifier
    integer(hid_t)      :: attr_id   ! Attribute identifier
    integer(hid_t)      :: dset_id  ! dataset identifier
    integer(hsize_t)    :: adims(1)  ! Attribute dimension
    logical             :: exists

    ! convert input data for the attribute to the precision required by the HDF library
    dim = 1
    adims = int(dim, kind=hsize_t)

    ! open the dataset
    call h5dopen_f(file_id, dsetname, dset_id, error)

    ! check if attribute exists
    call h5aexists_f(dset_id, aname, exists, error)

    if (exists) then
        ! open attribute
        call h5aopen_f(dset_id, aname, attr_id, error)
        ! Get dataspace for attribute
        call h5aget_space_f(attr_id, aspace_id, error)
        ! read attribute data
        call h5aread_f( attr_id, H5T_NATIVE_INTEGER, attribute_vct, adims, error)
        ! close attribute
        call h5aclose_f(attr_id,error) ! Close the attribute.
        call h5sclose_f(aspace_id,error) ! Terminate access to the data space.

        attribute = attribute_vct(1)
    else
        if (present(default_value)) then
            attribute = default_value
        else
            attribute = 0
        endif
    endif

    call h5dclose_f(dset_id,error)

  end subroutine read_attrib_int2


  ! overwrite and initialize file
  subroutine init_empty_file( fname )

    implicit none
    character (len=*), intent(in) :: fname

    open (55, file=fname, status='replace')
    close(55)

  end subroutine init_empty_file

end module module_hdf5_wrapper
