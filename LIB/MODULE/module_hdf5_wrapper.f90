! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_hdf5_wrapper.f90
! version: 0.4
! author: engels, msr
!
! !-------------------------------------------------------------------------------
! ! Read a field from a file
! ! Note a single file can contain many arrays, if required
! ! INPUT
! !   filename        file to read from
! !   dsetname        datasetname, i.e. the name of the array in the file
! !   rared           lower bounds of memory portion hold by the CPU
! !   rbred           upper bounds of memory portion hold by the CPU
! !                   NOTE: rared and rbred are 1:3 arrays. if running on one proc
! !                   they are rared=(/0,0,0/) rbred=(/nx-1,ny-1,nz-1/). If the data
! !                   is distributed among procs, each proc has to indicate which
! !                   portion of the array it holds
! !   field           actual data (ALLOCATED!!!)
! ! OUTPUT:
! !   field           data read from file
! !-------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 04/11/16 - linked old module version to wabbit (use module_params)
! ********************************************************************************************

module module_hdf5_wrapper

!---------------------------------------------------------------------------------------------
! modules

    use hdf5
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! Precision of doubles
    character(len=80) :: field_precision = "double" !"single"

    ! interface for writing attributes. an attribute is an object which is attached
    ! to a dataset, in our case a array saved in the file. we put useful information
    ! in attributes, for example the time, resolution and domain size
    ! both routines take vector values input
    interface write_attribute
    module procedure write_attrib_dble, write_attrib_int
    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

!-------------------------------------------------------------------------------
! write array to HDF5 file
! Note a single file can contain many arrays, if required
! INPUT
!   filename        file to write to
!   dsetname        datasetname, i.e. the name of the array in the file
!   field           actual data
!   overwrite       (optional) if .false., an existing file will not be erased and
!                   instead we just add the array to this file. Default is .true.
!                   NOTE: an error occurs if the dataset already exists, as hdf5
!                   does not provide a possibility to erase something already there
! OUTPUT:
!   none
!-------------------------------------------------------------------------------
subroutine write_field_hdf5( filename, dsetname, field, overwrite)
  implicit none

  character(len=*), intent (in) :: filename, dsetname
  real(kind=rk),intent(in) :: field(1:,1:)
  logical, intent(in), optional :: overwrite

  integer, parameter :: rank = 2 ! data dimensionality (2D or 3D)
  integer(hid_t) :: file_id   ! file identifier
  integer(hid_t) :: dset_id   ! dataset identifier
  integer(hid_t) :: filespace ! dataspace identifier in file
  integer(hid_t) :: memspace  ! dataspace identifier in memory
  integer(hid_t) :: plist_id  ! property list identifier
  integer(hid_t) :: file_precision

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dims_global
  ! hyperslab dimensions
  integer(hsize_t), dimension(rank) :: dims_local
  ! chunk dimensions
  integer(hsize_t), dimension(rank) :: chunk_dims
  ! how many blocks to select from dataspace
  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset = 0
  ! stride is spacing between elements, this is one here. striding is done in the
  ! caller; here, we just write the entire (possibly downsampled) field to disk.
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags
  integer :: nx, ny

  ! HDF attribute variables
  integer, parameter :: arank = 1
  !integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension

  !integer :: i, mindim, maxdim, mpicode
  logical :: exist1, ovrwrte

  if (present(overwrite)) then
    ovrwrte = overwrite
  else
    ! default is erase file and re-create it
    ovrwrte = .true.
  endif

  nx = size(field,1)
  ny = size(field,2)

  dims_global = (/nx, ny/)
  dims_local  = (/nx, ny/)
  chunk_dims  = (/nx, ny/)

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)
  ! Setup file access property list with parallel I/O access.
  ! this sets up a property list (plist_id) with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! setup serial driver (no MPI)
  call H5Pset_fapl_stdio_f(plist_id, error)

  !-----------------------------------------------------------------------------
  ! open the file
  !-----------------------------------------------------------------------------
  ! check if the file already exists
  inquire ( file=filename, exist=exist1 )
  if ((exist1 .eqv. .false. ) .or. (ovrwrte .eqv. .true.) ) then
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

  !-----------------------------------------------------------------------------
  ! create dataspace "filespace" to write to
  !-----------------------------------------------------------------------------
  ! Create the data space for the  dataset.
  ! Dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dims_global, filespace, error)

  ! Create chunked dataset.
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunk_dims, error)

  ! determine what precision to use when writing to disk
  if (field_precision=="double") then
    ! Output files in double precision
    file_precision = H5T_NATIVE_DOUBLE
  else
    ! Output files in single precision
    file_precision = H5T_NATIVE_REAL
  endif

  ! check if the dataset already exists
  call h5lexists_f(file_id, dsetname, exist1, error)
  if (exist1) then
    write(*,'(80("_"))')
    write(*,*) "ERROR: You are trying to write to an existing dataset...this is not supported."
    stop
  endif

  ! create the dataset
  call h5dcreate_f(file_id, dsetname, file_precision, filespace,dset_id, error, plist_id)
  call h5sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
  error, stride, dims_local)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)

  !-----------------------------------------------------------------------------
  ! create dataspace "memspace" to be written
  !-----------------------------------------------------------------------------
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dims_local, memspace, error)

  !-----------------------------------------------------------------------------
  ! actual writing of heavy data
  !-----------------------------------------------------------------------------
  ! Write the dataset collectively, double precision in memory
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dims_global, &
  error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

  !!! Close dataspaces:
  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error) ! Close the dataset.
  call h5pclose_f(plist_id, error) ! Close the property list.
  call h5fclose_f(file_id, error) ! Close the file.
  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.

end subroutine write_field_hdf5

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
subroutine write_attrib_dble(filename,dsetname,aname,attribute)

  implicit none

  character(len=*), intent (in) :: filename, dsetname, aname
  real(kind=rk), DIMENSION(:), intent (in) :: attribute

  integer, parameter :: arank = 1
  integer :: dim
  integer :: error  ! error flags
  integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id   ! Attribute identifier
  integer(hid_t) :: file_id
  integer(hid_t) :: dset_id  ! dataset identifier
  integer(hsize_t) :: adims(1)  ! Attribute dimension
  logical :: exists

  ! convert input data for the attribute to the precision required by the HDF library
  dim = size(attribute)
  adims = int(dim, kind=hsize_t)

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! open the file (existing file)
  call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error )

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
  call h5fclose_f(file_id,error)
  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.

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
subroutine write_attrib_int(filename,dsetname,aname,attribute)

  implicit none

  character(len=*), intent (in) :: filename, dsetname, aname
  integer, DIMENSION(:), intent (in) :: attribute

  integer(hsize_t) :: adims(1)  ! Attribute dimension
  integer, parameter :: arank = 1
  integer :: dim
  integer :: error  ! error flags
  integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id   ! Attribute identifier
  integer(hid_t) :: file_id
  integer(hid_t) :: dset_id  ! dataset identifier
  logical :: exists


  ! convert input data for the attribute to the precision required by the HDF library
  dim = size(attribute)
  adims = int(dim, kind=hsize_t)

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! open the file (existing file)
  call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error )

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
  call h5fclose_f(file_id,error)
  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.

end subroutine write_attrib_int

! overwrite and initialize file
subroutine init_empty_file( fname )

  implicit none
  character (len=*), intent(in) :: fname

  open (15, file=fname, status='replace')
  close(15)

end subroutine init_empty_file

end module module_hdf5_wrapper
