!> \file
! WABBIT
!> \name read_field_flusi.f90
!> \version 0.5
!> \author sm
!
!> \brief reads a field from a .h5 file saved in flusi format
!
! = log ======================================================================
!> \date  9/3/2018 - create hashcode: commit 
!-----------------------------------------------------------------------------
subroutine read_field_flusi ( fname, hvy_block, lgt_block, hvy_n ,hvy_active, params)


  implicit none
  !> file name
  character(len=*),intent(in)            :: fname
  !> heavy data array - block data
  real(kind=rk), intent(inout)           :: hvy_block(:, :, :, :, :)
  !> user defined parameter structure
  type (type_params), intent(in)         :: params
  integer(kind=ik), intent(in)           :: hvy_active(:)
  integer(kind=ik), intent(in)           :: lgt_block(:, :)
  integer(kind=ik), intent(in)           :: hvy_n

  ! grid parameter
  integer(kind=ik)                    :: Bs
  integer(kind=ik)                    :: k, lgt_id, start_x, start_y, start_z
  ! offset variables
  integer(kind=ik), dimension(3)      :: ubounds3D, lbounds3D
  integer(kind=ik), dimension(2)      :: ubounds2D, lbounds2D
  real(kind=rk), dimension(3)         :: x0, dx
  ! file id integer
  integer(hid_t)                      :: file_id
!----------------------------------------------------------------------------
  Bs = params%number_block_nodes

  call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)
    ! print a message
  if (params%rank==0) then
      write(*,'(80("_"))')
      write(*,'("READING: Reading Flusi datafield from file ",A)') &
          trim(adjustl(fname))
  end if

  do k=1, hvy_n
      call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
      !dx for flusi and for wabbit are the same!
      start_x = floor(x0(1)/dx(1))
      start_y = floor(x0(2)/dx(2))
      write(*,*) start_x, start_y
      if (params%threeD_case) then
          start_z = floor(x0(3)/dx(3))
          lbounds3D = (/start_x, start_y , start_z/)
          ubounds3D = (/start_x+Bs-1,start_y+Bs-1, start_z+Bs-1/) -1
          call read_dset_mpi_hdf5_3D(file_id, get_dsetname(fname), lbounds3D, ubounds3D, &
            hvy_block(1:Bs-1, 1:Bs-1, 1:Bs-1, 1, hvy_active(k)))
      else
          lbounds2D = (/start_x, start_y /)
          ubounds2D = (/start_x+Bs-1,start_y+Bs-1/) -1
          call read_dset_mpi_hdf5_2D(file_id, get_dsetname(fname), lbounds2D, ubounds2D, &
            hvy_block(1:Bs-1, 1:Bs-1, 1, 1, hvy_active(k)))
      end if
  end do

    ! close file and HDF5 library
    call close_file_hdf5(file_id)


end subroutine read_field_flusi

subroutine get_attributes_flusi(fname, nxyz, time, domain)

    implicit none
    !> file name
    character(len=*), intent(in)                  :: fname
    !> number of active blocks
    integer(kind=ik), dimension(3), intent(out)   :: nxyz
    !> time (to be read from file)
    real(kind=rk), intent(out)                    :: time
    !> domain size
    real(kind=rk), dimension(3), intent(out)      :: domain
    real(kind=rk), dimension(1)                   :: ttime
    integer(hid_t)                                :: file_id


    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)
    ! read attributes
    call read_attribute(file_id, trim(get_dsetname(fname)), "domain_size", domain)
    call read_attribute(file_id, trim(get_dsetname(fname)), "time", ttime)
    call read_attribute(file_id, trim(get_dsetname(fname)), "nxyz", nxyz)
    time = ttime(1)
    ! close file and HDF5 library
    call close_file_hdf5(file_id)
end subroutine get_attributes_flusi

character(len=80)  function get_dsetname(fname)
    implicit none
    character(len=*), intent(in) :: fname
    ! extract dsetname (from "/" until "_", excluding both)
    get_dsetname  = fname  ( index(fname,'/',.true.)+1:index( fname, '_',.true. )-1 )
    return
end function get_dsetname