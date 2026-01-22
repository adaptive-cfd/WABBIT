!> \brief reads a field from a .h5 file saved in flusi format
!-----------------------------------------------------------------------------
subroutine read_field_flusi(fname, hvy_block, params, Bs_f, tree_ID)

  implicit none
  character(len=*),intent(in)            :: fname                       !> file name
  real(kind=rk), intent(inout)           :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
  type (type_params), intent(in)         :: params                      !> user defined parameter structure
  integer(kind=ik), intent(in)           :: tree_ID
  integer(kind=ik), dimension(3), intent(in) :: Bs_f

  ! grid parameter
  integer(kind=ik), dimension(3)      :: Bs
  integer(kind=ik)                    :: k, lgt_id, start_x, start_y, start_z
  ! offset variables
  integer(kind=ik), dimension(3)      :: ubounds3D, lbounds3D
  integer(kind=ik), dimension(3)      :: ubounds2D, lbounds2D
  real(kind=rk), dimension(3)         :: x0, dx
  ! file id integer
  integer(hid_t)                      :: file_id
  real(kind=rk), dimension(:,:,:), allocatable   :: blockbuffer

!----------------------------------------------------------------------------
  Bs = params%Bs
  call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)
    ! print a message
  if (params%rank==0) then
      write(*,'(80("─"))')
      write(*,'("READING: Reading Flusi datafield from file ",A)') &
          trim(adjustl(fname))
  end if

  if (params%dim == 3) then
      allocate( blockbuffer(Bs_f(1)+1,Bs_f(2)+1,Bs_f(3)+1))
      lbounds3D = (/0, 0, 0/)
      ubounds3D = (/Bs_f(1), Bs_f(2), Bs_f(3)/)-1
      call read_dset_mpi_hdf5(file_id, get_dsetname(fname), lbounds3D, ubounds3D, &
          blockbuffer(1:Bs_f(1),1:Bs_f(2), 1:Bs_f(3)))
  else
      allocate( blockbuffer(1,Bs_f(1)+1,Bs_f(2)+1))
      lbounds2D = (/0, 0, 0/)
      ubounds2D = (/1, Bs_f(1), Bs_f(2)/)-1
      call read_dset_mpi_hdf5(file_id, get_dsetname(fname), lbounds2D, ubounds2D, &
          blockbuffer(1,1:Bs_f(1),1:Bs_f(2)))
  end if

  blockbuffer(:,Bs_f(2)+1,:) = blockbuffer(:,1,:)
  blockbuffer(:,:,Bs_f(3)+1) = blockbuffer(:,:,1)
  if (params%dim == 3) blockbuffer(Bs_f(1)+1,:,:) = blockbuffer(1,:,:)

  do k=1, hvy_n(tree_ID)

      call hvy2lgt(lgt_id, hvy_active(k,tree_ID), params%rank, params%number_blocks)
      call get_block_spacing_origin( params, lgt_id, x0, dx )

      start_x = nint(x0(1)/dx(1)) + 1
      start_y = nint(x0(2)/dx(2)) + 1
      if (params%dim == 3) then
          start_z = nint(x0(3)/dx(3)) + 1
          hvy_block(1:Bs(1), 1:Bs(2), 1:Bs(3), 1, hvy_active(k,tree_ID)) = blockbuffer(start_x:start_x+Bs(1)-1,&
              start_y:start_y+Bs(2)-1,start_z:start_z+Bs(3)-1)
      else
          hvy_block(1:Bs(1), 1:Bs(2), 1, 1, hvy_active(k,tree_ID)) = blockbuffer(1,&
              start_x:start_x+Bs(1)-1,start_y:start_y+Bs(2)-1)
      end if
  end do

  ! close file and HDF5 library
  call close_file_hdf5(file_id)

end subroutine read_field_flusi

subroutine read_field_flusi_MPI(fname, hvy_block, params, tree_ID)

  implicit none
  !> file name
  character(len=*),intent(in)         :: fname
  !> heavy data array - block data
  real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
  !> user defined parameter structure
  type (type_params), intent(in)      :: params
  integer(kind=ik), intent(in)        :: tree_ID

  integer(kind=ik), dimension(3)      :: nxyz
  integer(kind=ik)                    :: g
  integer(kind=ik), dimension(3)      :: Bs
  integer(kind=ik)                    :: k, lgt_id, start_x, start_y, start_z
  ! offset variables
  integer(kind=ik), dimension(3)      :: ubounds, lbounds, num_Bs
  real(kind=rk), dimension(3)         :: x0, dx, domain
  real(kind=rk) :: time
  ! file id integer
  integer(hid_t)                      :: file_id
  real(kind=rk), dimension(:,:,:), allocatable   :: blockbuffer, data_flusi

  ! read attributes such as number of discretisation points, time, domain size
  call get_attributes_flusi(fname, nxyz, time, domain)

  Bs = params%Bs
  g  = params%g
  call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

  ! print a message
  if (params%rank==0) then
      write(*,'(80("─"))')
      write(*,'("READING: Reading Flusi datafield from file ",A)') trim(adjustl(fname))
  end if

  if (params%number_procs /= 1) then
      ! we read the entire flusi datafield to memory. This means the routine is
      ! monoprocessing only. It can be fixed, but I encountered a HDF5 problem and
      ! am too busy to fight it. -Thomas
      call abort(2023232999, "This routine works only with one CPU")
  endif

  allocate( data_flusi(0:nxyz(1)-1,0:nxyz(2)-1,0:nxyz(3)-1) )
  call read_dset_mpi_hdf5(file_id, get_dsetname(fname), (/0,0,0/), (/nxyz(1)-1,nxyz(2)-1,nxyz(3)-1/), data_flusi)

  ! would be fairly easy but I have no time.
  if (nxyz(1) /= 0) call abort(13738213, "Only 2D right now")

  !> \todo test for 3D
  do k = 1, hvy_n(tree_ID)
      call hvy2lgt(lgt_id, hvy_active(k,tree_ID), params%rank, params%number_blocks)
      call get_block_spacing_origin( params, lgt_id, x0, dx )

      ! from spacing and origin of the block, get position in flusi matrix
      start_x = nint( x0(1) / dx(1) )
      start_y = nint( x0(2) / dx(2) )

      ! ! if (params%dim == 3) then
      ! !     start_z = nint(x0(3)/dx(3))
      ! !
      ! !     lbounds = (/start_x, start_y, start_z/)
      ! !     ubounds = (/end_bound(start_x,Bs(1),nxyz(1)), end_bound(start_y,Bs(2),nxyz(2)),&
      ! !         end_bound(start_z,Bs(3),nxyz(3))/)
      ! !
      ! !     num_Bs = ubounds-lbounds+1
      ! !
      ! !     call read_dset_mpi_hdf5(file_id, get_dsetname(fname), lbounds, ubounds, &
      ! !     hvy_block(g+1:g+num_Bs(1), g+1:g+num_Bs(2), g+1:g+num_Bs(3), 1, hvy_active(k,tree_ID)))
      ! ! else
      !     lbounds = (/ 0, start_x        , start_y         /)
      !     ubounds = (/ 0, start_x+Bs(1)-1, start_y+Bs(2)-1 /)
      !     ! lbounds = (/ start_x        , start_y        , 0 /)
      !     ! ubounds = (/ start_x+Bs(1)-1, start_y+Bs(2)-1, 0 /)
      !     num_Bs = ubounds - lbounds + 1
      !
      !     write(*,*) "lbound", lbounds
      !     write(*,*) "ubound", ubounds
      !     write(*,*) "numbs", num_bs
      !
      !     call read_dset_mpi_hdf5(file_id, get_dsetname(fname), lbounds, ubounds, &
      !     blockbuffer(0:num_Bs(1)-1,0:num_Bs(2)-1,1))

          ! hvy_block(g+1:g+num_Bs(1),g+1:g+num_Bs(2), 1, 1, hvy_active(k,tree_ID)) = blockbuffer(0:num_Bs(1)-1,0:num_Bs(2)-1,1)
      hvy_block(g+1:g+Bs(1), g+1:g+Bs(2), 1, 1, hvy_active(k,tree_ID)) = &
      data_flusi(0, start_x:start_x+Bs(1)-1, start_y:start_y+Bs(2)-1 )
      ! end if
  end do

  ! close file and HDF5 library
  call close_file_hdf5(file_id)
  if ( allocated(blockbuffer) ) deallocate(blockbuffer)

end subroutine read_field_flusi_MPI


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
