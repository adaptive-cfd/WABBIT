!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_field.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief read data of a single datafield dF at iteration and time t
!
!>
!! input:
!!           - name of the file we want to read from
!!           - parameter array
!!           - heavy data array
!!           - number of active blocks (light and heavy)
!!           - 
!!
!! output:
!!           - heavy data array
!!
!!
!! = log ======================================================================================
!! \n
!! 22/09/17 - create
!
! ********************************************************************************************

subroutine read_field(fname, dF, params, hvy_block, lgt_n, hvy_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> file name
    character(len=*), intent(in)        :: fname
    !> datafield number
    integer(kind=ik), intent(in)        :: dF
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(in)        :: hvy_n, lgt_n

    ! block data buffer, need for compact data storage
    real(kind=rk), allocatable          :: myblockbuffer(:,:,:,:)

    ! file id integer
    integer(hid_t)                      :: file_id


    ! process rank
    integer(kind=ik)                    :: rank
    ! grid parameter
    integer(kind=ik)                    :: Bs, g

    ! offset variables
    integer(kind=ik), dimension(4)      :: ubounds3D, lbounds3D
    integer(kind=ik), dimension(3)      :: ubounds2D, lbounds2D

    real(kind=rk), dimension(3)         :: domain

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank = params%rank

    ! grid parameter
    Bs   = params%number_block_nodes
    g    = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

    hvy_block = 0.0_rk  ! oder 9.99e99_rk???
    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

    allocate(myblockbuffer( 1:Bs, 1:Bs, 1:Bs, 1:hvy_n ))

    if ( params%threeD_case ) then

        ! tell the hdf5 wrapper what part of the global [bs x bs x bs x n_active]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds3D = (/1,1,1,lgt_n+1/) - 1
        ubounds3D = (/Bs-1,Bs-1,Bs-1,lbounds3D(4)+hvy_n-1/)

    else

        ! tell the hdf5 wrapper what part of the global [bs x bs x bs x n_active]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds2D = (/1,1,lgt_n+1/) - 1
        ubounds2D = (/Bs-1,Bs-1,lbounds2D(3)+hvy_n-1/)

    endif

    ! DO WE NEED THIS?
    ! if (rank==0) then
    !     write(*,'(40("~"))')
    !     write(*,'("Reading from file ",A)') trim(adjustl(fname))
    !     write(*,'(" time=",g12.4') ,time(1)
    !     write(*,'("Lx=",g12.4," Ly=",g12.4," Lz=",g12.4)') domain

    !     ! if the domain size doesn't match, proceed, but yell.
    !     if ((params%Lx.ne.domain(1)).or.(params%Ly.ne.domain(2)).or.(params%Lz.ne.domain(3))) then
    !         write (*,'(A)') " WARNING! Domain size mismatch."
    !         write (*,'("in memory:   Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') params%Lx,params%Ly,params%Lz
    !         write (*,'("but in file: Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') domain
    !         write (*,'(A)') "proceed, with fingers crossed."
    !     end if
    ! end if


    ! actual reading of file
    if ( params%threeD_case ) then
        ! 3D data case
        call read_dset_mpi_hdf5_4D(file_id, "blocks", lbounds3D,ubounds3D, myblockbuffer)

    else
        ! 2D data case
        call read_dset_mpi_hdf5_3D(file_id, "blocks", lbounds2D, ubounds2D, myblockbuffer)
    end if
    
    ! close file and HDF5 library
    call close_file_hdf5(file_id)
    
    hvy_block(:,:,:,dF,1:hvy_n) = myblockbuffer(:,:,:,1:hvy_n)
    deallocate(myblockbuffer)
end subroutine read_field
