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
!!           - 
!!           - 
!!           - 
!!           - 
!!           - 
!!
!! output:
!!           -
!!
!!
!! = log ======================================================================================
!! \n
!! 22/09/17 - create
!
! ********************************************************************************************

subroutine read_field(fname, dF, params, )

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

    t1 = MPI_wtime()

    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

    ! call read_attribute(fname, "blocks", "domain_size", domain)
    ! call read_attribute(fname, "blocks", "time", ttime)
    ! call read_attribute(fname, "blocks", "iteration", iiteration)
    ! call read_attribute(fname, "blocks", "total_number_blocks", lgt_n)
    
    ! if (rank==0) then
    !     hvy_n = lgt_n/params%number_procs + modulo(lgt_n,params%number_procs)
    ! else
    !     hvy_n = lgt_n/params%number_procs
    ! end if

    allocate(myblockbuffer( 1:Bs, 1:Bs, 1:Bs, 1:hvy_n ))
    allocate(block_treecode(1:params%max_treelevel, 1:hvy_n))

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

    if (rank==0) then
        write(*,'(40("~"))')
        write(*,'("Reading from file ",A)') trim(adjustl(fname))
        write(*,'("nx=",i4," ny=",i4," nz=",i4," time=",g12.4') ,ttime(1)
        write(*,'("Lx=",g12.4," Ly=",g12.4," Lz=",g12.4)') domain

        ! if the domain size doesn't match, proceed, but yell.
        if ((params%Lx.ne.domain(1)).or.(params%Ly.ne.domain(2)).or.(params%Lz.ne.domain(3))) then
            write (*,'(A)') " WARNING! Domain size mismatch."
            write (*,'("in memory:   Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') params%Lx,params%Ly,params%Lz
            write (*,'("but in file: Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') domain
            write (*,'(A)') "proceed, with fingers crossed."
        end if

        ! ! if the resolutions do not match, yell and hang yourself
        ! if ((nx/=nxyz(1)).or.(ny/=nxyz(2)).or.(nz/=nxyz(3))) then
        !   write (*,'(A)') "ERROR! Resolution mismatch"
        !   write (*,'(A)') "This happens if ra(:) and rb(:) are not properly initialized."
        !   write (*,'("in memory:   nx=",i4," ny=",i4," nz=",i4)') nx,ny,nz
        !   write (*,'("but in file: nx=",i4," ny=",i4," nz=",i4)') nxyz
        !   stop
        ! endif
    end if


    ! actual reading of file
    if ( params%threeD_case ) then
        ! 3D data case
        call read_dset_mpi_hdf5_4D(file_id, "blocks", lbounds3D,ubounds3D, myblockbuffer)
        call read_dset_mpi_hdf5_2D(file_id, "block_treecode", (/0,lbounds3D(4)/), (/2,ubounds3D(4)/), block_treecode)

    else
        ! 2D data case
        call read_dset_mpi_hdf5_3D(file_id, "blocks", lbounds2D, ubounds2D, myblockbuffer)
        call read_dset_mpi_hdf5_2D(file_id, "block_treecode", (/0,lbounds2D(3)/), (/2,ubounds2D(3)/), block_treecode)
    end if
    
    ! close file and HDF5 library
    call close_file_hdf5(file_id)
    
    do k=1, hvy_n
        hvy_block(:,:,:,dF,k) = myblockbuffer(:,:,:,k)
        call hvy_id_to_lgt_id( lgt_id, k, rank, hvy_n )
        lgt_block(lgt_id,1:params%max_treelevel) = block_treecode(k,:)
        !lgt_block(rank*hvy_n+1:hvy_n+rank*hvy_n,1:params%max_treelevel) = block_treecode(:,:)
    end do

    !create hvy_active and lgt_active list
    call create_active_and_sorted_lists(params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, create_sorted_list)

    deallocate(myblockbuffer)
    deallocate(block_treecode)
end subroutine read_field
