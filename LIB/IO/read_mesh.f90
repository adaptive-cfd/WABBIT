!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_mesh.f90
!> \version 0.5
!> \author sm
!
!> \brief read mesh properties of saved mesh in hdf5-file at iteration and time t
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
!! 29/09/17 - create
!
! ********************************************************************************************

subroutine read_mesh(fname, dF, params, lgt_n, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> file name
    character(len=*), intent(in)        :: fname
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> number of active blocks (heavy and light data)
    integer(kind=ik), intent(in)        :: hvy_n, light_n

    ! file id integer
    integer(hid_t)                      :: file_id
    ! process rank
    integer(kind=ik)                    :: rank
    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! offset variables
    integer(kind=ik), dimension(2)      :: ubounds, lbounds

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

    allocate(block_treecode(1:params%max_treelevel, 1:hvy_n))

    ! tell the hdf5 wrapper what part of the global [ n_active x max_treelevel + 2]
    ! array we want to hold, so that all CPU can read from the same file simultaneously
    ! (note zero-based offset):
    lbounds = (/1, lgt_n + 1/) - 1
    ubounds = (/2, lgt_n + hvy_n - 1/)

    call read_dset_mpi_hdf5_2D(file_id, "block_treecode", lbounds, ubounds, block_treecode)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)
   
    do k=1, hvy_n
        call hvy_id_to_lgt_id( lgt_id, k, rank, hvy_n )
        lgt_block(lgt_id,1:params%max_treelevel) = block_treecode(k,:)
        !lgt_block(rank*hvy_n+1:hvy_n+rank*hvy_n,1:params%max_treelevel) = block_treecode(:,:)
    end do

    deallocate(block_treecode)
end subroutine read_mesh
