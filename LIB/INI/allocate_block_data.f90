! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: allocate_block_data.f90
! version: 0.4
! author: msr, engels
!
! allocate and reset heavy block data
!
! heavy data array  -> dim 4: block id  ( 1:number_blocks )
!                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 3: data type ( field_1, 2:number_data_fields+1, data_old, k1, k2, k3,
!                                       k4 [for runge kutta] )
!           field_1 (to save mixed data):   line 1: x coordinates
!                                           line 2: y coordinates
!
! input:    - maximal number of blocks per process
!           - grid parameter
!           - number of data fields
! output:   - empty heavy data array
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! ********************************************************************************************

subroutine  allocate_block_data(block_data, number_blocks, Bs, g, dF)

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! heavy data array
    real(kind=rk), allocatable, intent(out)     :: block_data(:, :, :, :)
    ! number of heavy data
    integer(kind=ik), intent(in)                :: number_blocks
    ! grid parameter, blocksize (Bs), ghostnodes (g), number of fields (F)
    integer(kind=ik), intent(in)                :: Bs, g, dF

    ! allocation error variable
    integer(kind=ik)                            :: allocate_error

    ! MPI error variable
    integer(kind=ik)                            :: ierr
    ! process rank
    integer(kind=ik)                            :: rank

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! allocate memory
    allocate( block_data( Bs+2*g, Bs+2*g, dF, number_blocks ), stat=allocate_error )

    ! reset data
    !
    block_data(:, :, :, :) = 0.0_rk

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') number_blocks, dF
    end if

end subroutine allocate_block_data
