! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: allocate_work_data.f90
! version: 0.4
! author: msr
!
! allocate and reset heavy work data
!
! heavy data array  -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 3: data type ( data_old, k1, k2, k3, k4 )
!                   -> dim 4: block id  ( 1:number_blocks )
!
! input:    - maximal number of blocks per process
!           - grid parameter
! output:   - empty heavy data array
!
! = log ======================================================================================
!
! 07/12/16 - create
! ********************************************************************************************

subroutine  allocate_work_data(hvy_work, number_blocks, Bs, g, dF)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! heavy work data array
    real(kind=rk), allocatable, intent(out)     :: hvy_work(:, :, :, :)
    ! number of heavy data
    integer(kind=ik), intent(in)                :: number_blocks
    ! grid parameter, blocksize (Bs), ghostnodes (g), number of datafields (F)
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
    allocate( hvy_work( Bs+2*g, Bs+2*g, dF*5, number_blocks ), stat=allocate_error )

    ! reset data
    !
    hvy_work(:, :, :, :) = 0.0_rk

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') number_blocks
    end if

end subroutine allocate_work_data
