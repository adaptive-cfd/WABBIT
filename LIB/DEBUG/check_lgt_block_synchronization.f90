!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name check_lgt_block_synchronization.f90
!> \version 0.4
!> \author msr
!
!> \brief debug lgt_block data
!
!>
!!  proc 0 gather all data and compare the data to his own light data \n
!!
!!
!! input:    - params, light data \n
!! output:   - status of lgt_block synchronzation \n
!!
!!
!! = log ======================================================================================
!! \n
!! 29/11/16 - create
! ********************************************************************************************

subroutine check_lgt_block_synchronization( params, lgt_block)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)

    ! local light data array
    integer(kind=ik)                    :: my_lgt_block( size(lgt_block,1) , size(lgt_block,2)), lgt_block_0( size(lgt_block,1) , size(lgt_block,2))

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! loop variables
    integer(kind=ik)                    :: k, l, lgt_start
    !> MPI communicator
    integer(kind=ik)                    :: WABBIT_COMM

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determine process rank
    rank        = params%rank
    WABBIT_COMM = params%WABBIT_COMM

    ! light data start id
    call proc_to_lgt_data_start_id( lgt_start, rank, params%number_blocks )

    ! reset local light data
    lgt_block_0  = 0
    my_lgt_block = 0

    my_lgt_block( lgt_start : lgt_start + params%number_blocks - 1 , :) = lgt_block( lgt_start : lgt_start + params%number_blocks - 1 , :)

!---------------------------------------------------------------------------------------------
! main body

    ! gather light data
    call MPI_Allreduce(my_lgt_block, lgt_block_0, size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_SUM, WABBIT_COMM, ierr)

    ! loop over all blocks
    do k = 1, size(lgt_block, 1)
        do l = 1, size(lgt_block, 2)

            ! compare light data
            if ( lgt_block_0(k,l) /= lgt_block(k,l) ) then
                ! error case
                write(*,'(80("!"))')
                write(*, '("ERROR: lgt_block is not synchron, lgt_block id: ",i5," position: ", i5)') k, l
                stop
            end if

        end do
    end do

end subroutine check_lgt_block_synchronization
