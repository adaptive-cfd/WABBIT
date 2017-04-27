!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_hvy_active_list.f90
!> \version 0.4 
!> \author msr
!
!> \brief  create a list of all active blocks (heavy data)
!> \details dim 1: heavy data id
!! \n
!! input:   
!!             - light data
!!
!! output:  
!!             - list of active blocks, heavy data id
!!             - number of active blocks
!!
!! = log ======================================================================================
!! \n
!! 24/11/16 - create subroutine
! ********************************************************************************************

subroutine create_hvy_active_list( lgt_block, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)

    !> list of active blocks (light data)
    integer(kind=ik), intent(out)       :: hvy_active(:)

    !> number of active blocks (light data)
    integer(kind=ik), intent(out)       :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: k, lgt_start, N, heavy_id

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization
    ! reset active list
    hvy_active = -1

    ! list index
    hvy_n = 0

    ! number of blocks per proc
    N = size(hvy_active, 1)
!---------------------------------------------------------------------------------------------
! main body

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! start index on light data
    call proc_to_lgt_data_start_id( lgt_start, rank, N )

    ! loop over all light data
    ! every proc loops only over the data corresponding to his rank
    do k = lgt_start, lgt_start + N - 1

        ! block is active
        if ( lgt_block(k, 1) /= -1 ) then
            ! convert light data id into heavy data id
            call lgt_id_to_hvy_id( heavy_id, k, rank, N)
            hvy_active( hvy_n + 1 ) = heavy_id
            hvy_n = hvy_n + 1
        end if

    end do

end subroutine create_hvy_active_list
