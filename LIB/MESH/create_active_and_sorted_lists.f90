!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_active_and_sorted_lists.f90
!> \version 0.5
!> \author msr
!
!> \brief create all active (lgt/hvy) lists, create also sorted list of active
!! light blocks with numerical treecodes
!! \n
!! input:    
!!           - light data
!!
!! output:
!!           - list of active blocks
!!           - number of active blocks
!!           - sorted light data with numerical treecodes
!!
!> = log ======================================================================================
!! \n
!! 14/06/17 - create subroutine
!!
! ********************************************************************************************

subroutine create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, create_sorted_list )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)

    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)

    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n

    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)

    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_n

    !> sorted light data with numerical treecodes
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)

    !> switch for sorted list creation
    logical, intent(in)                 :: create_sorted_list

    ! loop variables
    integer(kind=ik)                    :: k, N, heavy_id, block_rank

    ! process rank
    integer(kind=ik)                    :: rank

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset old lists, use old numbers of active blocks
    lgt_active(1:lgt_n)          = -1
    hvy_active(1:hvy_n)          = -1
    lgt_sortednumlist(1:lgt_n,:) = -1

    ! reset active block numbers
    lgt_n = 0
    hvy_n = 0

    rank = params%rank
    N    = params%number_blocks

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all light data
    do k = 1, size(lgt_block, 1)
        ! block is active
        if ( lgt_block(k, 1) /= -1 ) then

            ! save lgt id as active block
            lgt_active( lgt_n + 1 ) = k
            lgt_n                   = lgt_n + 1

            ! save heavy id, only if proc responsable for block
            call lgt_id_to_proc_rank( block_rank, k, N )
            if ( rank == block_rank ) then
                ! convert light data id into heavy data id
                call lgt_id_to_hvy_id( heavy_id, k, rank, N)
                hvy_active( hvy_n + 1 ) = heavy_id
                hvy_n                   = hvy_n + 1
            end if

            if (create_sorted_list) then
                ! sorted list
                ! first index stores the light id of the block
                lgt_sortednumlist(lgt_n, 1) = k
                ! second index stores the numerical treecode
                lgt_sortednumlist(lgt_n, 2) = treecode2int( lgt_block(k, 1:params%max_treelevel) )
            end if

        end if
    end do

    if (create_sorted_list) then
        ! sort list
        if (lgt_n > 1) then
            call quicksort(lgt_sortednumlist, 1, lgt_n, 2)
        end if
    end if

end subroutine create_active_and_sorted_lists
