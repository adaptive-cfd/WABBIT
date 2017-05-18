!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name update_neighbors.f90
!> \version 0.4
!> \author engels
!
!> \brief functional wrapper for the 2D and 3D version of update neighbors.
!
!>
!! input:
!!            - light data array
!!            - params struct
!!
!! output:
!!            - neighbor list array
!!
!! = log ======================================================================================
!! \n
!! 31/03/17 - create (tired of if-clause everywhere...)
!
! ********************************************************************************************

subroutine update_neighbors(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n)

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    if ( params%threeD_case ) then
        ! 3D:
        call update_neighbors_3D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n)
    else
        ! 2D:
        call update_neighbors_2D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n)
    end if

end subroutine update_neighbors
