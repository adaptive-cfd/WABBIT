!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name ensure_completeness_2D.f90
!> \version 0.4
!> \author msr
!
!> \brief sets refinement status to -2 for all sister blocks, if coarsening is possible
!
!> \details
!! input:    - light data array \n
!! output:   - light data array
!!
!!
!! = log ======================================================================================
!! \n
!! 10/11/16 - switch to v0.4 \n
!! 05/04/17 - works for 2D and 3D data and uses readable find_sisters routine.
! ********************************************************************************************
!> \image html completeness.svg "Ensure Completeness" width=400

subroutine ensure_completeness( params, lgt_block, lgt_id, sisters )

!---------------------------------------------------------------------------------------------
! modules


!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :), sisters(:)
    integer(kind=ik), intent(in)        :: lgt_id

    ! max treelevel
    integer(kind=ik)                    :: Jmax
    ! loop variables
    integer(kind=ik)                    :: k, l, N_sisters, status

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    Jmax = params%max_treelevel
    N_sisters = size(sisters)

!---------------------------------------------------------------------------------------------
! main body

    ! if all sisters exists, then the array should not contain values smaller
    ! zero (-1 would mean not found)
    if ( minval(sisters(1:N_sisters)) > 0 ) then

        ! now loop over all sisters, check if they also want to coarsen and have status -1
        ! only if all sisters agree to coarsen, they can all be merged into their mother block.
        status = -1
        do l = 1, N_sisters
            status = max( status, lgt_block(sisters(l), Jmax+idx_refine_sts) )
        end do

        ! if all agree and share the status -1, then we can indeed coarsen, keep -1 status
        if ( status == -1 ) then
            do l = 1, N_sisters
                lgt_block( sisters(l), Jmax+idx_refine_sts )  = -1
            end do
        else
            ! We found all sister blocks, but they do not all share the -1 status: none
            ! of them can be coarsened, remove the status.
            do l = 1, N_sisters
                lgt_block( sisters(l), Jmax+idx_refine_sts )  = 0
            end do
        end if
    else
        ! We did not even find all sisters, that means a part of the four blocks is already
        ! refined. Therefore, they cannot be coarsened in any case, and we remove the coarsen
        ! flag
        do l = 1, N_sisters
            ! change status only for the existing sisters
            if (sisters(l)>0) then
                lgt_block( sisters(l), Jmax+idx_refine_sts )  = 0
            endif
        end do
    end if

end subroutine ensure_completeness
