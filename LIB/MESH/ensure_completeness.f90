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

subroutine ensure_completeness( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )

!---------------------------------------------------------------------------------------------
! modules


!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)

    ! max treelevel
    integer(kind=ik)                    :: max_treelevel
    ! loop variables
    integer(kind=ik)                    :: k, l, N_sisters, status
    ! sister ids
    integer(kind=ik)                    :: id(1:8)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    max_treelevel = params%max_treelevel
    if (params%threeD_case) then
      N_sisters = 8
    else
      N_sisters = 4
    end if
    id = -1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active blocks
    do k = 1, lgt_n

        ! if the block wants to coarsen, check refinement status of sister blocks
        if ( lgt_block( lgt_active(k), max_treelevel+2 ) == -1) then
            ! find sister IDs of the block we're looking at. If a sister is not found, -1
            ! is returned in the array id
            ! NOTE: the array "id" also contains the block itself
            call find_sisters( params, lgt_active(k), id(1:N_sisters), lgt_block, lgt_n, lgt_sortednumlist )

            ! if all sisters exists, then the array should not contain values smaller
            ! zero (-1 would mean not found)
            if ( minval(id(1:N_sisters)) > 0 ) then

                ! now loop over all sisters, check if they also want to coarsen and have status -1
                ! only if all sisters agree to coarsen, they can all be merged into their mother block.
                status = -1
                do l = 1, N_sisters
                  status = max( status, lgt_block(id(l), max_treelevel+2) )
                end do

                ! if all agree and share the status -1, then we can indeed coarsen, keep -1 status
                if ( status == -1 ) then
                    do l = 1, N_sisters
                        lgt_block( id(l), max_treelevel+2 )  = -1
                    end do
                else
                    ! We found all sister blocks, but they do not all share the -1 status: none
                    ! of them can be coarsened, remove the status.
                    do l = 1, N_sisters
                        lgt_block( id(l), max_treelevel+2 )  = 0
                    end do
                end if
            else
                ! We did not even find all sisters, that means a part of the four blocks is already
                ! refined. Therefore, they cannot be coarsened in any case, and we remove the coarsen
                ! flag
                do l = 1, N_sisters
                    ! change status only for the existing sisters
                    if (id(l)>0) then
                        lgt_block( id(l), max_treelevel+2 )  = 0
                    endif
                end do
            end if
        end if
    end do


    ! NOTE: this routine runs redundantly on all procs - they all do the same. That means
    ! by consequence that the light data here does not have to be synchronized.

end subroutine ensure_completeness
