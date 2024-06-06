!> \brief sets refinement status to -2 for all sister blocks, if coarsening is possible
!! input:    - light data array \n
!! output:   - light data array
! ********************************************************************************************
!> \image html completeness.svg "Ensure Completeness" width=400

subroutine ensure_completeness( params, lgt_id, sisters, mark_TMP_flag )

    implicit none

    type (type_params), intent(in)      :: params                         !< user defined parameter structure
    integer(kind=ik), intent(inout)     :: sisters(:)                     !< light data array
    integer(kind=ik), intent(in)        :: lgt_id                         !< Concerned Block
    logical, intent(in), optional       :: mark_TMP_flag                  !< Set refinement of completeness to 0 or temporary flag
    integer(kind=ik)                    :: Jmax, k, l                     ! max treelevel
    integer(kind=ik)                    :: N_sisters, status, markTMPflag ! loop variables

    Jmax = params%Jmax
    N_sisters = size(sisters)

    markTMPflag = 0
    if (present(mark_TMP_flag)) then
        if (mark_TMP_flag) markTMPflag = REF_TMP_TREATED_COARSEN
    endif

    ! if all sisters exists, then the array should not contain values smaller
    ! zero (-1 would mean not found)
    if ( minval(sisters(1:N_sisters)) > 0 ) then

        ! now loop over all sisters, check if they also want to coarsen and have status -1
        ! only if all sisters agree to coarsen, they can all be merged into their mother block.
        status = -1
        do l = 1, N_sisters
            status = max( status, lgt_block(sisters(l), IDX_REFINE_STS) )
        end do

        ! if all agree and share the status -1, then we can indeed coarsen, keep -1 status
        if ( status == -1 ) then
            do l = 1, N_sisters
                lgt_block( sisters(l), IDX_REFINE_STS )  = -1
            end do
        elseif ( status >= 0 ) then
            ! We found all sister blocks, but they do not all share the -1 status: none
            ! of them can be coarsened, remove the status.
            do l = 1, N_sisters
                lgt_block( sisters(l), IDX_REFINE_STS )  = 0
            end do
        else
            call abort(197005, "This is odd - all sisters have temporary flag?")
        end if
    else
        ! We did not even find all sisters, that means a part of the four blocks is already
        ! refined. Therefore, they cannot be coarsened in any case, and we remove the coarsen
        ! flag
        do l = 1, N_sisters
            ! change status only for the existing sisters, mark temporary for leaf-wise as maybe this sister will be created soon
            if (sisters(l)>0) then
                lgt_block( sisters(l), IDX_REFINE_STS )  = markTMPflag
            endif
        end do
    end if

end subroutine ensure_completeness
