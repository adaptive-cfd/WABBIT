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
    integer(kind=ik)                    :: N_sisters, markTMPflag, lgt_sisters(1:2**params%dim)  ! loop variables

    Jmax = params%Jmax
    N_sisters = size(sisters)

    markTMPflag = 0
    if (present(mark_TMP_flag)) then
        if (mark_TMP_flag) markTMPflag = REF_TMP_TREATED_COARSEN
    endif

    ! if all sisters exists, then the array should not contain -1
    if ( all(sisters(1:N_sisters) /= -1) ) then

        lgt_sisters = lgt_block(sisters(1:N_sisters), IDX_REFINE_STS )

        ! first check : any block wants to stay so all blocks have to stay
        if ( any(lgt_sisters(:) == 0) ) then
            lgt_block( sisters(1:N_sisters), IDX_REFINE_STS )  = 0

        ! second check: any block has temp_gradedness flag > 0 - all blocks that want to stay have to wait with temp flag
        elseif ( any(lgt_sisters(:) > 0) ) then
            do l = 1, N_sisters
                if (lgt_block( sisters(l), IDX_REFINE_STS )  == -1) lgt_block( sisters(l), IDX_REFINE_STS )  = markTMPflag
            end do
        ! third check: all blocks have either temp flag or want to coarsen and can do it so lets change all flags to coarsen (-1)
        elseif ( all(lgt_sisters(:) < 0)) then
            lgt_block( sisters(1:N_sisters), IDX_REFINE_STS )  = -1

        ! other cases should not occur currently
        else
            call abort(197, "How did you end up here?")
        endif
    else
        ! We did not even find all sisters, that means a part of the four blocks is already
        ! refined. Therefore, they cannot be coarsened in any case, and we remove the coarsen
        ! flag
        do l = 1, N_sisters
            ! change status only for the existing sisters, mark temporary for leaf-wise as maybe this sister will be created soon
            if (sisters(l) /= -1) then
                lgt_block( sisters(l), IDX_REFINE_STS )  = markTMPflag
            endif
        end do
    end if

end subroutine ensure_completeness
