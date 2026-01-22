!> \brief sets refinement status to -2 for all sister blocks, if coarsening is possible
!! input:    - light data array \n
!! output:   - light data array
! ********************************************************************************************
!> \image html completeness.svg "Ensure Completeness" width=400

subroutine ensure_completeness_block( params, lgt_id, sisters, stay_value, check_daughters )

    implicit none

    type (type_params), intent(in)  :: params                !< user defined parameter structure
    integer(kind=ik), intent(inout) :: sisters(:)            !< light data array
    integer(kind=ik), intent(in)    :: lgt_id                !< Concerned Block
    !> in case the newly set values for blocks that stay should be temporary flags, set_value is provided with this flag value
    integer(kind=ik), intent(in), optional  :: stay_value
    logical, intent(in), optional   :: check_daughters       !< For full tree CVS grids completeness has to check the mother as well

    integer(kind=ik)                :: Jmax, k, l, ld, hvy_id, lgt_daughter, tmp_stay, N_sisters, markTMPflag, lgt_sisters(1:2**params%dim)
    logical                         :: checkDaughters

    Jmax = params%Jmax
    N_sisters = size(sisters)
    call lgt2hvy(hvy_id, lgt_id, params%rank, params%number_blocks)  ! needed for mother check

    ! If we loop leaf-wise we do not want blocks to be set to 0 if not all sisters exist
    ! as we want to keep the wavelet-decomposed values and maybe just have to wait a step longer
    checkDaughters = .false.
    if (present(check_daughters)) checkDaughters = check_daughters

    ! Sometimes we want to do refinement flag magic - this is important to differentiate between significant blocks (that have a specific value until now)
    ! or blocks kept or refined due to completeness or gradedness. For this we give them a flag value which can be provided, that shows which block stay, that are not significant
    tmp_stay = 0
    if (present(stay_value)) tmp_stay = stay_value

    ! if all sisters exists, then the array should not contain -1
    if ( all(sisters(1:N_sisters) /= -1)  ) then

        lgt_sisters = lgt_block(sisters(1:N_sisters), IDX_REFINE_STS )
        ! zeroth check : any block exists but has flag REF_TMP_EMPTY, so it has no values and all have to wait
        if ( any(lgt_sisters(:) == REF_TMP_EMPTY) ) then
            do l = 1, N_sisters
                if (lgt_block( sisters(l), IDX_REFINE_STS )  < 0) lgt_block( sisters(l), IDX_REFINE_STS )  = tmp_stay
            enddo

        ! first check : any block wants to stay so all blocks have to stay
        elseif ( any(lgt_sisters(:) == 0) .or. any(lgt_sisters(:) == tmp_stay) ) then
            do l = 1, N_sisters
                if (lgt_block( sisters(l), IDX_REFINE_STS )  < 0) lgt_block( sisters(l), IDX_REFINE_STS )  = tmp_stay
            enddo

        ! second check: all blocks want to coarsen and can do it so lets change all flags to coarsen (-1)
        elseif ( all(lgt_sisters(:) < 0)) then
            lgt_block( sisters(1:N_sisters), IDX_REFINE_STS )  = -1

        ! other cases should not occur currently
        else
            call abort(197, "How did you end up here?")
        endif

        ! for CVS grids we need to check if we have a daughter and if so, if she'd like to coarsen as well
        ! all sisters share the same mother, so this hvy operation will be done correctly on all processors
        if (checkDaughters .and. .not. block_is_leaf(params, hvy_id)) then
            ! family takes care of one another - if any daughter wants to stay then me and all my sisters have to stay as well
            do ld = 2+2**params%dim, 1+2**(params%dim+1)
                lgt_daughter = hvy_family(hvy_ID, ld)
                if (lgt_daughter /= -1) then
                    if (lgt_block(lgt_daughter, IDX_REFINE_STS) == 0 .or. lgt_block(lgt_daughter, IDX_REFINE_STS) == tmp_stay) then
                        do l = 1, N_sisters
                            if (lgt_block( sisters(l), IDX_REFINE_STS )  < 0) lgt_block( sisters(l), IDX_REFINE_STS )  = tmp_stay
                        enddo
                    endif
                endif
            enddo
        endif
    else
        ! We did not even find all sisters, that means a part of the four blocks is already
        ! refined. Therefore, they cannot be coarsened in any case, and we remove the coarsen
        ! flag
        do l = 1, N_sisters
            ! change status only for the existing sisters, mark temporary for leaf-wise as maybe this sister will be created soon
            if (sisters(l) /= -1) then
                if (lgt_block( sisters(l), IDX_REFINE_STS ) < -1) lgt_block( sisters(l), IDX_REFINE_STS )  = tmp_stay
            endif
        end do
    end if

end subroutine ensure_completeness_block