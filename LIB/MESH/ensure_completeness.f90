!> \brief sets refinement status to -2 for all sister blocks, if coarsening is possible
!! input:    - light data array \n
!! output:   - light data array
! ********************************************************************************************
!> \image html completeness.svg "Ensure Completeness" width=400

subroutine ensure_completeness( params, lgt_id, sisters, mark_TMP_flag, check_daughters )

    implicit none

    type (type_params), intent(in)  :: params                !< user defined parameter structure
    integer(kind=ik), intent(inout) :: sisters(:)            !< light data array
    integer(kind=ik), intent(in)    :: lgt_id                !< Concerned Block
    logical, intent(in), optional   :: mark_TMP_flag         !< Set refinement of completeness to 0 or temporary flag
    logical, intent(in), optional   :: check_daughters          !< For overfull CVS grids completeness has to check the mother as well

    integer(kind=ik)                :: Jmax, k, l, hvy_id, lgt_daughter
    integer(kind=ik)                :: N_sisters, markTMPflag, lgt_sisters(1:2**params%dim)  ! loop variables
    logical                         :: checkDaughters

    Jmax = params%Jmax
    N_sisters = size(sisters)
    call lgt2hvy(hvy_id, lgt_id, params%rank, params%number_blocks)  ! needed for mother check

    markTMPflag = 0
    if (present(mark_TMP_flag)) then
        if (mark_TMP_flag) markTMPflag = REF_TMP_TREATED_COARSEN
    endif
    checkDaughters = .false.
    if (present(check_daughters)) checkDaughters = check_daughters

    ! if all sisters exists, then the array should not contain -1
    if ( all(sisters(1:N_sisters) /= -1) ) then

        lgt_sisters = lgt_block(sisters(1:N_sisters), IDX_REFINE_STS )

        ! first check : any block wants to stay so all blocks have to stay
        if ( any(lgt_sisters(:) == 0) ) then
            lgt_block( sisters(1:N_sisters), IDX_REFINE_STS )  = 0

        ! second check: any block has temp_gradedness flag > 1 - all blocks that want to stay have to wait with temp flag
        elseif ( any(lgt_sisters(:) > 1)) then
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

        ! for CVS grids we need to check if we have a daughter and if so, if she'd like to coarsen as well
        ! all sisters share the same mother, so this hvy operation will be done correctly on all processors
        if (checkDaughters .and. any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) then
            ! family takes care of one another - if any wants to stay then me and all my sisters have to stay as well
            do l = 2+2**params%dim, 1+2**(params%dim+1)
                lgt_daughter = hvy_family(hvy_ID, l)
                if (lgt_daughter /= -1) then
                    if (lgt_block(lgt_daughter, IDX_REFINE_STS) == 0) then
                        lgt_block( sisters(1:N_sisters), IDX_REFINE_STS ) = 0
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
                lgt_block( sisters(l), IDX_REFINE_STS )  = markTMPflag
            endif
        end do
    end if

end subroutine ensure_completeness



!> \brief For CVS, non-leaf blocks are created from finest to coarsest level. Meaning that every non-leaf block will always have a medium-lvl neighbor if they have a coarse-lvl neighbor.
!! Leaf blocks can always spawn their mother, however for non-leaf blocks to spawn their mother with neighbors for all patches, they need to have all possible
!! medium-lvl neighbors. If that is not the case, the mother would not have all neighbors to sync with for all her patches, so this block has to wait one iteration.
subroutine ensure_mother_neighbors( params, hvy_id, mark_TMP_flag)

    implicit none

    type (type_params), intent(in)  :: params                !< user defined parameter structure
    integer(kind=ik), intent(in)    :: hvy_id                !< Concerned Block
    logical, intent(in), optional   :: mark_TMP_flag         !< Set refinement of completeness to 0 or temporary flag

    integer(kind=ik)                :: i_n, lgt_id, markTMPflag
    logical                         :: mother_finds_all_neighbors

    markTMPflag = 0
    if (present(mark_TMP_flag)) then
        ! if this is REF_TMP_TREATED_COARSEN, then ensure_completeness would just overwrite it so we set it to a different value
        if (mark_TMP_flag) markTMPflag = REF_TMP_GRADED_STAY
    endif

    call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)

    mother_finds_all_neighbors = .true.
    ! loop over all possible patches that a block can have
    do i_n = 1,56
        ! skip patches not available for 2D
        if (params%dim == 2 .and. is_3D_neighbor(i_n)) cycle

        ! check if block is non-leaf and if there is a medium-lvl neighbor in those patches
        if (all(hvy_neighbor(hvy_id, np_l(i_n, 0):np_u(i_n, 0)) == -1) .and. any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) then
            mother_finds_all_neighbors = .false.
        endif
    enddo

    ! non-leaf blocks need all medium-lvl neighbors in order to spawn their mother, so if thats not the case it needs to wait
    if (.not. mother_finds_all_neighbors) then
        lgt_block( lgt_id, IDX_REFINE_STS ) = markTMPflag
    endif

end subroutine ensure_mother_neighbors
