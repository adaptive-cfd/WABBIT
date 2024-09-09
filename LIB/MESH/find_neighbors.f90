!> \brief Given a direction XYZ as +1 / 0 / -1 for each coordinate,
!! find the neighbor and assign it to the correct index in hvy_neighbor
!
!> neighbor codes: \n
!  ---------------
!>   1- 56 : lvl_diff =  0  (same level)
!>  57-112 : lvl_diff = +1  (coarser neighbor)
!> 113-168 : lvl_diff = -1  (finer   neighbor)
!> For each range, the different 56 entries are:
!> 01-08 : X side (4-,4+)
!> 09-16 : Y-side (4-,4+)
!> 17-24 : Z-side (4-,4+)
!> 25-32 : X-Y edge (2--, 2+-, 2-+, 2++)
!> 33-40 : X-Z edge (2--, 2+-, 2-+, 2++)
!> 41-48 : Y-Z edge (2--, 2+-, 2-+, 2++)
!> 49-56 : corners (---, +--, -+-, ++-, --+, +-+, -++, +++)
! ********************************************************************************************
subroutine find_neighbor(params, hvyID_block, lgtID_block, dir, error, n_domain, search_overlapping, verbose)

    implicit none
    type (type_params), intent(in)      :: params                   !< user defined parameter structure
    integer(kind=ik), intent(in)        :: hvyID_block
    integer(kind=ik), intent(in)        :: lgtID_block
    !> direction for neighbor search - number where each digit represents a cardinal direction
    !> 652 -> first 6 (bottom, z-1), then 5 (north, x-1) then 2 (front, y-1) 
    integer(kind=ik), intent(in)        :: dir                      
    logical, intent(inout)              :: error
    integer(kind=2), intent(in)         :: n_domain(1:3)
    logical, intent(in)                 :: search_overlapping  !< for CVS multiple neighbors can coexist, so we search all of them
    logical, intent(in), optional       :: verbose  ! no matter the value, if it is present we print additional output
    
    integer(kind=ik)                    :: neighborDirCode_sameLevel
    integer(kind=ik)                    :: neighborDirCode_coarserLevel, tcFinerAppendDigit(4)
    integer(kind=ik)                    :: neighborDirCode_finerLevel(4)
    integer(kind=ik)                    :: k, lgtID_neighbor, level, tc_last, tree_ID
    integer(kind=tsize)                 :: tcb_Block, tcb_Neighbor, tcb_Virtual
    logical                             :: exists

    ! new variables
    integer(kind=ik)                    :: dir_dim(1:3), dir_free, i_dim, i_dig, apply_free, vary_tc(3)

    level      = lgt_block( lgtID_block, IDX_MESH_LVL )
    tree_ID    = lgt_block( lgtID_block, IDX_TREE_ID )
    neighborDirCode_sameLevel     = -1
    neighborDirCode_coarserLevel  = -1
    neighborDirCode_finerLevel    = -1
    tcFinerAppendDigit            = -1

    ! we have to init tcBlock before we set it elsewise fortran doesnt like it
    tcb_Block    = get_tc(lgt_block(lgtID_block, IDX_TC_1 : IDX_TC_2))
    ! last digit is used very often so we only extract it once
    tc_last = tc_get_digit_at_level_b(tcb_Block, dim=params%dim, level=level, max_level=params%Jmax)

    ! extract the direction for each dimension from direction for neighbor search - number where each digit represents a cardinal direction
    ! XYZ where each digit can be 0 (no change), 1 (for + direction) or 9 (for - direction)
    dir_dim(1) = mod(dir/100, 10)
    dir_dim(2) = mod(dir/10, 10)
    dir_dim(3) = mod(dir/1, 10)

    ! compute if we are looking at a face (2), edge (1) or corner (0) by how many directions are 0
    ! in total, we have 2**(dir_free - (params%dim==3)) possible finer neighbors or configurations for coarser neighbors
    dir_free = 2**count(dir_dim == 0)
    if (params%dim == 2) dir_free = dir_free / 2

    ! compute the digits that are free, do this by looping over the free directions
    tcFinerAppendDigit(1:4) = 0
    apply_free = 1  ! needed to switch between if a variable is the first or second "free" variable
    vary_tc = (/ 2, 1, 4/)  ! x change varies treecode by 2 and y by 1, this is a dumb convention
    do i_dim = 1,params%dim
        ! if this direction can vary (0), we let the treecodes vary by it
        if (dir_dim(i_dim) == 0) then
            do i_dig = 1,4
                tcFinerAppendDigit(i_dig) = tcFinerAppendDigit(i_dig) + vary_tc(i_dim)*mod((i_dig-1)/apply_free,2)
            enddo
            apply_free = apply_free +1
        ! if this direction is fixed and goes positive (1), we shift the treecodes
        elseif (dir_dim(i_dim) == 1) then
            tcFinerAppendDigit(1:4) = tcFinerAppendDigit(1:4) + vary_tc(i_dim)
        endif
    enddo

    ! now we need to find out the corresponding indices where we wanna start
    if (count(dir_dim == 0) == 2) then
        ! for sides, indices start at 1
        neighborDirCode_sameLevel = 1
        ! shift according to entries
        do i_dim = 1,3
            ! first x is free, then y, then z giving this shift
            if (dir_dim(i_dim) /= 0) neighborDirCode_sameLevel = neighborDirCode_sameLevel + 8*(i_dim-1)
            ! shift indices for positive direction
            if (dir_dim(i_dim) == 1) neighborDirCode_sameLevel = neighborDirCode_sameLevel + 4
        enddo
    elseif (count(dir_dim == 0) == 1) then
        ! for edges, indices start at 24+1
        neighborDirCode_sameLevel = 25
        ! shift according to entries
        apply_free = 1  ! needed to switch between if a variable is the first or second "fixed" variable
        do i_dim = 1,3
            if (dir_dim(i_dim) == 0) then
                ! first z varies, then y, then x giving this shift
                neighborDirCode_sameLevel = neighborDirCode_sameLevel + 8*(3-i_dim)
            ! shift indices for positive direction
            elseif (dir_dim(i_dim) == 1) then
                neighborDirCode_sameLevel = neighborDirCode_sameLevel + apply_free*2
                apply_free = apply_free + 1
            else
                apply_free = apply_free + 1
            endif
        enddo
    elseif (count(dir_dim == 0) == 0) then
        ! for corners, indices start at 48+1
        neighborDirCode_sameLevel = 49
        ! shift indices for positive direction
        do i_dim = 1,3
            if (dir_dim(i_dim) == 1) neighborDirCode_sameLevel = neighborDirCode_sameLevel + 2**(i_dim-1)
        enddo
    endif

    ! lets find the index for level-down neighbor, offset is not yet added
    do i_dig = 1, dir_free
        if (tc_last == tcFinerAppendDigit(i_dig)) then
            neighborDirCode_coarserLevel = neighborDirCode_sameLevel + i_dig - 1
        endif
    enddo

    ! debug
    ! write(*, '(A, i3, A, 4(i2), A, i2)') "Dir= ", dir, " append= ", tcFinerAppendDigit(1:4), " start_id= ", neighborDirCode_sameLevel


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 1) Check if we find a neighbor on the SAME LEVEL
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! calculate treecode for neighbor on same level
    call adjacent_wrapper_b(tcb_Block, tcb_Neighbor, dir, level=level, dim=params%dim, max_level=params%Jmax)

    ! check if (hypothetical) neighbor exists and if so find its lgtID
    call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level, tree_id=tree_ID, max_level=params%Jmax)

    if (exists) then
        ! we found the neighbor on the same level.
        hvy_neighbor( hvyID_block, neighborDirCode_sameLevel ) = lgtID_neighbor
        if (.not. search_overlapping) return
    endif


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2) Check if we find a neighbor on the FINER LEVEL
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Note: if those exists, then we always need to find all 4 (for 3D) or 2 of them
    ! We might be interested to search for coarser neighbors first because their is only one, however for CVS we might have an full tree grid
    ! but still only want to sync ANY values. Then coarsest neighbors should be chosen last as their values will be interpolated with WC=0
    ! However, if we did not find any same-lvl neighbor we should not find finer neighbors, it is therefore not too important but for clarity
    ! in this order
    if (level < params%Jmax) then
        do i_dig = 1, dir_free
            ! first neighbor virtual treecode, one level up
            tcb_Virtual = tc_set_digit_at_level_b(tcb_Block, tcFinerAppendDigit(i_dig), level=level+1, max_level=params%Jmax, dim=params%dim)

            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_wrapper_b(tcb_Virtual, tcb_Neighbor, dir, level=level+1, max_level=params%Jmax, dim=params%dim)
            ! check if (hypothetical) neighbor exists and if so find its lgtID
            call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level+1, tree_id=tree_ID, max_level=params%Jmax)

            if (exists) then
                hvy_neighbor( hvyID_block, neighborDirCode_sameLevel + i_dig-1  + 2*56) = lgtID_neighbor
            else
                exit  ! no need to serch for other ones if one already is not found
            endif

            ! we can only return here if we found all blocks
            if (.not. search_overlapping .and. i_dig == dir_free) return
        end do
    endif


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3) Check if we find a neighbor on the COARSER LEVEL (if that is possible at all)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! We did not find the neighbor on the same level, and now check on the coarser level.
    ! Just clear all levels before the lower one
    tcb_Neighbor = tc_clear_until_level_b(tcb_Neighbor, dim=params%dim, level=level-1, max_level=params%Jmax)
    ! only continue if coarser neighbor can at all exist, consider:
    ! a c E E
    ! b d E E
    ! Then to the right, block b cannot have a coarser neighbor
    if (neighborDirCode_coarserLevel /= -1 .and. level > params%Jmin) then
        ! check if (hypothetical) neighbor exists and if so find its lgtID
        call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level-1, tree_id=tree_ID, max_level=params%Jmax)

        if ( exists ) then
            ! neighbor is one level down (coarser)
            hvy_neighbor( hvyID_block, neighborDirCode_coarserLevel + 56 ) = lgtID_neighbor
            if (.not. search_overlapping) return
        endif
    endif

    ! we did not find a neighbor. that may be a bad grid error, or simply, there is none
    ! we have to find a neighbor for faces, edges and corners might have no neighbors if a coarser block is also at a face
    if (count(dir_dim == 0) == 2 .and. .not. exists .and. .not. search_overlapping&
        .and. ( ALL(params%periodic_BC) .or. maxval(abs(n_domain))==0)) then
        write(*, '("Rank: ", i0, ", found no neighbor in direction: ", i3, ", lgtID-", i6, " lvl-", i2, " TC-", i21, "-", b64.64)') &
            params%rank, dir, lgtID_block, level, tcb_Block, tcb_Block
        error = .true.
    endif


end subroutine