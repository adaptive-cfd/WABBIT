!> \brief Given a direction XYZ as +1 / 0 / -1 for each coordinate,
!! find the neighbor and assign it to the correct index in hvy_neighbor
!
!> neighbor codes: \n
!  ---------------
!>   1- 56 : lvl_diff =  0  (same level)
!>  57-112 : lvl_diff = +1  (coarser neighbor)
!> 113-168 : lvl_diff = -1  (finer   neighbor)
!> For each range, the different 56 entries are:
!> 01-08 : X side
!> 09-16 : Y-side
!> 17-24 : Z-side
!> 25-32 : X-Y edge
!> 33-40 : X-Z edge
!> 41-48 : Y-Z edge
!> 49-56 : corners
! ********************************************************************************************
subroutine find_neighbor(params, hvyID_block, lgtID_block, Jmax, dir, error, n_domain, verbose)

    implicit none
    type (type_params), intent(in)      :: params                   !< user defined parameter structure
    integer(kind=ik), intent(in)        :: hvyID_block
    integer(kind=ik), intent(in)        :: lgtID_block
    integer(kind=ik), intent(in)        :: Jmax
    !> direction for neighbor search - number where each digit represents a cardinal direction
    !> 652 -> first 6 (bottom, z-1), then 5 (north, x-1) then 2 (front, y-1) 
    integer(kind=ik), intent(in)        :: dir                      
    logical, intent(inout)              :: error
    integer(kind=2), intent(in)         :: n_domain(1:3)
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
    tc_last = tc_get_digit_at_level_b(tcb_Block, dim=params%dim, level=level, max_level=Jmax)

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
    ! call adjacent(tcBlock, tcNeighbor, dir, level, Jmax, params%dim)
    call adjacent_wrapper_b(tcb_Block, tcb_Neighbor, dir, level=level, dim=params%dim, max_level=Jmax)

    ! check if (hypothetical) neighbor exists and if so find its lgtID
    call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level, tree_id=tree_ID, max_level=params%Jmax)

    if (exists) then
        ! we found the neighbor on the same level.
        hvy_neighbor( hvyID_block, neighborDirCode_sameLevel ) = lgtID_neighbor
        return
    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2) Check if we find a neighbor on the COARSER LEVEL (if that is possible at all)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! We did not find the neighbor on the same level, and now check on the coarser level.
    ! Just clear all levels before the lower one
    tcb_Neighbor = tc_clear_until_level_b(tcb_Neighbor, dim=params%dim, level=level-1, max_level=Jmax)
    ! only continue if coarser neighbor can at all exist, consider:
    ! a c E E
    ! b d E E
    ! Then to the right, block b cannot have a coarser neighbor
    if (neighborDirCode_coarserLevel /= -1) then
        ! check if (hypothetical) neighbor exists and if so find its lgtID
        call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level-1, tree_id=tree_ID, max_level=params%Jmax)

        if ( exists ) then
            ! neighbor is one level down (coarser)
            hvy_neighbor( hvyID_block, neighborDirCode_coarserLevel + 56 ) = lgtID_neighbor
            return
        endif
    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3) Check if we find a neighbor on the FINER LEVEL
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Note: there are several neighbors possible on finer levels (up to four in 3D!)
    ! loop over all 4 possible neighbors
    if (level < Jmax) then
        do i_dig = 1, dir_free
            ! first neighbor virtual treecode, one level up
            tcb_Virtual = tc_set_digit_at_level_b(tcb_Block, tcFinerAppendDigit(i_dig), level=level+1, max_level=Jmax, dim=params%dim)

            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_wrapper_b(tcb_Virtual, tcb_Neighbor, dir, level=level+1, max_level=Jmax, dim=params%dim)
            ! check if (hypothetical) neighbor exists and if so find its lgtID
            call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level+1, tree_id=tree_ID, max_level=params%Jmax)

            if (exists) then
                hvy_neighbor( hvyID_block, neighborDirCode_sameLevel + i_dig-1  + 2*56) = lgtID_neighbor
            end if

            ! we did not find a neighbor. that may be a bad grid error, or simply, there is none
            ! we have to find a neighbor for faces, edges and corners might have no neighbors if a coarser block is also at a face
            ! because symmetry conditions are used, check for .not. error to only print it once
            if (count(dir_dim == 0) == 2 .and. .not. error) then
                if (.not. exists .and. ( ALL(params%periodic_BC) .or. maxval(abs(n_domain))==0)) then
                    call adjacent_wrapper_b(tcb_Block, tcb_Virtual, dir, level=level, dim=params%dim, max_level=Jmax)
                    write(*, '("Rank: ", i0, ", found no neighbor in direction: ", i0, ", lgtID-", i0, " lvl-", i0, " TC-", i0, "-", b64.64, A, "Checked same-lvl TC-", i0, "-", b64.64, " and lower-lvl TC-", i0, "-", b64.64)') &
                        params%rank, dir, lgtID_block, level, tcb_Block, tcb_Block, NEW_LINE('a'), tcb_Virtual, tcb_Virtual, tcb_Neighbor, tcb_Neighbor
                    error = .true.
                endif
            endif
        end do
    endif


end subroutine