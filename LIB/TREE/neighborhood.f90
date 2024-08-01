!> \brief From a neighbourhood, compute the indices of the sender patch of thickness N_s (left end) or N_e (right end)
!> This function is kind of like the inverse of getting the ghost node patches, used for the sender patches inside the domain.
!> Additionally it has the option to extend those patches but does not do so correctly in the ghost point layer
!> Every call that also edits ghost patches (looking at you, CE!) should therefore also make sure to treat those points with ghost patches as well
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
subroutine get_indices_of_modify_patch(g, dim, relation, idx, N_xyz, N_s, N_e, g_m, g_p, lvl_diff)
    implicit none

    integer(kind=ik), intent(in)             :: g                   !> params%g
    integer(kind=ik), intent(in)             :: dim                 !> params%dim
    integer(kind=ik), intent(in)             :: relation            !> Which relation to apply manipulation
    integer(kind=ik), intent(out)            :: idx(1:2, 1:3)       !> Output of indices, first entry is l/r and second is dim
    integer(kind=ik), intent(in)             :: N_xyz(1:3)          !> Size of blocks including ghost points, should be size(block, i_dim)
    integer(kind=ik), intent(in)             :: N_s(1:3)             !> Number of points to modify at start, vec with entry for each dimension
    integer(kind=ik), intent(in)             :: N_e(1:3)             !> Number of points to modify at end, vec with entry for each dimension
    integer(kind=ik), intent(in), optional   :: g_p(1:3)             !> Ghost points left side, will be skipped
    integer(kind=ik), intent(in), optional   :: g_m(1:3)             !> Ghost points right side, will be skipped
    integer(kind=ik), intent(in), optional   :: lvl_diff             !> lvl_diff, special treatment for -1

    integer(kind=ik) :: i_dim, lvlDiff, relation_temp
    ! Indexes where to start or finish, short name because elsewise this list gets looong
    integer(kind=ik) :: I_s(1:3), I_e(1:3), gp(1:3), gm(1:3)

    ! for neighborhoods, +56 and +112 describe the level differences, as for CVS multiple neighbors can exist
    ! for patches we provide lvl_diff anyways so let's project it out of the relation
    if (relation > 56) then
        relation_temp = mod(relation-1, 56) +1
    else    
        relation_temp = relation
    endif

    lvlDiff = 0
    if (present(lvl_diff)) lvlDiff = lvl_diff

    ! set g_p g_m from params%g if not provided
    gp(1:3) = g
    gm = g
    if (present(g_p)) gp(1:3) = g_p(1:3)
    if (present(g_m)) gm(1:3) = g_m(1:3)

    ! Only select interior domain, example with g=1
    ! g g g g g           - - - - -
    ! g i i i g           - s s s -
    ! g i i i g           - s s s -
    ! g i i i g           - s s s -
    ! g g g g g           - - - - -
    do i_dim = 1, dim
        idx(1, i_dim) = 1 + gm(i_dim)
        idx(2, i_dim) = N_xyz(i_dim) - gp(i_dim)
    enddo

    ! now lets get to the select beast, this further limits / restricts the patch for the specific neighborhood
    ! example with g=1, N_s=N_e=1 and selecting neighborhood 9 or -x edge for 2D
    ! g g g g g           - - - - -
    ! g i i i g           - s - - -
    ! g i i i g           - s - - -
    ! g i i i g           - s - - -
    ! g g g g g           - - - - -
    if (any((/  1, 2, 3, 4,  25,26,29,30,33,34,37,38,  49,51,53,55 /) == relation_temp)) &
        idx(2, 1) = N_s(1)+gm(1)  ! -x
    if (any((/  5, 6, 7, 8,  27,28,31,32,35,36,39,40,  50,52,54,56 /) == relation_temp)) &
        idx(1, 1) = N_xyz(1) - N_e(1)-gp(1) + 1  ! +x
    if (any((/  9,10,11,12,  25,26,27,28,41,42,45,46,  49,50,53,54 /) == relation_temp)) &
        idx(2, 2) = N_s(2)+gm(2)  ! -y
    if (any((/ 13,14,15,16,  29,30,31,32,43,44,47,48,  51,52,55,56 /) == relation_temp)) &
        idx(1, 2) = N_xyz(2) - N_e(2)-gp(2) + 1  ! +y
    if (any((/ 17,18,19,20,  33,34,35,36,41,42,43,44,  49,50,51,52 /) == relation_temp)) &
        idx(2, 3) = N_s(3)+gm(3)  ! -z
    if (any((/ 21,22,23,24,  37,38,39,40,45,46,47,48,  53,54,55,56 /) == relation_temp)) &
        idx(1, 3) = N_xyz(3) - N_e(3)-gp(3) + 1  ! +z

    ! for leveldiff = -1 we have to restrict some patches to half the length (at edges)
    ! Additionally, the patch stretches further for the receiver corner patch which we send as well
    if (lvlDiff == -1) then
        ! +x border
        if (any(relation_temp == (/ 10,12,14,16,18,20,22,24,  42,44,46,48 /))) then
            idx(1,1) = N_xyz(1)/2 - N_e(1)+2 ! +x border
        endif
        ! -x border
        if (any(relation_temp == (/  9,11,13,15,17,19,21,23,  41,43,45,47 /))) then
            idx(2,1) = N_xyz(1)/2 + N_s(1)  ! -x border
        endif
        ! +y border
        if (any(relation_temp == (/  2, 4, 6, 8,19,20,23,24,  34,36,38,40 /))) then
            idx(1,2) = N_xyz(2)/2 - N_e(2)+2 ! +y border
        endif
        ! -y border
        if (any(relation_temp == (/  1, 3, 5, 7,17,18,21,22,  33,35,37,39 /))) then
            idx(2,2) = N_xyz(2)/2 + N_s(2)  ! -y border
        endif
        if (dim == 3) then
            ! +z border
            if (any(relation_temp == (/  3, 4, 7, 8,11,12,15,16,  26,28,30,32 /))) then
                idx(1,3) = N_xyz(3)/2 - N_e(3)+2 ! +z border
            endif
            ! -z border
            if (any(relation_temp == (/  1, 2, 5, 6, 9,10,13,14,  25,27,29,31 /))) then
                idx(2,3) = N_xyz(3)/2 + N_s(3)  ! -z border
            endif
        endif
    endif

end subroutine get_indices_of_modify_patch



!> \brief From a neighbourhood, compute the indices of the receiver ghost node patches of thickness gminus (left end) or gplus (right end)
!> For every relation and lvl_diff the necessary patch-sizes are computed. They start at the interior domain and then extend by gminus or gplus into the ghost node patch
!
! g g g g g g g g        g g g g g g g g
! g g g g g g g g        g s s s s s s g
! g g u u u u g g        g s u u u u s g
! g g u u u u g g        g s u u u u s g
! g g u u u u g g        g s u u u u s g
! g g u u u u g g        g s u u u u s g
! g g g g g g g g        g s s s s s s g
! g g g g g g g g        g g g g g g g g
!
! g: ghost u: interior s: selected
! here, Bs=4, params%g=2 and gminus=gplus=1
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
subroutine get_indices_of_ghost_patch( Bs, g, dim, relation, idx, gminus, gplus, lvl_diff)
    implicit none

    integer(kind=ik), intent(in)                    :: Bs(1:3)  !< params%Bs
    integer(kind=ik), intent(in)                    :: g        !< params%g 
    integer(kind=ik), intent(in)                    :: dim      !< params%dim
    !> idx
    integer(kind=ik), intent(inout)                 :: idx(2,3)
    !> neighborhood or family relation, id from dirs
    !! -8:-1 is mother/daughter relation
    !! 0 is full block relation, level
    !! 1:56*3 is neighborhood relation
    integer(kind=ik), intent(in)                    :: relation
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: lvl_diff, gminus, gplus

    integer(kind=ik) :: i_dim, relation_temp

    ! for neighborhoods, +56 and +112 describe the level differences, as for CVS multiple neighbors can exist
    ! for patches we provide lvl_diff anyways so let's project it out of the relation
    if (relation > 56) then
        relation_temp = mod(relation-1, 56) +1
    else    
        relation_temp = relation
    endif

    ! set 1 and not -1 (or anything else), because 2D bounds ignore 3rd dimension
    ! and thus cycle from 1:1
    idx(:,:) = 1

    ! preset to full block and then limit depending on relation, this is also for relation 0
    do i_dim = 1, dim
        idx(1,i_dim) = g+1
        idx(2,i_dim) = g+Bs(i_dim)
    enddo

    ! limit sides to specific edges
    ! every entry is first lvl_same faces, lvl_same edges, corners, lvl_diff faces, lvl_diff edges
    ! -x border
    if (any(relation_temp == (/  1, 2, 3, 4,  25,26,29,30,33,34,37,38,  49,51,53,55 /))) then
        idx(1,1) = g-gMinus+1
        idx(2,1) = g
    endif
    ! +x border
    if (any(relation_temp == (/  5, 6, 7, 8,  27,28,31,32,35,36,39,40,  50,52,54,56 /))) then
        idx(1,1) = Bs(1)+g+1
        idx(2,1) = Bs(1)+g+gplus
    endif
    ! -y border
    if (any(relation_temp == (/  9,10,11,12,  25,26,27,28,41,42,45,46,  49,50,53,54 /))) then
        idx(1,2) = g-gMinus+1
        idx(2,2) = g
    endif
    ! +y border
    if (any(relation_temp == (/ 13,14,15,16,  29,30,31,32,43,44,47,48,  51,52,55,56 /))) then
        idx(1,2) = Bs(2)+g+1
        idx(2,2) = Bs(2)+g+gplus
    endif
    ! -z border
    if (any(relation_temp == (/ 17,18,19,20,  33,34,35,36,41,42,43,44,  49,50,51,52 /))) then
        idx(1,3) = g-gMinus+1
        idx(2,3) = g
    endif
    ! +z border
    if (any(relation_temp == (/ 21,22,23,24,  37,38,39,40,45,46,47,48,  53,54,55,56 /))) then
        idx(1,3) = Bs(3)+g+1
        idx(2,3) = Bs(3)+g+gplus
    endif

    ! now for lvl_diff=1 and lvl_diff=-1 set the free directions
    ! for lvl_diff=-1 this limits it to 1/2 of the length of that direction
    ! for lvl_diff=+1 this extends it into the edges of that direction
    ! -x border change, variational +x
    if (any(relation_temp == (/ 10,12,14,16,18,20,22,24,  42,44,46,48 /))) then
        if (lvl_diff == +1) then
            idx(1,1) = g-gMinus+1
        elseif (lvl_diff == -1) then
            idx(1,1) = g+(Bs(1))/2 + 1
        endif
        idx(2,1) = Bs(1)+g
    endif
    ! +x border change, variational -x
    if (any(relation_temp == (/  9,11,13,15,17,19,21,23,  41,43,45,47 /))) then
        idx(1,1) = g+1
        if (lvl_diff == +1) then
            idx(2,1) = Bs(1)+g+gplus
        elseif (lvl_diff == -1) then
            idx(2,1) = g+(Bs(1))/2
        endif
    endif
    ! -y border change, variational +y
    if (any(relation_temp == (/  2, 4, 6, 8,19,20,23,24,  34,36,38,40 /))) then
        if (lvl_diff == +1) then
            idx(1,2) = g-gMinus+1
        elseif (lvl_diff == -1) then
            idx(1,2) = g+(Bs(2))/2 + 1
        endif
        idx(2,2) = Bs(2)+g
    endif
    ! +y border change, variational -y
    if (any(relation_temp == (/  1, 3, 5, 7,17,18,21,22,  33,35,37,39 /))) then
        idx(1,2) = g+1
        if (lvl_diff == +1) then
            idx(2,2) = Bs(2)+g+gplus
        elseif (lvl_diff == -1) then
            idx(2,2) = g+(Bs(2))/2
        endif
    endif
    if (dim == 3) then
        ! -z border change, variational +z
        if (any(relation_temp == (/  3, 4, 7, 8,11,12,15,16,  26,28,30,32 /))) then
            if (lvl_diff == +1) then
                idx(1,3) = g-gMinus+1
            elseif (lvl_diff == -1) then
                idx(1,3) = g+(Bs(3))/2 + 1
            endif
            idx(2,3) = Bs(3)+g
        endif
        ! +z border change, variational -z
        if (any(relation_temp == (/  1, 2, 5, 6, 9,10,13,14,  25,27,29,31 /))) then
            idx(1,3) = g+1
            if (lvl_diff == +1) then
                idx(2,3) = Bs(3)+g+gplus
            elseif (lvl_diff == -1) then
                idx(2,3) = g+(Bs(3))/2
            endif
        endif
    endif

    !---2D--3D--family relation, assume values on finer side are already WDed in mallat-ordering
    if (any(relation_temp == (/ -1, -2, -3, -4, -5, -6, -7, -8 /))) then   
        if (lvl_diff == 0) then
            return
        ! lvl_diff=+1 - family is mother and I am the daughter
        ! only transfer SC if this is the finer block, from Mallat ordering
        ! take into account that this depends on if g is even or odd
        elseif (lvl_diff == +1) then
            do i_dim = 1, dim
                idx(1,i_dim) = ceiling(g/2.0)+1
                idx(2,i_dim) = ceiling(g/2.0)+Bs(i_dim)/2
            enddo
        ! lvl_diff=-1 - family is daughter and I am the mother
        ! Transfer directly from domain wich is not decomposed so we only select the quarter or eigth part
        ! keep in mind: 0,1 varies in y- and 0,2 varies in x-direction
        else
            if (modulo(-relation_temp-1, 2) == 0) then
                idx(2,2) = g+Bs(1)/2
            else
                idx(1,2) = g+1+Bs(1)/2
            endif
            if (modulo((-relation_temp-1)/2, 2) == 0) then
                idx(2,1) = g+Bs(2)/2
            else
                idx(1,1) = g+1+Bs(2)/2
            endif
            if (dim == 3) then
                if (modulo((-relation_temp-1)/4, 2) == 0) then
                    idx(2,3) = g+Bs(3)/2
                else
                    idx(1,3) = g+1+Bs(3)/2
                endif
            endif
        endif
    endif

end subroutine get_indices_of_ghost_patch



!> Invert the relation between sender and receiver. With that the directional variables are inverted and free directions stay the same.
!! So a +x-y-- goes to -x+y-- where -- as z is still free and can be +z or -z variation for finer blocks or as configuration for coarser block
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
subroutine inverse_relation(relation, relation_inverse)
    implicit none

    integer(kind=ik), intent(in)                    :: relation
    integer(kind=ik), intent(out)                   :: relation_inverse
    integer(kind=ik) :: i_dim

    ! full block and family relations invert to themself
    if (relation <= 0) then
        relation_inverse = relation
        return
    endif

    !  0-24: faces with one directional change, every 4
    ! 25-48: edges with two directional changes, every 4 and 2
    ! 49-56: corners with three directional changes, every 4, 2 and 1
    relation_inverse = mod(relation-1, 56)+1
    ! invert lvl_diff being +56 or +112
    if (relation > 56) relation_inverse = relation_inverse + (3- (relation-1) / 56)*56
    do i_dim = 3, 3-mod(relation-1,56)/24, -1
        if (mod(relation-1, 2**i_dim) >= 2**(i_dim-1)) then
            relation_inverse = relation_inverse - 2**(i_dim-1)
        else
            relation_inverse = relation_inverse + 2**(i_dim-1)
        endif
    enddo
end subroutine inverse_relation



!> \brief Convert neighborhood to corresponding patch
!> The selection of ID is taken from the neighborhood relation, the structure goes from -2+, x2z, face2edge2corner.
!! Which numbers are taken is actually not important as long as it is consistent and unique
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
subroutine neighborhood2patchID(neighborhood, patchID)
    implicit none

    integer(kind=ik), intent(in)   :: neighborhood  !< neighborhood of significant patch
    integer(kind=ik), intent(out)  :: patchID       !< location of neighbor, where is the neighbor physically?
    integer(kind=ik) :: n_temp

    n_temp = mod(neighborhood-1, 56)+1  ! project level diffs out of the neighborhood

    ! faces - divide by 4 and add ofset
    if (n_temp <= 24) then
        patchID = (n_temp-1) / 4 + 1
    ! edges, divide by 2 and add offset
    elseif (n_temp <= 48) then
        patchID = (n_temp-25) / 2 + 1+6
    ! corners, leave them as they are
    elseif (n_temp <= 56) then
        patchID = (n_temp-49) + 1+6+12
    else
        call abort(20240723, "You ended up in the wrong neighborhood")
    endif
end 



!> \brief print to output the neighborhood array in a more readable fashion
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
subroutine write_neighborhood_info(hvy_neighbor, dim)
    implicit none

    integer(kind=ik), intent(in)  :: hvy_neighbor(1:3*56)  !< hvy_neighbor(hvy_id, :)
    integer(kind=ik), intent(in)  :: dim  !< params%dim
    integer(kind=ik) i_n, i_l
    character(len=200)  :: write_s
    logical :: found_n = .false.

    ! sides - indices 01-24
    write_s = ""
    write(write_s, '(A8)') "Sides:"
    do i_n = 1,24  ! every side has 4 possible finer blocks or coarser configurations
        do i_l = 0,2  ! for lvl_diff = 0, -1, 1
            if (hvy_neighbor(i_n + i_l*56) /= -1) then
                write(write_s(len_trim(write_s)+1:), '(A, i0)') " ", hvy_neighbor(i_n + i_l*56)
                found_n = .true.
            endif
        enddo
        ! print info which sides
        if (i_n == 4 .and. found_n) then
            write(write_s(len_trim(write_s)+1:), '(A)') " (-x)"
        elseif (i_n == 8 .and. found_n) then
            write(write_s(len_trim(write_s)+1:), '(A)') " (+x)"
        elseif (i_n == 12 .and. found_n) then
            write(write_s(len_trim(write_s)+1:), '(A)') " (-y)"
        elseif (i_n == 16 .and. found_n) then
            write(write_s(len_trim(write_s)+1:), '(A)') " (+y)"
        elseif (i_n == 20 .and. found_n) then
            write(write_s(len_trim(write_s)+1:), '(A)') " (-z)"
        elseif (i_n == 24 .and. found_n) then
            write(write_s(len_trim(write_s)+1:), '(A)') " (+z)"
        endif
        if (mod(i_n, 4) == 0) found_n = .false.
    enddo
    ! now print output
    write(*, '(A)') write_s(1:len_trim(write_s))

    ! edges - indices 25-48
    if (any(hvy_neighbor(25:48) /= -1) .or. any(hvy_neighbor(25+56:48+56) /= -1) .or. any(hvy_neighbor(25+56:48+56) /= -1)) then
        write_s = ""
        write(write_s, '(A8)') "Edges:"
        do i_n = 25,48  ! every side has 2 possible finer blocks or coarser configurations
            do i_l = 0,2  ! for lvl_diff = 0, -1, 1
                if (hvy_neighbor(i_n + i_l*56) /= -1) then
                    write(write_s(len_trim(write_s)+1:), '(A, i0)') " ", hvy_neighbor(i_n + i_l*56)
                    found_n = .true.
                endif
            enddo
            ! print info which sides
            if (i_n == 26 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x-y)"
            elseif (i_n == 28 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x-y)"
            elseif (i_n == 30 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x+y)"
            elseif (i_n == 32 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x+y)"
            elseif (i_n == 34 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x-z)"
            elseif (i_n == 36 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x-z)"
            elseif (i_n == 38 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x+z)"
            elseif (i_n == 40 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x+z)"
            elseif (i_n == 42 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-y-z)"
            elseif (i_n == 44 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+y-z)"
            elseif (i_n == 46 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-y+z)"
            elseif (i_n == 48 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+y+z)"
            endif
            if (mod(i_n, 2) == 0) found_n = .false.
        enddo
        ! now print output
        write(*, '(A)') write_s(1:len_trim(write_s))
    endif

    ! corners - indices 49-56
    if (any(hvy_neighbor(49:56) /= -1) .or. any(hvy_neighbor(49+56:56+56) /= -1) .or. any(hvy_neighbor(49+56:56+56) /= -1)) then
        write_s = ""
        write(write_s, '(A8)') "Corners:"
        do i_n = 49,56  ! every side has 2 possible finer blocks or coarser configurations
            do i_l = 0,2  ! for lvl_diff = 0, -1, 1
                if (hvy_neighbor(i_n + i_l*56) /= -1) then
                    write(write_s(len_trim(write_s)+1:), '(A, i0)') " ", hvy_neighbor(i_n + i_l*56)
                    found_n = .true.
                endif
            enddo
            ! print info which sides
            if (i_n == 49 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x-y-z)"
            elseif (i_n == 50 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x-y-z)"
            elseif (i_n == 51 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x+y-z)"
            elseif (i_n == 52 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x+y-z)"
            elseif (i_n == 53 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x-y+z)"
            elseif (i_n == 54 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x-y+z)"
            elseif (i_n == 55 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (-x+y+z)"
            elseif (i_n == 56 .and. found_n) then
                write(write_s(len_trim(write_s)+1:), '(A)') " (+x+y+z)"
            endif
            if (mod(i_n, 2) == 0) found_n = .false.
        enddo
        ! now print output
        write(*, '(A)') write_s(1:len_trim(write_s))
    endif
end subroutine