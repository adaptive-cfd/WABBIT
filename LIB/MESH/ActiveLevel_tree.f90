!> \brief return coarsest level in active block list
! ********************************************************************************************

function minActiveLevel_tree( tree_ID, use_active_list )
    implicit none

    integer(kind=ik), intent(in)           :: tree_ID
    logical, intent(in), optional          :: use_active_list

    integer(kind=ik)                       :: minActiveLevel_tree          ! return value
    integer(kind=ik)                       :: k, Jmin, max_treelevel    ! loop variables
    logical :: use_active_list2

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! max_treelevel = size( lgt_block, 2) - EXTRA_LGT_FIELDS
    ! this is not adapted to params%Jmax but it only matters when all blocks are inactive
    max_treelevel = maxdigits
    Jmin = max_treelevel

    use_active_list2 = .true.
    if (present(use_active_list)) use_active_list2 = use_active_list

    if (use_active_list2) then
        ! call with active lists (to be preferred, much faster)
        ! loop over all active blocks
        do k = 1, lgt_n(tree_ID)
            Jmin = min(Jmin, lgt_block( lgt_active(k, tree_ID), IDX_MESH_LVL))
        end do

    else
        ! call without active lists
        do k = 1, size(lgt_block, 1)
            if (lgt_block( k, IDX_TREE_ID) == tree_ID) then
                if (lgt_block( k, IDX_MESH_LVL)>=0) then
                    ! skip blocks marked with -1 (they are inactive)
                    Jmin = min(Jmin, lgt_block( k, IDX_MESH_LVL))
                endif
            endif
        end do
    end if

    minActiveLevel_tree = Jmin
end function



function maxActiveLevel_tree( tree_ID, use_active_list )

    implicit none

    integer(kind=ik), intent(in)           :: tree_ID
    logical, intent(in), optional          :: use_active_list

    integer(kind=ik)                    :: maxActiveLevel_tree           ! return value
    integer(kind=ik)                    :: k, Jmax                       ! loop variables
    logical :: use_active_list2

    Jmax = 0

    use_active_list2 = .true.
    if (present(use_active_list)) use_active_list2 = use_active_list

    if (use_active_list2) then
        ! call with active lists (to be preferred, much faster)
        ! loop over all active blocks
        do k = 1, lgt_n(tree_ID)
            Jmax = max(Jmax, lgt_block( lgt_active(k, tree_ID), IDX_MESH_LVL))
        end do

    else
        ! call without active lists
        do k = 1, size(lgt_block, 1)
            if (lgt_block( k, IDX_TREE_ID) == tree_ID) then
                Jmax = max(Jmax, lgt_block( k, IDX_MESH_LVL))
            endif
        end do
    end if

    maxActiveLevel_tree = Jmax
end function
