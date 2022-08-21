!> \brief return coarsest level in active block list
! ********************************************************************************************

function minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )
    implicit none

    integer(kind=ik), intent(in)           :: lgt_block(:, :)           !> light data array
    integer(kind=ik), intent(in)           :: tree_ID
    integer(kind=ik), intent(in), optional :: lgt_active(:,:)             !> list of active blocks (light data)
    integer(kind=ik), intent(in), optional :: lgt_n(:)                     !> number of active blocks (light data)

    integer(kind=ik)                       :: minActiveLevel_tree          ! return value
    integer(kind=ik)                       :: k, Jmin, max_treelevel    ! loop variables

    max_treelevel = size( lgt_block, 2) - EXTRA_LGT_FIELDS
    Jmin = max_treelevel

    if (present(lgt_active) .and. present(lgt_n)) then
        ! call with active lists (to be preferred, much faster)
        ! loop over all active blocks
        do k = 1, lgt_n(tree_ID)
            Jmin = min(Jmin, lgt_block( lgt_active(k, tree_ID), max_treelevel + IDX_MESH_LVL))
        end do

    else
        ! call without active lists
        do k = 1, size(lgt_block, 1)
            if (lgt_block( k, max_treelevel+IDX_TREE_ID) == tree_ID) then
                if (lgt_block( k, max_treelevel + IDX_MESH_LVL)>=0) then
                    ! skip blocks marked with -1 (they are inactive)
                    Jmin = min(Jmin, lgt_block( k, max_treelevel+IDX_MESH_LVL))
                endif
            endif
        end do
    end if

    minActiveLevel_tree = Jmin
end function



function maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )

    implicit none

    integer(kind=ik), intent(in)           :: lgt_block(:, :)           !> light data array
    integer(kind=ik), intent(in)           :: tree_ID
    integer(kind=ik), intent(in), optional :: lgt_active(:,:)             !> list of active blocks (light data)
    integer(kind=ik), intent(in), optional :: lgt_n(:)                     !> number of active blocks (light data)

    integer(kind=ik)                    :: maxActiveLevel_tree           ! return value
    integer(kind=ik)                    :: k, Jmax, max_treelevel     ! loop variables

    max_treelevel = size( lgt_block, 2) - EXTRA_LGT_FIELDS
    Jmax = 0

    if (present(lgt_active) .and. present(lgt_n)) then
        ! call with active lists (to be preferred, much faster)
        ! loop over all active blocks
        do k = 1, lgt_n(tree_ID)
            Jmax = max(Jmax, lgt_block( lgt_active(k, tree_ID), max_treelevel+IDX_MESH_LVL))
        end do

    else
        ! call without active lists
        do k = 1, size(lgt_block, 1)
            if (lgt_block( k, max_treelevel+IDX_TREE_ID) == tree_ID) then
                Jmax = max(Jmax, lgt_block( k, max_treelevel+IDX_MESH_LVL))
            endif
        end do
    end if

    maxActiveLevel_tree = Jmax
end function
