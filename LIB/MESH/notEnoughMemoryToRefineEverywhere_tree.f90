
! check if we still have enough memory left: for very large simulations
! we cannot affort to have them fail without the possibility to resume them
logical function notEnoughMemoryToRefineEverywhere_tree(params, tree_ID)
    implicit none
    type (type_params), intent(inout) :: params
    integer(kind=ik), intent(in)      :: tree_ID

    integer :: lgt_n_max, k, lgt_n_afterRefinement

    notEnoughMemoryToRefineEverywhere_tree = .false.
    params%out_of_memory = .false.

    ! without adaptivity, this routine makes no sense, as the memory is constant
    ! the run either crashes right away or never
    if (params%adapt_tree .eqv. .false.) then
        notEnoughMemoryToRefineEverywhere_tree = .false.
        return
    endif

    ! this is the available maximum number of active blocks.
    lgt_n_max = params%number_blocks*params%number_procs

    ! remove blocks already used for mask etc
    lgt_n_max = lgt_n_max - sum(lgt_n(2:size(lgt_n)))

    ! safety margin. inhomogenoeus distribution
    ! of blocks can make trouble. Therefore, we raise the alarm earlier.
    lgt_n_max = lgt_n_max * 9 / 10


    ! however, this routine is called BEFORE refinement. figure out how many blocks
    ! that we be after refinement
    if (params%force_maxlevel_dealiasing) then
        ! with dealiasing, the relation is very simple.
        lgt_n_afterRefinement = lgt_n(tree_ID) * (2**params%dim)
    else
        lgt_n_afterRefinement = 0
        do k = 1, lgt_n(tree_ID)
            if (lgt_block( lgt_active(k, tree_ID), params%Jmax + IDX_MESH_LVL ) < params%Jmax) then
                ! this block can be refined (and will be) (it is refined to it has 8 children
                ! but disappears, so 7 new blocks)
                lgt_n_afterRefinement = lgt_n_afterRefinement + (2**params%dim)-1
            else
                ! this block is at maximum refinement and will not be refined, but not be deleted neither
                lgt_n_afterRefinement = lgt_n_afterRefinement + 1
            endif
        enddo

    endif


    if ( lgt_n_afterRefinement > lgt_n_max) then
        ! oh-oh.
        notEnoughMemoryToRefineEverywhere_tree = .true.

        if (params%rank==0) then
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("Nblocks_total=",i7," Nblocks_fluid=",i7)') sum(lgt_n), lgt_n(tree_ID)
            write(*,'("Nblocks_allocated=",i7," Nblocks_fluid_available=",i7)') params%number_blocks*params%number_procs, lgt_n_max
            write(*,'("Nblocks_after_refinement=",i7)') lgt_n_afterRefinement
            write(*,'("Nblocks_after_refinement=",i7)') lgt_n_afterRefinement
            write(*,'("-------------_> oh-oh.")')
            write(*,'("-----------OUT OF MEMORY--------------")')
            write(*,'("-----------OUT OF MEMORY--------------")')
        endif

        params%out_of_memory = .true.
    endif


end function
