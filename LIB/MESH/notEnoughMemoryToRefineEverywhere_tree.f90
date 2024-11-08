
! check if we still have enough memory left: for very large simulations
! we cannot affort to have them fail without the possibility to resume them
logical function notEnoughMemoryToRefineEverywhere_tree(params, tree_ID, time)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none
    type (type_params), intent(inout) :: params
    integer(kind=ik), intent(in)      :: tree_ID
    real(kind=rk), intent(in)         :: time  ! used to differentiate for force_maxlevel_dealiasing

    integer :: hvy_n_max, k, hvy_n_afterRefinement, hvyID, lgtID

    notEnoughMemoryToRefineEverywhere_tree = .false.
    params%out_of_memory = .false.

    ! without adaptivity, this routine makes no sense, as the memory is constant
    ! the run either crashes right away or never
    if (params%adapt_tree .eqv. .false. .or. params%Jmax == params%Jmin) then
        return
    endif

    ! this is the available maximum number of active blocks on a single mpirank
    hvy_n_max = nint(0.95_rk * real(params%number_blocks, kind=rk) )

    ! remove blocks already used for mask etc
    hvy_n_max = hvy_n_max - sum(hvy_n(2:size(hvy_n)))

    ! Routine is called BEFORE refinement. figure out how many blocks
    ! that will be after refinement
    if (params%force_maxlevel_dealiasing .and. time > 0.0) then
        ! with dealiasing, the relation is very simple, as all blocks have to be on level JMax-1 or smaller
        ! however, this might not be true when a file is read in and has not been treated yet
        hvy_n_afterRefinement = hvy_n(tree_ID) * (2**params%dim)

    else
        hvy_n_afterRefinement = 0
        do k = 1, hvy_n(tree_ID)
            hvyID = hvy_active(k, tree_ID)
            call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

            if (lgt_block( lgtID, IDX_MESH_LVL ) < params%Jmax) then
                ! this block can be refined (and will be) (it is refined to it has 8 children
                ! but disappears, so 7 new blocks)
                hvy_n_afterRefinement = hvy_n_afterRefinement + (2**params%dim) - 1
            else
                ! this block is at maximum refinement and will not be refined, but not deleted neither
                hvy_n_afterRefinement = hvy_n_afterRefinement + 1
            endif
        enddo
    endif

    ! the factor 8/7 or 4/3 comes from the full wavelet transformation being done in adapt_tree
    hvy_n_afterRefinement = ceiling(hvy_n_afterRefinement * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk))

    params%out_of_memory = (hvy_n_afterRefinement > hvy_n_max)

    if (params%out_of_memory ) then
        write(*,'("OUT OF MEMORY: rank=",i5," Nb_before_refinement=",i7," Nb_after_refinement=",i7,"> Nb_limit=",i7, " lgt=",i7)') &
        params%rank, hvy_n(tree_ID), hvy_n_afterRefinement, hvy_n_max, sum(lgt_n)
    endif

    call MPI_Allreduce(MPI_IN_PLACE, params%out_of_memory, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, k )


    notEnoughMemoryToRefineEverywhere_tree = params%out_of_memory

end function
