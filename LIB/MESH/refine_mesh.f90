!> \brief Refine the mash (tag / interpolate / update lists)
!
!> \details This routine first sets the refinement flag for all blocks to +1
!! and then executes the refinement directly. Blocks that cannot be refined because they
!! are already on the finest allowed level are unaltered.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!!
!! Note we assume, on call, that active lists / neighbor are up-to-date
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
! ********************************************************************************************

subroutine refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, indicator, tree_ID  )

    use module_indicators

    implicit none

    type (type_params), intent(in)         :: params                      !> user defined parameter structure
    integer(kind=ik), intent(inout)        :: lgt_block(:, :)             !> light data array
    real(kind=rk), intent(inout)           :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    integer(kind=ik), intent(inout)        :: hvy_neighbor(:,:)           !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)        :: lgt_active(:)               !> list of active blocks (light data)
    integer(kind=ik), intent(inout)        :: lgt_n                       !> number of active blocks (light data)
    integer(kind=tsize), intent(inout)     :: lgt_sortednumlist(:,:)      !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(inout)        :: hvy_active(:)               !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_n                       !> number of active blocks (heavy data)
    character(len=*), intent(in)           :: indicator                   !> how to choose blocks for refinement
    integer(kind=ik), intent(in)           :: tree_ID

    ! cpu time variables for running time calculation
    real(kind=rk)                          :: t0, t1, t2, t_misc
    integer(kind=ik)                       :: k, hvy_n_afterRefinement, lgt_id
    ! start time
    t0 = MPI_Wtime()
    t_misc = 0.0_rk


    !> (a) loop over the blocks and set their refinement status.
    t1 = MPI_Wtime()
    call refinement_indicator( params, lgt_block, lgt_active, lgt_n, hvy_block, hvy_active, hvy_n, indicator )
    call toc( "refine_mesh (refinement_indicator)", MPI_Wtime()-t1 )


    !> (b) check if block has reached maximal level, if so, remove refinement flags
    t1 = MPI_Wtime()
    call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
    call toc( "refine_mesh (respect_min_max_treelevel)", MPI_Wtime()-t1 )


    !> (c) ensure gradedness of mesh. If the refinement is done everywhere, there is
    !! no way gradedness can be damaged, so we skip the call in this case. However,
    !! in all other indicators, this step is very important.
    t1 = MPI_Wtime()
    if ( indicator /= "everywhere" ) then
      call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
      lgt_sortednumlist, hvy_active, hvy_n )
    endif
    call toc( "refine_mesh (ensure_gradedness)", MPI_Wtime()-t1 )


    !---------------------------------------------------------------------------
    ! check if the refinement step will suceed or if we will run out of memory
    !---------------------------------------------------------------------------
    ! NOTE: in the simulation part, a second check is performed BEFORE this routine
    ! is called. It then aborts the time loop and dumps a backup for late resuming.
    ! If the following check fails, then this will be in postprocessing, and it will kill the code.
    hvy_n_afterRefinement = 0
    do k = 1, hvy_n
        ! loop over hvy is safer, because we can tell for each mpirank if it fails or not
        ! so we detect also problems arising from load-imbalancing.
        call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)

        if ( lgt_block(lgt_id, params%max_treelevel+IDX_REFINE_STS) == +1) then
            ! this block will be refined, so it creates 2**D new blocks but is deleted
            hvy_n_afterRefinement = hvy_n_afterRefinement + (2**params%dim - 1)
        else
            ! this block will not be refined, so it remains.
            hvy_n_afterRefinement = hvy_n_afterRefinement + 1
        endif
    enddo

    ! oh-oh case.
    if ( hvy_n_afterRefinement > params%number_blocks ) then
        write(*,'("On rank:",i5," hvy_n=",i5," will refine to ",i5," but limit is ",i5)') &
        params%rank, hvy_n, hvy_n_afterRefinement, params%number_blocks

        call abort (1909181827,"[refine_mesh.f90]: The refinement step will fail (not enough memory).")
    endif
    !---------------------------------------------------------------------------


    !> (d) execute refinement, interpolate the new mesh. All blocks go one level up
    !! except if they are already on the highest level. Note that those blocks have
    !! the status +11
    t1 = MPI_Wtime()
    if ( params%dim == 3 ) then
        call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
    else
        call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_active, hvy_n )
    end if
    call toc( "refine_mesh (refinement_execute)", MPI_Wtime()-t1 )


    !> (e) as the grid changed now with the refinement, we have to update the list of
    !! active blocks so other routines can loop just over these active blocks
    !! and do not have to ensure that the active list is up-to-date
    ! update list of sorted nunmerical treecodes, used for finding blocks
    t2 = MPI_wtime()
    call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n, tree_ID)
    t_misc = MPI_wtime() - t2


    !> (f) At this point the refinement is done. Since not all blocks are refined, namely only those
    !! that were not on Jmax, Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now. EXCEPTION: if the flag force_maxlevel_dealiasing is set
    !! then we force blocks on Jmax to coarsen, even if their details are large. Hence, each refinement
    !! step is a true "everywhere". Then, there is no need for balancing, as all mpiranks automatically
    !! hold the same number of blocks, if they started with a balanced distribution (which is the
    !! output of adapt_mesh.)
    !! This of course implies any other indicator than "everywhere" requires balancing here.
    if ((params%force_maxlevel_dealiasing .eqv. .false.) .or. (indicator/="everywhere")) then
        t1 = MPI_Wtime()
        call balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_active, hvy_n, tree_ID )
        call toc( "refine_mesh (balance_load)", MPI_Wtime()-t1 )
    endif


    call toc( "refine_mesh (lists+neighbors)", t_misc )
    call toc( "refine_mesh (TOTAL)", MPI_wtime()-t0 )

end subroutine refine_mesh
