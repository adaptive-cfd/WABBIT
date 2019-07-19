!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name refine_mesh.f90
!> \version 0.5
!> \author msr, engels
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
!
!>
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
!! \n
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 03/02/17 - 3D heavy data \n
!! 04/04/17 - include active lists, neighbor relations and load balancing, symmetrical to adapt_mesh \n
!! 05/04/17 - Provide an interface to use different criteria for refinement, rename routines
! ********************************************************************************************

subroutine refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, indicator, tree_ID  )

    use module_indicators


    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)         :: params
    !> light data array
    integer(kind=ik), intent(inout)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)           :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)        :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)     :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)           :: indicator
    integer(kind=ik), intent(in)           :: tree_ID

    ! cpu time variables for running time calculation
    real(kind=rk)                          :: t0, t1, t2, t_misc
    integer(kind=ik)                       :: k
    ! start time
    t0 = MPI_Wtime()
    t_misc = 0.0_rk

    ! if the refinement is "everwhere" and we use filtering on the max level
    ! (force_maxlvel_dealiasing=1), then we end up with exactly 8*hvy_n blocks (on this rank)
    ! and we can say right now if the memory is sufficient or not. Advantage: you
    ! get a cleaner error message than the one issued by get_free_local_light_id
    if ( indicator == "everywhere" .and. params%force_maxlevel_dealiasing ) then
        if ( (2**params%dim)*hvy_n > params%number_blocks ) then
            write(*,'("On rank:",i5," hvy_n=",i6," which will give ",i1,"*hvy_n=",i7," but limit is ",i6)') &
            params%rank, hvy_n, 2**params%dim, (2**params%dim)*hvy_n, params%number_blocks

            call abort (1909181827,"[refine_mesh.f90]: The refinement step will fail because we do not&
            & have enough memory. Try increasing --memory.")
        endif
    endif

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
