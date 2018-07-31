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

subroutine refine_mesh( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, indicator  )

!---------------------------------------------------------------------------------------------
! modules
    use module_indicators

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
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

    ! cpu time variables for running time calculation
    real(kind=rk)                          :: t0, t1, t2, t_misc
    integer(kind=ik)                       :: k

    ! start time
    t0 = MPI_Wtime()
    t_misc = 0.0_rk

    !> (a) loop over the blocks and set their refinement status.
    t1 = MPI_Wtime()
    call refinement_indicator( params, lgt_block, lgt_active, lgt_n, indicator )
    call toc( params, "refine_mesh (refinement_indicator)", MPI_Wtime()-t1 )


    !> (b) check if block has reached maximal level, if so, remove refinement flags
    t1 = MPI_Wtime()
    call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
    call toc( params, "refine_mesh (respect_min_max_treelevel)", MPI_Wtime()-t1 )


    !> (c) ensure gradedness of mesh. If the refinement is done everywhere, there is
    !! no way gradedness can be damaged, so we skip the call in this case. However,
    !! in all other indicators, this step is very important.
    t1 = MPI_Wtime()
    if ( indicator /= "everywhere" ) then
      call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
      lgt_sortednumlist, hvy_active, hvy_n )
    endif
    call toc( params, "refine_mesh (ensure_gradedness)", MPI_Wtime()-t1 )


    !> (d) execute refinement, interpolate the new mesh. All blocks go one level up
    !! except if they are already on the highest level. Note that those blocks have
    !! the status +11
    t1 = MPI_Wtime()
    if ( params%threeD_case ) then
        call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
    else
        call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_active, hvy_n )
    end if
    call toc( params, "refine_mesh (refinement_execute)", MPI_Wtime()-t1 )


    !> (e) as the grid changed now with the refinement, we have to update the list of
    !! active blocks so other routines can loop just over these active blocks
    !! and do not have to ensure that the active list is up-to-date
    ! update list of sorted nunmerical treecodes, used for finding blocks
    t2 = MPI_wtime()
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
         hvy_active, hvy_n, lgt_sortednumlist, .true. )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
         lgt_sortednumlist, hvy_active, hvy_n )
    t_misc = MPI_wtime() - t2


    !> (f) At this point the refinement is done. Since not all blocks are refined, namely only those
    !! that were not on Jmax, Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now. EXCEPTION: if the flag force_maxlevel_dealiasing is set
    !! then we force blocks on Jmax to coarsen, even if their details are large. Hence, each refinement
    !! step is a true "everywhere". Then, there is no need for balancing, as all mpiranks automatically
    !! hold the same number of blocks, if they started with a balanced distribution (which is the
    !! output of adapt_mesh.)
    if (params%force_maxlevel_dealiasing .eqv. .false.) then
        t1 = MPI_Wtime()
        call balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
        hvy_active, hvy_n, hvy_work )
        call toc( params, "refine_mesh (balance_load)", MPI_Wtime()-t1 )

        t2 = MPI_wtime()
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
             hvy_active, hvy_n, lgt_sortednumlist, .true. )

        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
             lgt_sortednumlist, hvy_active, hvy_n )
        t_misc = t_misc + MPI_wtime() - t2
    endif


    call toc( params, "refine_mesh (lists+neighbors)", t_misc )
    call toc( params, "refine_mesh (TOTAL)", MPI_wtime()-t0 )

end subroutine refine_mesh
