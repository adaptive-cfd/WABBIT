!> \brief This routine initializes the block data, i.e. it evaluates the initial condition on the grid
!! input:    - parameter array
!!           - light data array
!!           - heavy data array
!!           - neighbor data array
!!           - light and heavy active block list
!! output:   - filled user defined data structure for global params
!!           - initialized light and heavy data arrays
! ********************************************************************************************

subroutine set_initial_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_n, hvy_n, lgt_sortednumlist, tree_ID, adapt, time, iteration, hvy_mask, hvy_tmp)

    implicit none

    type (type_params), intent(inout)    :: params                      !> user defined parameter structure
    integer(kind=ik), intent(inout)      :: lgt_block(:, :)             !> light data array
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)   !> heavy work data array - block data.
    real(kind=rk), intent(inout)         :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(inout)      :: hvy_neighbor(:,:)           !> neighbor array (heavy data)
    integer(kind=ik), intent(inout)      :: lgt_active(:,:)             !> list of active blocks light data)
    integer(kind=ik), intent(inout)      :: hvy_active(:,:)             !> list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_n(:), lgt_n(:)          !> number of heavy and light active blocks
    integer(kind=ik), intent(in)         :: tree_ID
    integer(kind=tsize), intent(inout)   :: lgt_sortednumlist(:,:,:)    !> sorted list of numerical treecodes, used for block finding
    real(kind=rk), intent(inout)         :: time                        !> time loop variables
    integer(kind=ik), intent(inout)      :: iteration

    !> if .false. the code initializes on the coarsest grid, if .true. iterations
    !> are performed and the mesh is refined to gurantee the error eps
    logical, intent(in) :: adapt

    logical :: tmp
    integer(kind=ik) :: lgt_n_old, k, iter, lgt_n_tmp

    lgt_n_old = 9999999
    iter = 0
    time = 0.0_rk

    if (params%rank==0) then
        write(*,*) "(((((((((((((((((((inicond)))))))))))))))))))"
        write(*,*) "Setting initial condition on all blocks."
        write(*,*) "Adaptive initial condition is: ", adapt
    endif

    ! this is a HACK
    if (params%physics_type == 'ACM-new') then
        tmp = params%threshold_mask
        params%threshold_mask = .true.
        if (.not. params%penalization) params%threshold_mask = .false.
    endif

    ! choose between reading from files and creating datafields analytically
    if (params%read_from_files) then
        if (params%rank==0) write(*,*) "Initial condition is read from file!"
        !-----------------------------------------------------------------------
        ! read initial condition from file
        !-----------------------------------------------------------------------
        ! Note that reading from file is something all physics modules share - it
        ! is a wabbit routine and not affiliated with a specific physics module
        ! therefore, there is still a grid-level (=wabbit) parameter read_from_files
        call get_inicond_from_file(params, lgt_block, hvy_block, hvy_n, lgt_n, tree_ID, time, iteration)

        ! create lists of active blocks (light and heavy data)
        call updateMetadata_tree(params, lgt_block, hvy_neighbor, lgt_active, &
        lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID)

        ! synching is required for the adaptation step
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

        ! even if we read the initial condition from file, we can still adapt it immediately
        ! so the max error is eps. This is useful if we read a dense field where we already know that
        ! it contains many zeros. Can be used for wavelet compression tests as well.
        iter = 0
        if (adapt) then
            ! now, evaluate the coarsening criterion on each block, and coarsen the grid where possible.
            ! adapt-mesh also performs neighbor and active lists updates
            if (present(hvy_mask) .and. params%threshold_mask) then
                lgt_n_tmp = -99
                ! what happens on very coarse grids is that the first coarsening interation (internal to adapt_tree) removes
                ! the mask completely... This is because we do not currently reconstruct the mask inside this iteration loop (to avoid overhead)
                ! we therefore outsource the iteration loop here. (argument external_loop to adapt_tree)

                do while ( lgt_n_tmp /= lgt_n(tree_ID) )
                    lgt_n_tmp = lgt_n(tree_ID)

                    call create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
                    hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist, all_parts=.true. )

                    call adapt_tree( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
                    lgt_sortednumlist, hvy_active, hvy_n, tree_ID, params%coarsening_indicator, hvy_tmp, hvy_mask=hvy_mask, external_loop=.true. )
                enddo
            else
                call adapt_tree( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
                lgt_sortednumlist, hvy_active, hvy_n, tree_ID, params%coarsening_indicator, hvy_tmp )
            endif

            iter = iter + 1
            if (params%rank == 0) then
                write(*,'(" did ",i2," mesh adaptation for the initial condition. Nblocks=",i6, " Jmin=",i2, " Jmax=",i2)') &
                iter, lgt_n(tree_ID), &
                minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ), &
                maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )
            endif
        endif

        ! HACK
        ! Navier stokes has different statevector (sqrt(rho),sqrt(rho)u,sqrt(rho)v,p) than
        ! the statevector saved to file (rho,u,v,p)
        ! we therefore convert it once here
        if (params%physics_type == 'navier_stokes') then
            call set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n)
        end if

        !-----------------------------------------------------------------------
        ! in some situations, it is necessary to create the intial grid, and then refine it for a couple of times.
        ! for example if one does non-adaptive non-equidistant spatial convergence tests
        if (params%inicond_refinements > 0) then
            do k = 1, params%inicond_refinements
                ! refine entire mesh.
                call refine_tree( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
                lgt_sortednumlist, hvy_active, hvy_n, "everywhere", tree_ID)

                if (params%rank == 0) then
                    write(*,'(" did ",i2," refinement stage (beyond what is required for the &
                    &prescribed precision eps) Nblocks=",i6, " Jmin=",i2, " Jmax=",i2)') k, lgt_n, &
                    minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ),&
                    maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )
                endif
            enddo
        endif

    else
        if (params%inicond_grid_from_file == "no") then
            if (params%rank==0) write(*,*) "Initial condition is defined by physics modules!"
            if (params%rank==0) write(*,*) "Grid is auto-generated!"

            !---------------------------------------------------------------------------
            ! Create the first mesh on the coarsest treelevel
            !---------------------------------------------------------------------------
            call createEquidistantGrid_tree( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, params%min_treelevel, .true., tree_ID )

            !---------------------------------------------------------------------------
            ! on the grid, evaluate the initial condition
            !---------------------------------------------------------------------------
            call set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n)

            !---------------------------------------------------------------------------
            ! grid adaptation
            !---------------------------------------------------------------------------
            if (adapt) then
                if (params%rank==0) write(*,'("INIT: initial grid is adaptive")')

                ! we have to repeat the adapation process until the grid has reached a final
                ! state. Since we start on the coarsest level, in each iteration we cannot loose
                ! blocks, but only gain or no change. Therefore, iterate until lgt_n is constant.
                do while ( lgt_n(tree_ID) /= lgt_n_old  .and. iter<(params%max_treelevel-params%min_treelevel))
                    lgt_n_old = lgt_n(tree_ID)

                    ! push up the entire grid one level.
                    !> \todo It would be better to selectively
                    !! go up one level where a refinement indicator tells us to do so, but in the current code
                    !! versions it is easier to use everywhere. NOTE: you actually should call sync_ghosts before
                    !! but it shouldnt be necessary as the inicond is set also in the ghost nodes layer.
                    call refine_tree( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
                    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, "everywhere", tree_ID  )

                    ! It may seem surprising, but we now have to re-set the inicond on the blocks. if
                    ! not, the detail coefficients for all blocks are zero. In the time stepper, this
                    ! corresponds to advancing the solution in time, it's just that here we know the exact
                    ! solution (the inicond)
                    call set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n)

                    ! now, evaluate the refinement criterion on each block, and coarsen the grid where possible.
                    ! adapt-mesh also performs neighbor and active lists updates
                    ! NOTE: the grid adaptation can take the mask function into account (such that the fluid/solid
                    ! interface is on the finest level).
                    if (present(hvy_mask) .and. params%threshold_mask) then
                        lgt_n_tmp = -99
                        ! what happens on very coarse grids is that the first coarsening interation removes
                        ! the mask completely...
                        ! we therefore outsource the iteration loop here. (argument external_loop to
                        ! adapt_tree). This loop iterates until the grid does no longer change, then lgt_n_tmp = lgt_n(tree_ID)
                        do while ( lgt_n_tmp /= lgt_n(tree_ID) )
                            lgt_n_tmp = lgt_n(tree_ID)

                            call create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
                            hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist, all_parts=.true. )

                            call adapt_tree( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
                            lgt_sortednumlist, hvy_active, hvy_n, tree_ID, params%coarsening_indicator_inicond, hvy_tmp, hvy_mask=hvy_mask, external_loop=.true. )
                        enddo
                    else
                        call adapt_tree( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
                        lgt_sortednumlist, hvy_active, hvy_n, tree_ID, params%coarsening_indicator_inicond, hvy_tmp )
                    endif

                    iter = iter + 1
                    if (params%rank == 0) then
                        write(*,'("INIT: did ",i2," mesh adaptation for the initial condition. Nblocks=",i7, " Jmin=",i2, " Jmax=",i2)') iter, lgt_n(tree_ID), &
                        minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ), &
                        maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )
                    endif
                enddo
            endif

            !-----------------------------------------------------------------------
            ! in some situations, it is necessary to create the intial grid, and then refine it for a couple of times.
            ! for example if one does non-adaptive non-equidistant spatial convergence tests
            if (params%inicond_refinements > 0) then
                do k = 1, params%inicond_refinements
                    ! refine entire mesh.
                    call refine_tree( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
                    lgt_sortednumlist, hvy_active, hvy_n, "everywhere", tree_ID)

                    ! set initial condition
                    call set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n)

                    if (params%rank == 0) then
                        write(*,'(" did ",i2," refinement stage (beyond what is required for the &
                        &prescribed precision eps) Nblocks=",i6, " Jmin=",i2, " Jmax=",i2)') k, lgt_n(tree_ID), &
                        minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ),&
                        maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )
                    endif
                enddo
            endif
        else
            if (params%rank==0) write(*,*) "Initial condition is defined by physics modules!"
            if (params%rank==0) write(*,*) "Grid is read from file!"

            ! read grid from file
            call read_mesh(params%inicond_grid_from_file, params, lgt_n, hvy_n, lgt_block, tree_ID)

            ! create lists of active blocks (light and heavy data)
            call updateMetadata_tree(params, lgt_block, hvy_neighbor, lgt_active, &
            lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID)

            !---------------------------------------------------------------------------
            ! on the grid, evaluate the initial condition
            !---------------------------------------------------------------------------
            call set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n)
        endif
    end if

    !---------------------------------------------------------------------------
    ! Finalization. Note this routine has modified the grid, added some blocks, removed some blocks
    ! so the neighboring information and active lists are outdated. It is good practice that a routine
    ! that modifies the grid shall return a WORKING grid, with at least active / neighboring lists
    ! up to date. here, we also perform an initial load balancing step which is not overly important,
    ! since in the time loop it is performed as well, but its still good to have a neatly balanced initial
    ! condition. This is especially true if the initial condition is read from file.
    !
    ! NOTE: in fact, adapt_tree returns a working mesh (good practice!), so the update here may be
    ! redundant, but as we initialize only once we can live with the extra cost (and increased security)
    !---------------------------------------------------------------------------

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call updateMetadata_tree(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, tree_ID)

    ! balance the load
    call balanceLoad_tree(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, tree_ID)

    ! synchronize ghosts now, in order to start with a clean grid. NOTE this can actually be removed, but
    ! it is a safety issue. Better simply keep it.
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    ! footer...and done!
    if (params%rank == 0) then
        write(*,'("Resulting grid for initial condition: Nblocks=",i6, " Jmin=",i2, " Jmax=",i2)') lgt_n(tree_ID), &
        minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ), &
        maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )
        write(*,'("Initial grid and initial condition terminated.")')
        write(*,*) "(((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))"
    endif

    ! HACK
    if (params%physics_type == 'ACM-new') then
        params%threshold_mask = tmp
    endif
end subroutine set_initial_grid
