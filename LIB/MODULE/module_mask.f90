module module_mask

    use module_physics_metamodule
    use module_mesh
    use module_forest
    use module_params
    use module_precision
    use module_globals

    implicit none

contains

    subroutine create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
        hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist, all_parts)

        implicit none

        !> user defined parameter structure
        type (type_params), intent(in)      :: params
        real(kind=rk), intent(in)           :: time
        !> light data array
        integer(kind=ik), intent(inout)     :: lgt_block(:, :)
        real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
        real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
        !> heavy data array - neighbor data
        integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
        !> list of active blocks (heavy data)
        integer(kind=ik), intent(inout)        :: hvy_active(:,:)
        !> list of active blocks (light data)
        integer(kind=ik), intent(inout)        :: lgt_active(:,:)
        !> number of active blocks (heavy data)
        integer(kind=ik), intent(inout)        :: hvy_n(:)
        !> number of active blocks (light data)
        integer(kind=ik), intent(inout)        :: lgt_n(:)
        !> sorted list of numerical treecodes, used for block finding
        integer(kind=tsize), intent(inout)     :: lgt_sortednumlist(:,:,:)
        logical, intent(in), optional :: all_parts

        integer(kind=ik) :: k, lgt_id, Bs(1:3), g, tree_n, hvy_id, iter, Jactive, Jmax
        real(kind=rk) :: x0(1:3), dx(1:3)
        logical, save :: time_independent_part_ready = .false.
        logical :: force_all_parts

        Bs = params%Bs
        g  = params%n_ghosts
        Jactive = max_active_level(lgt_block,lgt_active(:,tree_ID_flow),lgt_n(tree_ID_flow))
        Jmax = params%max_treelevel

        ! without penalization, do nothing.
        if ( .not. params%penalization ) return


        force_all_parts = .false.
        if (present(all_parts)) force_all_parts = all_parts

        ! The advances pruned-tree technology works only if the mask is sufficiently
        ! fine, either on Jmax or Jmax-1. If this condition is not met, generate the
        ! entire mask.
! As a hint, we use the flag (params%threshold_mask) which forces the
        if ((Jactive < Jmax-1) .or. (params%threshold_mask .eqv. .false.) .or. (force_all_parts)) then
            do k = 1, hvy_n(tree_ID_flow)
                hvy_id = hvy_active(k, tree_ID_flow)

                ! convert given hvy_id to lgt_id for block spacing routine
                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

                ! get block spacing for RHS
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                hvy_mask(:,:,:,:,hvy_id), "all-parts" )
            enddo

            if (params%rank==0) then
                write(*,'("Generating mask without pruned trees.. Jactive=",i2," Jmax=",i2," threshold_mask=",L1," force_all=",L1)') &
                Jactive, Jmax, params%threshold_mask, force_all_parts
            endif

            ! we're done, leave routine now
            return
        endif


        ! * the mask coloring and overlapping parts
        ! * sponge and what to do with it



        ! initialization of time-independent part, if this is required.
        ! done only once. mask is generated starting from the coarsest to the
        ! finest level (refined only where interesting).
        ! At most, mask in generated (Jmax-Jmin) times.
        if ( (.not. time_independent_part_ready) .and. (params%mask_time_independent_part) ) then
            if (params%rank==0) then
                write(*,'(80("~"))')
                write(*,*) "creating time-independent part of mask function NOW"
                write(*,'(80("~"))')
            endif

            ! start with an equidistant grid on coarsest level.
            ! routine also deletes any existing mesh in the tree.
            call create_equidistant_grid( params, lgt_block, hvy_neighbor, &
            lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask), &
            lgt_sortednumlist(:,:,tree_ID_mask), hvy_active(:,tree_ID_mask), &
            hvy_n(tree_ID_mask), params%min_treelevel, .true., tree_ID_mask )

            do k = 1, hvy_n(tree_ID_mask)
                hvy_id = hvy_active(k, tree_ID_mask)
                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                hvy_mask(:,:,:,:,hvy_id), "time-independent-part" )
            enddo


            ! the constant part of the mask needs to be generated on Jmax (where RHS is computed)
            do iter = 1, params%max_treelevel - params%min_treelevel
                ! synchronization before refinement (because the interpolation takes place on the extended blocks
                ! including the ghost nodes)
                ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
                ! Snych'ing becomes much mor expensive one the grid is refined.
                ! sync possible only before pruning
                call sync_ghosts( params, lgt_block, hvy_mask, hvy_neighbor, hvy_active(:,tree_ID_mask), hvy_n(tree_ID_mask) )


                ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
                ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
                call refine_mesh( params, lgt_block, hvy_mask, hvy_neighbor, &
                lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask), &
                lgt_sortednumlist(:,:,tree_ID_mask), hvy_active(:,tree_ID_mask), &
                hvy_n(tree_ID_mask), "mask-threshold", tree_ID_mask )


                if (params%rank==0) then
                    write(*,'("Did one iteration for time-independent mask. Now: Jmax=",i2, " Nb=",i7,&
                    &" lgt_n=",(4(i6,1x)))') &
                    max_active_level(lgt_block,lgt_active(:,tree_ID_mask),lgt_n(tree_ID_mask)), lgt_n(tree_ID_mask), &
                    lgt_n
                endif

                ! the constant part needs to be generated on Jmax (where RHS is computed)
                do k = 1, hvy_n(tree_ID_mask)
                    hvy_id = hvy_active(k, tree_ID_mask)
                    call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
                    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                    call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                    hvy_mask(:,:,:,:,hvy_id), "time-independent-part" )
                enddo

            enddo

            ! we need the mask funcition both on Jmax (during the RHS) and Jmax-1
            ! (during saving after coarsening)
            call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_ID_mask_coarser, tree_ID_mask)

            ! coarsen by one level only.
            if (params%rank==0) write(*,*) "Coarsening the mask by one level (to Jmax-1)"
            call adapt_mesh( time, params, lgt_block, hvy_mask, hvy_neighbor, &
            lgt_active(:,tree_ID_mask_coarser), lgt_n(tree_ID_mask_coarser), &
            lgt_sortednumlist(:,:,tree_ID_mask_coarser), hvy_active(:,tree_ID_mask_coarser), &
            hvy_n(tree_ID_mask_coarser), tree_ID_mask_coarser, "everywhere", hvy_tmp, external_loop=.true.)

            ! prune both masks
            if (params%rank==0) write(*,*) "Pruning mask tree (on Jmax)"
            call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_ID_mask)

            if (params%rank==0) write(*,*) "Pruning mask tree (on Jmax-1)"
            call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_ID_mask_coarser)

            time_independent_part_ready = .true.
            if (params%rank==0) then
                write(*,'(80("~"))')
                write(*,*) "DONE creating time-independent part!"
                write(*,'(80("~"))')
            endif
        endif

        ! create "time-dependent-part" here, add the existing "time-independent-part"
        ! if it is available, return the complete mask incl. all parts

        ! read in time-indepent
        ! call read_field2tree(params, (/"chi_00.h5"/), 1, 2, tree_n, &
        ! lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        ! hvy_mask, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
        !

        if ( params%mask_time_dependent_part ) then
            do k = 1, hvy_n(tree_ID_flow)
                hvy_id = hvy_active(k, tree_ID_flow)

                ! convert given hvy_id to lgt_id for block spacing routine
                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

                ! get block spacing for RHS
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                ! note the meta-routine also resets to zero (the entire mask)
                call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                hvy_mask(:,:,:,:,hvy_id), "time-dependent-part" )
            enddo
        else
            ! at some point the mask requires resetting, before we add the time-independent part
            ! to it.
            do k = 1, hvy_n(tree_ID_flow)
                hvy_id = hvy_active(k, tree_ID_flow)
                hvy_mask(:,:,:,:,hvy_id) = 0.0_rk
                ! hvy_mask(:,:,:,6,hvy_id) = 0.0_rk ! sponge. as we keep it time-dependent, it is not required.
            enddo
        endif


        if ( params%mask_time_independent_part ) then
            if (Jactive == params%max_treelevel ) then
                ! flow is on finest level: add complete mask on finest level
                ! this is the case when computing the right hand side.
                call add_pruned_to_full_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_ID_mask, tree_ID_flow)

            elseif (Jactive == params%max_treelevel-1 ) then

                call add_pruned_to_full_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_ID_mask_coarser, tree_ID_flow)
            else
                write(*,*) "WARING: mask generation fails (flow grid neither on Jmax nor Jmax-1)"
            endif
        endif

    end subroutine
end module
