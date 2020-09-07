module module_mask

    use module_physics_metamodule
    use module_mesh
    use module_forest
    use module_params
    use module_precision
    use module_globals
    use module_MPI

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
        Jactive = max_active_level(lgt_block, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow))
        Jmax = params%max_treelevel
        tree_n = params%forest_size ! used only for resetting at this point

        ! without penalization, do nothing.
        if (.not. params%penalization) return

        ! HACK
        if (params%physics_type /= "ACM-new") return


        if (params%forest_size < 3) call abort(190719,"Forest size is too small (increase to at least 3 in parameter file)")

        ! default is false
        force_all_parts = .false.
        if (present(all_parts)) force_all_parts = all_parts
        if (params%dim == 2) force_all_parts = .true.
        if (params%dont_use_pruned_tree_mask) force_all_parts = .true.

        !-----------------------------------------------------------------------
        ! complete generation of all parts, if required
        !-----------------------------------------------------------------------
        ! The advanced pruned-tree technology works only if the mask is sufficiently
        ! fine, either on Jmax or Jmax-1. If this condition is not met, generate the
        ! entire mask.
        ! A closer look on the conditions:
        ! (1) Jactive < Jmax-1
        !     If the highest level Jactive on the current grid is not at least Jmax-1, we cannot use the
        !     pruned trees here, because we have the time-independent part prepared only on Jmax and Jmax-1.
        !     NOTE: unfortunately the inverse is not true: if Jactive==Jmax, it does in fact NOT necessarily mean that we can
        !     use pruned trees.
        ! (2) params%threshold_mask .eqv. .false.
        !     If the mask is not considered for thresholding, there is no guarantee that the fluid-solid interface
        !     is either on Jmax or Jmax-1. While it may still be the case, it will not always be like this. Therefore
        !     if mask thresholding is not used, always generate the entire mask (and be done with the routine)
        !     NOTE: expensive mask functions should always use mask_thresholding.
        ! (3) force_all
        !     Sometimes we know that the pruned-trees will not work, specifically during the initial condition.
        !     In this case, we force the code to generate all parts of the mask.
        !-----------------------------------------------------------------------

        if ((Jactive < Jmax-1) .or. (params%threshold_mask .eqv. .false.) .or. (force_all_parts)) then
            ! generate complete mask (may be expensive)
            call create_complete_mask(params, time, lgt_block, hvy_mask, hvy_active, hvy_n)

            ! we're done, all parts of mask function are created, leave routine now
            return
        endif


        !-----------------------------------------------------------------------
        ! advanced mask generation, if possible
        !-----------------------------------------------------------------------
        ! initialization of time-independent part, if this is required.
        ! done only once. mask is generated starting from the coarsest to the
        ! finest level (refined only where interesting).
        ! At most, mask in generated (Jmax-Jmin) times.
        if ( (.not. time_independent_part_ready) .and. (params%mask_time_independent_part) ) then
            call create_time_independent_mask(params, time, lgt_block, hvy_mask, hvy_tmp, &
            hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist)

            time_independent_part_ready = .true.
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




    subroutine create_complete_mask(params, time, lgt_block, hvy_mask, hvy_active, hvy_n)

        implicit none

        !> user defined parameter structure
        type (type_params), intent(in)      :: params
        real(kind=rk), intent(in)           :: time
        !> light data array
        integer(kind=ik), intent(inout)     :: lgt_block(:, :)
        real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
        !> list of active blocks (heavy data)
        integer(kind=ik), intent(inout)     :: hvy_active(:,:)
        !> number of active blocks (heavy data)
        integer(kind=ik), intent(inout)     :: hvy_n(:)
        integer :: k, hvy_id, Bs(1:3), g, lgt_id
        real(kind=rk) :: x0(1:3), dx(1:3)

        Bs = params%Bs
        g  = params%n_ghosts

        do k = 1, hvy_n(tree_ID_flow)
            hvy_id = hvy_active(k, tree_ID_flow)

            ! convert given hvy_id to lgt_id for block spacing routine
            call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! get block spacing for RHS
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
            hvy_mask(:,:,:,:,hvy_id), "all-parts" )
        enddo

    end subroutine




    subroutine create_time_independent_mask(params, time, lgt_block, hvy_mask, hvy_tmp, &
        hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist)

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
        integer(kind=ik), intent(inout)     :: hvy_active(:,:)
        !> list of active blocks (light data)
        integer(kind=ik), intent(inout)     :: lgt_active(:,:)
        !> number of active blocks (heavy data)
        integer(kind=ik), intent(inout)     :: hvy_n(:)
        !> number of active blocks (light data)
        integer(kind=ik), intent(inout)     :: lgt_n(:)
        !> sorted list of numerical treecodes, used for block finding
        integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:,:)
        integer :: k, hvy_id, Bs(1:3), g, Jactive, Jmax, tree_n, iter, lgt_id, Jmin
        real(kind=rk) :: x0(1:3), dx(1:3)

        if (params%rank==0) then
            write(*,'(80("~"))')
            write(*,*) "creating time-independent part of mask function NOW"
            write(*,'(80("~"))')
        endif

        Bs = params%Bs
        g  = params%n_ghosts
        Jactive = max_active_level(lgt_block,lgt_active(:,tree_ID_flow),lgt_n(tree_ID_flow))
        Jmax = params%max_treelevel
        Jmin = params%min_treelevel
        tree_n = params%forest_size ! used only for resetting at this point

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
        do iter = 1, Jmax - Jmin
            ! synchronization before refinement (because the interpolation takes place on the extended blocks
            ! including the ghost nodes)
            ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
            ! Snyc'ing becomes much more expensive once the grid is refined.
            ! sync possible only before pruning
            call sync_ghosts( params, lgt_block, hvy_mask, hvy_neighbor, hvy_active(:,tree_ID_mask), hvy_n(tree_ID_mask) )


            ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
            ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
            call refine_mesh( params, lgt_block, hvy_mask, hvy_neighbor, &
            lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask), &
            lgt_sortednumlist(:,:,tree_ID_mask), hvy_active(:,tree_ID_mask), &
            hvy_n(tree_ID_mask), "mask-threshold", tree_ID_mask )


            if (params%rank==0) then
                write(*,'("Did refinement for time-independent mask. Now: Jmax=",i2, " Nb=",i7," lgt_n=",(4(i6,1x)))') &
                max_active_level(lgt_block,lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask)), lgt_n(tree_ID_mask), lgt_n
            endif

            ! the constant part needs to be generated on Jmax (where RHS is computed)
            do k = 1, hvy_n(tree_ID_mask)
                hvy_id = hvy_active(k, tree_ID_mask)

                ! NOTE: if I am not mistaken, we could at this point also escape zero-valued blocks (Thomas, Yokohama, 23 Oct 2019)
                ! ==> you are mistaken. some blocks contain garbage and will not be removed
                ! probably you could set those blocks to zero but until we use the µCT really, we should not bother.
                ! if (maxval(hvy_mask(:,:,:,1,hvy_id)) > 1.0e-9_rk .and. iter>Jmax-Jmin-2   ) then
                    call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
                    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                    call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                    hvy_mask(:,:,:,:,hvy_id), "time-independent-part" )
                ! endif
            enddo

            ! we found that sometimes, we end up with more blocks than expected
            ! and some zero-blocks can be coarsened again. doing that turned out to be
            ! important for large-scale simulations
            call adapt_mesh( time, params, lgt_block, hvy_mask, hvy_neighbor, &
            lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask), &
            lgt_sortednumlist(:,:,tree_ID_mask), hvy_active(:,tree_ID_mask), &
            hvy_n(tree_ID_mask), tree_ID_mask, "mask-allzero-noghosts", hvy_tmp, external_loop=.false., ignore_maxlevel=.true.)


            if (params%rank==0) then
                write(*,'("Did coarsening for time-independent mask. Now: Jmax=",i2, " Nb=",i7," lgt_n=",(4(i6,1x)))') &
                max_active_level(lgt_block,lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask)), lgt_n(tree_ID_mask), lgt_n
            endif
        enddo

        ! required strictly speaking only if we intent to save TREE_ID_MASK separately to disk.
        ! is called only once so performance does not matter here.
        call sync_ghosts( params, lgt_block, hvy_mask, hvy_neighbor, hvy_active(:,tree_ID_mask), hvy_n(tree_ID_mask) )

        ! we need the mask function both on Jmax (during the RHS) and Jmax-1
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
        if (params%rank==0) write(*,'("Pruning mask tree (on Jmax = ",i3,")")') Jmax
        call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_ID_mask)

        if (params%rank==0) write(*,'("Pruning mask tree (on Jmax-1 = ",i3,")")') Jmax-1
        call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_ID_mask_coarser)

        if (params%rank==0) then
            write(*,'(80("~"))')
            write(*,*) "DONE creating time-independent part!"
            write(*,'(80("~"))')
        endif

    end subroutine
end module
