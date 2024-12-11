subroutine createMask_tree(params, time, hvy_mask, hvy_tmp, all_parts)
    ! to test update insect here (HACK)
    use module_ACM
    implicit none

    type (type_params), intent(inout)      :: params
    real(kind=rk), intent(in)              :: time
    real(kind=rk), intent(inout)           :: hvy_mask(:, :, :, :, :)
    real(kind=rk), intent(inout)           :: hvy_tmp(:, :, :, :, :)
    logical, intent(in), optional          :: all_parts
    integer(kind=ik)                       :: k, lgt_id, Bs(1:3), g, hvy_id, iter, Jactive, Jmax
    real(kind=rk)                          :: x0(1:3), dx(1:3)
    logical, save                          :: time_independent_part_ready = .false.
    logical                                :: force_all_parts

    Bs      = params%Bs
    g       = params%g
    Jactive = maxActiveLevel_tree(tree_ID_flow)
    Jmax    = params%Jmax
    tree_n  = params%forest_size ! used only for resetting at this point

    ! without penalization, do nothing.
    if (.not. params%penalization) return

    ! HACK
    if (params%physics_type /= "ACM-new") return

    ! HACK:
    ! some mask functions have initialization routines (insects) which are to be called once and not for each
    ! block (efficiency). Usually, this would be a staging concept as well, but as only Thomas uses it anyways, cleanup
    ! is left as FIXME
    if (params_acm%geometry=="Insect") then
        call Update_Insect_wrapper(time)
    endif

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
    !     In this case, we force the code to generate all parts of the mask. NB: as of 08/2022, even in the
    !     initial condition, pruned trees are used.
    !-----------------------------------------------------------------------
    if ((Jactive < Jmax-1) .or. (params%threshold_mask .eqv. .false.) .or. (force_all_parts)) then
        ! generate complete mask (may be expensive)
        call createCompleteMaskDirect_tree(params, time, hvy_mask)

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
    if ((.not. time_independent_part_ready) .and. (params%mask_time_independent_part)) then
        call createTimeIndependentMask_tree(params, time, hvy_mask, hvy_tmp)

        time_independent_part_ready = .true.
    endif

    ! create "time-dependent-part" here, add the existing "time-independent-part"
    ! if it is available, return the complete mask incl. all parts
    if ( params%mask_time_dependent_part ) then
        do k = 1, hvy_n(tree_ID_flow)
            hvy_id = hvy_active(k, tree_ID_flow)

            ! convert given hvy_id to lgt_id for block spacing routine
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! get block spacing for RHS
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            ! 11 Oct 2020: this seems to be required, as the mask array can rarely
            ! contain some garbage which is not deleted by the create_mask subroutines
            hvy_mask(:,:,:,:,hvy_id) = 0.0_rk

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
        if (Jactive == params%Jmax ) then
            ! flow is on finest level: add complete mask on finest level
            ! this is the case when computing the right hand side.
            call add_pruned_to_full_tree( params, hvy_mask, tree_ID_mask, tree_ID_flow)
            ! For adaption with the full tree, we also have the flow map present on JMax-1
            ! This is, as force-dealiasing overrides the mask importance on JMax, but on JMax-1 we want to respect it, so we have the mask on both levels
            ! Here we add both blocks therefore and in adaption process being in full tree formulation in will respect both at the same time
            ! For leaf-grid computations during the time-step we simply ignore the blocks on JMax-1
            if (params%force_maxlevel_dealiasing) then
                call add_pruned_to_full_tree( params, hvy_mask, tree_ID_mask_coarser, tree_ID_flow)
            endif

        elseif (Jactive == params%Jmax-1 .and. params%force_maxlevel_dealiasing) then

            call add_pruned_to_full_tree( params, hvy_mask, tree_ID_mask_coarser, tree_ID_flow)
        else
            write(*,*) "WARING: mask generation fails (flow grid neither on Jmax nor Jmax-1)"
        endif
    endif

end subroutine




subroutine createCompleteMaskDirect_tree(params, time, hvy_mask)
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(in)           :: time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    integer :: k, hvy_id, Bs(1:3), g, lgt_id
    real(kind=rk) :: x0(1:3), dx(1:3)

    Bs = params%Bs
    g  = params%g

    do k = 1, hvy_n(tree_ID_flow)
        hvy_id = hvy_active(k, tree_ID_flow)

        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
        hvy_mask(:,:,:,:,hvy_id), "all-parts" )
    enddo

end subroutine




subroutine createTimeIndependentMask_tree(params, time, hvy_mask, hvy_tmp)
    implicit none

    type (type_params), intent(inout)   :: params
    real(kind=rk), intent(in)           :: time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer :: k, hvy_id, Bs(1:3), g, Jactive, Jmax, iter, lgt_id, Jmin
    real(kind=rk)                       :: x0(1:3), dx(1:3)

    logical :: error_OOM
    if (params%rank==0) then
        write(*,'(80("~"))')
        write(*,*) "creating time-independent part of mask function NOW"
        write(*,'(80("~"))')
    endif

    Bs = params%Bs
    g  = params%g
    Jactive = maxActiveLevel_tree(tree_ID_flow)
    Jmax = params%Jmax
    Jmin = params%Jmin
    tree_n = params%forest_size ! used only for resetting at this point

    ! start with an equidistant grid on coarsest level.
    ! routine also deletes any existing mesh in the tree.
    call createEquidistantGrid_tree( params, hvy_mask, params%Jmin, .true., tree_ID_mask )


    do k = 1, hvy_n(tree_ID_mask)
        hvy_id = hvy_active(k, tree_ID_mask)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        call get_block_spacing_origin( params, lgt_id, x0, dx )
        call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
        hvy_mask(:,:,:,:,hvy_id), "time-independent-part" )
    enddo


    ! the constant part of the mask needs to be generated on Jmax (where RHS is computed)
    do iter = 1, Jmax - Jmin
        ! synchronization before refinement (because the interpolation takes place on the extended blocks
        ! including the ghost nodes)
        ! Note: at this point he grid is rather coarse (fewer blocks), and the sync step is rather cheap.
        ! Snyc'ing becomes much more expensive once the grid is refined.
        ! sync possible only before pruning
        call sync_ghosts_tree( params, hvy_mask, tree_ID_mask )


        ! refine the mesh
        call refine_tree( params, hvy_mask, hvy_tmp, "mask-threshold", tree_ID_mask, error_OOM)

        if (error_OOM) call abort(2512112,"Refinement failed, out of memory. Try with more memory.")

        ! if its mask-anynonzero, then the grid is refined inside the body. however, this is not what we assume
        ! in the add-pruned-tree: blocks on the actual grid are then coarser than in the pruned one, and still all 0
        ! call refine_tree( params, hvy_mask, "mask-anynonzero", tree_ID_mask )


        if (params%rank==0) then
            write(*,'("Did refinement for time-independent mask. Now: Jmax=",i2, " Nb=",i7," lgt_n=",(4(i6,1x)))') &
            maxActiveLevel_tree(tree_ID_mask), lgt_n(tree_ID_mask), lgt_n
        endif

        ! the constant part needs to be generated on Jmax (where RHS is computed)
        do k = 1, hvy_n(tree_ID_mask)
            hvy_id = hvy_active(k, tree_ID_mask)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            ! expensive:
            call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, hvy_mask(:,:,:,:,hvy_id), "time-independent-part" )
        enddo

        ! Deactivated following code (TE, 28jun2024), it seems unnecessary? We refine blocks that contain traces of the mask 
        ! function, then set the mask function again. How should we create blocks that can be coarsened again? -> skip adapt_tree part here

        ! ! we found that sometimes, we end up with more blocks than expected
        ! ! and some zero-blocks can be coarsened again. doing that turned out to be
        ! ! important for large-scale simulations
        ! call adapt_tree( time, params, hvy_mask, tree_ID_mask, "threshold-state-vector", hvy_tmp, ignore_maxlevel=.true.)

        ! ! TODO (TE, 28jun2024): It is not 100% clear to me why the constant mask looks like every 2nd point is deleted at this
        ! ! point. 

        ! ! re-create mask on the adapted grid
        ! do k = 1, hvy_n(tree_ID_mask)
        !     hvy_id = hvy_active(k, tree_ID_mask)

        !     ! As the adapt_tree routine coarsened the grid, it modified the data and we need to set it again.
        !     ! Note, however, that we can at least exploit the fact that we have already created it once, and use the
        !     ! modified mask data as an indicator where to re-create the mask now.
        !     if (maxval(hvy_mask(:,:,:,1,hvy_id)) > 1.0e-9_rk ) then
        !         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        !         call get_block_spacing_origin( params, lgt_id, x0, dx )

        !         ! expensive:
        !         call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, hvy_mask(:,:,:,:,hvy_id), "time-independent-part" )
        !     endif
        ! enddo

        ! if (params%rank==0) then
        !     write(*,'("Did coarsening for time-independent mask. Now: Jmax=",i2, " Nb=",i7," lgt_n=",(4(i6,1x)))') &
        !     maxActiveLevel_tree(tree_ID_mask), lgt_n(tree_ID_mask), lgt_n
        ! endif
    enddo

    ! call saveHDF5_tree('constmask_1.h5', time, 0_ik, 1, params, hvy_mask, tree_ID_mask)

    ! syncing now and pruning later keeps the ghost nodes of the time-independent mask function sync'ed (as they
    ! do not change)
    call sync_ghosts_tree( params, hvy_mask, tree_ID_mask )

    ! we need the mask function both on Jmax (during the RHS) and Jmax-1
    ! (during saving after coarsening)
    call copy_tree(params, hvy_mask, tree_ID_mask_coarser, tree_ID_mask)

    if (params%force_maxlevel_dealiasing) then
        ! coarsen by one level only - only needed when we force dealiasing
        if (params%rank==0) write(*,*) "Coarsening the mask by one level (to Jmax-1)"
        call adapt_tree( time, params, hvy_mask, tree_ID_mask_coarser, "everywhere", hvy_tmp)

        ! pruning (Jmax-1)
        if (params%rank==0) write(*,'("Pruning mask tree (on Jmax-1 = ",i3,")")') Jmax-1
        call prune_tree( params, hvy_mask, tree_ID_mask_coarser)
    endif

    ! pruning (Jmax)
    if (params%rank==0) write(*,'("Pruning mask tree (on Jmax = ",i3,")")') Jmax
    call prune_tree( params, hvy_mask, tree_ID_mask)

    if (params%rank==0) then
        write(*,'(80("~"))')
        write(*,*) "DONE creating time-independent part!"
        write(*,'(80("~"))')
    endif

end subroutine
