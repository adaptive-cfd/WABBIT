subroutine createMask_tree(params, time, hvy_mask, hvy_tmp, all_parts)
    use module_timing
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
    real(kind=rk)                          :: t_block, t_cycle

    t_cycle = MPI_wtime()
    Bs      = params%Bs
    g       = params%g
    x0       = 0.0_rk
    dx       = 0.0_rk
    Jactive = maxActiveLevel_tree(tree_ID_flow)
    Jmax    = params%Jmax
    tree_n  = params%forest_size ! used only for resetting at this point

    ! without penalization, do nothing.
    if (.not. params%penalization) return
    ! HACK
    if (params%physics_type /= "ACM-new" .and. params%physics_type /= "NSPP") return

    ! Initialize the mask-related geometries and update the insect kinematics.
    ! This stage is executed once per call, before processing individual blocks.
    call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, hvy_mask(:,:,:,:,1), "init_stage" )

    if (params%forest_size < 3) call abort(190719,"Forest size is too small (increase to at least 3 in parameter file)")

    ! By default, recompute all mask components.
    ! Selective updates must be explicitly requested (currently implemented only for insects 
    ! using RK4) and require a time-independent mask component and a 3D simulation.
    force_all_parts                                              = .true.
    if (present(all_parts)) force_all_parts                      = all_parts
    if (.not. params%mask_time_independent_part) force_all_parts = .true.
    if (params%dim == 2) force_all_parts                         = .true.


    !-----------------------------------------------------------------------
    ! complete generation of all parts, if required
    !-----------------------------------------------------------------------
    ! The advanced pruned-tree technology works only if the mask is sufficiently
    ! fine, either on Jmax or Jmax-1. If this condition is not met, generate the
    ! entire mask.
    ! A closer look on the conditions:
    ! (1) Jactive < Jmax or Jmax-1 (for force_maxlevel_dealiasing)
    !     If the highest level Jactive on the current grid is not at least Jmax(-1), we cannot use the
    !     pruned trees here, because we have the time-independent part prepared only on Jmax (and Jmax-1).
    !     NOTE: unfortunately the inverse is not true: if Jactive==Jmax, it does in fact NOT necessarily mean that we can
    !     use pruned trees.
    ! (2) params%threshold_mask .eqv. .false.
    !     If the mask is not considered for thresholding, there is no guarantee that the fluid-solid interface
    !     is either on Jmax (or Jmax-1). While it may still be the case, it will not always be like this. Therefore
    !     if mask thresholding is not used, always generate the entire mask (and be done with the routine)
    !     NOTE: expensive mask functions should always use mask_thresholding.
    ! (3) force_all
    !     Sometimes we know that the pruned-trees will not work, specifically during the initial condition.
    !     In this case, we force the code to generate all parts of the mask. NB: as of 08/2022, even in the
    !     initial condition, pruned trees are used.
    ! (4) time_independent_part_ready
    !     Once the mask has been fully generated on Jmax, the time-independent mask components are considered 
    !     initialized and can be preserved. Subsequent calls can then use selective mask updates instead of
    !     recomputing the complete mask.
    !-----------------------------------------------------------------------
    if (((Jactive < Jmax-1 .and. params%force_maxlevel_dealiasing) .or. Jactive < Jmax) & 
        .or. (params%threshold_mask .eqv. .false.) .or. (force_all_parts) .or. (.not. time_independent_part_ready)) then
        ! generate complete mask (may be expensive)
        t_block = MPI_Wtime()
        call createCompleteMaskDirect_tree(params, time, hvy_mask)

        call toc( "createMask_tree (createCompleteMaskDirect_tree)", 66, MPI_Wtime()-t_block )


        time_independent_part_ready = (Jactive == Jmax)

        ! we're done, all parts of mask function are created, leave routine now
        return
    endif

    !-----------------------------------------------------------------------
    ! Selective mask update
    !-----------------------------------------------------------------------
    ! Preserve the existing time-independent mask components and recompute
    ! only the time-dependent components.
    ! The required selective resets are handled by the physics routines,
    ! rather than resetting the complete mask here.

    if ( params%mask_time_dependent_part ) then
        t_block = MPI_Wtime()

        do k = 1, hvy_n(tree_ID_flow)
            hvy_id = hvy_active(k, tree_ID_flow)

            ! convert given hvy_id to lgt_id for block spacing routine
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            ! get block spacing for RHS
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                                    hvy_mask(:,:,:,:,hvy_id), "time-dependent-part" )

        enddo
  
        call toc( "createMask_tree (time-dependent-part)", 68, MPI_Wtime()-t_block )
    endif

    call toc( "createMask_tree (TOTAL)", 65, MPI_Wtime()-t_cycle )

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

