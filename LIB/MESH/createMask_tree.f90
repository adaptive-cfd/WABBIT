subroutine createMask_tree(params, time, hvy_mask, hvy_tmp)
    use module_timing
    implicit none

    type (type_params), intent(inout)      :: params
    real(kind=rk), intent(in)              :: time
    real(kind=rk), intent(inout)           :: hvy_mask(:, :, :, :, :)
    real(kind=rk), intent(inout)           :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                       :: k, lgt_id, Bs(1:3), g, hvy_id, iter, Jactive, Jmax
    real(kind=rk)                          :: x0(1:3), dx(1:3)
    logical                                :: force_all_parts
    real(kind=rk)                          :: t_block, t_cycle

    t_cycle = MPI_wtime()
    Bs      = params%Bs
    g       = params%g
    x0      = 0.0_rk
    dx      = 0.0_rk
    Jactive = maxActiveLevel_tree(tree_ID_flow)
    Jmax    = params%Jmax

    ! without penalization, do nothing.
    if (.not. params%penalization) return
    ! HACK
    if (params%physics_type /= "ACM-new" .and. params%physics_type /= "NSPP") return

    ! Initialize the mask-related geometries and update the insect kinematics.
    ! This stage is executed once per call, before processing individual blocks.
    call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, hvy_mask(:,:,:,:,1), "init_stage" )

    ! By default, recompute all mask components.
    force_all_parts                                              = .true.
    ! Selective updates must be explicitly requested (currently implemented only for insects 
    ! using RK4) and require a time-independent mask component and a 3D simulation.
    if (.not. params%selective_mask_update) force_all_parts      = .true.
    ! This flag is set in the calling routines and tells us if the selective mask update could now be used
    if (.not. params%selectiveMaskUpdateReady)  force_all_parts  = .true.
    ! If no time independent part is present (according to the user and the PARAMS file. no selective mask update)
    ! (actually done already in ini_file_to_params, redundant = kept for safety)
    if (.not. params%mask_time_independent_part) force_all_parts = .true.
    ! maxlevel dealiasing forces the blocks on Jmax to coarsen and then the selective mask update cannot work
    if (params%force_maxlevel_dealiasing) force_all_parts        = .true.
    ! if we try to generate the mask and the grid is not yet on Jmax, then we are sure that the selective mask
    ! update cannot work (this might happen during initial condition or dry run, but it should not happen, because selectiveMaskUpdateReady=.false.)
    if (Jactive < Jmax) force_all_parts                          = .true.
    ! Should we not use threshold_mask, the grid with the mask can change, and no selective update is possible
    ! (actually done already in ini_file_to_params, redundant = kept for safety)
    if (.not. params%threshold_mask)  force_all_parts            = .true.
    ! selective update not implemented for 2D runs (actually done already in ini_file_to_params, redundant = kept for safety)
    if (params%dim == 2) force_all_parts                         = .true.

    ! if (params%rank == 0) write(*,*) "selective mask", time, .not. force_all_parts, ":", params%selective_mask_update, params%selectiveMaskUpdateReady, params%mask_time_independent_part

    !-----------------------------------------------------------------------
    ! complete generation of all parts, if required
    !-----------------------------------------------------------------------
    ! The selective mask update does not delete time_independent parts of the 
    ! mask (the insect's body), it deletes other parts and leaves the body as is.
    ! This can only work if the conditions above are satisfied. 
    ! Otherwise, we always create the full mask.
    !-----------------------------------------------------------------------
    if (force_all_parts) then
        ! generate complete mask (may be expensive)
        t_block = MPI_Wtime()

        do k = 1, hvy_n(tree_ID_flow)
            hvy_id = hvy_active(k, tree_ID_flow)

            ! convert given hvy_id to lgt_id for block spacing routine
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            ! get block spacing for RHS
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, hvy_mask(:,:,:,:,hvy_id), "all-parts" )

        enddo
        call toc( "createMask_tree (createCompleteMaskDirect_tree)", 66, MPI_Wtime()-t_block )

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

            call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, hvy_mask(:,:,:,:,hvy_id), "time-dependent-part" )

        enddo
  
        call toc( "createMask_tree (time-dependent-part)", 68, MPI_Wtime()-t_block )
    endif

    call toc( "createMask_tree (TOTAL)", 65, MPI_Wtime()-t_cycle )

end subroutine
