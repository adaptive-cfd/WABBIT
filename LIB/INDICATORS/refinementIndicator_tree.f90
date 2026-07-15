!> \brief Set refinement status for active blocks, different methods possible
!> \details This routine sets the refinement flag for all blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very compley, but expected to grow
!! in the future.
!! ------------------
!! Refinement status:
!! ------------------
!! +1 refine
!! 0 do nothing
!! -1 block wants to refine (ignoring other constraints, such as gradedness)
!! -2 block will refine and be merged with her sisters
! ********************************************************************************************

subroutine refinementIndicator_tree(params, hvy_block, tree_ID, indicator, time)
    use module_physics_metamodule, only : geometry_indicator_meta  ! this is unfortunate
    implicit none
    type (type_params), intent(in)      :: params
    character(len=*), intent(in)        :: indicator                            !> how to choose blocks for refinement
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    real(kind=rk), intent(in), optional :: time  !> current simulation time, used for mask checking

    integer(kind=ik) :: k, Jmax, max_blocks, ierr                               ! local variables
    ! chance for block refinement, random number
    real(kind=rk) :: ref_chance, r, max_grid_density, current_grid_density
    integer(kind=ik) :: hvy_id, lgt_id, Bs(1:3), g, tags, level, ref_status
    real(kind=rk) :: mask_max, mask_min, x0(1:3), dx(1:3)

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas



    Jmax = params%Jmax
    Bs = params%Bs
    g = params%g

    !> loop over the blocks and set their refinement status.
    !! \note refinement is an absolute statement, that means once set, the block will be refined
    !! (which is not the case in block coarsening), it may even entrail other blocks in
    !! its vicinity to be refined as well.
    select case (indicator)
    case ("mask-threshold")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! this criterion uses the first component (and ONLY the first one) of the given data
        ! (useful only if this is the mask function), checks the details and tags the block
        ! for refinement if the details are significant.

        ! we check from g+1 until BS+g+1, so the first ghost layer is included, in order to check the bounds around all possible new points
        ! if not, the point at position BS+g+1/2 might be on the ghost layer but we would not detect that, even though it should be refined.

        ! reset refinement status to "stay"
        do k = 1, lgt_n(tree_ID)
           lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS) = 0
        enddo

        ! each CPU decides for its blocks if they're refined or not
        do k = 1, hvy_n(tree_ID)
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k, tree_ID)
            ! light id of this block
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            level = lgt_block( lgt_id, IDX_MESH_LVL)
            ! get block spacing and origin for the geometry indicator
            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
                params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)

            if (any(hvy_block(g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g+1, params%dim == 2), 1, hvy_id) > 1.0e-12_rk .and. &
                hvy_block(g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g+1, params%dim == 2), 1, hvy_id) < 1.0_rk-1.0e-12_rk)) then
                lgt_block(lgt_id, IDX_REFINE_STS) = +1
            endif
            mask_max = maxval(hvy_block(g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g+1, params%dim == 2), 1, hvy_id))
            mask_min = minval(hvy_block(g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g+1, params%dim == 2), 1, hvy_id))

            ! exclude blocks which are all zero or all one from refinement.
            ! they are boring.
            if ( abs(mask_max - mask_min) > 1.0e-12_rk ) then
                ! max and min value are different - the mask is not all constant on this block.
                ! set block status to REFINE.
                lgt_block(lgt_id, IDX_REFINE_STS) = +1
            endif

            ! In rare cases (mostly with large domains and small Jmin), the grid may be so coarse that (dx >> object size). In this case, no
            ! point on the grid may be assigned a nonzero mask value, even though the object (geometry) lies within this block. The "geometry_indicator"
            ! checks if the origin of the object lies within the blocks extend, and returns +1 if this is the case.
            call geometry_indicator_meta(params%physics_type, time, params%Bs, params%g, x0, dx, ref_status, "coarsening")
            if (ref_status == +1) then
                ! The origin of the object (geometry) is in fact inside this block.
                ! Block has to refine if geometry_indicator says so AND the whole mask is 0 - then 
                ! the rare exception did occur and we accidentally did not create the mask. 
                if (mask_max < 1.0e-12_rk) then
                    lgt_block(lgt_ID, IDX_REFINE_STS) = +1 ! REFINE
                endif
            endif
        enddo

        ! very important: CPU1 cannot decide if blocks on CPU0 have to be refined.
        ! therefore we have to sync the lgt data
        call synchronize_lgt_data( params, refinement_status_only=.true. )


    case ("mask-anynonzero")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! this criterion tags any block which has a value greater 0.0 in the first
        ! component to be refined.

        ! we check from g+1 until BS+g+1, so the first ghost layer is included, in order to check the bounds around all possible new points
        ! if not, the point at position BS+g+1/2 might be on the ghost layer but we would not detect that, even though it should be refined.

        ! reset refinement status to "stay"
        do k = 1, lgt_n(tree_ID)
           lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS) = 0
        enddo

        ! each CPU decides for its blocks if they're refined or not
        do k = 1, hvy_n(tree_ID)
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k, tree_ID)
            ! light id of this block
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! In rare cases (mostly with large domains and small Jmin), the grid may be so coarse that (dx >> object size). In this case, no
            ! point on the grid may be assigned a nonzero mask value, even though the object (geometry) lies within this block. The "geometry_indicator"
            ! checks if the origin of the object lies within the blocks extend, and returns +1 if this is the case.
            call geometry_indicator_meta(params%physics_type, time, params%Bs, params%g, x0, dx, ref_status, "coarsening")
            if (ref_status == +1) then
                ! The origin of the object (geometry) is in fact inside this block.
                ! Block has to refine if geometry_indicator says so AND the whole mask is 0 - then 
                ! the rare exception did occur and we accidentally did not create the mask. 
                if (mask_max < 1.0e-9_rk) then
                    lgt_block(lgt_ID, IDX_REFINE_STS) = +1 ! REFINE
                endif
            endif

            ! merge selects 2D or 3D bounds depending on params%dim
            if (any(hvy_block(g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g+1, params%dim == 2), 1, hvy_id) > 0.0_rk)) then
                lgt_block(lgt_id, IDX_REFINE_STS) = +1
            endif
        enddo

        ! very important: CPU1 cannot decide if blocks on CPU0 have to be refined.
        ! therefore we have to sync the lgt data
        call synchronize_lgt_data( params, refinement_status_only=.true. )

    case ("significant")
        ! we assume, that the following refinement flags have been set in the last adapt_tree call and still persist:
        !                      0 : this block is significant, it definetly needs to refine -> +1
        ! REF_UNSIGNIFICANT_STAY : this block was only kept to ensure completeness or a graded mesh, we do not need to refine it -> 0
        !                     -1 : this should not exist, as all blocks with that have been deleted

        ! coarseningIndicator and Securityzone work for coarsening, so -1, we now have to translate this to refinement range
        do k = 1, lgt_n(tree_ID)
            lgt_ID = lgt_active(k, tree_ID)

            ! 1 and -1 should not exist here
            if (any(lgt_block( lgt_id, IDX_REFINE_STS ) == (/ -1, 1/))) call abort(241119, "I am very confused by what is going on here and do not like it!")
            ! set 0 (significant) to +1 (refine), -9 (non-significant) to 0 (stay as it is)
            if (lgt_block( lgt_id, IDX_REFINE_STS ) == 0) lgt_block( lgt_id, IDX_REFINE_STS ) = +1
            if (lgt_block( lgt_id, IDX_REFINE_STS ) == REF_UNSIGNIFICANT_STAY) lgt_block( lgt_id, IDX_REFINE_STS ) = 0
        enddo

    case ("everywhere")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! set status "refine" for all active blocks, which is just setting the
        ! last index in the light data block list to +1. This indicator is used
        ! to refine the entire mesh at the beginning of a time step, if error
        ! control is desired.
        do k = 1, lgt_n(tree_ID)
            lgt_block( lgt_active(k, tree_ID), IDX_REFINE_STS ) = +1
        end do
    
    case ("nowhere")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! set status "stay" for all active blocks, which is just setting the
        ! last index in the light data block list to 0. This indicator is used
        ! to keep the same grid. Might sound stupid but might be useful for debugging purposes
        do k = 1, lgt_n(tree_ID)
            lgt_block( lgt_active(k, tree_ID), IDX_REFINE_STS ) = 0
        end do

    case ("random")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! set random seed
        call init_random_seed()

        ! random refinement can set at most this many blocks to refine (avoid errors
        ! due to insufficient memory)
        max_grid_density = params%max_grid_density
        current_grid_density = dble(lgt_n(tree_ID)) / dble(size(lgt_block,1))
        ! we allow at most this number of blocks to be active:
        max_blocks = floor( max_grid_density * dble(size(lgt_block,1)) )
        ! we already have lgt_n blocks we can set the status
        ! at most for Nmax-lgt_n blocks
        max_blocks = max_blocks - lgt_n(tree_ID)
        ! each block flagged for refinement creates (2**d -1) new blocks
        max_blocks = max_blocks / (2**params%dim - 1)
        ! chance for randomized refinement

        ! the idea is to generate a grid as close as possible to the desired max_grid_density
        ! so if the desired density is high, and the grid sparse, the probability to refine is high
        ! but we still cap it at some value
        ref_chance = 1.0_rk - current_grid_density / max_grid_density

        ! unset all refinement flags
        lgt_block( :, IDX_REFINE_STS ) = 0

        ! only root sets the flag, then we sync. It is messy if all procs set a
        ! random value which is not sync'ed
        if (params%rank == 0) then
            ! First pass: check which blocks would like to refine based on random chance
            tags = 0
            do k = 1, lgt_n(tree_ID)
                call random_number(r)
                ! temporarily mark blocks that pass the random test
                if (r<=ref_chance .and. lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS)==0) then
                    lgt_block( lgt_active(k, tree_ID), IDX_REFINE_STS ) = 1
                    tags = tags + 1
                end if
            end do
            
            ! Second pass: if too many blocks want to refine, randomly deselect some
            ! This ensures spatial uniformity rather than bias toward early blocks
            if (tags > max_blocks) then
                do k = 1, lgt_n(tree_ID)
                    if (lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS) == 1) then
                        call random_number(r)
                        ! keep only a fraction of the tagged blocks
                        if (r > dble(max_blocks)/dble(tags)) then
                            lgt_block( lgt_active(k, tree_ID), IDX_REFINE_STS ) = 0
                        end if
                    end if
                end do
            end if
        endif

        ! sync light data, as only root sets random refinement
        call MPI_BCAST( lgt_block(:, IDX_REFINE_STS), size(lgt_block,1), MPI_INTEGER4, 0, WABBIT_COMM, ierr )

    ! In some cases the refinement_status is set up by a routine other than this one. This is the case
    ! in some forest processing (=handling multiple trees). In such a case, we do nothing here (in particular
    ! we do not reset the refinement status)
    case ('nothing (external)')
        return

    case default
        call abort("ERROR: refine_tree: the refinement indicator is unkown")

    end select

end subroutine refinementIndicator_tree
