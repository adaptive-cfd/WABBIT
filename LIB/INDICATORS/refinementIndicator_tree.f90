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

subroutine refinementIndicator_tree(params, hvy_block, tree_ID, indicator)
    implicit none
    type (type_params), intent(in)      :: params
    character(len=*), intent(in)        :: indicator                            !> how to choose blocks for refinement
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik) :: k, Jmax, max_blocks, ierr                               ! local variables
    ! chance for block refinement, random number
    real(kind=rk) :: ref_chance, r, nnorm(1:size(hvy_block,4)), max_grid_density, current_grid_density
    integer(kind=ik) :: hvy_id, lgt_id, Bs(1:3), g, tags, level
    real(kind=rk) :: a, b

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

            ! do not use normalizaiton (mask is inherently normalized to 0...1)
            nnorm = 1.0_rk

            if (params%dim == 3) then
                a = minval(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id))
                b = maxval(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id))
            else
                a = minval(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1, hvy_id))
                b = maxval(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1, hvy_id))
            endif

            ! exclude blocks which are all zero or all one from refinement.
            ! they are boring.
            if ( abs(a - b)>1.0e-7_rk ) then
                lgt_block(lgt_id, IDX_REFINE_STS) = +1
            endif
        enddo

        ! very important: CPU1 cannot decide if blocks on CPU0 have to be refined.
        ! therefore we have to sync the lgt data
        call synchronize_lgt_data( params, refinement_status_only=.true. )


    case ("mask-anynonzero")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! this criterion tags any block which has a value greater 0.0 in the first
        ! component to be refined.

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

            ! merge selects 2D or 3D bounds depending on params%dim
            if (any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g, params%dim == 2), 1, hvy_id) > 0.0_rk)) then
                lgt_block(lgt_id, IDX_REFINE_STS) = +1
            endif
            ! if (any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g, params%dim == 2), 6, hvy_id) > 0.0_rk)) then
            !     lgt_block(lgt_id, IDX_REFINE_STS) = +1
            ! endif
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
            tags = 0
            do k = 1, lgt_n(tree_ID)
                call random_number(r)
                ! set refinement status to refine
                if (r<=ref_chance .and. tags<max_blocks .and. lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS)==0) then
                    tags = tags + 1
                    lgt_block( lgt_active(k, tree_ID), IDX_REFINE_STS ) = 1
                end if
            end do
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
