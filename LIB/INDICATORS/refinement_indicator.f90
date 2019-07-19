!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name refinement_indicator.f90
!> \version 0.5
!> \author engels
!> \brief Set refinement status for active blocks, different methods possible
!
!> \details This routine sets the refinement flag for all blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very compley, but expected to grow
!! in the future.
!! \n
!! ------------------
!! Refinement status:
!! ------------------
!! +1 refine
!! 0 do nothing
!! -1 block wants to refine (ignoring other constraints, such as gradedness)
!! -2 block will refine and be merged with her sisters
!! ------------------
!! \n
!! = log ======================================================================================
!! \n
!! 23/05/2017 create
! ********************************************************************************************

subroutine refinement_indicator( params, lgt_block, lgt_active, lgt_n, hvy_block, hvy_active, hvy_n, indicator )
    use module_helpers
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    !> heavy data array - block data
    real(kind=rk), intent(inout)           :: hvy_block(:, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)        :: hvy_n
    ! local variables
    integer(kind=ik) :: k, Jmax, max_blocks, ierr
    ! chance for block refinement, random number
    real(kind=rk) :: ref_chance, r, nnorm(1:size(hvy_block,4))
    logical :: used(1:size(hvy_block,4))
    integer(kind=ik) :: hvy_id, lgt_id, Bs(1:3), g


    Jmax = params%max_treelevel
    Bs = params%Bs
    g = params%n_ghosts

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
        do k = 1, lgt_n
           lgt_block(lgt_active(k), Jmax + IDX_REFINE_STS) = 0
        enddo

        ! each CPU decides for its blocks if they're refined or not
        do k = 1, hvy_n
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k)

            ! light id of this block
            call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! do not use normalizaiton (mask is inherently normalized to 0...1)
            nnorm = 1.0_rk
            ! threshold only the first component of the vector (which is the mask, the others
            ! are forcing values, colors, sponges, etc)
            used = .false.
            used(1) = .true.

            call threshold_block( params, hvy_block(:,:,:,:, hvy_id), used, &
            lgt_block(lgt_id, Jmax + IDX_REFINE_STS), nnorm, eps=1.0e-6_rk )

            ! hack: currently, threshold_block assigns only -1 or 0
            lgt_block(lgt_id, Jmax + IDX_REFINE_STS) = lgt_block(lgt_id, Jmax + IDX_REFINE_STS)+1
        enddo

        ! very important: CPU1 cannot decide if blocks on CPU0 have to be refined.
        ! therefore we have to sync the lgt data
        call synchronize_lgt_data( params, lgt_block, refinement_status_only=.true. )


    case ("mask-anynonzero")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! this criterion tags any block which has a value greater 0.0 in the first
        ! component to be refined.

        ! reset refinement status to "stay"
        do k = 1, lgt_n
           lgt_block(lgt_active(k), Jmax + IDX_REFINE_STS) = 0
        enddo

        ! each CPU decides for its blocks if they're refined or not
        if (params%dim==3) then
            do k = 1, hvy_n
                ! hvy_id of the block we're looking at
                hvy_id = hvy_active(k)

                ! light id of this block
                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

                if (any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id) > 0.0_rk)) then
                    lgt_block(lgt_id, Jmax + IDX_REFINE_STS) = +1
                endif

                ! if (any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 6, hvy_id) > 0.0_rk)) then
                !     lgt_block(lgt_id, Jmax + IDX_REFINE_STS) = +1
                ! endif
            enddo
        else
            do k = 1, hvy_n
                ! hvy_id of the block we're looking at
                hvy_id = hvy_active(k)

                ! light id of this block
                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

                if (any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1, hvy_id) > 0.0_rk)) then
                    lgt_block(lgt_id, Jmax + IDX_REFINE_STS) = +1
                endif
                ! if (any(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 6, hvy_id) > 0.0_rk)) then
                !     lgt_block(lgt_id, Jmax + IDX_REFINE_STS) = +1
                ! endif
            enddo
        endif

        ! very important: CPU1 cannot decide if blocks on CPU0 have to be refined.
        ! therefore we have to sync the lgt data
        call synchronize_lgt_data( params, lgt_block, refinement_status_only=.true. )

    case ("everywhere")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! set status "refine" for all active blocks, which is just setting the
        ! last index in the light data block list to +1. This indicator is used
        ! to refine the entire mesh at the beginning of a time step, if error
        ! control is desired.
        do k = 1, lgt_n
            ! do not refine blocks with +11 status, as they are on the maxlevel
            ! already (so no refinement allowed)
            if ( lgt_block( lgt_active(k), Jmax + IDX_REFINE_STS ) /= 11 ) then
                lgt_block( lgt_active(k), Jmax + IDX_REFINE_STS ) = +1
            end if
        end do

    case ("random")
        !(((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))
        ! randomized refinement. This can be used to generate debug meshes for
        ! testing purposes. For example the unit tests use this.
        ref_chance = 0.10_rk
        ! random refinement can set at most this many blocks to refine (avoid errors
        ! due to insufficient memory) (since we already have lgt_n blocks we can set the status
        ! at most for Nmax-lgt_n blocks)
        max_blocks = (size(lgt_block,1) - lgt_n ) ! / 4 ! 4: safety
        ! each block flagged for refinement creates (2**d-1) new blocks
        max_blocks = ( max_blocks / (2**params%dim-1) )
        ! safety
        max_blocks = max_blocks / 40

        ! set random seed
        call init_random_seed()

        ! unset all refinement flags
        lgt_block( :, Jmax + IDX_REFINE_STS ) = 0

        ! only root rank sets the flag, then we sync. It is messy if all procs set a
        ! random value which is not sync'ed
        if (params%rank == 0) then
            do k = 1, lgt_n
                ! random number
                call random_number(r)
                ! set refinement status to refine
                if ( r <= ref_chance .and. sum(lgt_block(:, Jmax+ IDX_REFINE_STS)) <= max_blocks) then
                    lgt_block( lgt_active(k), Jmax + IDX_REFINE_STS ) = 1
                else
                    lgt_block( lgt_active(k), Jmax + IDX_REFINE_STS ) = 0
                end if
            end do
        endif
        ! sync light data, as only root sets random refinement
        call MPI_BCAST( lgt_block(:, Jmax + IDX_REFINE_STS), size(lgt_block,1), MPI_INTEGER4, 0, WABBIT_COMM, ierr )

    case default
        call abort("ERROR: refine_mesh: the refinement indicator is unkown")

    end select

end subroutine refinement_indicator
