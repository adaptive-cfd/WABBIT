!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name coarsening_indicator.f90
!> \version 0.5
!> \author engels
!> \brief Set coarsening status for a single blocks, different methods possible
!
!> \details This routine sets the coarsening flag for a single blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very complex, but expected to grow
!! in the future.
!! \n
!! ------------------ \n
!! Refinement status: \n
!! ------------------ \n
!! +1 refine \n
!! 0 do nothing \n
!! -1 block wants to coarsen \n
!! ------------------ \n
!! \n
!! = log ======================================================================================
!! \n
!! 23/05/2017 create
! ********************************************************************************************

subroutine block_coarsening_indicator( params, block_data, block_work, dx, x0, indicator, &
    iteration, refinement_status, norm, block_mask)

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array.
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    ! NOTE: Here, the mask is required only if grid adaptation is also done on the mask.
    real(kind=rk), intent(inout), optional :: block_mask(:, :, :, :)
    !> heavy work data array (expected to hold the VORTICITY if thresholding is applied to vorticity)
    real(kind=rk), intent(inout)        :: block_work(:, :, :, :)
    !> block spacing and origin
    real(kind=rk), intent(in)           :: dx(1:3), x0(1:3)
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    !> coarsening iteration index. coarsening is done until the grid has reached
    !! the steady state; therefore, this routine is called several times for
    integer(kind=ik), intent(in)        :: iteration
    !> output is the refinement_status
    integer(kind=ik), intent(out)       :: refinement_status
    !
    real(kind=rk), intent(inout)        :: norm(1:size(block_data,4))

    ! local variables
    integer(kind=ik) :: k, Jmax, d, j, hvy_id, g, refinement_status_mask
    integer(kind=ik), dimension(3) :: Bs
    ! chance for block refinement, random number
    real(kind=rk) :: crsn_chance, r, nnorm(1)
    logical :: thresholding_component(1:size(block_data,4))
    logical :: tmp_threshold(1:20) ! just take a larger one...lazy tommy
    real(kind=rk) :: nnorm2(1:20) ! just take a larger one...lazy tommy

!---------------------------------------------------------------------------------------------
! variables initialization

    Jmax = params%max_treelevel
    Bs = params%Bs
    g = params%n_ghosts


!---------------------------------------------------------------------------------------------
! main body

    !> This routine sets the -1 coarsening flat on a block. it uses different methods to
    !! decide where to coarsen, each act on one block. Note due to gradedness and completeness
    !! this status may be revoked later in the computation.
    select case (indicator)
    case ("everywhere")
        ! simply coarsen the entire grid. Note that this means that adapt_mesh will coarsen down to Jmin
        ! if the iteration loop is on (you can bypass that behavior using external_loop=.true.)
        refinement_status = -1

    case ("threshold-vorticity")
        !! use thresholding, but on the vorticity rather than the state vector.
        !! this should allow to consider only the rotational part, not the divergent one.

        thresholding_component = .false.
        if (params%dim == 3) then
            thresholding_component(1:3) = .true.
        else
            thresholding_component(1) = .true.
        endif

        !! note we assume block_work contains the vorticity
        call threshold_block( params, block_work, thresholding_component, refinement_status, norm )

    case ("threshold-state-vector","primary-variables")
        !! use wavelet indicator to check where to coarsen. Note here, active components are considered
        !! and the max over all active components results in the coarsening state -1. The components
        !! to be used can be specified in the PARAMS file. default is all componants.

        if (.not. allocated(params%threshold_state_vector_component)) then
            call abort(7363823, "params%threshold_state_vector_component not allocated....")
        endif

        thresholding_component = params%threshold_state_vector_component
        call threshold_block( params, block_data, thresholding_component, refinement_status, norm )

    case ("random")
        !! randomly coarse some blocks. used for testing. note we tag for coarsening
        !! only once in the first iteration. this is important: as adapt_mesh is an iterative
        !! routine that calls the coarsening until the grid does not change anymore. without
        !! the iteration==0 check, it will always keep on coarsening.

        if (iteration == 0) then
            ! call init_random_seed()
            ! the chance for coarsening:
            crsn_chance = 0.25_rk
            ! random number
            call random_number(r)
            ! set refinement status to coarsen based on random numbers.
            if ( r <= crsn_chance ) refinement_status = -1
        endif

    case default
        call abort(151413,"ERROR: unknown coarsening operator: "//trim(adjustl(indicator)))

    end select


    ! mask thresholding on top of regular thresholding?
    ! it can be useful to also use the mask function (if penalization is used) for grid adaptation.
    ! i.e. the grid is always at the finest level on mask interfaces. Careful though: the Penalization
    ! is implemented on physics-module level, i.e. it is not available for all modules.  If it is
    ! not available, the option is useless but can cause errors.
    if (params%threshold_mask .and. present(block_mask)) then
        ! assuming block_mask holds mask function
        nnorm2 = 1.0_rk
        tmp_threshold = .false.
        tmp_threshold(1) = .true.

        ! even if the global eps is very large, we want the mask to be on the finest
        ! level. hence, here we set a small value (just for this call) to be sure that the
        ! mask interface is on Jmax
        call threshold_block( params, block_mask, tmp_threshold(1:size(block_mask,4)), &
        refinement_status_mask, nnorm2(1:size(block_mask,4)), eps=1.0e-4_rk )

        ! refinement_status_state: -1 refinemet_status_mask: -1 ==>  -1
        ! refinement_status_state: 0  refinemet_status_mask: -1 ==>   0
        ! refinement_status_state: 0  refinemet_status_mask: 0  ==>   0

        refinement_status = max(refinement_status, refinement_status_mask)
    endif

end subroutine block_coarsening_indicator
