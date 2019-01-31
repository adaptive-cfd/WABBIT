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
!! -1 block wants to coarsen (ignoring other constraints, such as gradedness) \n
!! -2 block will coarsen and be merged with her sisters \n
!! ------------------ \n
!! \n
!! = log ======================================================================================
!! \n
!! 23/05/2017 create
! ********************************************************************************************

subroutine block_coarsening_indicator( params, block_data, block_work, dx, x0, indicator, &
    iteration, refinement_status, norm)

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array.
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
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
    real(kind=rk), intent(inout)        :: norm(1:params%n_eqn)

    ! local variables
    integer(kind=ik) :: k, Jmax, d, j, hvy_id, g
    integer(kind=ik), dimension(3) :: Bs
    ! chance for block refinement, random number
    real(kind=rk) :: crsn_chance, r
    logical :: thresholding_component(1:params%n_eqn)

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
    case ("threshold-vorticity")
        !! use thresholding, but on the vorticity rather than the state vector.
        !! this should allow to consider only the rotational part, not the divergent one.

        thresholding_component = .false.
        if (params%threeD_case) then
            thresholding_component(1:3) = .true.
        else
            thresholding_component(1) = .true.
        endif

        !! note we assume hvy_work contains the vorticity
        call threshold_block( params, block_work, thresholding_component, refinement_status, norm )

    case ("threshold-state-vector")
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

end subroutine block_coarsening_indicator
