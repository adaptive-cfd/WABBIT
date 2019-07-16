!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_random_grid.f90
!> \version 0.5
!> \author engels
subroutine create_random_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, Jmin, verbosity, iterations )

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)   :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: 00, k1, k2 etc)
    real(kind=rk), intent(out)          :: hvy_work(:, :, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)          :: hvy_tmp(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> what level to initialize?
    integer(kind=ik), intent(in)        :: Jmin, iterations
    !> write output
    logical, intent(in)                 :: verbosity

    integer :: l

    ! set all blocks to free (since if we call inicond twice, all blocks are used in the second call)
    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )

    ! setup the coarsest grid level with some data (we don't care what data, we'll erase it)
    ! Note that active lists + neighbor relations are updated inside this routine as well, as
    ! the grid is modified
    call create_equidistant_grid( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, 2, .true. )

    !---------------------------------------------------------------------------------------------
    ! second: refine some blocks (random), coarsen some blocks (random)
    do l = 1, iterations
        if (params%rank==0 .and. verbosity) then
            write(*,'("RANDOM GRID GENERATION: iteration ",i1," active=",i4," Jmin=",i2," Jmax=",i2)') &
            l, lgt_n, &
            min_active_level( lgt_block, lgt_active, lgt_n ), &
            max_active_level( lgt_block, lgt_active, lgt_n )
        endif

        ! randomly refine some blocks
        call refine_mesh( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, lgt_active, lgt_n, &
             lgt_sortednumlist, hvy_active, hvy_n, "random" )

        ! randomly coarsen some blocks
        call adapt_mesh( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
             lgt_n, lgt_sortednumlist, hvy_active, hvy_n, "random", hvy_tmp )
    end do
end subroutine create_random_grid
