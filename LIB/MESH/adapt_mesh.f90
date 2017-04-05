! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: adapt_mesh.f90
! version: 0.4
! author: msr, engels
!
! This routine performs the coarsing of the mesh, where possible. For the given mesh
! we compute the details-coefficients on all blocks. If four sister blocks have maximum
! details below the specified tolerance, (so they are insignificant), they are merged to
! one coarser block one level below. This process is repeated until the grid does not change
! anymore.
!
! As the grid changes, active lists and neighbor relations are updated, and load balancing
! is applied.
!
! NOTE: The block thresholding is done with the restriction/prediction operators acting on the
! entire block, INCLUDING GHOST NODES. Ghost node syncing is performed in threshold_block.
!
! NOTE: It is well possible to start with a very fine mesh and end up with only one active
! block after this routine. You do *NOT* have to call it several times.
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    ! heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    ! heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    ! list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: k, lgt_id, j
    real(kind=rk) :: xx0(1:3), ddx(1:3)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

do k = 1, hvy_n
    ! light id of this block
    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
    ! compute blocks' spacing from treecode
    call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

    ! HACK repair first datafield, as we're about to remove it
    hvy_block(:,:,:,1,hvy_active(k)) = 0.0_rk
    hvy_block(1,1,:,1,hvy_active(k)) = xx0(1)
    hvy_block(2,1,:,1,hvy_active(k)) = xx0(2)
    hvy_block(3,1,:,1,hvy_active(k)) = xx0(3)

    hvy_block(1,2,:,1,hvy_active(k)) = xx0(1)+ddx(1)
    hvy_block(2,2,:,1,hvy_active(k)) = xx0(2)+ddx(2)
    hvy_block(3,2,:,1,hvy_active(k)) = xx0(3)+ddx(3)
end do

    ! maximal number of loops to coarsen the mesh == one block go down from max_treelevel to min_treelevel
    do k = 1, (params%max_treelevel - params%min_treelevel)

        ! check where to coarsen (refinement done with safety zone)
        call threshold_block( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

        ! update lists of active blocks (light and heavy data)
        call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
        call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

        ! unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )

        ! ensure completeness
        call ensure_completeness( params, lgt_block, lgt_active, lgt_n )

        ! adapt the mesh
        if ( params%threeD_case ) then
            ! 3D:
            call coarse_mesh_3D( params, lgt_block, hvy_block, lgt_active, lgt_n )

        else
            ! 2D:
            call coarse_mesh_2D( params, lgt_block, hvy_block(:,:,1,:,:), lgt_active, lgt_n )

        end if

        ! update lists of active blocks (light and heavy data)
        call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
        call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

    end do

    ! balance load
    if ( params%threeD_case ) then
        ! 3D:
        call balance_load_3D( params, lgt_block, hvy_block, lgt_active, lgt_n )
    else
        ! 2D:
        call balance_load_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )
    end if

    ! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

end subroutine adapt_mesh
