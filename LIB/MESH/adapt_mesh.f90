! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: adapt_mesh.f90
! version: 0.4
! author: msr
!
! mesh adapting main function
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
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :)
    ! heavy data array - neifghbor data
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
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization


!---------------------------------------------------------------------------------------------
! main body

    ! maximal number of loops to coarsen the mesh == one block go down from max_treelevel to min_treelevel
    do k = 1, (params%max_treelevel - params%min_treelevel)

        ! check where to coarsen (refinement done with safety zone)
        call threshold_block( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

        ! unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )

        ! ensure completeness
        call ensure_completeness( params, lgt_block, lgt_active, lgt_n )

        ! adapt the mesh
        call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n )

        ! update lists of active blocks (light and heavy data)
        call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
        call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

    end do
!
!    ! balance load
!    call balance_load( params, block_list, block_data, neighbor_list )

    ! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

end subroutine adapt_mesh
