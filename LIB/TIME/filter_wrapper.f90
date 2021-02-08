
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \file
!> \callgraph
!> \brief wrapper for filter
!> \version 0.5
!> \author Pkrah
!! \date 30/04/18 - create
!!
!
!**********************************************************************************************

subroutine filter_wrapper(time, params, hvy_block, hvy_tmp, lgt_block, hvy_active, hvy_n)
   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure, hvy_active
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> global integral
    real(kind=rk), dimension(3)         :: volume_int

    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, neqn, lgt_id
    ! grid parameter, error variable
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    !  surface normal
    integer(kind=2) :: n_domain(1:3)
    integer         :: level

    Bs    = params%Bs
    g     = params%n_ghosts
    n_domain = 0


    do k = 1, hvy_n
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

        ! level of the block:
        level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)

        if ((params%filter_only_maxlevel .and. level==params%max_treelevel) .or. .not. params%filter_only_maxlevel) then
            ! get block spacing for RHS
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            if ( .not. All(params%periodic_BC) ) then
                ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
                params%domain_size, params%Bs, params%dim, n_domain )
            endif

            call filter_meta(params%physics_type, time, hvy_block(:,:,:,:, hvy_active(k)), g, x0, dx,&
            hvy_tmp(:,:,:,:,hvy_active(k)), n_domain)
        endif
    enddo

end subroutine filter_wrapper
