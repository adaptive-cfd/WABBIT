
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

subroutine filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, hvy_active, hvy_n, hvy_neighbor)
   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure, hvy_active
    type (type_params), intent(inout)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> hvy_mask are qty that depend on the grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:, :)

    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, i
    ! grid parameter, error variable
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    ! surface normal
    integer(kind=2) :: n_domain(1:3)
    integer(kind=ik) :: level

    Bs = params%Bs
    g  = params%n_ghosts
    n_domain = 0

    if (params%filter_type == "wavelet_filter") then

        ! this filter removes all details (either on all blocks, only maxlevel or all except maxlevel)
        call wavelet_filter(time, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, hvy_active, hvy_n, hvy_neighbor)

    else

        do k = 1, hvy_n
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
            level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)

            if ((params%filter_only_maxlevel .and. level==params%max_treelevel) .or. &
                (.not.params%filter_only_maxlevel .and. .not.params%filter_all_except_maxlevel) .or. &
                (params%filter_all_except_maxlevel .and. level/=params%max_treelevel) ) then

                ! get block spacing for RHS
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                if ( .not. All(params%periodic_BC) ) then
                    ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                    call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
                    params%domain_size, params%Bs, params%dim, n_domain )
                endif

                call filter_meta(params%physics_type, time, hvy_block(:,:,:,:, hvy_active(k)), g, x0, dx,&
                hvy_tmp(:,:,:,:,hvy_active(k)), hvy_mask(:,:,:,:,hvy_active(k)), n_domain)
            endif
        enddo

    endif

end subroutine filter_wrapper


! In its standard setting, this filter does EXACTLY the same as coarsen everywhere followed by refine everywhere
! (.not.params%filter_only_maxlevel .and. .not.params%filter_all_except_maxlevel)
! ATTENTION: This is not entirely true, because in fact, we should sync_ghosts after coarsening, but that is not
! possible in the block based fashion. Thus, the mesh-based 
!   call adapt_mesh(..everywhere..)
!   call refine_mesh(..everywhere..)
! does something slightly different. You can also see that from the simple fact that coarseing of an arbitrary
! grid is not always possible for all blocks (completeness condition)
subroutine wavelet_filter(time, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, hvy_active, hvy_n, hvy_neighbor)
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> hvy_mask are qty that depend on the grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:, :)

    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, i
    ! grid parameter, error variable
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=ik) :: level
    real(kind=rk), allocatable, save :: u3(:,:,:)
    real(kind=rk), allocatable, save :: u2(:,:,:)

    Bs = params%Bs
    g  = params%n_ghosts

    !---------------------------------------------------------------------------
    ! lifted wavelet case (applies a low-pass filter first)
    !---------------------------------------------------------------------------
    if (params%wavelet_transform_type == 'biorthogonal') then
        do k = 1, hvy_n
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
            level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)

            if ((params%filter_only_maxlevel .and. level==params%max_treelevel) .or. &
                (.not.params%filter_only_maxlevel .and. .not.params%filter_all_except_maxlevel) .or. &
                (params%filter_all_except_maxlevel .and. level/=params%max_treelevel) ) then

                do i = 1, params%n_eqn
                    call restriction_prefilter_2D(hvy_block(:,:,1,i,hvy_active(k)), hvy_tmp(:,:,1,i,hvy_active(k)), params%wavelet)
                    hvy_block(:,:,1,i,hvy_active(k)) = hvy_tmp(:,:,1,i,hvy_active(k))
                enddo
            endif
        enddo

        ! required...was a bug. the first interpolated point uses one ghost node
        ! and ghost nodes are NOT correct after the filtering
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    endif


    !---------------------------------------------------------------------------
    ! Cycle of restriction/prediction to remove details. For non-lifted
    ! wavelets (CDF40, CDF20: harten-multiresolution) this is the only thing we do here
    !---------------------------------------------------------------------------
    if (params%dim == 3) then
        ! ********** 3D **********
        if (.not.allocated(u2)) allocate(u2( 1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g ))
        if (.not.allocated(u3)) allocate(u3((Bs(1)+1)/2+g, (Bs(2)+1)/2+g, (Bs(3)+1)/2+g))

        do k = 1, hvy_n
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
            level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)

            if ((params%filter_only_maxlevel .and. level==params%max_treelevel) .or. &
                (.not.params%filter_only_maxlevel .and. .not.params%filter_all_except_maxlevel) .or. &
                (params%filter_all_except_maxlevel .and. level/=params%max_treelevel) ) then

                ! now, coarsen array u1 (restriction)
                do i = 1, params%n_eqn
                    call restriction_3D( hvy_block(:,:,:,i,hvy_active(k)), u3 )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediction)
                    call prediction_3D( u3, hvy_block(:,:,:,i,hvy_active(k)), params%order_predictor )  ! coarse, fine
                enddo
            endif
        enddo
    else
        ! ********** 2D **********
        if (.not.allocated(u2)) allocate(u2( 1:Bs(1)+2*g, 1:Bs(2)+2*g, 1 ))
        if (.not.allocated(u3)) allocate(u3((Bs(1)+1)/2+g, (Bs(2)+1)/2+g, 1))

        do k = 1, hvy_n
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
            level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)

            if ((params%filter_only_maxlevel .and. level==params%max_treelevel) .or. &
                (.not.params%filter_only_maxlevel .and. .not.params%filter_all_except_maxlevel) .or. &
                (params%filter_all_except_maxlevel .and. level/=params%max_treelevel) ) then

                do i = 1, params%n_eqn
                    ! now, coarsen array u1 (restriction)
                    call restriction_2D( hvy_block(:,:,1,i,hvy_active(k)), u3(:,:,1) )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediction)
                    call prediction_2D( u3(:,:,1), hvy_block(:,:,1,i,hvy_active(k)), params%order_predictor )  ! coarse, fine
                enddo
            endif
        enddo
    endif
end subroutine wavelet_filter
