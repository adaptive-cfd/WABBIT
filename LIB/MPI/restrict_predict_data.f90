subroutine restrict_predict_data( params, res_pre_data, ijk, neighborhood, &
    level_diff, hvy_block, num_eqn, hvy_id, ignore_Filter )

    implicit none

    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(out)                      :: res_pre_data(:,:,:,:)
    !> indices in x,y,z direction of the ghost node patch
    integer(kind=ik), intent(in)                    :: ijk(2,3)
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)                    :: num_eqn      !< How many components? Needed as in between we use hvy_tmp
    integer(kind=ik), intent(in)                    :: hvy_id
    logical, intent(in)                             :: ignore_Filter !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets

    integer(kind=ik) :: iy  ! debug variable

    ! some neighborhoods are intrinsically on the same level (level_diff=0)
    ! and thus it makes no sense to call the up/downsampling routine for those
#ifdef DEV
    if ( params%dim == 3 .and. (neighborhood<=18) ) call abort(323223,"this case shouldnt appear")
    if ( params%dim == 2 .and. (neighborhood<=4) ) call abort(323223,"this case shouldnt appear")
#endif

    if ( level_diff == -1 ) then
        ! The neighbor is finer: we have to predict the data
        call predict_data( params, res_pre_data, ijk, hvy_block, num_eqn, hvy_id )

    elseif ( level_diff == +1) then
        ! The neighbor is coarser: we have to downsample the data
        call restrict_data( params, res_pre_data, ijk, hvy_block, num_eqn, hvy_id, ignore_Filter )

    else
        call abort(123005, "Lord Vader, restrict_predict_data is called with leveldiff /= -+1")

    end if

end subroutine restrict_predict_data

subroutine restrict_data( params, res_data, ijk, hvy_block, num_eqn, hvy_id, ignore_Filter )
    implicit none

    type (type_params), intent(in)  :: params
    !> data buffer
    real(kind=rk), intent(out)      :: res_data(:,:,:,:)
    !> ijk
    integer(kind=ik), intent(in)    :: ijk(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)    :: num_eqn      !< How many components? Needed as in between we use hvy_tmp
    integer(kind=ik), intent(in)    :: hvy_id
    logical, intent(in)             :: ignore_Filter !< If set, coarsening will be done only with loose downsampling, not applying HD filter even in the case of lifted wavelets

    integer(kind=ik)                :: iy, nc, nx, ny, nz

    nx = size(hvy_block,1)
    ny = size(hvy_block,2)
    nz = size(hvy_block,3)
    nc = num_eqn

#ifdef DEV
    if (.not. allocated(params%HD)) call abort(230301051, "Pirates! Maybe setup_wavelet was not called?")
#endif

    if (.not. ignore_Filter) then
        ! applying the filter is expensive, and we therefore apply it only once to the entire
        ! block. 
        ! Pro: we maybe save a bit of CPU time, as we usually do not need the entire block filtered
        ! Con: more work and maybe we compute some values twice (if patches overlap)
        ! The if-clause checks if the blocks hvy_ID is equal to the hvy_ID of the block currently held
        ! (filtered) in the block-array "hvy_restricted". This way, we hold only one filtered block 
        ! in memory, but still apply the filter only once to every block. Filter is applied to entire
        ! block.
        if (hvy_ID /= restricted_hvy_ID) then
            call restrict_copy_at_CE(params, hvy_block, hvy_ID, nc)
            ! block work array now contains filtered version of hvy_ID:
            restricted_hvy_ID = hvy_ID
        endif

        ! ! apply the filter but only on the patch ijk on which we are interested, as especially for 3D computing it on the whole block is expensive
        ! call restrict_copy_at_CE(params, hvy_block, hvy_ID, nc, ijk)

        res_data( 1:(ijk(2,1)-ijk(1,1))/2+1, 1:(ijk(2,2)-ijk(1,2))/2+1, 1:(ijk(2,3)-ijk(1,3))/2+1, 1:nc) = &
            hvy_restricted( ijk(1,1):ijk(2,1):2, ijk(1,2):ijk(2,2):2, ijk(1,3):ijk(2,3):2, 1:nc)
    else
        ! ignoring the filter - copy from the hvy_block array directly (corresponds to using CDFX0 wavelets, even if globally lifted ones are used.)
        ! This functionality is used in the RHS evaluation, because we do not perform a coarseExtension there. It saves on CPU time to then bypass
        ! the filter, which should not matter anyways, as the first few ghost nodes on a coarse block that are sync'ed with its fine neighbor
        ! are copied anyways (this is part of the coarseExtension; the filter HD cannot be applied near a coarse/fine interface). The FD stencils are 
        ! shorter than the SC copied in soarseExtension.
        res_data( 1:(ijk(2,1)-ijk(1,1))/2+1, 1:(ijk(2,2)-ijk(1,2))/2+1, 1:(ijk(2,3)-ijk(1,3))/2+1, 1:nc) = &
            hvy_block( ijk(1,1):ijk(2,1):2, ijk(1,2):ijk(2,2):2, ijk(1,3):ijk(2,3):2, 1:nc, hvy_ID)
    endif
end subroutine restrict_data



!> \brief Uses HD filter at all interfaces and copies values everywhere where it would use from coarser or finer neighbor
!> Used for sync in order to match values with coarse extension. Coarser AND finer neighbors cannot be filtered as there values
!! are not present in the ghost layers, so for all those patches where filters need those points the values will be copied.
!
!                      g g o o o o o o         - - - - - - - -
!                      g g o o o o o o         - - - - - - - -
!      m c c c         g g i i i i g g         - - i i i i - -
!      m b b m         g g i i i i g g         - - f f f f - -
!      m b b m         g g i i i i g g         - - f f f f - -
!      f m m m         g g i i i i g g         - - i f f f - -
!                      o o g g g g g g         - - - - - - - -
!                      o o g g g g g g         - - - - - - - -
!   Block relations       Block in                Block out
!
! b = block, m = medium neighbor (same level), c = coarse neighbor, f = fine neighbor
! g = ghost point, o = outdated ghost point, i = interior point (input data), f = filtered point, - = value should not be used
! Note: The bottom left value is copied (i) due to the finer neighbor in the bottom left corner, the top values due to the coarser top neighbor.
! Ghost points of output should not be used, as they may have been filtered correctly, wrongly, partially or left unfiltered arbitrarily
subroutine restrict_copy_at_CE(params, hvy_data, hvy_ID, num_eqn, ijk)

    implicit none
    type (type_params), intent(in) :: params  !< good ol params
    real(kind=rk), intent(inout)   :: hvy_data(:, :, :, :, :)  !< heavy data array - block data
    integer(kind=ik), intent(in)   :: hvy_id  !< which block to look at
    integer(kind=ik), intent(in)   :: num_eqn !< How many components? Needed as it could vary
    integer(kind=ik), intent(in), optional   :: ijk(2,3)  !< ijk, if this is present we only filter this part

    integer(kind=ik) :: k_n, lvl_me, lvl_diff, lgt_ID, lgt_ID_n, HD_l, HD_r

    HD_l = lbound(params%HD, dim=1)
    HD_r = ubound(params%HD, dim=1)

    if (.not. present(ijk)) then
        ! filter whole block and don't care about neighbors first
        call blockFilterXYZ_vct( params, hvy_data(:,:,:,1:num_eqn, hvy_id), hvy_restricted(:,:,:,1:num_eqn), params%HD, HD_l, HD_r, do_restriction=.true.)
    else
        ! only filter parts of the block that we need
        ! we pass do_restriction=.true. in order to not filter parts we will later restrict
        if (params%dim == 2) then
            call blockFilterXYZ_wherepossible_vct( params, &
                hvy_data(      ijk(1,1)+HD_l:ijk(2,1)+HD_r, ijk(1,2)+HD_l:ijk(2,2)+HD_r, ijk(1,3):ijk(2,3), 1:num_eqn, hvy_id), &
                hvy_restricted(ijk(1,1)+HD_l:ijk(2,1)+HD_r, ijk(1,2)+HD_l:ijk(2,2)+HD_r, ijk(1,3):ijk(2,3), 1:num_eqn), params%HD, HD_l, HD_r, do_restriction=.true.)
        else
            call blockFilterXYZ_wherepossible_vct( params, &
                hvy_data(      ijk(1,1)+HD_l:ijk(2,1)+HD_r, ijk(1,2)+HD_l:ijk(2,2)+HD_r, ijk(1,3)+HD_l:ijk(2,3)+HD_r, 1:num_eqn, hvy_id), &
                hvy_restricted(ijk(1,1)+HD_l:ijk(2,1)+HD_r, ijk(1,2)+HD_l:ijk(2,2)+HD_r, ijk(1,3)+HD_l:ijk(2,3)+HD_r, 1:num_eqn), params%HD, HD_l, HD_r, do_restriction=.true.)
        endif
        ! call abort(197, "Not finished yet but also how did you end up here?")
    endif


    call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
    lvl_me = lgt_block(lgt_ID, IDX_MESH_LVL)

    ! apply copying, do it for fine and coarse neighbors as we assume for coarse n independency and for fine n that they are not synched
    do k_n = 1, size(hvy_neighbor, 2)
        ! neighbor exists?
        lgt_ID_n = hvy_neighbor( hvy_ID, k_n )
        if ( lgt_ID_n /= -1 ) then
            lvl_diff = lvl_me - lgt_block(lgt_ID_n, IDX_MESH_LVL)
            if (lvl_diff == -1 .or. lvl_diff == +1) then
                call coarseExtensionManipulateSC_block(params, hvy_restricted(:,:,:,1:num_eqn), hvy_data(:,:,:,1:num_eqn, hvy_id), k_n, skip_ghosts=.true., ijk=ijk)
            endif
        endif
    end do
end subroutine



subroutine predict_data( params, pre_data, ijk, hvy_block, num_eqn, hvy_id )
    implicit none

    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(out)                      :: pre_data(:,:,:,:)
    !> ijk
    integer(kind=ik), intent(in)                    :: ijk(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)                    :: num_eqn      !< How many components? Needed as in between we use hvy_tmp
    integer(kind=ik), intent(in)                    :: hvy_id

    integer(kind=ik) :: dF, nx, ny, nz, nc

    ! data size
    nx = ijk(2,1) - ijk(1,1) + 1
    ny = ijk(2,2) - ijk(1,2) + 1
    nz = ijk(2,3) - ijk(1,3) + 1
    nc = num_eqn

    ! The neighbor is finer: we have to interpolate the data
    ! Notice how the indices are now in the beginning of the array, this was once decided whysoever (I hope I change it at some point)
    if ( params%dim == 3 ) then
    ! 3D
        do dF = 1, nc
            call prediction_3D( hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), dF, hvy_id ), &
                pre_data( 1:2*nx-1, 1:2*ny-1, 1:2*nz-1, dF), &
            params%order_predictor)
        end do
    else
    ! 2D
        do dF = 1, nc
            call prediction_2D( hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), 1, dF, hvy_id ), &
                pre_data( 1:2*nx-1, 1:2*ny-1, 1, dF),  params%order_predictor)
        end do
    end if
end subroutine predict_data
