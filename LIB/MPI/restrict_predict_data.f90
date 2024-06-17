subroutine restrict_predict_data( params, res_pre_data, ijk, neighborhood, &
    level_diff, hvy_block, num_eqn, hvy_id )

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

    integer(kind=ik) :: iy  ! debug variable

    ! some neighborhoods are intrinsically on the same level (level_diff=0)
    ! and thus it makes no sense to call the up/downsampling routine for those
    if ( params%dim == 3 .and. (neighborhood<=18) ) call abort(323223,"this case shouldnt appear")
    if ( params%dim == 2 .and. (neighborhood<=4) ) call abort(323223,"this case shouldnt appear")


    if ( level_diff == -1 ) then
        ! The neighbor is finer: we have to predict the data
        call predict_data( params, res_pre_data, ijk, hvy_block, num_eqn, hvy_id )

    elseif ( level_diff == +1) then
        ! The neighbor is coarser: we have to downsample the data
        call restrict_data( params, res_pre_data, ijk, hvy_block, num_eqn, hvy_id )

    else
        call abort(123005, "Lord Vader, restrict_predict_data is called with leveldiff /= -+1")

    end if

end subroutine restrict_predict_data

subroutine restrict_data( params, res_data, ijk, hvy_block, num_eqn, hvy_id )
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

    integer(kind=ik)                :: ix, iy, iz, dF, nc, nx, ny, nz

    nx = size(hvy_block,1)
    ny = size(hvy_block,2)
    nz = size(hvy_block,3)
    nc = num_eqn

#ifdef DEV
    if (.not. allocated(params%HD)) call abort(230301051, "Pirates! Maybe setup_wavelet was not called?")
#endif

    ! applying the filter is expensive, and we therefor apply it only once to the entire
    ! block. Result is stored in a large work array. Note we could try to optimize this
    ! by applying the filter only in the required patch.
    ! Pro: we maybe save a bi of CPU time, as we usually do not need the entire block filtered
   - - - - - - - -
!                      g g g g g g g g         - - - - - - - -
!   Block relations       Block in                Block out
! b = block, m = medium neighbor, c = coarse n, f = fine n, g = ghost point, i = interior point, f = filtered point, - = value should not be used
! Note: The bottom left value is copied (i) due to the finer neighbor in the bottom left corner, the top values due to the coarser top neighbor.
! Ghost points of output should not be used, as they may have been filtered correctly, partially or left unfiltered arbitrarily    ! Con: more work and maybe we compute some values twice (if patches overlap)
    if (hvy_ID /= restricted_hvy_ID) then
        call restrict_copy_at_CE(params, hvy_block, hvy_ID, nc)

        ! call blockFilterXYZ_interior_vct( params, hvy_block(:,:,:,1:nc,hvy_id), hvy_restricted(:,:,:,1:nc), params%HD, &
        ! lbound(params%HD, dim=1), ubound(params%HD, dim=1), params%g)

        ! if (hvy_id == 17 .and. params%rank == 1) THEN
        !     write(*, '("Restrict block R-0 BH-7 IJK-", 3(i0, "-"), i0)') ijk(:,1:2)
        !     do iy = 2,25, 2
        !         write(*, '(25(es9.2))') hvy_restricted(2:25:2, iy, 1, 1)
        !     enddo
        ! endif

        restricted_hvy_ID = hvy_ID
    endif

    do dF = 1, nc
        do iz = ijk(1,3), ijk(2,3), 2
            do iy = ijk(1,2), ijk(2,2), 2
                do ix = ijk(1,1), ijk(2,1), 2

                    ! write restricted (downsampled) data
                    res_data( (ix-ijk(1,1))/2+1, (iy-ijk(1,2))/2+1, (iz-ijk(1,3))/2+1, dF) &
                    = hvy_restricted( ix, iy, iz, dF)
                    ! res_data( (ix-ijk(1,1))/2+1, (iy-ijk(1,2))/2+1, (iz-ijk(1,3))/2+1, dF) &
                    ! = hvy_block( ix, iy, iz, dF, hvy_id )
                end do
            end do
        end do
    end do
end subroutine restrict_data



!> \brief Uses HD filter at all interfaces and copies values everywhere where it would use from coarser or finer neighbor
!> Used for sync in order to match values with coarse extension. Coarser AND finer neighbors cannot be filtered as there values
!! are not present in the ghost layers, so for all those patches where filters need those points the values will be copied.
!
!                      g g g g g g g g         - - - - - - - -
!                      g g g g g g g g         - - - - - - - -
!      m c c c         g g i i i i g g         - - i i i i - -
!      m b b m         g g i i i i g g         - - f f f f - -
!      m b b m         g g i i i i g g         - - f f f f - -
!      f m m m         g g i i i i g g         - - i f f f - -
!                      g g g g g g g g      
subroutine restrict_copy_at_CE(params, hvy_data, hvy_ID, num_eqn)

    implicit none
    type (type_params), intent(in) :: params  !< good ol params
    real(kind=rk), intent(inout)   :: hvy_data(:, :, :, :, :)  !< heavy data array - block data
    integer(kind=ik), intent(in)   :: hvy_id  !< which block to look at
    integer(kind=ik), intent(in)   :: num_eqn !< How many components? Needed as it could vary

    integer(kind=ik) :: k_n, lvl_me, lvl_diff, lgt_ID, lgt_ID_n, HD_l, HD_r

    HD_l = lbound(params%HD, dim=1)
    HD_r = ubound(params%HD, dim=1)

    ! filter whole block and don't care about neighbors first
    call blockFilterXYZ_vct( params, hvy_data(:,:,:,1:num_eqn, hvy_id), hvy_restricted(:,:,:,1:num_eqn), params%HD, HD_l, HD_r)


    call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
    lvl_me = lgt_block(lgt_ID, IDX_MESH_LVL)

    ! apply copying, do it for fine and coarse neighbors as we assume for coarse n independency and for fine n that they are not synched
    do k_n = 1, size(hvy_neighbor, 2)
        ! neighbor exists?
        lgt_ID_n = hvy_neighbor( hvy_ID, k_n )
        if ( lgt_ID_n /= -1 ) then
            lvl_diff = lvl_me - lgt_block(lgt_ID_n, IDX_MESH_LVL)
            if (lvl_diff == -1 .or. lvl_diff == +1) then
                call coarseExtensionManipulateSC_block(params, hvy_restricted(:,:,:,1:num_eqn), hvy_data(:,:,:,1:num_eqn, hvy_id), k_n, skip_ghosts=.true.)
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



! ! note this is a tree-level routine (as is everything in the ghost nodes module)
! ! hence, HVY_N is a scalar, not an array, and no tree_ID is passed
! subroutine fixInterpolatedPoints_postSync( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
!     implicit none
!     type (type_params), intent(in) :: params
!     !> light data array
!     integer(kind=ik), intent(in)   :: lgt_block(:, :)
!     !> heavy data array - block data
!     real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)
!     !> heavy data array - neighbor data
!     integer(kind=ik), intent(in)   :: hvy_neighbor(:,:)
!     !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)   :: hvy_active(:)
!     !> number of active blocks (heavy data)
!     integer(kind=ik), intent(in)   :: hvy_n
!     real(kind=rk), allocatable, save :: u_tmp(:,:,:,:)
!
!     integer(kind=ik) :: k, hvyID, lgtID, neighborhood, ix, iy, iz, g, Bs(1:3), level_diff
!     integer(kind=ik) :: n_minus=3, n_plus=2, ix2, iy2, iz2, level_me, level_neighbor, lgtID_neighbor
!     integer(kind=ik) :: ix_low, ix_high, iy_low, iy_high, iz_low, iz_high
!     integer(kind=ik) :: ny,nx,nz,nc, x1patch,y1patch,z1patch,x2patch,y2patch,z2patch, &
!     x1reinterpolate,y1reinterpolate,z1reinterpolate,x2reinterpolate,y2reinterpolate,z2reinterpolate
!
!     g = params%g
!     Bs = params%Bs
!     nx = size(hvy_block, 1)
!     ny = size(hvy_block, 2)
!     nz = size(hvy_block, 3)
!     nc = size(hvy_block, 4)
!
!     if (allocated(u_tmp)) then
!         if (.not. areArraysSameSize(u_tmp, hvy_block(:,:,:,:,1))) deallocate(u_tmp)
!     endif
!     if (.not. allocated(u_tmp)) then
!         allocate(u_tmp(1:size(hvy_block,1), 1:size(hvy_block,2), 1:size(hvy_block,3), 1:size(hvy_block,4)))
!     endif
!
!     do k = 1, hvy_n
!         hvyID = hvy_active(k)
!         call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! FIXME: this can be done more efficiently then Predict(Restrict()) the entire block
! u_tmp = 0.0_rk
! ! loose downsampling ( because we deal with the ghost nodes here, there is no HD filter)
! if (modulo(g, 2) == 0) then
!     u_tmp(1:nx:2,1:ny:2,1:nz:2, 1:nc) = hvy_block(1:nx:2,1:ny:2,1:nz:2, 1:nc, hvyID)
! else
!     u_tmp(2:nx-1:2,2:ny-1:2,2:nz-1:2, 1:nc) = hvy_block(2:nx-1:2,2:ny-1:2,2:nz-1:2, 1:nc, hvyID)
! endif
! ! interpolation
! call blockFilterXYZ_wherePossible_vct(params, u_tmp, u_tmp, params%HR, lbound(params%HR, dim=1), ubound(params%HR, dim=1))
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!         do neighborhood = 1, size(hvy_neighbor,2)
!             ! neighbor exists ?
!             if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
!                 ! neighbor light data id
!                 lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
!                 level_me       = lgt_block( lgtID,          IDX_MESH_LVL )
!                 level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )
!
!                 ! this routine corrects for mistakes made on the coarser block when it interpolates
!                 ! the values for its finer neighbor, on the finer one. Hence, it is active if the
!                 ! neighbor is coarser.
!                 if (level_neighbor < level_me) then
!
!                     level_diff = level_neighbor - level_me
!
!                     ! ijkPatches([start,end], [dir (x,y,z)], [neighborhood], [level-diff], [sender/receiver/up-downsampled])
!                     ! leveldiff = J_sender - J_recver
!                     ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
!                     ! leveldiff =  0 : sender is same level as recver
!                     ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
! ! because the ghost nodes module loops over sender, and here we want the receiver, we need to take the inverse neighborhood
!                     x1patch = ijkPatches(1, 1, inverse_neighbor(neighborhood,dim), -1, 2)
!                     y1patch = ijkPatches(1, 2, inverse_neighbor(neighborhood,dim), -1, 2)
!                     z1patch = ijkPatches(1, 3, inverse_neighbor(neighborhood,dim), -1, 2)
!                     x2patch = ijkPatches(2, 1, inverse_neighbor(neighborhood,dim), -1, 2)
!                     y2patch = ijkPatches(2, 2, inverse_neighbor(neighborhood,dim), -1, 2)
!                     z2patch = ijkPatches(2, 3, inverse_neighbor(neighborhood,dim), -1, 2)
!
! ! hvy_block( x1patch:x2patch, y1patch:y2patch, z1patch:z2patch ,:, hvyID) = 0.0_rk
! ! endif
!                     x1reinterpolate = (g+1) - n_minus
!                     y1reinterpolate = (g+1) - n_minus
!                     z1reinterpolate = (g+1) - n_minus
!                     x2reinterpolate = Bs(1)+g+n_plus
!                     y2reinterpolate = Bs(2)+g+n_plus
!                     z2reinterpolate = Bs(3)+g+n_plus
!
!                     ! if ((x1reinterpolate > x1patch) .and. (x1reinterpolate < x2patch)) then
!                     !     ix_low = x1reinterpolate
!                     ! else
!                     !     ix_low = x1patch
!                     ! endif
!                     !
!                     ! if ((x2reinterpolate > x1patch) .and. (x2reinterpolate < x2patch)) then
!                     !     ix_high = x2reinterpolate
!                     ! else
!                     !     ix_high = x2patch
!                     ! endif
!                     !
!                     ! if ((y1reinterpolate > y1patch) .and. (y1reinterpolate < y2patch)) then
!                     !     iy_low = y1reinterpolate
!                     ! else
!                     !     iy_low = y1patch
!                     ! endif
!                     !
!                     ! if ((y2reinterpolate > y1patch) .and. (y2reinterpolate < y2patch)) then
!                     !     iy_high = y2reinterpolate
!                     ! else
!                     !     iy_high = y2patch
!                     ! endif
!                     !
!                     ! if ((z1reinterpolate > z1patch) .and. (z1reinterpolate < z2patch)) then
!                     !     iz_low = z1reinterpolate
!                     ! else
!                     !     iz_low = z1patch
!                     ! endif
!                     !
!                     ! if ((z2reinterpolate > z1patch) .and. (z2reinterpolate < z2patch)) then
!                     !     iz_high = z2reinterpolate
!                     ! else
!                     !     iz_high = z2patch
!                     ! endif
!
!                     ! the loop runs over the intersection of the two rectangles: the ghost nodes patch
!                     ! and the area for re-interpolation of the points (which have been wrongly interpolated
!                     ! because the coarse blocks ghost nodes had not been filled yet)
!                     ! 1st pass: interpolate in x-direction
!                     ! if (params%dim == 3) then
!                     !     do iz = iz_low, iz_high
!                     !         do iy = iy_low, iy_high
!                     !             do ix = ix_low, ix_high
!                     !                 ! set origin to first interior point
!                     !                 ix2 = ix -(g+1)
!                     !                 iy2 = iy -(g+1)
!                     !                 iz2 = iz -(g+1)
!                     !
!                     !                 ! now even points lie on the coarse grid
!                     !                 !     odd  points lie on the fine grid (and have thus been interpolated on the coarse block)
!                     !                 if ( .not. ( (modulo(ix2,2)==0).and.(modulo(iy2,2)==0).and.(modulo(iz2,2)==0) )) then
!                     !                     ! copy interpolated values
!                     !                     hvy_block(ix,iy,iz,1:nc,hvyID) = 0.0_rk!u_tmp(ix,iy,iz,1:nc)
!                     !                 endif
!                     !             enddo
!                     !         enddo
!                     !     enddo
!                     ! else
!
!
!
!                     do iy = 1, Bs(2)+2*g
!                         do ix = 1, Bs(2)+2*g
!                             if ((ix>=x1reinterpolate).and.(ix<=x2reinterpolate).and.(iy>=y1reinterpolate).and.(iy<=y2reinterpolate)) then
!                                 if ((ix>=x1patch).and.(ix<=x2patch).and.(iy>=y1patch).and.(iy<=y2patch)) then
!                                     ! set origin to first interior point
!                                     ix2 = ix -(g+1)
!                                     iy2 = iy -(g+1)
!
!                                     ! now even points lie on the coarse grid
!                                     !     odd  points lie on the fine grid (and have thus been interpolated on the coarse block)
!                                     ! if ( .not. ( (modulo(ix2,2)==0).and.(modulo(iy2,2)==0) )) then
!                                         ! copy interpolated values
!                                         hvy_block(ix,iy,:,1:nc,hvyID) = 2.0_rk!u_tmp(ix,iy,1,1:nc)
!                                         ! write(*,*) "yes"
!                                     ! endif
!                                 endif
!                             endif
!                         enddo
!                     enddo
!                     do iy = 1, Bs(2)+2*g
!                         do ix = 1, Bs(2)+2*g
!                             ! set origin to first interior point
!                             ix2 = ix -(g+1)
!                             iy2 = iy -(g+1)
!
!                             ! now even points lie on the coarse grid
!                             !     odd  points lie on the fine grid (and have thus been interpolated on the coarse block)
!                             if ( ( (modulo(ix2,2)==0).and.(modulo(iy2,2)==0) )) then
!                                 hvy_block(ix,iy,:,1:nc,hvyID) = 5.55_rk!
!                             endif
!                         enddo
!                     enddo
!
!                     !     do iy = iy_low, iy_high
!                     !         do ix = ix_low, ix_high
!                     !             ! set origin to first interior point
!                     !             ix2 = ix -(g+1)
!                     !             iy2 = iy -(g+1)
!                     !
!                     !             ! now even points lie on the coarse grid
!                     !             !     odd  points lie on the fine grid (and have thus been interpolated on the coarse block)
!                     !             if ( .not. ( (modulo(ix2,2)==0).and.(modulo(iy2,2)==0) )) then
!                     !                 ! copy interpolated values
!                     !                 hvy_block(ix,iy,1,1:nc,hvyID) = 0.0_rk!u_tmp(ix,iy,1,1:nc)
!                     !             endif
!                     !         enddo
!                     !     enddo
!                     ! endif
!
!
!                 endif
!             endif
!         enddo
!     enddo
!
! end subroutine
