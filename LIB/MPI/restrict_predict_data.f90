subroutine restrict_predict_data( params, res_pre_data, ijk, neighborhood, &
    level_diff, hvy_block, hvy_id )

    implicit none

    type (type_params), intent(in)                  :: params                   !> user defined parameter structure
    real(kind=rk), intent(out)                      :: res_pre_data(:,:,:,:)    !> data buffer
    integer(kind=ik), intent(in)                    :: ijk(2,3)                 !> indices in x,y,z direction of the ghost node patch
    integer(kind=ik), intent(in)                    :: neighborhood             !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: level_diff               !> difference between block levels
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :) !> heavy data array - block data
    integer(kind=ik), intent(in)                    :: hvy_id

    ! some neighborhoods are intrinsically on the same level (level_diff=0)
    ! and thus it makes no sense to call the up/downsampling routine for those
    if ( params%dim == 3 .and. (neighborhood<=18) ) call abort(323223,"this case shouldnt appear")
    if ( params%dim == 2 .and. (neighborhood<=4) ) call abort(323223,"this case shouldnt appear")

    if ( level_diff == -1 ) then
        ! The neighbor is finer: we have to predict the data
        call predict_data( params, res_pre_data, ijk, hvy_block, hvy_id )

    elseif ( level_diff == +1) then
        ! The neighbor is coarser: we have to downsample the data
        call restrict_data( params, res_pre_data, ijk, hvy_block, hvy_id )

    else
        call abort(123005, "Lord Vader, restrict_predict_data is called with leveldiff /= -+1")

    end if

end subroutine restrict_predict_data


subroutine restrict_data( params, res_data, ijk, hvy_block, hvy_id )
    implicit none

    type (type_params), intent(in)                  :: params                   !> user defined parameter structure
    real(kind=rk), intent(out)                      :: res_data(:,:,:,:)        !> data buffer
    integer(kind=ik), intent(in)                    :: ijk(2,3)
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :) !> heavy data array - block data
    integer(kind=ik), intent(in)                    :: hvy_id
    integer(kind=ik)                                :: ix, iy, iz, dF           ! local variables

    real(kind=rk), allocatable:: tmp_block(:,:,:)

    if (filter) then
        allocate(tmp_block(size(hvy_block,1), size(hvy_block,2), size(hvy_block,3)))

        do dF = 1, size(hvy_block,4)
            if (params%dim==2) then
                call restriction_prefilter_2D(hvy_block(:,:,1,dF,hvy_id), tmp_block(:,:,1), params%wavelet)
            else
                call restriction_prefilter_3D(hvy_block(:,:,:,dF,hvy_id), tmp_block(:,:,:), params%wavelet)
            endif

            do iz = ijk(1,3), ijk(2,3), 2
                do iy = ijk(1,2), ijk(2,2), 2
                    do ix = ijk(1,1), ijk(2,1), 2

                        ! write restricted (downsampled) data
                        res_data( (ix-ijk(1,1))/2+1, (iy-ijk(1,2))/2+1, (iz-ijk(1,3))/2+1, dF) &
                        = tmp_block( ix, iy, iz )

                    end do
                end do
            end do
        end do

        deallocate(tmp_block)
    else
        do dF = 1, size(hvy_block,4)
            do iz = ijk(1,3), ijk(2,3), 2
                do iy = ijk(1,2), ijk(2,2), 2
                    do ix = ijk(1,1), ijk(2,1), 2

                        ! write restricted (downsampled) data
                        res_data( (ix-ijk(1,1))/2+1, (iy-ijk(1,2))/2+1, (iz-ijk(1,3))/2+1, dF) &
                        = hvy_block( ix, iy, iz, dF, hvy_id )

                    end do
                end do
            end do
        end do
    endif

end subroutine restrict_data


subroutine predict_data( params, pre_data, ijk, hvy_block, hvy_id )
    implicit none

    type (type_params), intent(in)                  :: params                   !> user defined parameter structure
    real(kind=rk), intent(out)                      :: pre_data(:,:,:,:)        !> data buffer
    integer(kind=ik), intent(in)                    :: ijk(2,3)
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :) !> heavy data array - block data
    integer(kind=ik), intent(in)                    :: hvy_id
    integer(kind=ik)                                :: dF, nx, ny, nz           ! local variables
    integer(kind=ik)                                :: NdF


    NdF = size(hvy_block,4)

    ! data size
    nx = ijk(2,1) - ijk(1,1) + 1
    ny = ijk(2,2) - ijk(1,2) + 1
    nz = ijk(2,3) - ijk(1,3) + 1

    ! The neighbor is finer: we have to interpolate the data

    if ( params%dim == 3 ) then      ! 3D

        do dF = 1, NdF
            call prediction_3D( hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), &
            ijk(1,3):ijk(2,3), dF, hvy_id ), pre_data( 1:2*nx-1, 1:2*ny-1, 1:2*nz-1, dF), &
            params%order_predictor)
        end do

    else      ! 2D
        do dF = 1, NdF
            call prediction_2D( hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2),&
            1, dF, hvy_id ), pre_data( 1:2*nx-1, 1:2*ny-1, 1, dF),  params%order_predictor)
        end do

    end if !3d/2D
end subroutine predict_data
