subroutine restrict_predict_data( params, res_pre_data, data_bounds, neighborhood, &
    level_diff, hvy_block, hvy_id )

    !---------------------------------------------------------------------------------------------
    ! modules

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(out)                      :: res_pre_data(:,:,:,:)
    !> data_bounds
    integer(kind=ik), intent(in)                    :: data_bounds(2,3)
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id

    ! local variables
    integer(kind=ik)                                :: i, j, k, dF, iN, jN, kN
    integer(kind=ik)                                :: NdF

!---------------------------------------------------------------------------------------------
! variables initialization

    NdF = params%number_data_fields

    ! data size
    iN = data_bounds(2,1) - data_bounds(1,1) + 1
    jN = data_bounds(2,2) - data_bounds(1,2) + 1
    kN = data_bounds(2,3) - data_bounds(1,3) + 1

!---------------------------------------------------------------------------------------------
! main body

    if ( params%threeD_case ) then
        !-----------------------------------------------------------------------
        ! 3D
        !-----------------------------------------------------------------------
        select case(neighborhood)
            case(1:18)
                ! nothing to do: those cases cannot occur!
                ! '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___'
                ! '_12/___', '_13/___', '_14/___', '_15/___'
                ! '_62/___', '_63/___', '_64/___', '_65/___'
                ! '_23/___', '_25/___', '_43/___', '_45/___'

            case(19:74)
                ! '123/___', '134/___', '145/___', '152/___'
                ! '623/___', '634/___', '645/___', '652/___'
                ! '__1/123', '__1/134', '__1/145', '__1/152',
                ! '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
                ! '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
                ! '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
                ! '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
                ! '_43/134', '_43/634', '_45/145', '_45/645'
                if ( level_diff == -1 ) then
                    ! The neighbor is finer: we have to interpolate the data
                    do dF = 1, NdF
                        ! interpolate data
                        call prediction_3D( hvy_block( data_bounds(1,1):data_bounds(2,1), &
                        data_bounds(1,2):data_bounds(2,2), &
                        data_bounds(1,3):data_bounds(2,3), dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1:kN*2-1, dF), params%order_predictor)
                    end do

                elseif ( level_diff == 1) then
                    ! The neighbor is coarser: we have to downsample the data
                    do dF = 1, NdF
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            do j = data_bounds(1,2), data_bounds(2,2), 2
                                do k = data_bounds(1,3), data_bounds(2,3), 2

                                    ! write restricted data
                                    res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, (k-data_bounds(1,3))/2+1, dF) &
                                    = hvy_block( i, j, k, dF, hvy_id )

                                end do
                            end do
                        end do
                    end do

                end if
            end select

    else
        !-----------------------------------------------------------------------
        ! 2D
        !-----------------------------------------------------------------------
        select case(neighborhood)
            case(1,2,3,4)
                ! nothing to do: those cases cannot occur!
                ! '__N' '__E' '__S' '__W'

            case(5,6,7,8)
                ! '_NE' '_NW' '_SE' '_SW'
                if ( level_diff == -1 ) then
                    ! The neighbor is finer: we have to interpolate the data
                    do dF = 1, NdF
                        ! interpolate data
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1, dF), params%order_predictor)
                    end do

                elseif ( level_diff == +1) then
                    ! The neighbor is coarser: we have to downsample the data
                    do dF = 1, NdF
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            do j = data_bounds(1,2), data_bounds(2,2), 2
                                ! write restricted data
                                res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, 1, dF) &
                                = hvy_block( i, j, 1, dF, hvy_id )
                            end do
                        end do
                    end do

                end if

            case(9,10, 11,12, 13,14, 15,16)
                ! 'NNE' 'NNW' 'SSE' 'SSW' ENE' 'ESE' 'WNW' 'WSW'
                if ( level_diff == -1 ) then
                    ! The neighbor is finer: we have to interpolate the data
                    do dF = 1, NdF
                        ! interpolate data
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1, dF), params%order_predictor)
                    end do

                elseif ( level_diff == 1 ) then
                    ! The neighbor is coarser: we have to downsample the data
                    do dF = 1, NdF
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            do j = data_bounds(1,2), data_bounds(2,2), 2
                                ! write restricted data
                                res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, 1, dF) &
                                = hvy_block( i, j, 1, dF, hvy_id )
                            end do
                        end do
                    end do
                end if
        end select ! neighborhood
    end if !3d/2D
end subroutine restrict_predict_data
