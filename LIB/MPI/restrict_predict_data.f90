subroutine restrict_predict_data( params, res_pre_data, data_bounds2, neighborhood, &
    level_diff, data_bounds_type, hvy_block, hvy_id, data_bounds )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(out)                 :: res_pre_data(:,:,:,:)
    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds2(2,3)
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff
    ! data_bounds_type
    integer(kind=ik), intent(in)                   :: data_bounds_type
    !> heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, iN, jN, kN

    ! grid parameter
    integer(kind=ik)                                :: Bs, g, NdF
    integer(kind=ik)                                :: sh_start, sh_end

!---------------------------------------------------------------------------------------------
! interfaces
! HACK
 data_bounds = data_bounds2
!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs  = params%number_block_nodes
    g   = params%number_ghost_nodes
    ndf = params%number_data_fields

    ! data size
    iN = data_bounds(2,1) - data_bounds(1,1) + 1
    jN = data_bounds(2,2) - data_bounds(1,2) + 1
    kN = data_bounds(2,3) - data_bounds(1,3) + 1

    sh_start = 0
    sh_end   = 0

    if ( data_bounds_type == exclude_redundant ) then
        sh_start = 1
    end if
    if ( data_bounds_type == only_redundant ) then
        sh_end = -g
    end if

!---------------------------------------------------------------------------------------------
! main body

    if ( params%threeD_case ) then
        ! 3D
        select case(neighborhood)
            ! nothing to do
            ! '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___'
            ! '_12/___', '_13/___', '_14/___', '_15/___'
            ! '_62/___', '_63/___', '_64/___', '_65/___'
            ! '_23/___', '_25/___', '_43/___', '_45/___'
            case(1:18)

            ! '123/___', '134/___', '145/___', '152/___'
            ! '623/___', '634/___', '645/___', '652/___'
            ! '__1/123', '__1/134', '__1/145', '__1/152',
            ! '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
            ! '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
            ! '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
            ! '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
            ! '_43/134', '_43/634', '_45/145', '_45/645'
            case(19:74)
                if ( level_diff == -1 ) then
                    ! loop over all data fields
                    do dF = 1, NdF
                        ! interpolate data
                        call prediction_3D( hvy_block( data_bounds(1,1):data_bounds(2,1), &
                                                       data_bounds(1,2):data_bounds(2,2), &
                                                       data_bounds(1,3):data_bounds(2,3), dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1:kN*2-1, dF), params%order_predictor)
                    end do
                    ! reset data bounds
                    select case(neighborhood)
                        case(19) ! '123/___'
                            data_bounds(1,1) = g-1-sh_end
                            data_bounds(2,1) = 2*g-1-sh_start
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = g-1-sh_end
                            data_bounds(2,3) = 2*g-1-sh_start

                        case(20) ! '134/___'
                            data_bounds(1,1) = g-1-sh_end
                            data_bounds(2,1) = 2*g-1-sh_start
                            data_bounds(1,2) = g-1-sh_end
                            data_bounds(2,2) = 2*g-1-sh_start
                            data_bounds(1,3) = g-1-sh_end
                            data_bounds(2,3) = 2*g-1-sh_start

                        case(21) ! '145/___'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = g-1-sh_end
                            data_bounds(2,2) = 2*g-1-sh_start
                            data_bounds(1,3) = g-1-sh_end
                            data_bounds(2,3) = 2*g-1-sh_start

                        case(22) ! '152/___'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = g-1-sh_end
                            data_bounds(2,3) = 2*g-1-sh_start

                        case(23) ! '623/___'
                            data_bounds(1,1) = g-1-sh_end
                            data_bounds(2,1) = 2*g-1-sh_start
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(24) ! '634/___'
                            data_bounds(1,1) = g-1-sh_end
                            data_bounds(2,1) = 2*g-1-sh_start
                            data_bounds(1,2) = g-1-sh_end
                            data_bounds(2,2) = 2*g-1-sh_start
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(25) ! '645/___'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = g-1-sh_end
                            data_bounds(2,2) = 2*g-1-sh_start
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(26) ! '652/___'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(27) ! '__1/123'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(28) ! '__1/134'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(29) ! '__1/145'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(30) ! '__1/152'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(31) ! '__2/123'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(32) ! '___2/623'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(33) ! '__2/152'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(34) ! '__2/652'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(35) ! '__3/123'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(37) ! '__3/134'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(38) ! '__3/634'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(36) ! '__3/623'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(40) ! '__4/634'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(39) ! '__4/134'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(41) ! '__4/145'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(42) ! '__4/645'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(45) ! '__5/152'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(43) ! '__5/145'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(44) ! '__5/645'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(46) ! '__5/652'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(47) ! '__6/623'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(48) ! '__6/634'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(49) ! '__6/645'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(50) ! '__6/652'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(51) ! '_12/123'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(52) ! '_12/152'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(54) ! '_13/134'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(53) ! '_13/123'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(55) ! '_14/134'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(56) ! '_14/145'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(57) ! '_15/145'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(58) ! '_15/152'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = Bs+g-sh_end
                            data_bounds(2,3) = Bs+2*g-sh_start

                        case(59) ! '_62/623'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(60) ! '_62/652'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(62) ! '_63/634'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(61) ! '_63/623'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(63) ! '_64/634'
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+2*g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(64) ! '_64/645'
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(65) ! '_65/645'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+2*g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(66) ! '_65/652'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g
                            data_bounds(1,3) = 1+sh_start
                            data_bounds(2,3) = g+1+sh_end

                        case(67) ! '_23/123'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(68) ! '_23/236''
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(69) ! '_25/152'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(70) ! '_25/652'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = 1+sh_start
                            data_bounds(2,2) = g+1+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(71) ! '_43/134'
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(72) ! '_43/634''
                            data_bounds(1,1) = Bs+g-sh_end
                            data_bounds(2,1) = Bs+2*g-sh_start
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                        case(73) ! '_45/145'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+2*g

                        case(74) ! '_45/645'
                            data_bounds(1,1) = 1+sh_start
                            data_bounds(2,1) = g+1+sh_end
                            data_bounds(1,2) = Bs+g-sh_end
                            data_bounds(2,2) = Bs+2*g-sh_start
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = Bs+g

                    end select

                elseif ( level_diff == 1) then
                    ! loop over all data fields
                    do dF = 1, NdF
                        ! first dimension
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            ! second dimension
                            do j = data_bounds(1,2), data_bounds(2,2), 2
                                ! third dimension
                                do k = data_bounds(1,3), data_bounds(2,3), 2

                                    ! write restricted data
                                    res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, (k-data_bounds(1,3))/2+1, dF) &
                                    = hvy_block( i, j, k, dF, hvy_id )

                                end do
                            end do
                        end do
                    end do

                    select case(neighborhood)

                        case(19:26)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(27:30)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(31:34)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = (Bs+1)/2

                        case(35:38)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = (Bs+1)/2

                        case(39:42)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = (Bs+1)/2

                        case(43:46)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = (Bs+1)/2

                        case(47:50)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(51:52)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(53:54)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(55:56)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(57:58)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(59:60)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(61:62)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(63:64)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = (Bs+1)/2
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(65:66)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = (Bs+1)/2
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = g+1-sh_start+sh_end

                        case(67:68)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = (Bs+1)/2

                        case(69:74)
                            ! reset data bounds
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = g+1-sh_start+sh_end
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = g+1-sh_start+sh_end
                            data_bounds(1,3) = 1
                            data_bounds(2,3) = (Bs+1)/2

                    end select

                end if
        end select

    else
        ! 2D
        select case(neighborhood)

            ! nothing to do
            ! '__N' '__E' '__S' '__W'
            case(1,2,3,4)

            ! '_NE' '_NW' '_SE' '_SW'
            case(5,6,7,8)
                if ( level_diff == -1 ) then
                    ! loop over all data fields
                    do dF = 1, NdF
                        ! interpolate data
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1, dF), params%order_predictor)
                    end do
                    ! reset data bounds
                    select case(neighborhood)
                        ! '_NE'
                        case(5)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-2

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 2*g-1
                            end select
                        ! '_NW'
                        case(6)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_SE'
                        case(7)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-2
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-2

                                case(include_redundant)
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-1
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 2*g-1
                                    data_bounds(1:2,2) = 2*g-1
                            end select
                        ! '_SW'
                        case(8)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-2
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case(include_redundant)
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 2*g-1
                                    data_bounds(1:2,2) = 1
                            end select
                    end select

                elseif ( level_diff == 1) then
                    ! loop over all data fields
                    do dF = 1, NdF
                        ! first dimension
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            ! second dimension
                            do j = data_bounds(1,2), data_bounds(2,2), 2

                                ! write restricted data
                                res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, 1, dF) &
                                = hvy_block( i, j, 1, dF, hvy_id )

                            end do
                        end do
                    end do
                    ! reset data bounds
                    select case(neighborhood)
                        ! '_NE'
                        case(5)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_NW'
                        case(6)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_SE'
                        case(7)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_SW'
                        case(8)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                    end select

                end if

            ! 'NNE' 'NNW' 'SSE' 'SSW' ENE' 'ESE' 'WNW' 'WSW'
            case(9,10,11,12,13,14,15,16)
                if ( level_diff == -1 ) then
                    ! loop over all data fields
                    do dF = 1, NdF
                        ! interpolate data
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1, dF), params%order_predictor)
                    end do
                    ! reset data bounds
                    select case(neighborhood)
                        ! 'NNE'
                        case(9)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        ! 'NNW'
                        case(10)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(only_redundant)
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                            end select

                        ! 'SSE'
                        case(11)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case(include_redundant)
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case(only_redundant)
                                    data_bounds(1:2,1) = Bs+2*g
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        ! 'SSW'
                        case(12)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(include_redundant)
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(only_redundant)
                                    data_bounds(1:2,1) = Bs+2*g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                            end select

                        ! 'ENE'
                        case(13)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g-1

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g

                                case(only_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1:2,2) = Bs+2*g

                            end select

                        ! 'ESE'
                        case(14)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g-1

                                case(include_redundant)
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g

                                case(only_redundant)
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1:2,2) = Bs+2*g

                            end select

                        ! 'WNW'
                        case(15)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case(include_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1:2,2) = 1

                            end select

                        ! 'WSW'
                        case(16)
                            select case(data_bounds_type)
                                case(exclude_redundant)
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case(include_redundant)
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case(only_redundant)
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1:2,2) = 1

                            end select

                        end select

                elseif ( level_diff == 1 ) then
                    ! loop over all data fields
                    do dF = 1, NdF
                        ! first dimension
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            ! second dimension
                            do j = data_bounds(1,2), data_bounds(2,2), 2

                                ! write restricted data
                                res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, 1, dF) &
                                = hvy_block( i, j, 1, dF, hvy_id )

                            end do
                        end do
                    end do
                    ! reset data bounds
                    data_bounds(1,1:2) = 1
                    data_bounds(2,1)   = (iN+1)/2
                    data_bounds(2,2)   = (jN+1)/2

                end if

        end select
    end if
end subroutine restrict_predict_data
