subroutine calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, sender_or_receiver)

    !---------------------------------------------------------------------------------------------
    ! modules

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff

    ! data_bounds_type
    integer(kind=ik), intent(in)                   :: data_bounds_type
    ! sender or reciver
    character(len=*), intent(in)                   :: sender_or_receiver

    integer(kind=ik) :: Bs, g, N, M
    integer(kind=ik) :: sender_interp_short_down1, sender_interp_short_down2, sender_interp_short_up1, sender_interp_short_up2
    integer(kind=ik) :: sender_interp_long_up1, sender_interp_long_up2, sender_interp_long_down1, sender_interp_long_down2
    integer(kind=ik) :: sh_start, sh_end

    !---------------------------------------------------------------------------------------------

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes
    ! midpoint of a block:
    M     = (Bs-1)/2 + g + 1

    ! the number S is how many extra coarse points on the sender side you use to
    ! avoid one-sided interpolation. The actual formula is S = (order-2 )/2 so for
    ! 6th order you would require 2 extra ones.
    ! NOTE: you can set it to S=0 and then the code uses one-sided interpolation stencils.
    ! NOTE: S is symmetric, we add the layer on all sides.
    if (params%order_predictor == "multiresolution_4th" ) then
        S  = 1
    elseif (params%order_predictor == "multiresolution_2nd" ) then
        S  = 0
    else
        call abort(2875490, "The predictor method is unknown")
    endif

    ! asymmetric shift. If the interpolation domain becomes too small for the stencil,
    ! then we can also extend it to the interior of the block, while still
    ! including only the redudant lines. The shift A (asymmetric shift) does that
    ! COARSE grid points (sender side, red area in python figs)
    A = 0

    ! these are derived indices useful for selection of interpolation slab on
    ! sender side. This means they are COARSE. NOTE: I derived the bounds for one additional safety point (S=1)
    ! but forgot to introduce the shift S. So I did that afterwards, hence the ...1 data_bounds
    ! have +1-S and the uper ones ...2 have -1+S
    sender_interp_short_down1 = g +1 -S
    sender_interp_short_down2 = g + (g+6)/2 -1 -1 +S + A
    sender_interp_short_up1   = Bs+g+1-(g+6)/2+1  +1 -S -A
    sender_interp_short_up2   = Bs+g+1-1 +S
    sender_interp_long_up1    = M - (g/2+1)  +1 -S -A
    sender_interp_long_up2    = Bs+g+1-1 +S
    sender_interp_long_down1  = g  +1 -S
    sender_interp_long_down2  = M + g/2 + 1-1 +S +A

    if (sender_interp_short_down2-sender_interp_short_down1+1 < 4 .and. params%order_predictor == "multiresolution_4th" ) then
        ! it can happen that for g=4 and S=0, the 4th order predictor selects a domain which is so small
        ! that the one-sided stencil does not have enough points.
        write(*,*) sender_interp_short_down2-sender_interp_short_down1+1
        call abort(2875491,"The ghost nodes patch for interpolation is too small for 4th order stencil.")
    endif

    sh_start = 0
    sh_end   = 0

    if ( data_bounds_type == exclude_redundant ) then
        sh_start = 1
    end if
    if ( data_bounds_type == only_redundant ) then
        sh_end = -g
    end if

    ! set 1 and not -1 (or anything else), because 2D bounds ignore 3rd dimension
    ! and thus cycle from 1:1
    data_bounds(:,:) = 1

    !---------------------------------------------------------------------------------------------

    select case(sender_or_receiver)
    case('sender')
        if ( params%threeD_case ) then
            !---3D------3D------3D------3D------3D------3D------3D------3D---
            select case(neighborhood)
            case(1)
                ! '__1/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = Bs-sh_end
                data_bounds(2,3) = Bs+g-sh_start

            case(2)
                ! '__2/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1+sh_start
                data_bounds(2,2) = g+1+g+sh_end
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(3)
                ! '__3/___'
                data_bounds(1,1) = Bs-sh_end
                data_bounds(2,1) = Bs+g-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(4)
                ! '__4/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs-sh_end
                data_bounds(2,2) = Bs+g-sh_start
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(5)
                ! '__5/___'
                data_bounds(1,1) = g+1+sh_start
                data_bounds(2,1) = g+1+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(6)
                ! '__6/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = g+1+sh_start
                data_bounds(2,3) = g+1+g+sh_end

            case(7)
                ! '_12/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1+sh_start
                data_bounds(2,2) = g+1+g+sh_end
                data_bounds(1,3) = Bs-sh_end
                data_bounds(2,3) = Bs+g-sh_start

            case(8)
                ! '_13/___'
                data_bounds(1,1) = Bs-sh_end
                data_bounds(2,1) = Bs+g-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = Bs-sh_end
                data_bounds(2,3) = Bs+g-sh_start

            case(9)
                ! '_14/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs-sh_end
                data_bounds(2,2) = Bs+g-sh_start
                data_bounds(1,3) = Bs-sh_end
                data_bounds(2,3) = Bs+g-sh_start

            case(10)
                ! '_15/___'
                data_bounds(1,1) = g+1+sh_start
                data_bounds(2,1) = g+1+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = Bs-sh_end
                data_bounds(2,3) = Bs+g-sh_start

            case(11)
                ! '_62/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1+sh_start
                data_bounds(2,2) = g+1+g+sh_end
                data_bounds(1,3) = g+1+sh_start
                data_bounds(2,3) = g+1+g+sh_end

            case(12)
                ! '_63/___'
                data_bounds(1,1) = Bs-sh_end
                data_bounds(2,1) = Bs+g-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = g+1+sh_start
                data_bounds(2,3) = g+1+g+sh_end

            case(13)
                ! '_64/___'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs-sh_end
                data_bounds(2,2) = Bs+g-sh_start
                data_bounds(1,3) = g+1+sh_start
                data_bounds(2,3) = g+1+g+sh_end

            case(14)
                ! '_65/___'
                data_bounds(1,1) = g+1+sh_start
                data_bounds(2,1) = g+1+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = g+1+sh_start
                data_bounds(2,3) = g+1+g+sh_end

            case(15)
                ! '_23/___'
                data_bounds(1,1) = Bs-sh_end
                data_bounds(2,1) = Bs+g-sh_start
                data_bounds(1,2) = g+1+sh_start
                data_bounds(2,2) = g+1+g+sh_end
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(16)
                ! '_25/___'
                data_bounds(1,1) = g+1+sh_start
                data_bounds(2,1) = g+1+g+sh_end
                data_bounds(1,2) = g+1+sh_start
                data_bounds(2,2) = g+1+g+sh_end
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(17)
                ! '_43/___'
                data_bounds(1,1) = Bs-sh_end
                data_bounds(2,1) = Bs+g-sh_start
                data_bounds(1,2) = Bs-sh_end
                data_bounds(2,2) = Bs+g-sh_start
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(18)
                ! '_45/___'
                data_bounds(1,1) = g+1+sh_start
                data_bounds(2,1) = g+1+g+sh_end
                data_bounds(1,2) = Bs-sh_end
                data_bounds(2,2) = Bs+g-sh_start
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(19,20,21,22)
                if ( level_diff == 0 ) then
                    data_bounds(1,3) = Bs-sh_end
                    data_bounds(2,3) = Bs+g-sh_start
                    select case(neighborhood)
                    case(19) ! '123/___'
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end

                    case(20) ! '134/___'
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start

                    case(21) ! '145/___'
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start

                    case(22) ! '152/___'
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end

                    end select

                elseif ( level_diff == -1 ) then
                    data_bounds(1,3) = Bs+1
                    data_bounds(2,3) = Bs+g
                    select case(neighborhood)
                    case(19) ! '123/___'
                        data_bounds(1,1) = Bs+1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+g

                    case(20) ! '134/___'
                        data_bounds(1,1) = Bs+1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = Bs+1
                        data_bounds(2,2) = Bs+g

                    case(21) ! '145/___'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+g
                        data_bounds(1,2) = Bs+1
                        data_bounds(2,2) = Bs+g

                    case(22) ! '152/___'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,3) = Bs-g-sh_end*2
                    data_bounds(2,3) = Bs+g-sh_start*2
                    select case(neighborhood)
                    case(19) ! '123/___'
                        data_bounds(1,1) = Bs-g-sh_end*2
                        data_bounds(2,1) = Bs+g-sh_start*2
                        data_bounds(1,2) = g+1+sh_start*2
                        data_bounds(2,2) = g+1+g+g+sh_end*2

                    case(20) ! '134/___'
                        data_bounds(1,1) = Bs-g-sh_end*2
                        data_bounds(2,1) = Bs+g-sh_start*2
                        data_bounds(1,2) = Bs-g-sh_end*2
                        data_bounds(2,2) = Bs+g-sh_start*2

                    case(21) ! '145/___'
                        data_bounds(1,1) = g+1+sh_start*2
                        data_bounds(2,1) = g+1+g+g+sh_end*2
                        data_bounds(1,2) = Bs-g-sh_end*2
                        data_bounds(2,2) = Bs+g-sh_start*2

                    case(22) ! '152/___'
                        data_bounds(1,1) = g+1+sh_start*2
                        data_bounds(2,1) = g+1+g+g+sh_end*2
                        data_bounds(1,2) = g+1+sh_start*2
                        data_bounds(2,2) = g+1+g+g+sh_end*2

                    end select
                end if

            case(23,24,25,26)
                if ( level_diff == 0 ) then
                    data_bounds(1,3) = g+1+sh_start
                    data_bounds(2,3) = g+1+g+sh_end
                    select case(neighborhood)
                    case(23) ! '623/___'
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end

                    case(24) ! '634/___'
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start

                    case(25) ! '645/___'
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start

                    case(26) ! '652/___'
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end

                    end select

                elseif ( level_diff == -1 ) then
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+g
                    select case(neighborhood)
                    case(23) ! '623/___'
                        data_bounds(1,1) = Bs+1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+g

                    case(24) ! '634/___'
                        data_bounds(1,1) = Bs+1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = Bs+1
                        data_bounds(2,2) = Bs+g

                    case(25) ! '645/___'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+g
                        data_bounds(1,2) = Bs+1
                        data_bounds(2,2) = Bs+g

                    case(26) ! '652/___'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,3) = g+1+sh_start*2
                    data_bounds(2,3) = g+1+g+g+sh_end*2
                    select case(neighborhood)
                    case(23) ! '623/___'
                        data_bounds(1,1) = Bs-g-sh_end*2
                        data_bounds(2,1) = Bs+g-sh_start*2
                        data_bounds(1,2) = g+1+sh_start*2
                        data_bounds(2,2) = g+1+g+g+sh_end*2

                    case(24) ! '634/___'
                        data_bounds(1,1) = Bs-g-sh_end*2
                        data_bounds(2,1) = Bs+g-sh_start*2
                        data_bounds(1,2) = Bs-g-sh_end*2
                        data_bounds(2,2) = Bs+g-sh_start*2

                    case(25) ! '645/___'
                        data_bounds(1,1) = g+1+sh_start*2
                        data_bounds(2,1) = g+1+g+g+sh_end*2
                        data_bounds(1,2) = Bs-g-sh_end*2
                        data_bounds(2,2) = Bs+g-sh_start*2

                    case(26) ! '652/___'
                        data_bounds(1,1) = g+1+sh_start*2
                        data_bounds(2,1) = g+1+g+g+sh_end*2
                        data_bounds(1,2) = g+1+sh_start*2
                        data_bounds(2,2) = g+1+g+g+sh_end*2

                    end select
                end if

            case(27,28,29,30)
                if ( level_diff == -1 ) then
                    data_bounds(1,3) = (Bs+1)/2
                    data_bounds(2,3) = Bs+g
                    select case(neighborhood)
                    case(27) ! '__1/123'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g

                    case(28) ! '__1/134'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(29) ! '__1/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(30) ! '__1/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,3) = Bs-g-sh_end*2
                    data_bounds(2,3) = Bs+g-sh_start*2
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g

                end if

            case(31,32,33,34)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = (Bs+1)/2+g+g
                    ! first, third dimension
                    select case(neighborhood)
                    case(32) ! '__2/623'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    case(31) ! '__2/123'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(33) ! '__2/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(34) ! '__2/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            case(35,36,37,38)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = (Bs+1)/2
                    data_bounds(2,1) = Bs+g
                    ! second, third dimension
                    select case(neighborhood)
                    case(35) ! '__3/123'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(37) ! '__3/134'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(38) ! '__3/634'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    case(36) ! '__3/623'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            case(39,40,41,42)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = (Bs+1)/2
                    data_bounds(2,2) = Bs+g
                    ! first, third dimension
                    select case(neighborhood)
                    case(40) ! '__4/634'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    case(39) ! '__4/134'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(41) ! '__4/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(42) ! '__4/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            case(43,44,45,46)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = (Bs+1)/2+g+g
                    ! second, third dimension
                    select case(neighborhood)
                    case(45) ! '__5/152'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(43) ! '__5/145'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(44) ! '__5/645'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    case(46) ! '__5/652'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            case(47,48,49,50)
                if ( level_diff == -1 ) then
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = (Bs+1)/2+g+g
                    select case(neighborhood)
                    case(47) ! '__6/623'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g

                    case(48) ! '__6/634'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(49) ! '__6/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(50) ! '__6/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,3) = g+1+sh_start*2
                    data_bounds(2,3) = g+1+g+g+sh_end*2
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g

                end if

            case(51,52)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = (Bs+1)/2+g+g
                    data_bounds(1,3) = (Bs+1)/2
                    data_bounds(2,3) = Bs+g
                    select case(neighborhood)
                    case(51) ! '_12/123'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(52) ! '_12/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2
                    data_bounds(1,3) = Bs-g-sh_end*2
                    data_bounds(2,3) = Bs+g-sh_start*2

                end if

            case(53,54)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = (Bs+1)/2
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,3) = (Bs+1)/2
                    data_bounds(2,3) = Bs+g
                    select case(neighborhood)
                    case(54) ! '_13/134'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(53) ! '_13/123'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = Bs-g-sh_end*2
                    data_bounds(2,3) = Bs+g-sh_start*2

                end if

            case(55,56)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = (Bs+1)/2
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = (Bs+1)/2
                    data_bounds(2,3) = Bs+g
                    select case(neighborhood)
                    case(55) ! '_14/134'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(56) ! '_14/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2
                    data_bounds(1,3) = Bs-g-sh_end*2
                    data_bounds(2,3) = Bs+g-sh_start*2

                end if

            case(57,58)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = (Bs+1)/2+g+g
                    data_bounds(1,3) = (Bs+1)/2
                    data_bounds(2,3) = Bs+g
                    select case(neighborhood)
                    case(57) ! '_15/145'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(58) ! '_15/152''
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = Bs-g-sh_end*2
                    data_bounds(2,3) = Bs+g-sh_start*2

                end if

            case(59,60)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = (Bs+1)/2+g+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = (Bs+1)/2+g+g
                    select case(neighborhood)
                    case(59) ! '_62/623'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(60) ! '_62/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2
                    data_bounds(1,3) = g+1+sh_start*2
                    data_bounds(2,3) = g+1+g+g+sh_end*2

                end if

            case(61,62)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = (Bs+1)/2
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = (Bs+1)/2+g+g
                    select case(neighborhood)
                    case(62) ! '_63/634'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(61) ! '_63/623'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = g+1+sh_start*2
                    data_bounds(2,3) = g+1+g+g+sh_end*2

                end if

            case(63,64)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = (Bs+1)/2
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = (Bs+1)/2+g+g
                    select case(neighborhood)
                    case(63) ! '_64/634'
                        data_bounds(1,1) = (Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(64) ! '_64/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2
                    data_bounds(1,3) = g+1+sh_start*2
                    data_bounds(2,3) = g+1+g+g+sh_end*2

                end if

            case(65,66)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = (Bs+1)/2+g+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = (Bs+1)/2+g+g
                    select case(neighborhood)
                    case(65) ! '_65/645'
                        data_bounds(1,2) = (Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(66) ! '_65/652'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = (Bs+1)/2+g+g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g
                    data_bounds(1,3) = g+1+sh_start*2
                    data_bounds(2,3) = g+1+g+g+sh_end*2

                end if

            case(67,68)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = (Bs+1)/2
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = (Bs+1)/2+g+g
                    select case(neighborhood)
                    case(67) ! '_23/123'
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(68) ! '_23/236''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            case(69,70)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = (Bs+1)/2+g+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = (Bs+1)/2+g+g
                    select case(neighborhood)
                    case(69) ! '_25/152'
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(70) ! '_25/652''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            case(71,72)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = (Bs+1)/2
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = (Bs+1)/2
                    data_bounds(2,2) = Bs+g
                    select case(neighborhood)
                    case(71) ! '_43/134'
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(72) ! '_43/634''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            case(73,74)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = (Bs+1)/2+g+g
                    data_bounds(1,2) = (Bs+1)/2
                    data_bounds(2,2) = Bs+g
                    select case(neighborhood)
                    case(73) ! '_45/145'
                        data_bounds(1,3) = (Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(74) ! '_45/645'
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = (Bs+1)/2+g+g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs+g

                end if

            end select

        else
            !---2D------2D------2D------2D------2D------2D------2D------2D---
            select case(neighborhood)
            case(1)
                ! '__N'
                data_bounds(1,1) = g+1+sh_start
                data_bounds(2,1) = g+1+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g

            case(2)
                ! '__E'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs-sh_end
                data_bounds(2,2) = Bs+g-sh_start

            case(3)
                ! '__S'
                data_bounds(1,1) = Bs-sh_end
                data_bounds(2,1) = Bs+g-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g

            case(4)
                ! '__W'
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1+sh_start
                data_bounds(2,2) = g+1+g+sh_end

            case(5)
                ! '_NE'
                if ( level_diff == 0 ) then
                    data_bounds(1,1) = g+1+sh_start
                    data_bounds(2,1) = g+1+g+sh_end
                    data_bounds(1,2) = Bs-sh_end
                    data_bounds(2,2) = Bs+g-sh_start

                elseif ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_down1, sender_interp_short_down2 /)
                    data_bounds(1:2,2) = (/ sender_interp_short_up1, sender_interp_short_up2 /)

                elseif ( level_diff == 1) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2

                end if

            case(6)
                ! '_NW'
                if ( level_diff == 0 ) then
                    data_bounds(1,1) = g+1+sh_start
                    data_bounds(2,1) = g+1+g+sh_end
                    data_bounds(1,2) = g+1+sh_start
                    data_bounds(2,2) = g+1+g+sh_end

                elseif ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_down1, sender_interp_short_down2 /)
                    data_bounds(1:2,2) = (/ sender_interp_short_down1, sender_interp_short_down2 /)

                elseif ( level_diff == 1) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2

                end if

            case(7)
                ! '_SE'
                if ( level_diff == 0 ) then
                    data_bounds(1,1) = Bs-sh_end
                    data_bounds(2,1) = Bs+g-sh_start
                    data_bounds(1,2) = Bs-sh_end
                    data_bounds(2,2) = Bs+g-sh_start

                elseif ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_up1, sender_interp_short_up2 /)
                    data_bounds(1:2,2) = (/ sender_interp_short_up1, sender_interp_short_up2 /)

                elseif ( level_diff == 1) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2

                end if

            case(8)
                ! '_SW'
                if ( level_diff == 0 ) then
                    data_bounds(1,1) = Bs-sh_end
                    data_bounds(2,1) = Bs+g-sh_start
                    data_bounds(1,2) = g+1+sh_start
                    data_bounds(2,2) = g+1+g+sh_end

                elseif ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_up1, sender_interp_short_up2 /)
                    data_bounds(1:2,2) = (/ sender_interp_short_down1, sender_interp_short_down2 /)

                elseif ( level_diff == 1) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2

                end if

            case(9)
                ! 'NNE'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_down1, sender_interp_short_down2 /)
                    data_bounds(1:2,2) = (/ sender_interp_long_up1, sender_interp_long_up2 /)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g

                end if

            case(10)
                ! 'NNW'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_down1, sender_interp_short_down2 /)
                    data_bounds(1:2,2) = (/ sender_interp_long_down1, sender_interp_long_down2 /)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1+sh_start*2
                    data_bounds(2,1) = g+1+g+g+sh_end*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g

                end if

            case(11)
                ! 'SSE'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_up1, sender_interp_short_up2 /)
                    data_bounds(1:2,2) = (/ sender_interp_long_up1, sender_interp_long_up2/)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g

                end if

            case(12)
                ! 'SSW'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_short_up1, sender_interp_short_up2 /)
                    data_bounds(1:2,2) = (/ sender_interp_long_down1, sender_interp_long_down2 /)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs-g-sh_end*2
                    data_bounds(2,1) = Bs+g-sh_start*2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g

                end if

            case(13)
                ! 'ENE'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_long_down1, sender_interp_long_down2 /)
                    data_bounds(1:2,2) = (/ sender_interp_short_up1, sender_interp_short_up2/)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2

                end if

            case(14)
                ! 'ESE'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/sender_interp_long_up1, sender_interp_long_up2/)
                    data_bounds(1:2,2) = (/sender_interp_short_up1, sender_interp_short_up2/)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = Bs-g-sh_end*2
                    data_bounds(2,2) = Bs+g-sh_start*2

                end if

            case(15)
                ! 'WNW'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_long_down1, sender_interp_long_down2/)
                    data_bounds(1:2,2) = (/ sender_interp_short_down1, sender_interp_short_down2/)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2

                end if

            case(16)
                ! 'WSW'
                if ( level_diff == -1 ) then
                    data_bounds(1:2,1) = (/ sender_interp_long_up1, sender_interp_long_up2 /)
                    data_bounds(1:2,2) = (/ sender_interp_short_down1, sender_interp_short_down2/)

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = g+1+sh_start*2
                    data_bounds(2,2) = g+1+g+g+sh_end*2

                end if

            end select
        end if

    case('receiver')
        if ( params%threeD_case ) then
            !---3D------3D------3D------3D------3D------3D------3D------3D---
            select case(neighborhood)
                ! '__1/___'
            case(1)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start

                ! '__2/___'
            case(2)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

                ! '__3/___'
            case(3)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

                ! '__4/___'
            case(4)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

                ! '__5/___'
            case(5)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

                ! '__6/___'
            case(6)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = Bs+g+sh_start
                data_bounds(2,3) = Bs+g+g+sh_end

                ! '_12/___'
            case(7)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start

                ! '_13/___'
            case(8)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start

                ! '_14/___'
            case(9)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start

                ! '_15/___'
            case(10)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start

                ! '_62/___'
            case(11)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end
                data_bounds(1,3) = Bs+g+sh_start
                data_bounds(2,3) = Bs+g+g+sh_end

                ! '_63/___'
            case(12)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = Bs+g+sh_start
                data_bounds(2,3) = Bs+g+g+sh_end

                ! '_64/___'
            case(13)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = Bs+g+sh_start
                data_bounds(2,3) = Bs+g+g+sh_end

                ! '_65/___'
            case(14)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g
                data_bounds(1,3) = Bs+g+sh_start
                data_bounds(2,3) = Bs+g+g+sh_end

                ! '_23/___'
            case(15)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

                ! '_25/___'
            case(16)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

                ! '_43/___'
            case(17)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

                ! '_45/___'
            case(18)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = g+1
                data_bounds(2,3) = Bs+g

            case(19,20,21,22)
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(19) ! '123/___'
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                case(20) ! '134/___'
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                case(21) ! '145/___'
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                case(22) ! '152/___'
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                end select

            case(23,24,25,26)
                data_bounds(1,3) = Bs+g+sh_start
                data_bounds(2,3) = Bs+g+g+sh_end
                select case(neighborhood)
                case(23) ! '623/___'
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                case(24) ! '634/___'
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                case(25) ! '645/___'
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                case(26) ! '652/___'
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                end select

            case(27,28,29,30)
                if ( level_diff == -1 ) then
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(27) ! '__1/123'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g

                    case(28) ! '__1/134'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(29) ! '__1/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(30) ! '__1/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(27) ! '__1/123'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2

                    case(28) ! '__1/134'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(29) ! '__1/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(30) ! '__1/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2
                    end select

                end if

            case(31,32,33,34)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    ! first, third dimension
                    select case(neighborhood)
                    case(31) ! '__2/123'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(32) ! '__2/623'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    case(33) ! '__2/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(34) ! '__2/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    ! first, third dimension
                    select case(neighborhood)
                    case(31) ! '__2/123'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(32) ! '__2/623'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    case(33) ! '__2/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(34) ! '__2/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    end select

                end if

            case(35,36,37,38)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    ! second, third dimension
                    select case(neighborhood)
                    case(35) ! '__3/123'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(37) ! '__3/134'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(38) ! '__3/634'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    case(36) ! '__3/623'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    ! second, third dimension
                    select case(neighborhood)
                    case(35) ! '__3/123'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(37) ! '__3/134'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(38) ! '__3/634'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    case(36) ! '__3/623'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    end select

                end if

            case(39,40,41,42)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    ! first, third dimension
                    select case(neighborhood)
                    case(40) ! '__4/634'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    case(39) ! '__4/134'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(41) ! '__4/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(42) ! '__4/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    ! first, third dimension
                    select case(neighborhood)
                    case(40) ! '__4/634'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    case(39) ! '__4/134'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(41) ! '__4/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(42) ! '__4/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    end select

                end if

            case(43,44,45,46)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    ! second, third dimension
                    select case(neighborhood)
                    case(45) ! '__5/152'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(43) ! '__5/145'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(44) ! '__5/645'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    case(46) ! '__5/652'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    ! second, third dimension
                    select case(neighborhood)
                    case(45) ! '__5/152'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(43) ! '__5/145'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(44) ! '__5/645'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    case(46) ! '__5/652'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2

                    end select

                end if

            case(47,48,49,50)
                if ( level_diff == -1 ) then
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(47) ! '__6/623'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g

                    case(48) ! '__6/634'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(49) ! '__6/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(50) ! '__6/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(47) ! '__6/623'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2

                    case(48) ! '__6/634'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(49) ! '__6/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(50) ! '__6/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2

                    end select

                end if

            case(51,52)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(51) ! '_12/123'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g

                    case(52) ! '_12/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(51) ! '_12/123'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(52) ! '_12/152'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                    end select

                end if

            case(53,54)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(54) ! '_13/134'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(53) ! '_13/123'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(54) ! '_13/134'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(53) ! '_13/123'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2
                    end select

                end if

            case(55,56)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(55) ! '_14/134'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g

                    case(56) ! '_14/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(55) ! '_14/134'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(56) ! '_14/145'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2

                    end select

                end if

            case(57,58)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(57) ! '_15/145'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(58) ! '_15/152''
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,3) = 1-sh_end
                    data_bounds(2,3) = g+1-sh_start
                    select case(neighborhood)
                    case(57) ! '_15/145'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(58) ! '_15/152''
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2

                    end select

                end if

            case(59,60)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(59) ! '_62/623'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g

                    case(60) ! '_62/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(59) ! '_62/623'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(60) ! '_62/652'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2
                    end select

                end if

            case(61,62)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(62) ! '_63/634'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(61) ! '_63/623'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(62) ! '_63/634'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(61) ! '_63/623'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2
                    end select

                end if

            case(63,64)
                if ( level_diff == -1 ) then
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(63) ! '_64/634'
                        data_bounds(1,1) = 1
                        data_bounds(2,1) = Bs+g

                    case(64) ! '_64/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(63) ! '_64/634'
                        data_bounds(1,1) = g+(Bs+1)/2
                        data_bounds(2,1) = Bs+g

                    case(64) ! '_64/645'
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = g+(Bs+1)/2

                    end select

                end if

            case(65,66)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(65) ! '_65/645'
                        data_bounds(1,2) = 1
                        data_bounds(2,2) = Bs+g

                    case(66) ! '_65/652'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+2*g

                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,3) = Bs+g+sh_start
                    data_bounds(2,3) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(65) ! '_65/645'
                        data_bounds(1,2) = g+(Bs+1)/2
                        data_bounds(2,2) = Bs+g

                    case(66) ! '_65/652'
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = g+(Bs+1)/2

                    end select

                end if

            case(67,68)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(67) ! '_23/123'
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(68) ! '_23/236''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(67) ! '_23/123'
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(68) ! '_23/236''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2
                    end select

                end if

            case(69,70)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(69) ! '_25/152'
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(70) ! '_25/652''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end
                    select case(neighborhood)
                    case(69) ! '_25/152'
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(70) ! '_25/652''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2
                    end select

                end if

            case(71,72)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    select case(neighborhood)
                    case(71) ! '_43/134'
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(72) ! '_43/634''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    select case(neighborhood)
                    case(71) ! '_43/134'
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(72) ! '_43/634''
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2
                    end select

                end if

            case(73,74)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    select case(neighborhood)
                    case(73) ! '_45/145'
                        data_bounds(1,3) = 1
                        data_bounds(2,3) = Bs+g

                    case(74) ! '_45/645'
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+2*g
                    end select

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start
                    select case(neighborhood)
                    case(73) ! '_45/145'
                        data_bounds(1,3) = g+(Bs+1)/2
                        data_bounds(2,3) = Bs+g

                    case(74) ! '_45/645'
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = g+(Bs+1)/2
                    end select

                end if

            end select
        else
            !---2D------2D------2D------2D------2D------2D------2D------2D---
            select case(neighborhood)
                ! '__N'
            case(1)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g

                ! '__E'
            case(2)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

                ! '__S'
            case(3)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs+g

                ! '__W'
            case(4)
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs+g
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end

                ! '_NE'
            case(5)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

                ! '_NW'
            case(6)
                data_bounds(1,1) = Bs+g+sh_start
                data_bounds(2,1) = Bs+g+g+sh_end
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end

                ! '_SE'
            case(7)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

                ! '_SW'
            case(8)
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = Bs+g+sh_start
                data_bounds(2,2) = Bs+g+g+sh_end

                ! 'NNE'
            case(9)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs+g

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = g+(Bs+1)/2
                    data_bounds(2,2) = Bs+g

                end if

                ! 'NNW'
            case(10)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g+g

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = Bs+g+sh_start
                    data_bounds(2,1) = Bs+g+g+sh_end
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs+1)/2

                end if

                ! 'SSE'
            case(11)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs+g

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = g+(Bs+1)/2
                    data_bounds(2,2) = Bs+g

                end if

                ! 'SSW'
            case(12)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs+g+g

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = 1-sh_end
                    data_bounds(2,1) = g+1-sh_start
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs+1)/2

                end if

                ! 'ENE'
            case(13)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g+g
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs+1)/2
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                end if

                ! 'ESE'
            case(14)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+(Bs+1)/2
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = 1-sh_end
                    data_bounds(2,2) = g+1-sh_start

                end if

                ! 'WNW'
            case(15)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs+g+g
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs+1)/2
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                end if

                ! 'WSW'
            case(16)
                if ( level_diff == -1 ) then
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                elseif ( level_diff == 1 ) then
                    data_bounds(1,1) = g+(Bs+1)/2
                    data_bounds(2,1) = Bs+g
                    data_bounds(1,2) = Bs+g+sh_start
                    data_bounds(2,2) = Bs+g+g+sh_end

                end if

            end select
        end if

    case('restricted-predicted')
        if ( params%threeD_case ) then
            !---3D------3D------3D------3D------3D------3D------3D------3D---
            if (level_diff == -1) then
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
            elseif (level_diff == 1) then
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
            endif

        else
            !---2D------2D------2D------2D------2D------2D------2D------2D---
            if (level_diff == -1) then

                ! we set only the origin (in the interpolated, fine domain), the upper bounds are dictated by what the
                ! recver expects
                select case(neighborhood)
                case(5)
                    data_bounds(1, 1:2) = (/1+2*S, 1+2*S+2*A/)
                case(6)
                    data_bounds(1, 1:2) = (/1+2*S, 1+2*S/)
                case(7)
                    data_bounds(1, 1:2) = (/1+2*S+2*A, 1+2*S+2*A/)
                case(8)
                    data_bounds(1, 1:2) = (/1+2*S+2*A, 1+2*S/)
                case(9)
                    data_bounds(1, 1:2) = (/1+2*S, 1+2*S+2*A/)
                case(10)
                    data_bounds(1, 1:2) = (/1+2*S, 1+2*S/)
                case(11)
                    data_bounds(1, 1:2) = (/1+2*S+2*A, 1+2*S+2*A/)
                case(12)
                    data_bounds(1, 1:2) = (/1+2*S+2*A, 1+2*S/)
                case(13)
                    data_bounds(1, 1:2) = (/1+2*S, 1+2*S+2*A/)
                case(14)
                    data_bounds(1, 1:2) = (/1+2*S+2*A, 1+2*S+2*A/)
                case(15)
                    data_bounds(1, 1:2) = (/1+2*S, 1+2*S/)
                case(16)
                    data_bounds(1, 1:2) = (/1+2*S+2*A, 1+2*S/)
                end select

                ! if the redundant point shall be excluded, there is just a small
                ! shift in the origin again
                if (data_bounds_type == exclude_redundant) then
                    select case(neighborhood)
                    case(5)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/1, 0/)
                    case(6)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/1, 1/)
                    case(7)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/0, 0/)
                    case(8)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/0, 1/)
                    case(9)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/1, -1/)
                    case(10)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/1, 1/)
                    case(11)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/0, -1/)
                    case(12)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/0, 1/)
                    case(13)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/1, 0/)
                    case(14)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/-1, 0/)
                    case(15)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/1, 1/)
                    case(16)
                        data_bounds(1, 1:2) = data_bounds(1, 1:2) + (/-1, 1/)
                    end select
                endif

                ! note for ONLY_REDUNDANT, the corners 5-8 degenerate to a single
                ! point, for the half-lines 9-16 it is a line
                if (data_bounds_type == ONLY_REDUNDANT) then
                    ! the point were interested in is
                    N = 2*g+1 - 3 -1
                    select case(neighborhood)
                    case(5)
                        data_bounds(1, 1:2) = (/1+2*S, N+2*A+2*S/)
                    case(6)
                        data_bounds(1, 1:2) = (/1+2*S, 1+2*S/)
                    case(7)
                        data_bounds(1, 1:2) = (/N+2*A+2*S, N+2*A+2*S/)
                    case(8)
                        data_bounds(1, 1:2) = (/N+2*A+2*S, 1+2*S/)
                    case(9)
                        ! works
                    case(10)
                        ! works
                    case(11)
                        data_bounds(1, 1) = N + 2*A + 2*S
                    case(12)
                        data_bounds(1, 1) = N + 2*A + 2*S
                    case(13)
                        data_bounds(1, 2) = N + 2*A + 2*S
                    case(14)
                        data_bounds(1, 2) = N + 2*A + 2*S
                    case(15)
                        ! works
                    case(16)
                        ! works
                    end select
                endif

                ! set upper bounds according to what recver expects
                N = ijkGhosts(2,1, neighborhood, level_diff, data_bounds_type, RECVER) &
                  - ijkGhosts(1,1, neighborhood, level_diff, data_bounds_type, RECVER) + 1
                data_bounds(2,1) = data_bounds(1,1) + N - 1

                N = ijkGhosts(2,2, neighborhood, level_diff, data_bounds_type, RECVER) &
                  - ijkGhosts(1,2, neighborhood, level_diff, data_bounds_type, RECVER) + 1
                data_bounds(2,2) = data_bounds(1,2) + N - 1

            elseif (level_diff == +1) then
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

                    ! the following bounds are special and particularily ugly.
                    ! the history is that mario had restrict_predict_data return array bounds
                    ! (which are the ones defined here in the 'restricted-predicted' fork).
                    ! For the 9:16 relations, he computed the return value from the incoming value,
                    ! which here is not possible (no input value). thus, I here explicitly write his
                    ! formula
                    ! NOTE: would well have been possible. Sender/recver bounds are filled before, and could thus be used.

                    ! 'NNE'
                case(9)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (g+1+g+g+sh_end*2 -(g+1+sh_start*2) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (Bs+g-(g+1) + 1 + 1)/2

                    ! 'NNW'
                case(10)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (g+1+g+g+sh_end*2-(g+1+sh_start*2) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (Bs+g-(g+1) + 1 + 1)/2

                    ! 'SSE'
                case(11)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (Bs+g-sh_start*2-(Bs-g-sh_end*2) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (Bs+g-(g+1) + 1 + 1)/2
                    ! 'SSW'
                case(12)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (Bs+g-sh_start*2-(Bs-g-sh_end*2) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (Bs+g-(g+1) + 1 + 1)/2

                    ! 'ENE'
                case(13)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (Bs+g-(g+1) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (Bs+g-sh_start*2-(Bs-g-sh_end*2) + 1 + 1)/2

                    ! 'ESE'
                case(14)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (Bs+g-(g+1) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (Bs+g-sh_start*2-(Bs-g-sh_end*2) + 1 + 1)/2

                    ! 'WNW'
                case(15)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (Bs+g-(g+1) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (g+1+g+g+sh_end*2-(g+1+sh_start*2) + 1 + 1)/2

                    ! 'WSW'
                case(16)
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = (Bs+g-(g+1) + 1 + 1)/2
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = (g+1+g+g+sh_end*2-(g+1+sh_start*2) + 1 + 1)/2

                end select
            endif
        end if

    case default
        call abort(06602338, "Calc-data-bounds: this is an unknown string.")
    end select
end subroutine calc_data_bounds
