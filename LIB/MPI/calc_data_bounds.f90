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

    ! grid parameter
    integer(kind=ik)                                :: Bs, g

    ! start and edn shift values
    integer(kind=ik)                                :: sh_start, sh_end

!---------------------------------------------------------------------------------------------
! interfaces

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    sh_start = 0
    sh_end   = 0

    if ( data_bounds_type == exclude_redundant ) then
        sh_start = 1
    end if
    if ( data_bounds_type == only_redundant ) then
        sh_end = -g
    end if

    ! reset data bounds
    data_bounds = 1

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    select case(sender_or_receiver)

        case('sender')

            if ( params%threeD_case ) then
                ! 3D
                select case(neighborhood)
                    ! '__1/___'
                    case(1)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = Bs-sh_end
                        data_bounds(2,3) = Bs+g-sh_start

                    ! '__2/___'
                    case(2)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__3/___'
                    case(3)
                        ! first dimension
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__4/___'
                    case(4)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__5/___'
                    case(5)
                        ! first dimension
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__6/___'
                    case(6)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = g+1+sh_start
                        data_bounds(2,3) = g+1+g+sh_end

                    ! '_12/___'
                    case(7)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = Bs-sh_end
                        data_bounds(2,3) = Bs+g-sh_start

                    ! '_13/___'
                    case(8)
                        ! first dimension
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = Bs-sh_end
                        data_bounds(2,3) = Bs+g-sh_start

                    ! '_14/___'
                    case(9)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start
                        ! third dimension
                        data_bounds(1,3) = Bs-sh_end
                        data_bounds(2,3) = Bs+g-sh_start

                    ! '_15/___'
                    case(10)
                        ! first dimension
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = Bs-sh_end
                        data_bounds(2,3) = Bs+g-sh_start

                      ! '_62/___'
                    case(11)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = g+1+sh_start
                        data_bounds(2,3) = g+1+g+sh_end

                    ! '_63/___'
                    case(12)
                        ! first dimension
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = g+1+sh_start
                        data_bounds(2,3) = g+1+g+sh_end

                    ! '_64/___'
                    case(13)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start
                        ! third dimension
                        data_bounds(1,3) = g+1+sh_start
                        data_bounds(2,3) = g+1+g+sh_end

                    ! '_65/___'
                    case(14)
                        ! first dimension
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = g+1+sh_start
                        data_bounds(2,3) = g+1+g+sh_end

                    ! '_23/___'
                    case(15)
                        ! first dimension
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '_25/___'
                    case(16)
                        ! first dimension
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = Bs+1+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '_43/___'
                    case(17)
                        ! first dimension
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        ! second dimension
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '_45/___'
                    case(18)
                        ! first dimension
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = Bs+1+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    case(19,20,21,22)
                        if ( level_diff == 0 ) then
                            ! third dimension
                            data_bounds(1,3) = Bs-sh_end
                            data_bounds(2,3) = Bs+g-sh_start
                            ! first, second dimension
                            select case(neighborhood)
                                case(19) ! '123/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-sh_end
                                    data_bounds(2,1) = Bs+g-sh_start
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start
                                    data_bounds(2,2) = g+1+g+sh_end

                                case(20) ! '134/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-sh_end
                                    data_bounds(2,1) = Bs+g-sh_start
                                    ! second dimension
                                    data_bounds(1,2) = Bs-sh_end
                                    data_bounds(2,2) = Bs+g-sh_start

                                case(21) ! '145/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start
                                    data_bounds(2,1) = g+1+g+sh_end
                                    ! second dimension
                                    data_bounds(1,2) = Bs-sh_end
                                    data_bounds(2,2) = Bs+g-sh_start

                                case(22) ! '152/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start
                                    data_bounds(2,1) = g+1+g+sh_end
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start
                                    data_bounds(2,2) = g+1+g+sh_end

                            end select

                        elseif ( level_diff == -1 ) then
                            ! third dimension
                            data_bounds(1,3) = Bs+1
                            data_bounds(2,3) = Bs+g
                            ! first, second dimension
                            select case(neighborhood)
                                case(19) ! '123/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs+1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+g

                                case(20) ! '134/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs+1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = Bs+1
                                    data_bounds(2,2) = Bs+g

                                case(21) ! '145/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+g
                                    ! second dimension
                                    data_bounds(1,2) = Bs+1
                                    data_bounds(2,2) = Bs+g

                                case(22) ! '152/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! third dimension
                            data_bounds(1,3) = Bs-g-sh_end*2
                            data_bounds(2,3) = Bs+g-sh_start*2
                            ! first, second dimension
                            select case(neighborhood)
                                case(19) ! '123/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-g-sh_end*2
                                    data_bounds(2,1) = Bs+g-sh_start*2
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start*2
                                    data_bounds(2,2) = g+1+g+g+sh_end*2

                                case(20) ! '134/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-g-sh_end*2
                                    data_bounds(2,1) = Bs+g-sh_start*2
                                    ! second dimension
                                    data_bounds(1,2) = Bs-g-sh_end*2
                                    data_bounds(2,2) = Bs+g-sh_start*2

                                case(21) ! '145/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start*2
                                    data_bounds(2,1) = g+1+g+g+sh_end*2
                                    ! second dimension
                                    data_bounds(1,2) = Bs-g-sh_end*2
                                    data_bounds(2,2) = Bs+g-sh_start*2

                                case(22) ! '152/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start*2
                                    data_bounds(2,1) = g+1+g+g+sh_end*2
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start*2
                                    data_bounds(2,2) = g+1+g+g+sh_end*2

                            end select
                        end if

                    case(23,24,25,26)
                        if ( level_diff == 0 ) then
                            ! third dimension
                            data_bounds(1,3) = g+1+sh_start
                            data_bounds(2,3) = g+1+g+sh_end
                            ! first, second dimension
                            select case(neighborhood)
                                case(23) ! '623/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-sh_end
                                    data_bounds(2,1) = Bs+g-sh_start
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start
                                    data_bounds(2,2) = g+1+g+sh_end

                                case(24) ! '634/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-sh_end
                                    data_bounds(2,1) = Bs+g-sh_start
                                    ! second dimension
                                    data_bounds(1,2) = Bs-sh_end
                                    data_bounds(2,2) = Bs+g-sh_start

                                case(25) ! '645/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start
                                    data_bounds(2,1) = g+1+g+sh_end
                                    ! second dimension
                                    data_bounds(1,2) = Bs-sh_end
                                    data_bounds(2,2) = Bs+g-sh_start

                                case(26) ! '652/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start
                                    data_bounds(2,1) = g+1+g+sh_end
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start
                                    data_bounds(2,2) = g+1+g+sh_end

                            end select

                        elseif ( level_diff == -1 ) then
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = g+g
                            ! first, second dimension
                            select case(neighborhood)
                                case(23) ! '623/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs+1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+g

                                case(24) ! '634/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs+1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = Bs+1
                                    data_bounds(2,2) = Bs+g

                                case(25) ! '645/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+g
                                    ! second dimension
                                    data_bounds(1,2) = Bs+1
                                    data_bounds(2,2) = Bs+g

                                case(26) ! '652/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! third dimension
                            data_bounds(1,3) = g+1+sh_start*2
                            data_bounds(2,3) = g+1+g+g+sh_end*2
                            ! first, second dimension
                            select case(neighborhood)
                                case(23) ! '623/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-g-sh_end*2
                                    data_bounds(2,1) = Bs+g-sh_start*2
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start*2
                                    data_bounds(2,2) = g+1+g+g+sh_end*2

                                case(24) ! '634/___'
                                    ! first dimension
                                    data_bounds(1,1) = Bs-g-sh_end*2
                                    data_bounds(2,1) = Bs+g-sh_start*2
                                    ! second dimension
                                    data_bounds(1,2) = Bs-g-sh_end*2
                                    data_bounds(2,2) = Bs+g-sh_start*2

                                case(25) ! '645/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start*2
                                    data_bounds(2,1) = g+1+g+g+sh_end*2
                                    ! second dimension
                                    data_bounds(1,2) = Bs-g-sh_end*2
                                    data_bounds(2,2) = Bs+g-sh_start*2

                                case(26) ! '652/___'
                                    ! first dimension
                                    data_bounds(1,1) = g+1+sh_start*2
                                    data_bounds(2,1) = g+1+g+g+sh_end*2
                                    ! second dimension
                                    data_bounds(1,2) = g+1+sh_start*2
                                    data_bounds(2,2) = g+1+g+g+sh_end*2

                            end select
                        end if

                    case(27,28,29,30)
                        if ( level_diff == -1 ) then
                            ! third dimension
                            data_bounds(1,3) = (Bs+1)/2
                            data_bounds(2,3) = Bs+g
                            ! first, second dimension
                            select case(neighborhood)
                                case(27) ! '__1/123'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g

                                case(28) ! '__1/134'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(29) ! '__1/145'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(30) ! '__1/152'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! third dimension
                            data_bounds(1,3) = Bs-g-sh_end*2
                            data_bounds(2,3) = Bs+g-sh_start*2
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g

                        end if

                    case(31,32,33,34)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g
                            ! first, third dimension
                            select case(neighborhood)
                                case(32) ! '__2/623'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                                case(31) ! '__2/123'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(33) ! '__2/152'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(34) ! '__2/652'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                    case(35,36,37,38)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second, third dimension
                            select case(neighborhood)
                                case(35) ! '__3/123'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(37) ! '__3/134'
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(38) ! '__3/634'
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                                case(36) ! '__3/623'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                    case(39,40,41,42)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g
                            ! first, third dimension
                            select case(neighborhood)
                                case(40) ! '__4/634'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                                case(39) ! '__4/134'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(41) ! '__4/145'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(42) ! '__4/645'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                    case(43,44,45,46)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! second, third dimension
                            select case(neighborhood)
                                case(45) ! '__5/152'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(43) ! '__5/145'
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(44) ! '__5/645'
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                                case(46) ! '__5/652'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                    case(47,48,49,50)
                        if ( level_diff == -1 ) then
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = (Bs+1)/2+g+g
                            ! first, second dimension
                            select case(neighborhood)
                                case(47) ! '__6/623'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g

                                case(48) ! '__6/634'
                                    ! first dimension
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(49) ! '__6/645'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! second dimension
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(50) ! '__6/652'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! third dimension
                            data_bounds(1,3) = g+1+sh_start*2
                            data_bounds(2,3) = g+1+g+g+sh_end*2
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g

                        end if

                    case(51,52)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g
                            ! third dimension
                            data_bounds(1,3) = (Bs+1)/2
                            data_bounds(2,3) = Bs+g
                            ! first dimension
                            select case(neighborhood)
                                case(51) ! '_12/123'
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g

                                case(52) ! '_12/152'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2
                            ! third dimension
                            data_bounds(1,3) = Bs-g-sh_end*2
                            data_bounds(2,3) = Bs+g-sh_start*2

                        end if

                    case(53,54)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = (Bs+1)/2
                            data_bounds(2,3) = Bs+g
                            ! second dimension
                            select case(neighborhood)
                                case(54) ! '_13/134'
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(53) ! '_13/123'
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = Bs-g-sh_end*2
                            data_bounds(2,3) = Bs+g-sh_start*2

                        end if

                    case(55,56)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = (Bs+1)/2
                            data_bounds(2,3) = Bs+g
                            ! first dimension
                            select case(neighborhood)
                                case(55) ! '_14/134'
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g

                                case(56) ! '_14/145'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2
                            ! third dimension
                            data_bounds(1,3) = Bs-g-sh_end*2
                            data_bounds(2,3) = Bs+g-sh_start*2

                        end if

                    case(57,58)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! third dimension
                            data_bounds(1,3) = (Bs+1)/2
                            data_bounds(2,3) = Bs+g
                            ! second dimension
                            select case(neighborhood)
                                case(57) ! '_15/145'
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(58) ! '_15/152''
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = Bs-g-sh_end*2
                            data_bounds(2,3) = Bs+g-sh_start*2

                        end if

                    case(59,60)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = (Bs+1)/2+g+g
                            ! first dimension
                            select case(neighborhood)
                                case(59) ! '_62/623'
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g

                                case(60) ! '_62/652'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2
                            ! third dimension
                            data_bounds(1,3) = g+1+sh_start*2
                            data_bounds(2,3) = g+1+g+g+sh_end*2

                        end if

                     case(61,62)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = (Bs+1)/2+g+g
                            ! second dimension
                            select case(neighborhood)
                                case(62) ! '_63/634'
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(61) ! '_63/623'
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1+sh_start*2
                            data_bounds(2,3) = g+1+g+g+sh_end*2

                        end if

                    case(63,64)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = (Bs+1)/2+g+g
                            ! first dimension
                            select case(neighborhood)
                                case(63) ! '_64/634'
                                    data_bounds(1,1) = (Bs+1)/2
                                    data_bounds(2,1) = Bs+g

                                case(64) ! '_64/645'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2
                            ! third dimension
                            data_bounds(1,3) = g+1+sh_start*2
                            data_bounds(2,3) = g+1+g+g+sh_end*2

                        end if

                    case(65,66)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = (Bs+1)/2+g+g
                            ! second dimension
                            select case(neighborhood)
                                case(65) ! '_65/645'
                                    data_bounds(1,2) = (Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(66) ! '_65/652'
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = (Bs+1)/2+g+g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            data_bounds(1,3) = g+1+sh_start*2
                            data_bounds(2,3) = g+1+g+g+sh_end*2

                        end if

                     case(67,68)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g
                            ! third dimension
                            select case(neighborhood)
                                case(67) ! '_23/123'
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(68) ! '_23/236''
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                     case(69,70)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g
                            ! third dimension
                            select case(neighborhood)
                                case(69) ! '_25/152'
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(70) ! '_25/652''
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                     case(71,72)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            select case(neighborhood)
                                case(71) ! '_43/134'
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(72) ! '_43/634''
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                     case(73,74)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g
                            ! third dimension
                            select case(neighborhood)
                                case(73) ! '_45/145'
                                    data_bounds(1,3) = (Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(74) ! '_45/645'
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = (Bs+1)/2+g+g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2
                            ! third dimension
                            data_bounds(1,3) = g+1
                            data_bounds(2,3) = Bs+g

                        end if

                end select

            else
                ! 2D
                select case(neighborhood)
                    ! '__N'
                    case(1)
                        ! first dimension
                        data_bounds(1,1) = g+1+sh_start
                        data_bounds(2,1) = g+1+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g

                    ! '__E'
                    case(2)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs-sh_end
                        data_bounds(2,2) = Bs+g-sh_start

                    ! '__S'
                    case(3)
                        ! first dimension
                        data_bounds(1,1) = Bs-sh_end
                        data_bounds(2,1) = Bs+g-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g

                    ! '__W'
                    case(4)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1+sh_start
                        data_bounds(2,2) = g+1+g+sh_end

                    ! '_NE'
                    case(5)
                        if ( level_diff == 0 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start
                            data_bounds(2,1) = g+1+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = Bs-sh_end
                            data_bounds(2,2) = Bs+g-sh_start

                        elseif ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = g+g
                            ! second dimension
                            data_bounds(1,2) = Bs+1
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2

                        end if

                    ! '_NW'
                    case(6)
                        if ( level_diff == 0 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start
                            data_bounds(2,1) = g+1+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start
                            data_bounds(2,2) = g+1+g+sh_end

                        elseif ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = g+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = g+g

                        elseif ( level_diff == 1) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2

                        end if

                    ! '_SE'
                    case(7)
                        if ( level_diff == 0 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-sh_end
                            data_bounds(2,1) = Bs+g-sh_start
                            ! second dimension
                            data_bounds(1,2) = Bs-sh_end
                            data_bounds(2,2) = Bs+g-sh_start

                        elseif ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = Bs+1
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2

                        end if

                    ! '_SW'
                    case(8)
                        if ( level_diff == 0 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-sh_end
                            data_bounds(2,1) = Bs+g-sh_start
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start
                            data_bounds(2,2) = g+1+g+sh_end

                        elseif ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = g+g

                        elseif ( level_diff == 1) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2

                        end if

                    ! 'NNE'
                    case(9)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g

                        end if

                    ! 'NNW'
                    case(10)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1+sh_start*2
                            data_bounds(2,1) = g+1+g+g+sh_end*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g

                        end if

                    ! 'SSE'
                    case(11)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g

                        end if

                    ! 'SSW'
                    case(12)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs-g-sh_end*2
                            data_bounds(2,1) = Bs+g-sh_start*2
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g

                        end if

                    ! 'ENE'
                    case(13)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2

                        end if

                    ! 'ESE'
                    case(14)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = (Bs+1)/2
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = Bs-g-sh_end*2
                            data_bounds(2,2) = Bs+g-sh_start*2

                        end if

                    ! 'WNW'
                    case(15)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = (Bs+1)/2+g+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2

                        end if

                    ! 'WSW'
                    case(16)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = (Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = (Bs+1)/2+g+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = g+1+sh_start*2
                            data_bounds(2,2) = g+1+g+g+sh_end*2

                        end if

                end select
            end if

        case('receiver')

            if ( params%threeD_case ) then
                ! 3D
                select case(neighborhood)
                    ! '__1/___'
                    case(1)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = 1-sh_end
                        data_bounds(2,3) = g+1-sh_start

                    ! '__2/___'
                    case(2)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__3/___'
                    case(3)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__4/___'
                    case(4)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__5/___'
                    case(5)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '__6/___'
                    case(6)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = Bs+g+sh_start
                        data_bounds(2,3) = Bs+g+g+sh_end

                    ! '_12/___'
                    case(7)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = 1-sh_end
                        data_bounds(2,3) = g+1-sh_start

                    ! '_13/___'
                    case(8)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = 1-sh_end
                        data_bounds(2,3) = g+1-sh_start

                    ! '_14/___'
                    case(9)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start
                        ! third dimension
                        data_bounds(1,3) = 1-sh_end
                        data_bounds(2,3) = g+1-sh_start

                    ! '_15/___'
                    case(10)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = 1-sh_end
                        data_bounds(2,3) = g+1-sh_start

                      ! '_62/___'
                    case(11)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = Bs+g+sh_start
                        data_bounds(2,3) = Bs+g+g+sh_end

                    ! '_63/___'
                    case(12)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = Bs+g+sh_start
                        data_bounds(2,3) = Bs+g+g+sh_end

                    ! '_64/___'
                    case(13)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start
                        ! third dimension
                        data_bounds(1,3) = Bs+g+sh_start
                        data_bounds(2,3) = Bs+g+g+sh_end

                    ! '_65/___'
                    case(14)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g
                        ! third dimension
                        data_bounds(1,3) = Bs+g+sh_start
                        data_bounds(2,3) = Bs+g+g+sh_end

                    ! '_23/___'
                    case(15)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '_25/___'
                    case(16)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '_43/___'
                    case(17)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    ! '_45/___'
                    case(18)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start
                        ! third dimension
                        data_bounds(1,3) = g+1
                        data_bounds(2,3) = Bs+g

                    case(19,20,21,22)
                        ! third dimension
                        data_bounds(1,3) = 1-sh_end
                        data_bounds(2,3) = g+1-sh_start
                        ! first, second dimension
                        select case(neighborhood)
                            case(19) ! '123/___'
                                ! first dimension
                                data_bounds(1,1) = 1-sh_end
                                data_bounds(2,1) = g+1-sh_start
                                ! second dimension
                                data_bounds(1,2) = Bs+g+sh_start
                                data_bounds(2,2) = Bs+g+g+sh_end

                            case(20) ! '134/___'
                                ! first dimension
                                data_bounds(1,1) = 1-sh_end
                                data_bounds(2,1) = g+1-sh_start
                                ! second dimension
                                data_bounds(1,2) = 1-sh_end
                                data_bounds(2,2) = g+1-sh_start

                            case(21) ! '145/___'
                                ! first dimension
                                data_bounds(1,1) = Bs+g+sh_start
                                data_bounds(2,1) = Bs+g+g+sh_end
                                ! second dimension
                                data_bounds(1,2) = 1-sh_end
                                data_bounds(2,2) = g+1-sh_start

                            case(22) ! '152/___'
                                ! first dimension
                                data_bounds(1,1) = Bs+g+sh_start
                                data_bounds(2,1) = Bs+g+g+sh_end
                                ! second dimension
                                data_bounds(1,2) = Bs+g+sh_start
                                data_bounds(2,2) = Bs+g+g+sh_end

                        end select

                    case(23,24,25,26)
                        ! third dimension
                        data_bounds(1,3) = Bs+g+sh_start
                        data_bounds(2,3) = Bs+g+g+sh_end
                        ! first, second dimension
                        select case(neighborhood)
                            case(23) ! '623/___'
                                ! first dimension
                                data_bounds(1,1) = 1-sh_end
                                data_bounds(2,1) = g+1-sh_start
                                ! second dimension
                                data_bounds(1,2) = Bs+g+sh_start
                                data_bounds(2,2) = Bs+g+g+sh_end

                            case(24) ! '634/___'
                                ! first dimension
                                data_bounds(1,1) = 1-sh_end
                                data_bounds(2,1) = g+1-sh_start
                                ! second dimension
                                data_bounds(1,2) = 1-sh_end
                                data_bounds(2,2) = g+1-sh_start

                            case(25) ! '645/___'
                                ! first dimension
                                data_bounds(1,1) = Bs+g+sh_start
                                data_bounds(2,1) = Bs+g+g+sh_end
                                ! second dimension
                                data_bounds(1,2) = 1-sh_end
                                data_bounds(2,2) = g+1-sh_start

                            case(26) ! '652/___'
                                ! first dimension
                                data_bounds(1,1) = Bs+g+sh_start
                                data_bounds(2,1) = Bs+g+g+sh_end
                                ! second dimension
                                data_bounds(1,2) = Bs+g+sh_start
                                data_bounds(2,2) = Bs+g+g+sh_end

                        end select

                    case(27,28,29,30)
                        if ( level_diff == -1 ) then
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! first, second dimension
                            select case(neighborhood)
                                case(27) ! '__1/123'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case(28) ! '__1/134'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(29) ! '__1/145'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(30) ! '__1/152'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! first, second dimension
                            select case(neighborhood)
                                case(27) ! '__1/123'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2

                                case(28) ! '__1/134'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(29) ! '__1/145'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(30) ! '__1/152'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2
                            end select

                        end if

                    case(31,32,33,34)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! first, third dimension
                            select case(neighborhood)
                                case(31) ! '__2/123'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(32) ! '__2/623'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                                case(33) ! '__2/152'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(34) ! '__2/652'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! first, third dimension
                            select case(neighborhood)
                                case(31) ! '__2/123'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(32) ! '__2/623'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                                case(33) ! '__2/152'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(34) ! '__2/652'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                            end select

                        end if

                    case(35,36,37,38)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second, third dimension
                            select case(neighborhood)
                                case(35) ! '__3/123'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(37) ! '__3/134'
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(38) ! '__3/634'
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                                case(36) ! '__3/623'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second, third dimension
                            select case(neighborhood)
                                case(35) ! '__3/123'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(37) ! '__3/134'
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(38) ! '__3/634'
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                                case(36) ! '__3/623'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                            end select

                        end if

                    case(39,40,41,42)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! first, third dimension
                            select case(neighborhood)
                                case(40) ! '__4/634'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                                case(39) ! '__4/134'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(41) ! '__4/145'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(42) ! '__4/645'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! first, third dimension
                            select case(neighborhood)
                                case(40) ! '__4/634'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                                case(39) ! '__4/134'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(41) ! '__4/145'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(42) ! '__4/645'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                            end select

                        end if

                    case(43,44,45,46)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second, third dimension
                            select case(neighborhood)
                                case(45) ! '__5/152'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(43) ! '__5/145'
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(44) ! '__5/645'
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                                case(46) ! '__5/652'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second, third dimension
                            select case(neighborhood)
                                case(45) ! '__5/152'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(43) ! '__5/145'
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+(Bs+1)/2
                                    data_bounds(2,3) = Bs+g

                                case(44) ! '__5/645'
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                                case(46) ! '__5/652'
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2
                                    ! third dimension
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = g+(Bs+1)/2

                            end select

                        end if

                    case(47,48,49,50)
                        if ( level_diff == -1 ) then
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! first, second dimension
                            select case(neighborhood)
                                case(47) ! '__6/623'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case(48) ! '__6/634'
                                    ! first dimension
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(49) ! '__6/645'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! second dimension
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(50) ! '__6/652'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! first, second dimension
                            select case(neighborhood)
                                case(47) ! '__6/623'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2

                                case(48) ! '__6/634'
                                    ! first dimension
                                    data_bounds(1,1) = g+(Bs+1)/2
                                    data_bounds(2,1) = Bs+g
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(49) ! '__6/645'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! second dimension
                                    data_bounds(1,2) = g+(Bs+1)/2
                                    data_bounds(2,2) = Bs+g

                                case(50) ! '__6/652'
                                    ! first dimension
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = g+(Bs+1)/2
                                    ! second dimension
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = g+(Bs+1)/2

                            end select

                        end if

                    case(51,52)
                        if ( level_diff == -1 ) then
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! first dimension
                            select case(neighborhood)
                                case(51) ! '_12/123'
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g

                                case(52) ! '_12/152'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! first dimension
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
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! second dimension
                            select case(neighborhood)
                                case(54) ! '_13/134'
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(53) ! '_13/123'
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! second dimension
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
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! first dimension
                            select case(neighborhood)
                                case(55) ! '_14/134'
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g

                                case(56) ! '_14/145'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! first dimension
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
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! second dimension
                            select case(neighborhood)
                                case(57) ! '_15/145'
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(58) ! '_15/152''
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = 1-sh_end
                            data_bounds(2,3) = g+1-sh_start
                            ! second dimension
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
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! first dimension
                            select case(neighborhood)
                                case(59) ! '_62/623'
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g

                                case(60) ! '_62/652'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! first dimension
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
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! second dimension
                            select case(neighborhood)
                                case(62) ! '_63/634'
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(61) ! '_63/623'
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! second dimension
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
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! first dimension
                            select case(neighborhood)
                                case(63) ! '_64/634'
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g

                                case(64) ! '_64/645'
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! first dimension
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
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! second dimension
                            select case(neighborhood)
                                case(65) ! '_65/645'
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case(66) ! '_65/652'
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! third dimension
                            data_bounds(1,3) = Bs+g+sh_start
                            data_bounds(2,3) = Bs+g+g+sh_end
                            ! second dimension
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
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
                            select case(neighborhood)
                                case(67) ! '_23/123'
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(68) ! '_23/236''
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
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
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
                            select case(neighborhood)
                                case(69) ! '_25/152'
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(70) ! '_25/652''
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end
                            ! third dimension
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
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
                            select case(neighborhood)
                                case(71) ! '_43/134'
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(72) ! '_43/634''
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
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
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
                            select case(neighborhood)
                                case(73) ! '_45/145'
                                    data_bounds(1,3) = 1
                                    data_bounds(2,3) = Bs+g

                                case(74) ! '_45/645'
                                    data_bounds(1,3) = g+1
                                    data_bounds(2,3) = Bs+2*g
                            end select

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start
                            ! third dimension
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
                ! 2D
                select case(neighborhood)
                    ! '__N'
                    case(1)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g

                    ! '__E'
                    case(2)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start

                    ! '__S'
                    case(3)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = g+1
                        data_bounds(2,2) = Bs+g

                    ! '__W'
                    case(4)
                        ! first dimension
                        data_bounds(1,1) = g+1
                        data_bounds(2,1) = Bs+g
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end

                    ! '_NE'
                    case(5)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start

                    ! '_NW'
                    case(6)
                        ! first dimension
                        data_bounds(1,1) = Bs+g+sh_start
                        data_bounds(2,1) = Bs+g+g+sh_end
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end

                    ! '_SE'
                    case(7)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = 1-sh_end
                        data_bounds(2,2) = g+1-sh_start

                    ! '_SW'
                    case(8)
                        ! first dimension
                        data_bounds(1,1) = 1-sh_end
                        data_bounds(2,1) = g+1-sh_start
                        ! second dimension
                        data_bounds(1,2) = Bs+g+sh_start
                        data_bounds(2,2) = Bs+g+g+sh_end

                    ! 'NNE'
                    case(9)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = g+(Bs+1)/2
                            data_bounds(2,2) = Bs+g

                        end if

                    ! 'NNW'
                    case(10)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = Bs+g+sh_start
                            data_bounds(2,1) = Bs+g+g+sh_end
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = g+(Bs+1)/2

                        end if

                    ! 'SSE'
                    case(11)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = 1
                            data_bounds(2,2) = Bs+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = g+(Bs+1)/2
                            data_bounds(2,2) = Bs+g

                        end if

                    ! 'SSW'
                    case(12)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = Bs+g+g

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1-sh_end
                            data_bounds(2,1) = g+1-sh_start
                            ! second dimension
                            data_bounds(1,2) = g+1
                            data_bounds(2,2) = g+(Bs+1)/2

                        end if

                    ! 'ENE'
                    case(13)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g+g
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = g+(Bs+1)/2
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start

                        end if

                    ! 'ESE'
                    case(14)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+(Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = 1-sh_end
                            data_bounds(2,2) = g+1-sh_start

                        end if

                    ! 'WNW'
                    case(15)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = Bs+g+g
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+1
                            data_bounds(2,1) = g+(Bs+1)/2
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end

                        end if

                    ! 'WSW'
                    case(16)
                        if ( level_diff == -1 ) then
                            ! first dimension
                            data_bounds(1,1) = 1
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end

                        elseif ( level_diff == 1 ) then
                            ! first dimension
                            data_bounds(1,1) = g+(Bs+1)/2
                            data_bounds(2,1) = Bs+g
                            ! second dimension
                            data_bounds(1,2) = Bs+g+sh_start
                            data_bounds(2,2) = Bs+g+g+sh_end

                        end if

                end select
            end if

    end select

end subroutine calc_data_bounds
