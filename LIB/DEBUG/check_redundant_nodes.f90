!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name check_redundant_nodes.f90
!> \version 0.5
!> \author msr
!
!> \brief check all redundant nodes, if there is a difference between nodes: stop
!> wabbit and write current state to file
!>
!> subroutine structure:
!> ---------------------
!>
!
!>
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
!!
!> \details
!! = log ======================================================================================
!! \n
!! 09/05/17 - create
!! 13/03/18 - add check for redundant ghost nodes, rework subroutine structure
!
! ********************************************************************************************

subroutine check_redundant_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer, stop_status )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)

    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)

    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)     :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)        :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

    ! MPI parameter
    integer(kind=ik)                    :: rank

    ! loop variables
    integer(kind=ik)                    :: N, k, neighborhood, neighbor_num, level_diff

    ! id integers
    integer(kind=ik)                    :: lgt_id, neighbor_light_id, neighbor_rank

    ! type of data bounds
    ! 'exclude_redundant', 'include_redundant', 'only_redundant'
    character(len=25)                   :: data_bounds_type
    integer(kind=ik), dimension(2,3)    :: data_bounds

    ! local send buffer, note: max size is (blocksize)*(ghost nodes size + 1)*(number of datafields)
    ! restricted/predicted data buffer
    real(kind=rk), allocatable :: data_buffer(:), res_pre_data(:,:,:,:)
    ! data buffer size
    integer(kind=ik)                        :: buffer_size

    ! grid parameter
    integer(kind=ik)                                :: Bs, g
    ! number of datafields
    integer(kind=ik)                                :: NdF


!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes
    ! number of datafields
    NdF   = params%number_data_fields

    ! set number of blocks
    N = params%number_blocks

    ! set MPI parameter
    rank = params%rank

    ! reset status
    stop_status = .false.

    ! set loop number for 2D/3D case
    neighbor_num = size(hvy_neighbor, 2)

    ! 'exclude_redundant', 'include_redundant', 'only_redundant'
    data_bounds_type = 'include_redundant'

    ! 2D only!
    allocate( data_buffer( (Bs+g)*(g+1)*NdF ), res_pre_data( Bs+2*g, Bs+2*g, Bs+2*g, NdF) )

!---------------------------------------------------------------------------------------------
! main body

    ! loop over active heavy data
    do k = 1, hvy_n
        ! loop over all neighbors
        do neighborhood = 1, neighbor_num
            ! neighbor exists
            if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                ! 0. ids bestimmen
                ! neighbor light data id
                neighbor_light_id = hvy_neighbor( hvy_active(k), neighborhood )
                ! calculate light id
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )
                ! calculate the difference between block levels
                ! define leveldiff: sender - receiver, so +1 means sender on higher level
                ! sender is active block (me)
                level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_light_id, params%max_treelevel+1 )

                ! 1. ich (aktiver block) ist der sender für seinen nachbarn
                ! lese daten und sortiere diese in bufferform
                ! wird auch für interne nachbarn gemacht, um gleiche routine für intern/extern zu verwenden
                ! um diue lesbarkeit zu erhöhen werden zunächst die datengrenzen bestimmt
                ! diese dann benutzt um die daten zu lesen
                ! 2D/3D wird bei der datengrenzbestimmung unterschieden, so dass die tatsächliche leseroutine stark vereinfacht ist
                ! da die interpolation bei leveldiff -1 erst bei der leseroutine stattfindet, werden als datengrenzen die für die interpolation noitwendigen bereiche angegeben
                ! auch für restriction ist der datengrenzenbereich größer, da dann auch hier später erst die restriction stattfindet
                call calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, 'sender' )

                ! vor dem schreiben der daten muss ggf interpoliert werden
                ! hier werden die datengrenzen ebenfalls angepasst
                ! interpolierte daten stehen in einem extra array
                ! dessen größe richtet sich nach dem größten möglichen interpolationsgebiet: (Bs+2*g)^3
                ! auch die vergröberten daten werden in den interpolationbuffer geschrieben und die datengrenzen angepasst
                if ( level_diff /= 0 ) then
                    call restrict_predict_data( params, res_pre_data, data_bounds, neighborhood, level_diff, data_bounds_type, hvy_block, hvy_active(k) )
                end if

                ! lese nun mit den datengrenzen die daten selbst
                ! die gelesenen daten werden als buffervektor umsortiert
                ! so können diese danach entweder in den buffer geschrieben werden oder an die schreiberoutine weitergegeben werden
                call read_hvy_data( params, data_buffer, buffer_size, data_bounds, level_diff, hvy_block, hvy_active(k) )

                ! daten werden jetzt entweder in den speicher geschrieben -> schreiberoutine
                ! oder in den send buffer geschrieben
                ! schreiberoutine erhält die date grenzen
                ! diese werden vorher durch erneuten calc data bounds aufruf berechnet
                ! achtung: die nachbarschaftsbeziehung wird hier wie eine interner Kopieren ausgewertet
                ! invertierung der nachbarschaftsbeziehung findet beim füllen desn sendbuffer statt
                if ( rank == neighbor_rank ) then
                call calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, 'receiver' )

                !

            end if
        end do
    end do

    ! clean up
    deallocate( data_buffer, res_pre_data )

end subroutine check_redundant_nodes

!############################################################################################################

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
    character(len=25), intent(in)                   :: data_bounds_type
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

    if ( data_bounds_type == 'exclude_redundant' ) then
        sh_start = 1
    end if
    if ( data_bounds_type == 'only_redundant' ) then
        sh_end = -g
    end if

    ! reset data bounds
    data_bounds = 0

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    select case(sender_or_receiver)

        case('sender')

            if ( params%threeD_case ) then
                ! 3D

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

    end select

end subroutine calc_data_bounds

!############################################################################################################

subroutine restrict_predict_data( params, res_pre_data, data_bounds, neighborhood, level_diff, data_bounds_type, hvy_block, hvy_id )

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
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff
    ! data_bounds_type
    character(len=25), intent(in)                   :: data_bounds_type
    !> heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, iN, jN, kN

    ! grid parameter
    integer(kind=ik)                                :: Bs, g

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    ! data size
    iN = data_bounds(2,1) - data_bounds(1,1) + 1
    jN = data_bounds(2,2) - data_bounds(1,2) + 1
    kN = data_bounds(2,3) - data_bounds(1,3) + 1

!---------------------------------------------------------------------------------------------
! main body

    if ( params%threeD_case ) then
        ! 3D

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
                    do dF = 1, params%number_data_fields
                        ! interpolate data
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1, dF), params%order_predictor)
                    end do
                    ! reset data bounds
                    select case(data_bounds_type)
                        case('exclude_redundant')
                            data_bounds(1,1:2) = 2
                            data_bounds(2,1:2) = g+1

                        case('include_redundant')
                            data_bounds(1,1:2) = 1
                            data_bounds(2,1:2) = g+1

                        case('only_redundant')
                            data_bounds(1:2,1:2) = 1

                    end select

                elseif ( level_diff == 1) then
                    ! loop over all data fields
                    do dF = 1, params%number_data_fields
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

            ! 'NNE' 'NNW' 'SSE' 'SSW' ENE' 'ESE' 'WNW' 'WSW'
            case(9,10,11,12,13,14,15,16)
                if ( level_diff == -1 ) then
                    ! loop over all data fields
                    do dF = 1, params%number_data_fields
                        ! interpolate data
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF, hvy_id ), &
                        res_pre_data( 1:iN*2-1, 1:jN*2-1, 1, dF), params%order_predictor)
                    end do
                    ! reset data bounds
                    select case(neighborhood)
                        ! 'NNE'
                        case(9)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        ! 'NNW'
                        case(10)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                            end select

                        ! 'SSE'
                        case(11)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = Bs+g+1
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case('include_redundant')
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case('only_redundant')
                                    data_bounds(1:2,1) = Bs+g
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        ! 'SSW'
                        case(12)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = Bs+g+1
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case('include_redundant')
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case('only_redundant')
                                    data_bounds(1:2,1) = Bs+g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                            end select

                        ! 'ENE'
                        case(13)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = Bs+g+1
                                    data_bounds(2,2) = Bs+2*g-1

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g-1

                                case('only_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1:2,2) = Bs+g

                            end select

                        ! 'ESE'
                        case(14)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = Bs+g+1
                                    data_bounds(2,2) = Bs+2*g-1

                                case('include_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g-1

                                case('only_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1:2,2) = Bs+g

                            end select

                        ! 'WNW'
                        case(15)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1:2,2) = Bs+g

                            end select

                        ! 'WSW'
                        case(16)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case('include_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1:2,2) = Bs+g

                            end select

                        end select

                elseif ( level_diff == 1 ) then
                    ! loop over all data fields
                    do dF = 1, params%number_data_fields
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


!############################################################################################################

subroutine read_hvy_data( params, data_buffer, buffer_counter, data_bounds, level_diff, hvy_block, hvy_id )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(out)                 :: data_buffer(:)
    ! buffer size
    integer(kind=ik), intent(out)                 :: buffer_counter
    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds(2,3)
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff
    !> heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, iN, jN, kN

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset buffer size
    buffer_counter = 0

    ! data size
    iN = data_bounds(2,1) - data_bounds(1,1) + 1
    jN = data_bounds(2,2) - data_bounds(1,2) + 1
    kN = data_bounds(2,3) - data_bounds(1,3) + 1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all data fields
    do dF = 1, params%number_data_fields
        ! first dimension
        do i = 1, iN
            ! second dimension
            do j = 1, jN
                ! third dimension, note: for 2D cases kN is allways 1
                do k = 1, kN

                    ! increase buffer size
                    buffer_counter = buffer_counter + 1
                    ! write data buffer
                    data_buffer(buffer_counter)   = hvy_block( i, j, k, dF, hvy_id )

                end do
            end do
        end do
    end do

end subroutine read_hvy_data

