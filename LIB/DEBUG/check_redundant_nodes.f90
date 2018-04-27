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

subroutine check_redundant_nodes( params, lgt_block, hvy_block, hvy_synch, hvy_neighbor,&
     hvy_active, hvy_n, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer, stop_status )

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
    !> heavy synch array
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)

    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)     :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)        :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

    ! MPI parameter
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! loop variables
    integer(kind=ik)                    :: N, k, dF, neighborhood, invert_neighborhood, neighbor_num, level_diff, l

    ! id integers
    integer(kind=ik)                    :: lgt_id, neighbor_light_id, neighbor_rank, hvy_id

    ! type of data bounds
    ! 'exclude_redundant', 'include_redundant', 'only_redundant'
    character(len=25)                   :: data_bounds_type
    integer(kind=ik), dimension(2,3)    :: data_bounds

    ! local send buffer, note: max size is (blocksize)*(ghost nodes size + 1)*(number of datafields)
    ! restricted/predicted data buffer
    real(kind=rk), allocatable :: data_buffer(:), res_pre_data(:,:,:,:)
    ! data buffer size
    integer(kind=ik)                        :: buffer_size, buffer_position

    ! grid parameter
    integer(kind=ik)                                :: Bs, g
    ! number of datafields
    integer(kind=ik)                                :: NdF

    ! type of data writing
    character(len=25)                   :: data_writing_type

    ! communications matrix (only 1 line)
    ! note: todo: check performance without allocation?
    ! todo: remove dummy com matrix, needed for old MPI subroutines
    integer(kind=ik), allocatable     :: com_matrix(:), dummy_matrix(:,:)

    ! position in integer buffer, need for every neighboring process
    integer(kind=ik), allocatable                                :: int_pos(:)

    ! synch stage loop variables
    integer(kind=ik) :: synch_stage, stages
    ! synch status
    ! synch == .true. : active block sends data to neighboring block
    ! neighbor_synch == .true. : neighbor block send data to active block
    logical    :: synch, neighbor_synch, test2

!---------------------------------------------------------------------------------------------
! interfaces

 !---------------------------------------------------------------------------------------------
! variables initialization

    ! hack to use subroutine as redundant nodes test and for ghost nodes synchronization
    if (stop_status) then
        ! synchronization
        ! 'exclude_redundant', 'include_redundant', 'only_redundant'
        data_bounds_type = 'include_redundant'
        ! 'average', 'simple', 'staging', 'compare'
        data_writing_type = 'staging'

    else
        ! nodes test
        ! 'exclude_redundant', 'include_redundant', 'only_redundant'
        data_bounds_type = 'only_redundant'
        ! 'average', 'simple', 'staging', 'compare'
        data_writing_type = 'compare'
        ! reset status
        stop_status = .false.

    end if

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes
    ! number of datafields
    NdF   = params%number_data_fields

    ! set number of blocks
    N = params%number_blocks

    ! set MPI parameter
    rank = params%rank
    number_procs = params%number_procs

    ! set loop number for 2D/3D case
    neighbor_num = size(hvy_neighbor, 2)

    ! 2D only!
    allocate( data_buffer( (Bs+g)*(g+1)*NdF ), res_pre_data( Bs+2*g, Bs+2*g, Bs+2*g, NdF), &
    com_matrix(number_procs), int_pos(number_procs), dummy_matrix(number_procs, number_procs) )

    ! reset ghost nodes for all blocks - for debugging
    ! todo: use reseting subroutine from MPI module
    ! if ( (params%test_ghost_nodes_synch) .and. (data_bounds_type /= 'only_redundant') ) then
    !     !-- x-direction
    !     hvy_block(1:g, :, :, :, : )           = 9.0e9_rk
    !     hvy_block(Bs+g+1:Bs+2*g, :, :, :, : ) = 9.0e9_rk
    !     !-- y-direction
    !     hvy_block(:, 1:g, :, :, : )           = 9.0e9_rk
    !     hvy_block(:, Bs+g+1:Bs+2*g, :, :, : ) = 9.0e9_rk
    !     !-- z-direction
    !     if ( params%threeD_case ) then
    !         hvy_block(:, :, 1:g, :, : )           = 9.0e9_rk
    !         hvy_block(:, :, Bs+g+1:Bs+2*g, :, : ) = 9.0e9_rk
    !     end if
    ! end if

    ! reset com matrix
    com_matrix = 0
    dummy_matrix = 0

    ! reseting all ghost nodes to zero
    if ( (data_writing_type == 'average') .and. (data_bounds_type /= 'only_redundant') ) then
        do k = 1, hvy_n
            !-- x-direction
            hvy_block(1:g, :, :, :, hvy_active(k) )               = 0.0_rk
            hvy_block(Bs+g+1:Bs+2*g, :, :, :, hvy_active(k) )     = 0.0_rk
            !-- y-direction
            hvy_block(:, 1:g, :, :, hvy_active(k) )               = 0.0_rk
            hvy_block(:, Bs+g+1:Bs+2*g, :, :, hvy_active(k) )     = 0.0_rk
            !-- z-direction
            if ( params%threeD_case ) then
                hvy_block(:, :, 1:g, :, hvy_active(k) )           = 0.0_rk
                hvy_block(:, :, Bs+g+1:Bs+2*g, :, hvy_active(k) ) = 0.0_rk
            end if
        end do
    end if

    ! set number of synch stages
    if ( data_writing_type == 'staging' ) then
        ! all four stages
        ! NOTE: WE CAN ACTUALLY SKIP STAGE 4
        stages = 3
    else
        ! only one stage
        stages = 1
    end if

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all synch stages
    do synch_stage = 1, stages

        ! reset integer send buffer position
        int_pos = 2
        ! reset first in send buffer position
        int_send_buffer( 1, : ) = 0
        int_send_buffer( 2, : ) = -99

        ! loop over active heavy data
        if (data_writing_type=="average") then
        do k = 1, hvy_n

            ! reset synch array
            ! alles auf null, knoten im block auf 1
            ! jeder später gespeicherte knoten erhöht wert um 1
            ! am ende der routine wird der wert aus dem synch array ggf. für die durchschnittsberechnung benutzt
            ! synch array hat die maximale anzahl von blöcken pro prozess alloziiert, so dass die heavy id unverändert
            ! benutzt werden kann
            ! ghost nodes layer auf 1 setzen, wenn nur die redundanten Knoten bearbeitet werden
            if (data_bounds_type == 'only_redundant') then
                hvy_synch(:, :, :, hvy_active(k)) = 1
            else
                hvy_synch(:, :, :, hvy_active(k)) = 0
            end if
            ! alles knoten im block werden auf 1 gesetzt

            ! todo: ist erstmal einfacher als nur die redundaten zu setzen, aber unnötig
            ! so gibt es aber nach der synch keine nullen mehr, kann ggf. als synch test verwendet werden?
            if ( params%threeD_case ) then
                hvy_synch( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, hvy_active(k)) = 1
            else
                hvy_synch( g+1:Bs+g, g+1:Bs+g, 1, hvy_active(k)) = 1
            end if

        end do
        end if

        ! loop over active heavy data
        do k = 1, hvy_n
            ! loop over all neighbors
            do neighborhood = 1, neighbor_num
                ! neighbor exists
                if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                    ! 0. ids bestimmen
                    ! neighbor light data id
                    neighbor_light_id = hvy_neighbor( hvy_active(k), neighborhood )
                    ! calculate neighbor rank
                    call lgt_id_to_proc_rank( neighbor_rank, neighbor_light_id, N )
                    ! calculate light id
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )
                    ! neighbor heavy id
                    call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, neighbor_rank, N )
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
                    if ( level_diff == 0 ) then
                        ! lese nun mit den datengrenzen die daten selbst
                        ! die gelesenen daten werden als buffervektor umsortiert
                        ! so können diese danach entweder in den buffer geschrieben werden oder an die schreiberoutine weitergegeben werden
                        ! in die lese routine werden nur die relevanten Daten (data bounds) übergeben
                        call read_hvy_data( params, data_buffer, buffer_size, &
                        hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), data_bounds(1,3):data_bounds(2,3), :, hvy_active(k)) )

                    else
                        ! interpoliere daten
                        call restrict_predict_data( params, res_pre_data, data_bounds, neighborhood, level_diff, data_bounds_type, hvy_block, hvy_active(k) )
                        ! lese daten, verwende interpolierte daten
                        call read_hvy_data( params, data_buffer, buffer_size, res_pre_data( data_bounds(1,1):data_bounds(2,1), &
                                                                                            data_bounds(1,2):data_bounds(2,2), &
                                                                                            data_bounds(1,3):data_bounds(2,3), &
                                                                                            :) )

                    end if

                    ! daten werden jetzt entweder in den speicher geschrieben -> schreiberoutine
                    ! oder in den send buffer geschrieben
                    ! schreiberoutine erhält die date grenzen
                    ! diese werden vorher durch erneuten calc data bounds aufruf berechnet
                    ! achtung: die nachbarschaftsbeziehung wird hier wie eine interner Kopieren ausgewertet
                    ! invertierung der nachbarschaftsbeziehung findet beim füllen des sendbuffer statt
                    if ( (rank==neighbor_rank).and.(data_writing_type=='simple') ) then
                        ! internal neighbor and direct writing method: copy the ghost nodes as soon as possible, without passing
                        ! via the buffers first.
                        ! data bounds
                        call calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, 'receiver' )
                        ! simply write data. No care
                        call write_hvy_data( params, data_buffer, data_bounds, hvy_block, hvy_id )

                    else
                        ! synch status for staging method
                        synch = .true.
                        if (data_writing_type == 'staging') then
                            call set_synch_status( synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_active(k), neighborhood )
                        end if

                        ! first: fill com matrix, count number of communication to neighboring process, needed for int buffer length
                        com_matrix(neighbor_rank+1) = com_matrix(neighbor_rank+1) + 1

                        if (synch) then
                            ! active block send data to his neighbor block
                            ! fill int/real buffer
                            call write_buffers( int_send_buffer, real_send_buffer, buffer_size, neighbor_rank, data_buffer, int_pos(neighbor_rank+1), hvy_id, neighborhood, level_diff )

                        else
                            ! neighbor block send data to active block
                            ! write -1 to int_send buffer, placeholder
                            int_send_buffer( int_pos(neighbor_rank+1) : int_pos(neighbor_rank+1)+4  , neighbor_rank+1 ) = -1
                        end if

                        ! increase int buffer position
                        int_pos(neighbor_rank+1) = int_pos(neighbor_rank+1) + 5

                        ! markiere das aktuelle ende des buffers, falls weitere elemente dazu kommen, wird die -99 wieder überschrieben
                        int_send_buffer( int_pos(neighbor_rank+1)  , neighbor_rank+1 ) = -99

                    end if

                end if
            end do
        end do

        ! pretend that no communication with myself takes place, in order to skip the
        ! MPI transfer in the following routine. NOTE: you can also skip this step and just have isend_irecv_data_2
        ! transfer the data, in which case you should skip the copy part directly after isend_irecv_data_2
        com_matrix(rank+1) = 0

        !***********************************************************************
        ! transfer part (send/recv)
        !***********************************************************************
        ! send/receive data
        ! note: todo, remove dummy subroutine
        ! note: new dummy subroutine sets receive buffer position accordingly to process rank
        ! note: todo: use more than non-blocking send/receive
        call isend_irecv_data_2( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix  )

        ! fill receive buffer for internal neighbors for averaging writing type
        if ( (data_writing_type == 'average') .or. (data_writing_type == 'compare') .or. (data_writing_type == 'staging') ) then
            ! fill receive buffer
            int_receive_buffer( 1:int_pos(rank+1)  , rank+1 ) = int_send_buffer( 1:int_pos(rank+1)  , rank+1 )
            real_receive_buffer( 1:int_receive_buffer(1,rank+1), rank+1 ) = real_send_buffer( 1:int_receive_buffer(1,rank+1), rank+1 )
            ! change com matrix, need to sort in buffers in next step
            com_matrix(rank+1) = 1
        end if

        !***********************************************************************
        ! Unpack received data in the ghost node layers
        !***********************************************************************
        ! Daten einsortieren
        ! für simple, average, compare: einfach die buffer einsortieren, Reihenfolge ist egal
        ! staging: erneuter loop über alle blöcke und nachbarschaften, wenn daten notwendig, werden diese in den buffern gesucht
        if ( data_writing_type /= 'staging' ) then
            ! sortiere den real buffer ein
            ! loop over all procs
            do k = 1, number_procs
                if ( com_matrix(k) /= 0 ) then
                    ! neighboring proc
                    ! first element in int buffer is real buffer size
                    l = 2
                    ! -99 marks end of data
                    do while ( int_receive_buffer(l, k) /= -99 )

                        ! hvy id
                        hvy_id = int_receive_buffer(l, k)
                        ! neighborhood
                        neighborhood = int_receive_buffer(l+1, k)
                        ! level diff
                        level_diff = int_receive_buffer(l+2, k)
                        ! buffer position
                        buffer_position = int_receive_buffer(l+3, k)
                        ! buffer size
                        buffer_size = int_receive_buffer(l+4, k)

                        ! data buffer
                        data_buffer(1:buffer_size) = real_receive_buffer( buffer_position : buffer_position-1 + buffer_size, k )

                        ! data bounds
                        call calc_data_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, 'receiver' )
                        ! write data, hängt vom jeweiligen Fall ab
                        ! average: schreibe daten, merke Anzahl der geschriebenen Daten, Durchschnitt nach dem Einsortieren des receive buffers berechnet
                        ! simple: schreibe ghost nodes einfach in den speicher (zum Testen?!)
                        ! staging: wende staging konzept an
                        ! compare: vergleiche werte mit vorhandenen werten (nur für redundante knoten sinnvoll, als check routine)
                        select case(data_writing_type)
                            case('simple')
                                ! simply write data
                                call write_hvy_data( params, data_buffer, data_bounds, hvy_block, hvy_id )

                            case('average')
                                ! add data
                                call add_hvy_data( params, data_buffer, data_bounds, hvy_block, hvy_synch, hvy_id )

                            case('compare')
                                ! compare data
                                call compare_hvy_data( params, data_buffer, data_bounds, hvy_block, hvy_id, stop_status, level_diff )

                        end select

                        ! increase buffer postion marker
                        l = l + 5

                    end do
                end if
            end do

            ! last averaging step
            if ( data_writing_type == 'average' ) then
                ! loop over active heavy data
                do k = 1, hvy_n
                    do dF = 1, NdF

                        ! calculate average for all nodes, todo: proof performance?
                        hvy_block(:, :, :, dF, hvy_active(k)) = hvy_block(:, :, :, dF, hvy_active(k)) / real( hvy_synch(:, :, :, hvy_active(k)) , kind=rk)

                    end do
                end do
            end if

        else
            ! staging type
            ! loop over active heavy data
            do k = 1, hvy_n
                ! loop over all neighbors
                do neighborhood = 1, neighbor_num
                    ! neighbor exists
                    if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                        ! invert neighborhood, needed for in buffer searching, because sender proc has invert neighborhood relation
                        call calc_invert_neighborhood( params, neighborhood, invert_neighborhood)

                        ! 0. ids bestimmen
                        ! neighbor light data id
                        neighbor_light_id = hvy_neighbor( hvy_active(k), neighborhood )
                        ! calculate neighbor rank
                        call lgt_id_to_proc_rank( neighbor_rank, neighbor_light_id, N )
                        ! calculate light id
                        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )
                        ! calculate the difference between block levels
                        ! define leveldiff: sender - receiver, so +1 means sender on higher level
                        ! sender is active block (me)
                        level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_light_id, params%max_treelevel+1 )

                        ! set synch status
                        call set_synch_status( synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_active(k), neighborhood )
                        ! synch == .true. bedeutet, dass der aktive block seinem nachbarn daten gibt
                        ! hier sind wir aber auf der seite des empfängers, das bedeutet, neighbor_synch muss ausgewertet werden

                        if (neighbor_synch) then

                            ! search buffers for synchronized data
                            ! first element in int buffer is real buffer size
                            l = 2

                            ! -99 marks end of data
                            test2 = .false.
                            do while ( int_receive_buffer(l, neighbor_rank+1) /= -99 )

                                ! proof heavy id and neighborhood id
                                if (  (int_receive_buffer( l,   neighbor_rank+1 ) == hvy_active(k) ) &
                                .and. (int_receive_buffer( l+1, neighbor_rank+1 ) == invert_neighborhood) ) then

                                    ! set parameter
                                    ! level diff, read from buffer because calculated level_diff is not sender-receiver
                                    level_diff = int_receive_buffer(l+2, neighbor_rank+1)
                                    ! buffer position
                                    buffer_position = int_receive_buffer(l+3, neighbor_rank+1)
                                    ! buffer size
                                    buffer_size = int_receive_buffer(l+4, neighbor_rank+1)

                                    ! data buffer
                                    data_buffer(1:buffer_size) = real_receive_buffer( buffer_position : buffer_position-1 + buffer_size, neighbor_rank+1 )

                                    ! data bounds
                                    call calc_data_bounds( params, data_bounds, invert_neighborhood, level_diff, data_bounds_type, 'receiver' )

                                    ! write data
                                    call write_hvy_data( params, data_buffer(1:buffer_size), data_bounds, hvy_block, hvy_active(k) )
                                    ! done, exit the while loop?
                                    test2=.true.
                                    exit
                                end if

                                ! increase buffer postion marker
                                l = l + 5

                            end do
                            if (test2 .eqv. .false.) call abort(777771,"not found")

                        end if

                    end if
                end do
            end do

        end if

    end do ! loop over stages

    if ( data_writing_type=='compare' ) then
        test2 = stop_status
        call MPI_Allreduce(test2, stop_status, 1,MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, k )
    endif

    ! clean up
    deallocate( data_buffer, res_pre_data, com_matrix, int_pos, dummy_matrix )

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
    !! level_diff = ME -  NEIGHBOR from SENDER point of view.
    !! +1 : I'm finer, my neighbor is coarser
    !! -1 : My neighbor is finer, I am coarser
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
    data_bounds = 1

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
                    select case(neighborhood)
                        ! '_NE'
                        case(5)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-2

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 2*g-1
                            end select
                        ! '_NW'
                        case(6)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 2
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_SE'
                        case(7)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-2
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-2

                                case('include_redundant')
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-1
                                    data_bounds(1,2) = g-1
                                    data_bounds(2,2) = 2*g-1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 2*g-1
                                    data_bounds(1:2,2) = 2*g-1
                            end select
                        ! '_SW'
                        case(8)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-2
                                    data_bounds(1,2) = 2
                                    data_bounds(2,2) = g+1

                                case('include_redundant')
                                    data_bounds(1,1) = g-1
                                    data_bounds(2,1) = 2*g-1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 2*g-1
                                    data_bounds(1:2,2) = 1
                            end select
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
                    select case(neighborhood)
                        ! '_NE'
                        case(5)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_NW'
                        case(6)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_SE'
                        case(7)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                        ! '_SW'
                        case(8)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = g+1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = g+1

                                case('only_redundant')
                                    data_bounds(1:2,1) = 1
                                    data_bounds(1:2,2) = 1
                            end select
                    end select

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
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case('include_redundant')
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                                case('only_redundant')
                                    data_bounds(1:2,1) = Bs+2*g
                                    data_bounds(1,2) = g+1
                                    data_bounds(2,2) = Bs+2*g

                            end select

                        ! 'SSW'
                        case(12)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g-1
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case('include_redundant')
                                    data_bounds(1,1) = Bs+g
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                                case('only_redundant')
                                    data_bounds(1:2,1) = Bs+2*g
                                    data_bounds(1,2) = 1
                                    data_bounds(2,2) = Bs+g

                            end select

                        ! 'ENE'
                        case(13)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g-1

                                case('include_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g

                                case('only_redundant')
                                    data_bounds(1,1) = 1
                                    data_bounds(2,1) = Bs+g
                                    data_bounds(1:2,2) = Bs+2*g

                            end select

                        ! 'ESE'
                        case(14)
                            select case(data_bounds_type)
                                case('exclude_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g-1

                                case('include_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1,2) = Bs+g
                                    data_bounds(2,2) = Bs+2*g

                                case('only_redundant')
                                    data_bounds(1,1) = g+1
                                    data_bounds(2,1) = Bs+2*g
                                    data_bounds(1:2,2) = Bs+2*g

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
                                    data_bounds(1:2,2) = 1

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
                                    data_bounds(1:2,2) = 1

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

subroutine read_hvy_data( params, data_buffer, buffer_counter, hvy_data )

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
    !> heavy block data, all data fields
    real(kind=rk), intent(in)                       :: hvy_data(:, :, :, :)

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset buffer size
    buffer_counter = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all data fields
    do dF = 1, params%number_data_fields
        ! first dimension
        do i = 1, size(hvy_data, 1)
            ! second dimension
            do j = 1, size(hvy_data, 2)
                ! third dimension, note: for 2D cases kN is allways 1
                do k = 1, size(hvy_data, 3)

                    ! increase buffer size
                    buffer_counter = buffer_counter + 1
                    ! write data buffer
                    data_buffer(buffer_counter)   = hvy_data( i, j, k, dF )

                end do
            end do
        end do
    end do

end subroutine read_hvy_data

!############################################################################################################

subroutine write_hvy_data( params, data_buffer, data_bounds, hvy_block, hvy_id )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(in)                 :: data_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                       :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    buffer_i = 1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all data fields
    do dF = 1, params%number_data_fields
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases kN is allways 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    ! write data buffer
                    hvy_block( i, j, k, dF, hvy_id ) = data_buffer( buffer_i )
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    end do

end subroutine write_hvy_data

!############################################################################################################

subroutine add_hvy_data( params, data_buffer, data_bounds, hvy_block, hvy_synch, hvy_id )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(in)                 :: data_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                       :: hvy_block(:, :, :, :, :)
    !> heavy synch array
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    buffer_i = 1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all data fields
    do dF = 1, params%number_data_fields
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases kN is allways 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    ! write data buffer
                    hvy_block( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id ) + data_buffer( buffer_i )

                    ! count synchronized data
                    ! note: only for first datafield
                    if (dF==1) hvy_synch( i, j, k, hvy_id ) = hvy_synch( i, j, k, hvy_id ) + 1

                    ! increase buffer counter
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    end do

end subroutine add_hvy_data

!############################################################################################################

subroutine compare_hvy_data( params, data_buffer, data_bounds, hvy_block, hvy_id, stop_status, level_diff )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(in)                 :: data_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                       :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id, level_diff
    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i, oddeven

    ! error threshold
    real(kind=rk)                                   :: eps

    ! error norm
    real(kind=rk)       :: error_norm

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    buffer_i = 1

    ! set error threshold
    eps = 1e-9_rk

    ! reset error norm
    error_norm = 0.0_rk

    ! the first index of the redundant points is (g+1, g+1, g+1)
    ! so if g is even, then we must compare the odd indices i,j,k on the lines
    ! of the redundant points.
    ! if g is odd, then we must compare the even ones
    ! Further note that BS is odd (always), so as odd+even=odd and odd+odd=even
    ! we can simply study the parity of g
    oddeven = mod(params%number_ghost_nodes,2)

!---------------------------------------------------------------------------------------------
! main body
    ! loop over all data fields
    do dF = 1, params%number_data_fields
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases k is always 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    if (level_diff/=-1) then
                        ! pointwise error norm
                        error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - data_buffer( buffer_i )))
                    else
                        ! if the level diff is -1, I compare with interpolated (upsampled) data. that means every EVEN
                        ! point is the result of interpolation, and not truely redundant.
                        ! Note this routine ALWAYS just compares the redundant nodes, so it will mostly be called
                        ! with a line of points (i.e. one dimension is length one)
                        !
                        ! This routine has been tested:
                        !   - old method (working version): no error found (okay)
                        !   - old method, non_uniform_mesh_correction=0; in params file -> plenty of errors (okay)
                        !   - old method, sync stage 4 deactivated: finds all occurances of "3finer blocks on corner problem" (okay)
                        !   - new method, averaging, no error found (makes sense: okay)
                        if (oddeven==0) then
                            ! even number of ghost nodes
                            if (((data_bounds(2,1)- data_bounds(1,1) == 0).and.(mod(j,2)/=0)) .or. ((data_bounds(2,2)-data_bounds(1,2) == 0).and.(mod(i,2)/=0))) then
                                error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - data_buffer( buffer_i )))
                            endif
                        else
                            ! odd number of ghost nodes
                            if (((data_bounds(2,1)- data_bounds(1,1) == 0).and.(mod(j,2)==0)) .or. ((data_bounds(2,2)-data_bounds(1,2) == 0).and.(mod(i,2)==0))) then
                                error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - data_buffer( buffer_i )))
                            endif
                        endif
                    endif
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    end do

    if (error_norm > eps)  then
        ! error message
        write(*,*) "ERROR: difference in redundant nodes", error_norm, level_diff, hvy_id, params%rank
        ! stop program
        stop_status = .true.
        ! mark block by putting a dot in the middle. this way, we can identify all
        ! blocks that have problems. If we set the entire block to 100, then subsequent
        ! synchronizations fail because of that. If we set just a point in the middle, it
        ! is far away from the boundary and thus does not affect other sync steps.
        hvy_block( size(hvy_block,1)/2, size(hvy_block,2)/2, :, :, hvy_id ) = 100.0_rk
    end if

end subroutine compare_hvy_data

!############################################################################################################

subroutine isend_irecv_data_2( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> send/receive buffer, integer and real
    integer(kind=ik), intent(in)        :: int_send_buffer(:,:)
    integer(kind=ik), intent(out)       :: int_receive_buffer(:,:)

    real(kind=rk), intent(in)           :: real_send_buffer(:,:)
    real(kind=rk), intent(out)          :: real_receive_buffer(:,:)

    !> communications matrix: neighboring proc rank
    !> com matrix pos: position in send buffer
    integer(kind=ik), intent(in)        :: com_matrix(:)

    ! process rank
    integer(kind=ik)                    :: rank
    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! MPI status
    !integer                             :: status(MPI_status_size)

    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI request
    integer(kind=ik)                    :: send_request(size(com_matrix,1)), recv_request(size(com_matrix,1))

    ! column number of send buffer, column number of receive buffer, real data buffer length
    integer(kind=ik)                    :: real_pos, int_length

    ! loop variable
    integer(kind=ik)                    :: k, i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank            = params%rank
    number_procs    = params%number_procs

    ! set message tag
    tag = 0

!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! first: integer data

    ! reset communication counter
    i = 0

    ! reset request arrays
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! loop over com matrix
    do k = 1, number_procs

        ! communication between proc rank and proc k-1
        if ( com_matrix(k) > 0 ) then

            ! legth of integer buffer
            int_length = 5*com_matrix(k) + 3

            ! increase communication counter
            i = i + 1

            ! tag
            tag = rank+1+k

            ! receive data
            call MPI_Irecv( int_receive_buffer(1, k), int_length, MPI_INTEGER4, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)

            ! send data
            call MPI_Isend( int_send_buffer(1, k), int_length, MPI_INTEGER4, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)

        end if

    end do

    !> \todo Please check if waiting twice is really necessary
    ! synchronize non-blocking communications
    ! note: single status variable do not work with all compilers, so use MPI_STATUSES_IGNORE instead
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
    end if

    ! ----------------------------------------------------------------------------------------
    ! second: real data

    ! reset communication couter
    i = 0

    ! reset request arrays
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! loop over corresponding com matrix line
    do k = 1, number_procs

        ! communication between proc rank and proc k-1
        if ( com_matrix(k) > 0 ) then

            ! increase communication counter
            i = i + 1

            tag = number_procs*10*(rank+1+k)

            ! real buffer length
            real_pos = int_receive_buffer(1, k)

            ! receive data
            call MPI_Irecv( real_receive_buffer(1, k), real_pos, MPI_REAL8, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)

            ! real buffer length
            real_pos = int_send_buffer(1, k)

            ! send data
            call MPI_Isend( real_send_buffer(1, k), real_pos, MPI_REAL8, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)

        end if

    end do

    ! synchronize non-blocking communications
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
    end if

end subroutine isend_irecv_data_2

!############################################################################################################

subroutine set_synch_status( synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_id, neighborhood )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! synch stage
    integer(kind=ik), intent(in)        :: synch_stage

    ! synch status
    logical, intent(inout)    :: synch, neighbor_synch

    ! level diffenrence
    integer(kind=ik), intent(in)        :: level_diff

    ! heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    ! list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_id

    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! set synch stage
    ! stage 1: level +1
    ! stage 2: level 0
    ! stage 3: level -1
    ! stage 4: special
    synch = .false.
    neighbor_synch = .false.

    ! stage 1
    if ( (synch_stage == 1) .and. (level_diff == 1) ) then
        ! block send data
        synch = .true.
    elseif ( (synch_stage == 1) .and. (level_diff == -1) ) then
        ! neighbor send data
        neighbor_synch = .true.
    end if

    ! stage 2
    if ( (synch_stage == 2) .and. (level_diff == 0) ) then
        ! block send data
        synch = .true.
        ! neighbor send data
        neighbor_synch = .true.
    end if

    ! stage 3
    if ( (synch_stage == 3) .and. (level_diff == -1) ) then
        ! block send data
        synch = .true.
    elseif ( (synch_stage == 3) .and. (level_diff == 1) ) then
        ! neighbor send data
        neighbor_synch = .true.
    end if

    ! stage 4
    if ( (synch_stage == 4) .and. (level_diff == 0) ) then
        ! neighborhood NE
        if ( neighborhood == 5 ) then
            if ( (hvy_neighbor( hvy_id, 9) /= -1) .or. (hvy_neighbor( hvy_id, 13) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
        ! neighborhood NW
        if ( neighborhood == 6 ) then
            if ( (hvy_neighbor( hvy_id, 10) /= -1) .or. (hvy_neighbor( hvy_id, 15) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
        ! neighborhood SE
        if ( neighborhood == 7 ) then
            if ( (hvy_neighbor( hvy_id, 11) /= -1) .or. (hvy_neighbor( hvy_id, 14) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
        ! neighborhood SW
        if ( neighborhood == 8 ) then
            if ( (hvy_neighbor( hvy_id, 12) /= -1) .or. (hvy_neighbor( hvy_id, 16) /= -1) ) then
                synch = .true.
                neighbor_synch = .true.
            end if
        end if
    end if

end subroutine set_synch_status

!############################################################################################################

subroutine write_buffers( int_send_buffer, real_send_buffer, buffer_size, neighbor_rank, data_buffer, int_pos, hvy_id, neighborhood, level_diff )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> send buffers, integer and real
    integer(kind=ik), intent(inout)        :: int_send_buffer(:,:)
    real(kind=rk), intent(inout)           :: real_send_buffer(:,:)

    ! data buffer size
    integer(kind=ik), intent(in)           :: buffer_size

    ! id integer
    integer(kind=ik), intent(in)           :: neighbor_rank

    ! restricted/predicted data buffer
    real(kind=rk), intent(in)              :: data_buffer(:)

    ! integer buffer position
    integer(kind=ik), intent(in)           :: int_pos

    ! data buffer intergers, receiver heavy id, neighborhood id, level diffenrence
    integer(kind=ik), intent(in)           :: hvy_id, neighborhood, level_diff

    ! buffer position
    integer(kind=ik)                       :: buffer_position

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! fill real buffer
    ! position in real buffer is stored in int buffer
    buffer_position = int_send_buffer(1  , neighbor_rank+1 ) + 1

    ! real data
    real_send_buffer( buffer_position : buffer_position-1 + buffer_size, neighbor_rank+1 ) = data_buffer(1:buffer_size)

    ! fill int buffer
    ! sum size of single buffers on first element
    int_send_buffer(1  , neighbor_rank+1 ) = int_send_buffer(1  , neighbor_rank+1 ) + buffer_size

    ! save: neighbor id, neighborhood, level diffenrence, buffer size
    int_send_buffer( int_pos,   neighbor_rank+1 ) = hvy_id
    int_send_buffer( int_pos+1, neighbor_rank+1 ) = neighborhood
    int_send_buffer( int_pos+2, neighbor_rank+1 ) = level_diff
    int_send_buffer( int_pos+3, neighbor_rank+1 ) = buffer_position
    int_send_buffer( int_pos+4, neighbor_rank+1 ) = buffer_size

end subroutine write_buffers

!############################################################################################################

subroutine calc_invert_neighborhood( params, neighborhood, invert_neighborhood )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> neighborhood id
    integer(kind=ik), intent(in)                 :: neighborhood
    !> invert neighborhood id
    integer(kind=ik), intent(out)                       :: invert_neighborhood

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! neighborhood
    if ( params%threeD_case ) then

        select case (neighborhood)
            !'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___'
            case(1)
                invert_neighborhood = 6
            case(2)
                invert_neighborhood = 4
            case(3)
                invert_neighborhood = 5
            case(4)
                invert_neighborhood = 2
            case(5)
                invert_neighborhood = 3
            case(6)
                invert_neighborhood = 1
            !'_12/___', '_13/___', '_14/___', '_15/___'
            case(7)
                invert_neighborhood = 13
            case(8)
                invert_neighborhood = 14
            case(9)
                invert_neighborhood = 11
            case(10)
                invert_neighborhood = 12
            !'_62/___', '_63/___', '_64/___', '_65/___'
            case(11)
                invert_neighborhood = 9
            case(12)
                invert_neighborhood = 10
            case(13)
                invert_neighborhood = 7
            case(14)
                invert_neighborhood = 8
            !'_23/___', '_25/___'
            case(15)
                invert_neighborhood = 18
            case(16)
                invert_neighborhood = 17
            !'_43/___', '_45/___'
            case(17)
                invert_neighborhood = 16
            case(18)
                invert_neighborhood = 15
            !'123/___', '134/___', '145/___', '152/___'
            case(19)
                invert_neighborhood = 25
            case(20)
                invert_neighborhood = 26
            case(21)
                invert_neighborhood = 23
            case(22)
                invert_neighborhood = 24
            !'623/___', '634/___', '645/___', '652/___'
            case(23)
                invert_neighborhood = 21
            case(24)
                invert_neighborhood = 22
            case(25)
                invert_neighborhood = 19
            case(26)
                invert_neighborhood = 20
            !'__1/123', '__1/134', '__1/145', '__1/152'
            case(27)
                invert_neighborhood = 47
            case(28)
                invert_neighborhood = 48
            case(29)
                invert_neighborhood = 49
            case(30)
                invert_neighborhood = 50
            !'__2/123', '__2/623', '__2/152', '__2/652'
            case(31)
                invert_neighborhood = 39
            case(32)
                invert_neighborhood = 40
            case(33)
                invert_neighborhood = 41
            case(34)
                invert_neighborhood = 42
            !'__3/123', '__3/623', '__3/134', '__3/634'
            case(35)
                invert_neighborhood = 45
            case(36)
                invert_neighborhood = 46
            case(37)
                invert_neighborhood = 43
            case(38)
                invert_neighborhood = 44
            !'__4/134', '__4/634', '__4/145', '__4/645'
            case(39)
                invert_neighborhood = 31
            case(40)
                invert_neighborhood = 32
            case(41)
                invert_neighborhood = 33
            case(42)
                invert_neighborhood = 34
            !'__5/145', '__5/645', '__5/152', '__5/652'
            case(43)
                invert_neighborhood = 37
            case(44)
                invert_neighborhood = 38
            case(45)
                invert_neighborhood = 35
            case(46)
                invert_neighborhood = 36
            !'__6/623', '__6/634', '__6/645', '__6/652'
            case(47)
                invert_neighborhood = 27
            case(48)
                invert_neighborhood = 28
            case(49)
                invert_neighborhood = 29
            case(50)
                invert_neighborhood = 30
            !'_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152'
            case(51)
                invert_neighborhood = 63
            case(52)
                invert_neighborhood = 64
            case(53)
                invert_neighborhood = 66!65
            case(54)
                invert_neighborhood = 65!66
            case(55)
                invert_neighborhood = 59
            case(56)
                invert_neighborhood = 60
            case(57)
                invert_neighborhood = 62
            case(58)
                invert_neighborhood = 61
            !'_62/623', '_62/652', '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652'
            case(59)
                invert_neighborhood = 55
            case(60)
                invert_neighborhood = 56
            case(61)
                invert_neighborhood = 58
            case(62)
                invert_neighborhood = 57
            case(63)
                invert_neighborhood = 51
            case(64)
                invert_neighborhood = 52
            case(65)
                invert_neighborhood = 54!53
            case(66)
                invert_neighborhood = 53!54
            !'_23/123', '_23/623', '_25/152', '_25/652'
            case(67)
                invert_neighborhood = 73
            case(68)
                invert_neighborhood = 74
            case(69)
                invert_neighborhood = 71
            case(70)
                invert_neighborhood = 72
            !'_43/134', '_43/634', '_45/145', '_45/645'
            case(71)
                invert_neighborhood = 69
            case(72)
                invert_neighborhood = 70
            case(73)
                invert_neighborhood = 67
            case(74)
                invert_neighborhood = 68
        end select
    else

        select case (neighborhood)
            case(1)
                invert_neighborhood = 3
            case(2)
                invert_neighborhood = 4
            case(3)
                invert_neighborhood = 1
            case(4)
                invert_neighborhood = 2
            case(5)
                invert_neighborhood = 8
            case(6)
                invert_neighborhood = 7
            case(7)
                invert_neighborhood = 6
            case(8)
                invert_neighborhood = 5
            case(9)
                invert_neighborhood = 11
            case(10)
                invert_neighborhood = 12
            case(11)
                invert_neighborhood = 9
            case(12)
                invert_neighborhood = 10
            case(13)
                invert_neighborhood = 15
            case(14)
                invert_neighborhood = 16
            case(15)
                invert_neighborhood = 13
            case(16)
                invert_neighborhood = 14
        end select
    end if

end subroutine calc_invert_neighborhood

!############################################################################################################
