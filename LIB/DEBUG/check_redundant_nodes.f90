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
     hvy_active, hvy_n, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer, stop_status, synch_method )

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

    ! synch method
    integer(kind=ik), intent(in)                   :: synch_method

    ! MPI parameter
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! loop variables
    integer(kind=ik)                    :: ix, iy, iz, N, k, dF, neighborhood, invert_neighborhood, neighbor_num, level_diff, l

    ! id integers
    integer(kind=ik)                    :: lgt_id, neighbor_light_id, neighbor_rank, hvy_id

    ! type of data bounds
    ! 'exclude_redundant', 'include_redundant', 'only_redundant'
    character(len=25)                   :: data_bounds_type
    integer(kind=ik), dimension(2,3)    :: data_bounds
    integer(kind=ik)                    :: data_size_x, data_size_y, data_size_z

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
    logical    :: synch, neighbor_synch

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! timing
    sub_t0 = MPI_Wtime()
    ! init timing for proc without blocks
    call toc( params, "---synchronize ghost: fill buffers", MPI_wtime()-sub_t0, .true. )
    call toc( params, "---synchronize ghost: send data", MPI_wtime()-sub_t0, .true. )
    call toc( params, "---synchronize ghost: write data", MPI_wtime()-sub_t0, .true. )

    ! 'exclude_redundant', 'include_redundant', 'only_redundant'
    ! 'average', 'simple', 'staging', 'compare'
    select case(synch_method)
        ! check redundant nodes
        case(0)
            data_bounds_type = 'only_redundant'
            data_writing_type = 'compare'
            stop_status = .false.
            stages = 1

        ! old staging
        case(1)
            data_bounds_type = 'include_redundant'
            data_writing_type = 'staging'
            stages = 4

        ! new staging
        ! todo: remove stage 4, switch stage 2 with 3, but needs additional redundant nodes correction

        ! averaging
        case(3)
            data_bounds_type = 'include_redundant'
            data_writing_type = 'average'
            stages = 1

        ! avarage redundant nodes
        case(9)
            data_bounds_type = 'only_redundant'
            data_writing_type = 'average'
            stages = 1

        ! correct redundant nodes, if Jmax is reached (stage0)
        case(11)
            data_bounds_type = 'only_redundant'
            data_writing_type = 'staging'
            stages = 1

        ! error case
        case default
             write(*,'(80("_"))')
             write(*,*) "ERROR: ghost nodes synchronization method is unknown"
             stop

    end select

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
    allocate( data_buffer( (Bs+2*g)*(Bs+2*g)*(Bs+2*g)*NdF ), res_pre_data( Bs+2*g, Bs+2*g, Bs+2*g, NdF), &
    com_matrix(number_procs), int_pos(number_procs), dummy_matrix(number_procs, number_procs) )

    ! reset ghost nodes for all blocks - for debugging
    ! todo: use reseting subroutine from MPI module
    if ( (params%test_ghost_nodes_synch) .and. (data_bounds_type /= 'only_redundant') ) then
        !-- x-direction
        hvy_block(1:g, :, :, :, : )           = 55.0e9_rk
        hvy_block(Bs+g+1:Bs+2*g, :, :, :, : ) = 55.0e9_rk
        !-- y-direction
        hvy_block(:, 1:g, :, :, : )           = 55.0e9_rk
        hvy_block(:, Bs+g+1:Bs+2*g, :, :, : ) = 55.0e9_rk
        !-- z-direction
        if ( params%threeD_case ) then
            hvy_block(:, :, 1:g, :, : )           = 55.0e9_rk
            hvy_block(:, :, Bs+g+1:Bs+2*g, :, : ) = 55.0e9_rk
        end if
    end if

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

!---------------------------------------------------------------------------------------------
! main body

    ! only work, if proc has blocks
    if (hvy_n/=0) then

        ! loop over all synch stages
        do synch_stage = 1, stages

            ! timing
            sub_t0 = MPI_Wtime()

            ! old staging, all 4 stages
            if ( synch_method == 1 ) then
                ! stage 3: coarse to fine, stage 4: correction step
                if (synch_stage==3) then
                    data_bounds_type = 'exclude_redundant'
                end if
                if (synch_stage==4) then
                    data_bounds_type = 'include_redundant'
                end if

!            elseif ( data_writing_type == 'staging_new' ) then
!                ! stage 2: coarse to fine, stage 3: level_diff=0
!                if (synch_stage==2) then
!                    data_bounds_type = 'exclude_redundant'
!                end if
!                if (synch_stage==3) then
!                    data_bounds_type = 'include_redundant'
!                end if
            end if

            ! reset integer send buffer position
            int_pos = 2
            ! reset first in send buffer position
            int_send_buffer( 1, : ) = 0

            ! nur wenn 'average' gemacht wird ist der nächste schritt notwendig
            if ( data_writing_type == 'average' ) then
                ! loop over active heavy data
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
                            data_size_x = data_bounds(2,1)-data_bounds(1,1)+1
                            data_size_y = data_bounds(2,2)-data_bounds(1,2)+1
                            data_size_z = data_bounds(2,3)-data_bounds(1,3)+1

                            buffer_size = data_size_x*data_size_y*data_size_z

                            call read_hvy_data( NdF, data_size_x, data_size_y, data_size_z, &
                                                     data_buffer( 1:data_size_x*data_size_y*data_size_z ), &
                                                     hvy_block( data_bounds(1,1):data_bounds(2,1), &
                                                                data_bounds(1,2):data_bounds(2,2), &
                                                                data_bounds(1,3):data_bounds(2,3), &
                                                                :, hvy_active(k)) )

                        else
                            ! interpoliere daten
                            call restrict_predict_data( params, Bs, g, NdF, res_pre_data, data_bounds, neighborhood, level_diff, data_bounds_type, hvy_block(:,:,:,:,hvy_active(k)) )

                            ! lese daten, verwende interpolierte daten
                            data_size_x = data_bounds(2,1)-data_bounds(1,1)+1
                            data_size_y = data_bounds(2,2)-data_bounds(1,2)+1
                            data_size_z = data_bounds(2,3)-data_bounds(1,3)+1

                            buffer_size = data_size_x*data_size_y*data_size_z

                            call read_hvy_data( NdF, data_size_x, data_size_y, data_size_z, &
                                                     data_buffer( 1:data_size_x*data_size_y*data_size_z ), &
                                                     res_pre_data( data_bounds(1,1):data_bounds(2,1), &
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
                        if ( rank == neighbor_rank ) then

                            ! interner nachbar
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
                                    data_size_x = data_bounds(2,1)-data_bounds(1,1)+1
                                    data_size_y = data_bounds(2,2)-data_bounds(1,2)+1
                                    data_size_z = data_bounds(2,3)-data_bounds(1,3)+1

                                    buffer_size = data_size_x*data_size_y*data_size_z

                                    call write_hvy_data( Bs, g, NdF, buffer_size, data_buffer, data_bounds, hvy_block(:,:,:,:,hvy_id) )

                                case('average', 'compare')
                                    ! treat internal neighbor like external neighbor
                                    ! but do not MPI_send/receive the data

                                    ! fill int/real buffer
                                    ! note: much better performance do fill the buffers here and not inside a additional subroutine
                                    ! this is of course worse programing style, but we run faster
                                    ! ---------------------------------------------------------------------------------------------
                                    buffer_position = int_send_buffer(1, rank+1 )+1

                                    real_send_buffer( buffer_position : buffer_position-1 + buffer_size, rank+1 ) = data_buffer(1:buffer_size)

                                    int_send_buffer( int_pos(rank+1), rank+1 )   = hvy_id
                                    int_send_buffer( int_pos(rank+1)+1, rank+1 ) = neighborhood
                                    int_send_buffer( int_pos(rank+1)+2, rank+1 ) = level_diff
                                    int_send_buffer( int_pos(rank+1)+3, rank+1 ) = buffer_position
                                    int_send_buffer( int_pos(rank+1)+4, rank+1 ) = buffer_size

                                    int_send_buffer(1, rank+1 ) = int_send_buffer(1, rank+1 ) + buffer_size

                                    ! increase int buffer position
                                    int_pos(rank+1) = int_pos(rank+1) + 5

                                    ! markiere das aktuelle ende des buffers, falls weitere elemente dazu kommen, wird die -99 wieder überschrieben
                                    int_send_buffer( int_pos(rank+1)  , rank+1 ) = -99

                                case('staging')
                                    ! set synch status
                                    call set_synch_status( params, synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_active(k), &
                                         lgt_id, neighbor_light_id, neighborhood, synch_method, &
                                         lgt_block(lgt_id,params%max_treelevel+2), lgt_block(neighbor_light_id,params%max_treelevel+2) )

                                    ! data has to synchronize in current stage
                                    if (synch) then

                                        ! fill int/real buffer
                                        ! note: much better performance do fill the buffers here and not inside a additional subroutine
                                        ! this is of course worse programing style, but we run faster
                                        ! ---------------------------------------------------------------------------------------------
                                        buffer_position = int_send_buffer(1, rank+1 )+1

                                        real_send_buffer( buffer_position : buffer_position-1 + buffer_size, rank+1 ) = data_buffer(1:buffer_size)

                                        int_send_buffer( int_pos(rank+1), rank+1 )   = hvy_id
                                        int_send_buffer( int_pos(rank+1)+1, rank+1 ) = neighborhood
                                        int_send_buffer( int_pos(rank+1)+2, rank+1 ) = level_diff
                                        int_send_buffer( int_pos(rank+1)+3, rank+1 ) = buffer_position
                                        int_send_buffer( int_pos(rank+1)+4, rank+1 ) = buffer_size

                                        int_send_buffer(1, rank+1 ) = int_send_buffer(1, rank+1 ) + buffer_size

                                        ! increase int buffer position
                                        int_pos(rank+1) = int_pos(rank+1) + 5

                                        ! markiere das aktuelle ende des buffers, falls weitere elemente dazu kommen, wird die -99 wieder überschrieben
                                        int_send_buffer( int_pos(rank+1)  , rank+1 ) = -99

                                    end if

                            end select

                        else
                            ! externer nachbar
                            ! synch status for staging method
                            if ( data_writing_type == 'staging' ) then
                                ! set synch status
                                call set_synch_status( params, synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_active(k), &
                                lgt_id, neighbor_light_id, neighborhood, synch_method, &
                                lgt_block(lgt_id,params%max_treelevel+2), lgt_block(neighbor_light_id,params%max_treelevel+2) )

                            else
                                ! synch status is allways true
                                synch = .true.
                            end if

                            ! first: fill com matrix, count number of communication to neighboring process, needed for int buffer length
                            com_matrix(neighbor_rank+1) = com_matrix(neighbor_rank+1) + 1

                            if (synch) then
                                ! active block send data to his neighbor block

                                ! fill int/real buffer
                                ! note: much better performance do fill the buffers here and not inside a additional subroutine
                                ! this is of course worse programing style, but we run faster
                                ! ---------------------------------------------------------------------------------------------
                                buffer_position = int_send_buffer(1, neighbor_rank+1 )+1

                                real_send_buffer( buffer_position : buffer_position-1 + buffer_size, neighbor_rank+1 ) = data_buffer(1:buffer_size)

                                int_send_buffer( int_pos(neighbor_rank+1), neighbor_rank+1 )   = hvy_id
                                int_send_buffer( int_pos(neighbor_rank+1)+1, neighbor_rank+1 ) = neighborhood
                                int_send_buffer( int_pos(neighbor_rank+1)+2, neighbor_rank+1 ) = level_diff
                                int_send_buffer( int_pos(neighbor_rank+1)+3, neighbor_rank+1 ) = buffer_position
                                int_send_buffer( int_pos(neighbor_rank+1)+4, neighbor_rank+1 ) = buffer_size

                                int_send_buffer(1, neighbor_rank+1 ) = int_send_buffer(1, neighbor_rank+1 ) + buffer_size

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

            ! timing
            call toc( params, "---synchronize ghost: fill buffers", MPI_wtime()-sub_t0, .false. )
            sub_t0 = MPI_Wtime()

            ! send/receive data
            ! note: todo, remove dummy subroutine
            ! note: new dummy subroutine sets receive buffer position accordingly to process rank
            ! note: todo: use more than non-blocking send/receive
            call isend_irecv_data_2( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix  )

            ! fill receive buffer for internal neighbors for averaging writing type
            ! note: only work if send buffer has data
            if ( (data_writing_type == 'average')     .or. (data_writing_type == 'compare')       &
            .or. (data_writing_type == 'staging') ) then
                if ( int_send_buffer(1 , rank+1 ) > 0 ) then
                    ! fill receive buffer
                    int_receive_buffer( 1:int_pos(rank+1)  , rank+1 ) = int_send_buffer( 1:int_pos(rank+1)  , rank+1 )
                    real_receive_buffer( 1:int_receive_buffer(1,rank+1), rank+1 ) = real_send_buffer( 1:int_receive_buffer(1,rank+1), rank+1 )
                    ! change com matrix, need to sort in buffers in next step
                    com_matrix(rank+1) = 1
                end if
            end if

            ! timing
            call toc( params, "---synchronize ghost: send data", MPI_wtime()-sub_t0, .false. )
            sub_t0 = MPI_Wtime()

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
                                    data_size_x = data_bounds(2,1)-data_bounds(1,1)+1
                                    data_size_y = data_bounds(2,2)-data_bounds(1,2)+1
                                    data_size_z = data_bounds(2,3)-data_bounds(1,3)+1

                                    buffer_size = data_size_x*data_size_y*data_size_z

                                    call write_hvy_data( Bs, g, NdF, buffer_size, data_buffer, data_bounds, hvy_block(:,:,:,:,hvy_id) )

                                case('average')
                                    ! add data
                                    call add_hvy_data( Bs, g, NdF, buffer_size, data_buffer(1:buffer_size), data_bounds, hvy_block(:,:,:,:,hvy_id), hvy_synch(:,:,:,hvy_id), data_writing_type )

                                case('staging')
                                    ! nothing to do

                                case('compare')
                                    ! only compare if until now everything is fine
                                    if ( stop_status ) then
                                        ! do nothing
                                    else
                                        ! compare data
                                        call compare_hvy_data( params, lgt_block, data_buffer, data_bounds, hvy_block, hvy_id, stop_status, level_diff, neighborhood )
                                    end if

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

                            ! calculate average for all nodes
                            do ix = 1, size(hvy_block,1)
                                do iy = 1, size(hvy_block,2)
                                    do iz = 1, size(hvy_block,3)
                                        if ( hvy_synch(ix, iy, iz, hvy_active(k)) > 1 ) then

                                            hvy_block(ix, iy, iz, dF, hvy_active(k)) = hvy_block(ix, iy, iz, dF, hvy_active(k)) &
                                                                                     / real( hvy_synch(ix, iy, iz, hvy_active(k)) , kind=rk)

                                        end if
                                    end do
                                end do
                            end do

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
                            call set_synch_status( params, synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_active(k), &
                            lgt_id, neighbor_light_id, neighborhood, synch_method, &
                            lgt_block(lgt_id,params%max_treelevel+2), lgt_block(neighbor_light_id,params%max_treelevel+2) )
                            ! synch == .true. bedeutet, dass der aktive block seinem nachbarn daten gibt
                            ! hier sind wir aber auf der seite des empfängers, das bedeutet, neighbor_synch muss ausgewertet werden

                            if (neighbor_synch) then

                                ! search buffers for synchronized data
                                ! first element in int buffer is real buffer size
                                l = 2

                                ! -99 marks end of data
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
                                        call write_hvy_data( Bs, g, NdF, buffer_size, data_buffer, data_bounds, hvy_block(:,:,:,:,hvy_active(k)) )

                                        ! done, exit the while loop?
                                        exit

                                    end if

                                    ! increase buffer postion marker
                                    l = l + 5

                                end do

                            end if

                        end if
                    end do
                end do

            end if

            ! timing
            call toc( params, "---synchronize ghost: write data", MPI_wtime()-sub_t0, .false. )
            sub_t0 = MPI_Wtime()

        end do

    end if

    ! clean up
    deallocate( data_buffer, res_pre_data, com_matrix, int_pos, dummy_matrix )

    ! timing
    call toc( params, "---synchronize ghost: ...", MPI_wtime()-sub_t0, .true. )
    sub_t0 = MPI_Wtime()

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

!############################################################################################################

subroutine restrict_predict_data( params, Bs, g, NdF, res_pre_data, data_bounds, neighborhood, level_diff, data_bounds_type, hvy_block )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params

    !> grid parameter
    integer(kind=ik), intent(in)                    :: Bs, g
    !> number of datafields
    integer(kind=ik), intent(in)                    :: NdF

    !> data buffer
    real(kind=rk), intent(out)                      :: res_pre_data(Bs+2*g, Bs+2*g, Bs+2*g, NdF)
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff
    ! data_bounds_type
    character(len=25), intent(in)                   :: data_bounds_type
    !> heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(Bs+2*g, Bs+2*g, Bs+2*g, NdF)

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, iN, jN, kN

    ! start and edn shift values
    integer(kind=ik)                                :: sh_start, sh_end

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! data size
    iN = data_bounds(2,1) - data_bounds(1,1) + 1
    jN = data_bounds(2,2) - data_bounds(1,2) + 1
    kN = data_bounds(2,3) - data_bounds(1,3) + 1

    sh_start = 0
    sh_end   = 0

    if ( data_bounds_type == 'exclude_redundant' ) then
        sh_start = 1
    end if
    if ( data_bounds_type == 'only_redundant' ) then
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
                                                       data_bounds(1,3):data_bounds(2,3), dF ), &
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
                                    = hvy_block( i, j, k, dF )

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
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF ), &
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
                    do dF = 1, NdF
                        ! first dimension
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            ! second dimension
                            do j = data_bounds(1,2), data_bounds(2,2), 2

                                ! write restricted data
                                res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, 1, dF) &
                                = hvy_block( i, j, 1, dF )

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
                    do dF = 1, NdF
                        ! interpolate data
                        call prediction_2D( hvy_block( data_bounds(1,1):data_bounds(2,1), data_bounds(1,2):data_bounds(2,2), 1, dF ), &
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
                    do dF = 1, NdF
                        ! first dimension
                        do i = data_bounds(1,1), data_bounds(2,1), 2
                            ! second dimension
                            do j = data_bounds(1,2), data_bounds(2,2), 2

                                ! write restricted data
                                res_pre_data( (i-data_bounds(1,1))/2+1, (j-data_bounds(1,2))/2+1, 1, dF) &
                                = hvy_block( i, j, 1, dF )

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

!subroutine read_hvy_data( params, data_buffer, buffer_counter, hvy_data )
subroutine read_hvy_data( NdF, data_size_x, data_size_y, data_size_z, data_buffer, hvy_data )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> number of datafields
    integer(kind=ik), intent(in)                    :: NdF

    !> data buffer
    real(kind=rk), intent(out)                 :: data_buffer(data_size_x*data_size_y*data_size_z)

    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_size_x, data_size_y, data_size_z

    !> heavy block data, all data fields
    real(kind=rk), intent(in)                       :: hvy_data(data_size_x, data_size_y, data_size_z, NdF)

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_counter

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset buffer size
    buffer_counter = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all data fields
    do dF = 1, NdF
        ! first dimension
        do i = 1, data_size_x
            ! second dimension
            do j = 1, data_size_y
                ! third dimension, note: for 2D cases kN is allways 1
                do k = 1, data_size_z

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

subroutine write_hvy_data( Bs, g, NdF, buffer_size, data_buffer, data_bounds, hvy_block )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)                    :: Bs, g
    !> number of datafields
    integer(kind=ik), intent(in)                    :: NdF
    !> size of received data buffer
    integer(kind=ik), intent(in)                    :: buffer_size

    !> data buffer
    real(kind=rk), intent(in)                 :: data_buffer(buffer_size)
    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                       :: hvy_block(Bs+2*g, Bs+2*g, Bs+2*g, NdF)

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
    do dF = 1, NdF
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases kN is allways 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    ! write data buffer
                    hvy_block( i, j, k, dF ) = data_buffer( buffer_i )
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    end do

end subroutine write_hvy_data

!############################################################################################################

subroutine add_hvy_data( Bs, g, NdF, buffer_size, data_buffer, data_bounds, block_data, hvy_synch, data_writing_type )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in)                    :: Bs, g
    !> number of datafields
    integer(kind=ik), intent(in)                    :: NdF
    !> size of received data buffer
    integer(kind=ik), intent(in)                    :: buffer_size

    !> data buffer
    real(kind=rk), intent(in)                 :: data_buffer(buffer_size)

    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds(2,3)

    !> heavy data array - block data
    real(kind=rk), intent(inout)                       :: block_data(Bs+2*g, Bs+2*g, Bs+2*g, NdF)

    !> heavy synch array
    integer(kind=1), intent(inout)      :: hvy_synch(Bs+2*g, Bs+2*g, Bs+2*g)

    ! type of data writing
    character(len=25), intent(in)                   :: data_writing_type

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    buffer_i = 1

!---------------------------------------------------------------------------------------------
! main body

    ! without condition check inside loop
    ! datafield 1:
    ! ------------
    if ( data_writing_type == 'average' ) then
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases kN is allways 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    ! write data buffer
                    block_data( i, j, k, 1 ) = block_data( i, j, k, 1 ) + data_buffer( buffer_i )

                    ! count synchronized data
                    hvy_synch( i, j, k ) = hvy_synch( i, j, k ) + int(1,kind=1)

                    ! increase buffer counter
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    else
        ! without hvy_synch array
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases kN is allways 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    ! write data buffer
                    block_data( i, j, k, 1 ) = block_data( i, j, k, 1 ) + data_buffer( buffer_i )

                    ! increase buffer counter
                    buffer_i = buffer_i + 1

                end do
            end do
        end do

    end if

    ! other datafields:
    ! -----------------
    ! loop over all data fields
    !do dF = 2, params%number_data_fields
    do dF = 2, NdF
        ! first dimension
        do i = data_bounds(1,1), data_bounds(2,1)
            ! second dimension
            do j = data_bounds(1,2), data_bounds(2,2)
                ! third dimension, note: for 2D cases kN is allways 1
                do k = data_bounds(1,3), data_bounds(2,3)

                    ! write data buffer
                    block_data( i, j, k, dF ) = block_data( i, j, k, dF ) + data_buffer( buffer_i )

                    ! increase buffer counter
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    end do

end subroutine add_hvy_data

!############################################################################################################

subroutine compare_hvy_data( params, lgt_block, data_buffer, data_bounds, hvy_block, hvy_id, stop_status, level_diff, neighborhood )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> data buffer
    real(kind=rk), intent(in)                 :: data_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(in)                 :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                       :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id, level_diff
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i, lgt_id

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
    eps = 1e-12_rk!1e-9_rk

    ! reset error norm
    error_norm = 0.0_rk

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

                    if (level_diff/=-1) then
                        ! pointwise error norm
                        error_norm = max(error_norm, dabs(hvy_block( i, j, k, dF, hvy_id ) - data_buffer( buffer_i )))

                    else
                        ! if the level diff is -1, I compare with interpolated (upsampled) data. that means every EVEN
                        ! point is the result of interpolation, and not truely redundant.
                        ! Note this routine ALWAYS just compares the redundant nodes, so it will mostly be called
                        ! with a line of points (i.e. one dimension is length one)
                        ! \todo: check if number of ghost nodes is odd or even
                        !
                        ! This routine has been tested:
                        !   - old method (working version): no error found (okay)
                        !   - old method, non_uniform_mesh_correction=0; in params file -> plenty of errors (okay)
                        !   - old method, sync stage 4 deactivated: finds all occurances of "3finer blocks on corner problem" (okay)
                        !   - new method, averaging, no error found (makes sense: okay)
                        if ( (mod(i,2)/=0) .and. (mod(j,2)/=0) .and. (mod(k,2)/=0) ) then
                            ! pointwise error norm
                            error_norm = max(error_norm, dabs(hvy_block( i, j, k, dF, hvy_id ) - data_buffer( buffer_i )))

                        end if

                    end if
                    buffer_i = buffer_i + 1

                end do
            end do
        end do
    end do

    ! calculate light id
    call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

    if (error_norm > eps)  then
        ! error message
        write(*,*) "ERROR: difference in redundant nodes", error_norm, level_diff, neighborhood, lgt_block(lgt_id,:)
        ! stop program
        stop_status = .true.
        ! mark block by putting a dot in the middle. this way, we can identify all
        ! blocks that have problems. If we set the entire block to 100, then subsequent
        ! synchronizations fail because of that. If we set just a point in the middle, it
        ! is far away from the boundary and thus does not affect other sync steps.
        !hvy_block( size(hvy_block,1)/2, size(hvy_block,2)/2, :, :, hvy_id ) = 100.0_rk
        hvy_block( :, :, :, :, hvy_id ) = 100.0_rk
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

subroutine set_synch_status( params, synch_stage, synch, neighbor_synch, level_diff, hvy_neighbor, hvy_id, lgt_id, &
                             neighbor_lgt_id, neighborhood, synch_method, my_ref_status, neighbor_ref_status )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params

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

    ! level diffenrence
    integer(kind=ik), intent(in)        :: lgt_id, neighbor_lgt_id

    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood

    ! type of data writing
    integer(kind=ik), intent(in)                   :: synch_method

    ! refinement status
    ! need to check synch status in stage 0
    integer(kind=ik), intent(in)                   :: my_ref_status, neighbor_ref_status

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    synch = .false.
    neighbor_synch = .false.

    select case(synch_method)

        ! 'staging_old'
        case(1)
            if ( params%threeD_case ) then
                ! 3D
                ! set synch stage
                ! stage 1: level +1
                ! stage 2: level 0
                ! stage 3: level -1
                ! stage 4: special

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

                    ! '_12/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 7 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 27) /= -1) .or. & !'__1/123'
                             (hvy_neighbor( hvy_id, 30) /= -1) .or. & !'__1/152'
                             (hvy_neighbor( hvy_id, 31) /= -1) .or. & !'__2/123'
                             (hvy_neighbor( hvy_id, 33) /= -1) ) then !'__2/152'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_13/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 8 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 27) /= -1) .or. & !'__1/123'
                             (hvy_neighbor( hvy_id, 28) /= -1) .or. & !'__1/134'
                             (hvy_neighbor( hvy_id, 35) /= -1) .or. & !'__3/123'
                             (hvy_neighbor( hvy_id, 37) /= -1) ) then !'__3/134'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_14/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 9 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 28) /= -1) .or. & !'__1/134'
                             (hvy_neighbor( hvy_id, 29) /= -1) .or. & !'__1/145'
                             (hvy_neighbor( hvy_id, 39) /= -1) .or. & !'__4/134'
                             (hvy_neighbor( hvy_id, 41) /= -1) ) then !'__4/145'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_15/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 10 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 29) /= -1) .or. & !'__1/145'
                             (hvy_neighbor( hvy_id, 30) /= -1) .or. & !'__1/152'
                             (hvy_neighbor( hvy_id, 43) /= -1) .or. & !'__5/145'
                             (hvy_neighbor( hvy_id, 45) /= -1) ) then !'__5/152'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_62/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 11 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 32) /= -1) .or. & !'__2/623'
                             (hvy_neighbor( hvy_id, 34) /= -1) .or. & !'__2/652'
                             (hvy_neighbor( hvy_id, 47) /= -1) .or. & !'__6/623'
                             (hvy_neighbor( hvy_id, 50) /= -1) ) then !'__6/652'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_63/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 12 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 36) /= -1) .or. & !'__3/623'
                             (hvy_neighbor( hvy_id, 38) /= -1) .or. & !'__3/634'
                             (hvy_neighbor( hvy_id, 47) /= -1) .or. & !'__6/623'
                             (hvy_neighbor( hvy_id, 48) /= -1) ) then !'__6/634'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_64/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 13 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 40) /= -1) .or. & !'__4/634'
                             (hvy_neighbor( hvy_id, 42) /= -1) .or. & !'__4/645'
                             (hvy_neighbor( hvy_id, 48) /= -1) .or. & !'__6/634'
                             (hvy_neighbor( hvy_id, 49) /= -1) ) then !'__6/645'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_65/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 14 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 44) /= -1) .or. & !'__5/645'
                             (hvy_neighbor( hvy_id, 46) /= -1) .or. & !'__5/652'
                             (hvy_neighbor( hvy_id, 49) /= -1) .or. & !'__6/645'
                             (hvy_neighbor( hvy_id, 50) /= -1) ) then !'__6/652'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_23/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 15 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 31) /= -1) .or. & !'__2/123'
                             (hvy_neighbor( hvy_id, 32) /= -1) .or. & !'__2/623'
                             (hvy_neighbor( hvy_id, 35) /= -1) .or. & !'__3/123'
                             (hvy_neighbor( hvy_id, 36) /= -1) ) then !'__3/623'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_25/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 16 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 33) /= -1) .or. & !'__2/152'
                             (hvy_neighbor( hvy_id, 34) /= -1) .or. & !'__2/652'
                             (hvy_neighbor( hvy_id, 45) /= -1) .or. & !'__5/152'
                             (hvy_neighbor( hvy_id, 46) /= -1) ) then !'__5/652'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_43/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 17 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 37) /= -1) .or. & !'__3/134'
                             (hvy_neighbor( hvy_id, 38) /= -1) .or. & !'__3/634'
                             (hvy_neighbor( hvy_id, 39) /= -1) .or. & !'__4/134'
                             (hvy_neighbor( hvy_id, 40) /= -1) ) then !'__4/634'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '_45/___'
                    ! coarser neighbor at face overwrides nodes from edge neighbor at same level -> need correction
                    ! check edge neigborhood, send data again
                    if ( neighborhood == 18 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 41) /= -1) .or. & !'__4/145'
                             (hvy_neighbor( hvy_id, 42) /= -1) .or. & !'__4/645'
                             (hvy_neighbor( hvy_id, 43) /= -1) .or. & !'__5/145'
                             (hvy_neighbor( hvy_id, 44) /= -1) ) then !'__5/645'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '123/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 19 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 27) /= -1) .or. & !'__1/123'
                             (hvy_neighbor( hvy_id, 31) /= -1) .or. & !'__2/123'
                             (hvy_neighbor( hvy_id, 35) /= -1) ) then !'__3/123'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '134/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 20 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 28) /= -1) .or. & !'__1/134'
                             (hvy_neighbor( hvy_id, 37) /= -1) .or. & !'__3/134'
                             (hvy_neighbor( hvy_id, 39) /= -1) ) then !'__4/134'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '145/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 21 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 29) /= -1) .or. & !'__1/145'
                             (hvy_neighbor( hvy_id, 41) /= -1) .or. & !'__4/145'
                             (hvy_neighbor( hvy_id, 43) /= -1) ) then !'__5/145'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '152/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 22 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 30) /= -1) .or. & !'__1/152'
                             (hvy_neighbor( hvy_id, 33) /= -1) .or. & !'__2/152'
                             (hvy_neighbor( hvy_id, 45) /= -1) ) then !'__5/152'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '623/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 23 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 32) /= -1) .or. & !'__2/623'
                             (hvy_neighbor( hvy_id, 36) /= -1) .or. & !'__3/623'
                             (hvy_neighbor( hvy_id, 47) /= -1) ) then !'__6/623'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '634/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 24 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 38) /= -1) .or. & !'__3/634'
                             (hvy_neighbor( hvy_id, 40) /= -1) .or. & !'__4/634'
                             (hvy_neighbor( hvy_id, 48) /= -1) ) then !'__6/634'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '645/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 25 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 42) /= -1) .or. & !'__4/645'
                             (hvy_neighbor( hvy_id, 44) /= -1) .or. & !'__5/645'
                             (hvy_neighbor( hvy_id, 49) /= -1) ) then !'__6/645'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                    ! '652/___'
                    ! coarser neighbor at edge overwrides nodes from corner neighbor at same level -> need correction
                    ! check corner and face neigborhood, send data again
                    if ( neighborhood == 26 ) then
                        ! edge neighbor exists
                        if ( (hvy_neighbor( hvy_id, 34) /= -1) .or. & !'__2/652'
                             (hvy_neighbor( hvy_id, 46) /= -1) .or. & !'__5/652'
                             (hvy_neighbor( hvy_id, 50) /= -1) ) then !'__6/652'
                            synch = .true.
                            neighbor_synch = .true.
                        end if
                    end if

                end if

            else
                ! set synch stage
                ! stage 1: level +1
                ! stage 2: level 0
                ! stage 3: level -1
                ! stage 4: special

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

                ! 2D
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

            end if

!        ! staging new
!        case(?)
!            ! set synch stage
!            ! stage 1: level +1
!            ! stage 2: level -1
!            ! stage 3: level 0
!
!            ! stage 1
!            if ( (synch_stage == 1) .and. (level_diff == 1) ) then
!                ! block send data
!                synch = .true.
!            elseif ( (synch_stage == 1) .and. (level_diff == -1) ) then
!                ! neighbor send data
!                neighbor_synch = .true.
!            end if
!
!            ! stage 2
!            if ( (synch_stage == 2) .and. (level_diff == -1) ) then
!                ! block send data
!                synch = .true.
!            elseif ( (synch_stage == 2) .and. (level_diff == 1) ) then
!                ! neighbor send data
!                neighbor_synch = .true.
!            end if
!
!            ! stage 3
!            if ( (synch_stage == 3) .and. (level_diff == 0) ) then
!                ! block send data
!                synch = .true.
!                ! neighbor send data
!                neighbor_synch = .true.
!            end if

        ! stage 0
        case(11)
            ! this is the zeroth stage. it corrects blocks that are on the same level, but have a different history. one is on Jmax from
            ! before, one has just gotten to Jmax via interpolation. In those cases, the former block has the status +11
            ! which indicates that its redundant nodes must overwrite the ones on the other block (which has been interpolated)
            if ( level_diff == 0 ) then
                if ((my_ref_status==11) .and. (neighbor_ref_status/=11)) then
                    ! if a block has the +11 status, it must send data to the neighbor, if that is not +11
                    synch = .true.
                elseif ((my_ref_status/=11) .and. (neighbor_ref_status==11)) then
                    ! if a block is not +11 and its neighbor is, then unpack data
                    neighbor_synch = .true.
                endif
            endif

    end select

end subroutine set_synch_status

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
