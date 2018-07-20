!>\dir
!> MPI synchronization routines

!> \file
!> \brief File of MPI Module


!> \brief This module implements all MPI synchronization routines
!> \details
!>          * synchronization of ghost nodes
!>          * synchronization of light data
!>          * creation of rank block rank lists
!>          * copy, write send and recieve buffers
!> \version 0.6
!> \author engels, reiss, msr
!! 12/01/17 - create
!! 18/07/2018 - remove all old ghost node routines, build in JR's improved version of MSR's new ghost nodes
module module_MPI

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! interpolation routines
    use module_interpolation
    use module_treelib

    implicit none
!---------------------------------------------------------------------------------------------
! variables

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! everything is save by default
    SAVE

    ! Just because we have MPISIZE ranks does not mean that everybody talks to everybody.
    ! While this CAN happen, it is much more likely that an MPIRANK talks only to a limited number of
    ! "friends". This is a direct consequence of the space-filling curves. Therefore, we do not allocate
    ! buffer for all MPSIZE ranks, but only for N_FRIENDS. If the number happens to be too small, we increase
    ! it dynamically (via deallocate / reallocate)
    integer(kind=ik) :: N_friends, N_friends_used

    ! We require two stages: first, we fill all ghost nodes which are simple copy,
    ! then in the second stage we can use interpolation and fill the remaining ones.
    ! In order not to send ALL data in both stages, we allocate one buffer for each stage.
    integer(kind=ik) :: Nstages = 2

    ! send/receive buffer, integer and real
    ! allocate in init substep not in synchronize subroutine, to avoid slow down when using
    ! large numbers of processes and blocks per process, when allocating on every call to the routine
    integer(kind=ik), allocatable :: int_send_buffer(:,:,:), int_receive_buffer(:,:,:)
    real(kind=rk), allocatable    :: real_send_buffer(:,:,:), real_receive_buffer(:,:,:)

    ! this array is used only in AVERAGING submodule of (deprecated) MSR ghost nodes.
    ! TODO: remove, as averaging did not work.
    integer(kind=1), allocatable  :: hvy_synch(:, :, :, :)

    ! an array to count how many messages we send to the other mpiranks. NOTE: the
    ! arrays communication_counter(:), int_pos(:) have the size N_friends
    integer(kind=ik), allocatable :: communication_counter(:,:), int_pos(:,:), mpirank2friend(:), friend2mpirank(:)

    ! internally, we flatten the ghost nodes layers to a line. this is stored in
    ! this buffer (NOTE max size is (blocksize)*(ghost nodes size + 1)*(number of datafields))
    real(kind=rk), allocatable    :: line_buffer(:)

    ! restricted/predicted data buffer
    real(kind=rk), allocatable    :: res_pre_data(:,:,:,:), tmp_block(:,:,:,:)

    ! it is faster to use named consts than strings, although strings are nicer to read
    integer(kind=ik), PARAMETER   :: EXCLUDE_REDUNDANT = 1_ik, INCLUDE_REDUNDANT = 2_ik, ONLY_REDUNDANT = 3_ik
    integer(kind=ik), PARAMETER   :: SENDER=1_ik, RECVER=2_ik, RESPRE=3_ik

    ! we set up a table that gives us directly the inverse neighbor relations.
    ! it is filled (once) in init_ghost_nodes
    integer(kind=ik), dimension(1:74,2:3) :: inverse_neighbor

    ! these arrays are used in the compare_hvy_data routines for the ghost nodes origin test. They allow
    ! much nicer handling of errors. Only allocated in the "check_unique_origin" routine (i.e. unused in production)
    real(kind=rk), allocatable :: hvy_block_test_err(:,:,:,:,:)
    real(kind=rk), allocatable :: hvy_block_test_val(:,:,:,:,:)
    real(kind=rk), allocatable :: hvy_block_test_interpref(:,:,:,:,:)


    ! We frequently need to know the indices of a ghost nodes chunk. Thus we save them
    ! once in a large, module-global array (which is faster than computing it every time with tons
    ! of IF-THEN clauses).
    ! This arrays indices are:
    ! ijkGhosts([start,end], [dir (x,y,z)], [neighborhood], [level-diff], [data_bounds_type], [sender/receiver/up-downsampled])
    integer(kind=ik), dimension(1:2, 1:3, 1:74, -1:1, 1:3, 1:3) :: ijkGhosts

    ! it is useful to keep a named constant for the dimensionality here (we use
    ! it to access e.g. two/three D arrays in inverse_neighbor)
    integer(kind=ik) :: dim = 2_ik

    ! we use this flag to call the allocation routine only once.
    logical          :: ghost_nodes_module_ready = .false.

    ! two shift parameters (asymmetric and symmetric ) used for selection of interpolation
    ! bounds on sender side. used to avoid one-sided interpolation if desired. They're
    ! global so we can save them easily to the ghosts_bounds.dat file
    integer(kind=ik) :: A, S
!---------------------------------------------------------------------------------------------
! public parts of this module

    PUBLIC :: sync_ghosts, blocks_per_mpirank, synchronize_lgt_data, reset_ghost_nodes
    PUBLIC :: check_redundant_nodes, synchronize_ghosts_generic_sequence, init_ghost_nodes, check_unique_origin

!---------------------------------------------------------------------------------------------
! main body

contains

    include "synchronize_ghosts.f90"
    include "synchronize_ghosts_generic.f90"
    include "blocks_per_mpirank.f90"
    include "synchronize_lgt_data.f90"
    include "reset_ghost_nodes.f90"
    include "check_redundant_nodes.f90"
    include "calc_data_bounds.f90"
    include "restrict_predict_data.f90"


! Just because we have MPISIZE ranks does not mean that everybody talks to everybody.
! While this CAN happen, it is much more likely that an MPIRANK talks only to a limited number of
! "friends". This is a direct consequence of the space-filling curves. Therefore, we do not allocate
! buffer for all MPSIZE ranks, but only for N_FRIENDS. If the number happens to be too small, we increase
! it dynamically (via deallocate / reallocate)
subroutine reallocate_buffers(params)
    implicit none
    type (type_params), intent(in) :: params
    integer(kind=ik), allocatable :: int_buffer_tmp(:,:,:)
    real(kind=rk), allocatable    :: real_buffer_tmp(:,:,:)
    integer(kind=ik), allocatable :: communication_counter_tmp(:,:), int_pos_tmp(:,:), &
    mpirank2friend_tmp(:), friend2mpirank_tmp(:)

    write(*,'("GHOSTS-runtime: rank=",i5," is changing buffer size to N_friends=",i4)') params%rank, N_friends

    allocate( int_buffer_tmp( size(int_send_buffer,1), size(int_send_buffer,2), 1:Nstages) )
        int_buffer_tmp = int_send_buffer
        deallocate(int_send_buffer)
        allocate( int_send_buffer(size(int_buffer_tmp,1), N_friends, 1:Nstages) )
        int_send_buffer(:, 1:size(int_buffer_tmp,2), : ) = int_buffer_tmp

        int_buffer_tmp = int_receive_buffer
        deallocate(int_receive_buffer)
        allocate( int_receive_buffer(size(int_buffer_tmp,1), N_friends, 1:Nstages) )
        int_receive_buffer(:, 1:size(int_buffer_tmp,2), : ) = int_buffer_tmp

        ! new appended buffer requires initialization (see main routine for doc)
        int_send_buffer( 1, N_friends, : ) = 0
        int_send_buffer( 2, N_friends, : ) = -99
    deallocate(int_buffer_tmp)

    allocate( real_buffer_tmp( size(real_send_buffer,1), size(real_send_buffer,2), 1:Nstages) )
        ! very slow...
        real_buffer_tmp = real_send_buffer
        deallocate(real_send_buffer)
        allocate( real_send_buffer(size(real_buffer_tmp,1), N_friends, 1:Nstages) )
        real_send_buffer(:, 1:size(real_buffer_tmp,2), : ) = real_buffer_tmp

        ! very slow...
        real_buffer_tmp = real_receive_buffer
        deallocate(real_receive_buffer)
        allocate( real_receive_buffer(size(real_buffer_tmp,1), N_friends, 1:Nstages) )
        real_receive_buffer(:, 1:size(real_buffer_tmp,2), : ) = real_buffer_tmp
    deallocate(real_buffer_tmp)


    allocate( communication_counter_tmp(size(communication_counter), 1:Nstages) )
    communication_counter_tmp = communication_counter
    deallocate(communication_counter)
    allocate( communication_counter(1:N_friends, 1:Nstages) )
    communication_counter(1:size(communication_counter_tmp), :) = communication_counter_tmp
    deallocate( communication_counter_tmp )
    communication_counter( N_friends, : ) = 0


    allocate( int_pos_tmp(size(int_pos), 1:Nstages) )
    int_pos_tmp = int_pos
    deallocate(int_pos)
    allocate( int_pos(1:N_friends, 1:Nstages) )
    int_pos(1:size(int_pos_tmp),:) = int_pos_tmp
    deallocate( int_pos_tmp )
    ! new appended buffer requires initialization
    int_pos(N_friends, 1:Nstages) = 2


    allocate( friend2mpirank_tmp(size(friend2mpirank)) )
    friend2mpirank_tmp = friend2mpirank
    deallocate(friend2mpirank)
    allocate( friend2mpirank(1:N_friends) )
    friend2mpirank(1:size(friend2mpirank_tmp)) = friend2mpirank_tmp
    deallocate( friend2mpirank_tmp )

end subroutine


!! The friends concept avoids to reserve memory so that all procs can talk to all
!! other procs. There is a simple, unique, invertible relation between mpirank and
!! Friend ID established here. If a proc wants to add a Friend and the pre-allocated
!! array is full, then the buffers are increased. Note this process is not for free
!! but rather time consuming. Best is not to use the functionality, by allocating enough
!! Friends at the start.
subroutine get_friend_id_for_mpirank( params, neighbor_rank, id_Friend )
    implicit none
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in) :: neighbor_rank
    integer(kind=ik), intent(out) :: id_Friend

    ! did we already add this proc to the friends list?
    if (mpirank2friend(neighbor_rank+1) < 0) then
        ! no, we didn't
        if (N_friends_used < N_friends) then ! some free friends-slots left?
            N_friends_used = N_friends_used +1
            mpirank2friend(neighbor_rank+1) = N_friends_used ! one-based
            friend2mpirank(N_friends_used) = neighbor_rank+1 ! one-based
        else
            ! no space left for friends, re-allocate
            N_friends = N_friends + 1
            N_friends_used = N_friends
            call reallocate_buffers(params)
            mpirank2friend(neighbor_rank+1) = N_friends_used ! one-based
            friend2mpirank(N_friends_used) = neighbor_rank+1 ! one-based
        endif
    endif
    id_Friend = mpirank2friend(neighbor_rank+1)
end subroutine



!! initialize ghost nodes module. allocate buffers and create data bounds array,
!! which we use to rapidly identify a ghost nodes layer
subroutine init_ghost_nodes( params )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in) :: params
    ! local variables
    integer(kind=ik) :: buffer_N_int, buffer_N, Bs, g, Neqn, number_blocks, rank
    integer(kind=ik) :: ineighbor, Nneighbor, leveldiff, idata_bounds_type
    integer(kind=ik) ::  j, rx0, rx1, ry0, ry1, rz0, rz1, sx0, sx1, sy0, sy1, sz0, sz1
    integer(kind=ik) :: i, k, status(1:4)
    integer(kind=ik) :: ijkrecv(2,3)
    integer(kind=ik) :: ijkbuffer(2,3)
    integer(kind=ik) :: ijksend(2,3)

    ! on second call, nothing happens
    if (.not. ghost_nodes_module_ready) then
        number_blocks   = params%number_blocks
        Bs              = params%number_block_nodes
        g               = params%number_ghost_nodes
        Neqn            = params%number_data_fields
        rank            = params%rank

        if (rank==0) write(*,'("---------------------------------------------------------")')
        if (rank==0) write(*,'("                     GHOST-INIT ")')
        if (rank==0) write(*,'("---------------------------------------------------------")')

        if (g>=(Bs+1)/2) then
            call abort(921151369, "Young skywalker, you failed at set g>=(Bs+1)/2 which implies &
            & that the ghost nodes layer can span beyond an entire finer block. Either decrease &
            & number_ghost_nodes or increase number_block_nodes.")
        endif

        ! set default number of "friends", that is mpiranks we exchange data with.
        ! NOTE: their number can be increased if necessary
        N_friends = min( params%number_procs, params%N_friends )

        ! synchronize buffer length
        ! assume: all blocks are used, all blocks have external neighbors,
        ! max neighbor number: 2D = 12, 3D = 56
        ! max neighborhood size, 2D: (Bs+g+1)*(g+1)
        ! max neighborhood size, 3D: (Bs+g+1)*(g+1)*(g+1)
        if ( params%threeD_case ) then
            !---3d---3d---
            ! space dimensions: used in the static arrays as index
            dim = 3

            buffer_N = number_blocks * Neqn * ((Bs+2*g)**dim - Bs**dim) / N_friends
            ! buffer_N = number_blocks * 56 * (Bs+g+1)*(g+1)*(g+1) * Neqn
            buffer_N_int = number_blocks * 56 * 3
            ! how many possible neighbor relations are there?
            Nneighbor = 74

            allocate( tmp_block( Bs+2*g, Bs+2*g, Bs+2*g, Neqn) )
        else
            !---2d---2d---
            ! space dimensions: used in the static arrays as index
            dim = 2

            buffer_N = number_blocks * Neqn * ((Bs+2*g)**dim - Bs**dim) / N_friends
            ! buffer_N = number_blocks * 12 * (Bs+g+1)*(g+1) * Neqn
            buffer_N_int = number_blocks * 12 * 3
            ! how many possible neighbor relations are there?
            Nneighbor = 16

            allocate( tmp_block( Bs+2*g, Bs+2*g, 1, Neqn) )
        end if

        !-----------------------------------------------------------------------
        ! allocate auxiliary memory
        !-----------------------------------------------------------------------
        ! allocate synch buffer
        if (rank==0) then
            write(*,'("GHOSTS-INIT: Attempting to allocate the ghost-sync-buffer.")')

            write(*,'("GHOSTS-INIT: buffer_N_int=",i12," buffer_N=",i12," N_friends=",i3," Nstages=",i1)') &
            buffer_N_int, buffer_N, N_friends, Nstages

            write(*,'("GHOSTS-INIT: On each MPIRANK, Int  buffer:", f9.4, "GB")') &
                2.0*dble(buffer_N_int)*dble(N_friends)*dble(Nstages)*8e-9

            write(*,'("GHOSTS-INIT: On each MPIRANK, Real buffer:", f9.4, "GB")') &
                2.0*dble(buffer_N)*dble(N_friends)*dble(Nstages)*8e-9
            write(*,'("---------------- allocating now ----------------")')
        endif

        ! wait now so that if allocation fails, we get at least the above info
        call MPI_barrier( WABBIT_COMM, status(1))

        allocate( int_send_buffer( 1:buffer_N_int, 1:N_friends, 1:Nstages), stat=status(1) )
        allocate( int_receive_buffer( 1:buffer_N_int, 1:N_friends, 1:Nstages), stat=status(2) )
        allocate( real_send_buffer( 1:buffer_N, 1:N_friends, 1:Nstages), stat=status(3) )
        allocate( real_receive_buffer( 1:buffer_N, 1:N_friends, 1:Nstages), stat=status(4) )

        if (maxval(status) /= 0) call abort(999999, "Buffer allocation failed. Not enough memory?")

        ! synch array, use for ghost nodes synchronization
        if (params%threeD_case) then
            allocate( hvy_synch( Bs+2*g, Bs+2*g, Bs+2*g, number_blocks ) )
        else
            allocate( hvy_synch( Bs+2*g, Bs+2*g, 1, number_blocks ) )
        endif

        if (rank==0) then
            write(*,'("GHOSTS-INIT: initial N_friends=",i4)') N_friends
            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "hvy_synch", shape(hvy_synch)

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "real_receive_buffer", shape(real_receive_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "real_send_buffer", shape(real_send_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "int_send_buffer", shape(int_send_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "int_receive_buffer", shape(int_receive_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, Real buffer size is",f9.4," GB ")') &
             2.0*dble(buffer_N)*dble(N_friends)*dble(Nstages)*8e-9

            write(*,'("GHOSTS-INIT: on each mpirank, Int  buffer size is",f9.4," GB ")') &
             2.0*dble(buffer_N_int)*dble(N_friends)*dble(Nstages)*8e-9
        endif

        ! this is a list of communications with all other procs
        allocate( communication_counter(1:N_friends, 1:Nstages) )
        allocate( int_pos(1:N_friends, 1:Nstages) )
        ! this is the list friend <-> mpirank
        allocate( mpirank2friend(1:params%number_procs) )
        allocate( friend2mpirank(1:N_friends) )

        allocate( line_buffer( Neqn*(Bs+2*g)**(dim) ) )

        !-----------------------------------------------------------------------
        ! set up constant arrays
        !-----------------------------------------------------------------------
        ! We frequently need to know the indices of a ghost nodes chunk. Thus we save them
        ! once in a large, module-global array (which is faster than computing it every time with tons
        ! of IF-THEN clauses).
        ! This arrays indices are:
        ! ijkGhosts([start,end], [dir], [ineighbor], [leveldiff], [idata_bounds_type], [isendrecv])

        ijkGhosts = 1
        do ineighbor = 1, Nneighbor
            do leveldiff = -1, 1
                do idata_bounds_type = 1, 3
                    call set_recv_bounds( params, ijkrecv, ineighbor, leveldiff, idata_bounds_type, 'receiver')
                    ijkGhosts(1:2, 1:3, ineighbor, leveldiff, idata_bounds_type, RECVER) = ijkrecv

                    call compute_sender_buffer_bounds(params, ijkrecv, ijksend, ijkbuffer, ineighbor, leveldiff, idata_bounds_type)
                    ijkGhosts(1:2, 1:3, ineighbor, leveldiff, idata_bounds_type, SENDER) = ijksend
                    ijkGhosts(1:2, 1:3, ineighbor, leveldiff, idata_bounds_type, RESPRE) = ijkbuffer
                enddo
            enddo
        enddo

        ! now we know how large the patches are we'd like to store in the RESPRE buffer
        i = maxval( ijkGhosts(2,1,:,:,:,RESPRE) )*2
        j = maxval( ijkGhosts(2,2,:,:,:,RESPRE) )*2
        k = maxval( ijkGhosts(2,3,:,:,:,RESPRE) )*2

        allocate( res_pre_data( i, j, k, Neqn) )

        ! this output can be plotted using the python script
        if (params%rank==0) Then
            open(16,file='ghost_bounds.dat',status='replace')
            write(16,'(i3)') Bs
            write(16,'(i3)') g
            write(16,'(i3)') S
            write(16,'(i3)') A
            do ineighbor = 1, Nneighbor
                do leveldiff = -1, 1
                    do idata_bounds_type = 1, 3
                        do i = 1, 2
                            do j = 1, 3
                                do k = 1, 3
                                    write(16,'(i3)') ijkGhosts(i, j, ineighbor, leveldiff, idata_bounds_type, k)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            close(16)
        endif

        ! set up table with inverse neighbor relations
        inverse_neighbor = -1
        ! ---2D------2D------2D---
        inverse_neighbor(1,2) = 3
        inverse_neighbor(2,2) = 4
        inverse_neighbor(3,2) = 1
        inverse_neighbor(4,2) = 2
        inverse_neighbor(5,2) = 8
        inverse_neighbor(6,2) = 7
        inverse_neighbor(7,2) = 6
        inverse_neighbor(8,2) = 5
        inverse_neighbor(9,2) = 11
        inverse_neighbor(10,2) = 12
        inverse_neighbor(11,2) = 9
        inverse_neighbor(12,2) = 10
        inverse_neighbor(13,2) = 15
        inverse_neighbor(14,2) = 16
        inverse_neighbor(15,2) = 13
        inverse_neighbor(16,2) = 14
        ! ---3D------3D------3D------3D---
        !'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___'
        inverse_neighbor(1,3) = 6
        inverse_neighbor(2,3) = 4
        inverse_neighbor(3,3) = 5
        inverse_neighbor(4,3) = 2
        inverse_neighbor(5,3) = 3
        inverse_neighbor(6,3) = 1
        !'_12/___', '_13/___', '_14/___', '_15/___'
        inverse_neighbor(7,3) = 13
        inverse_neighbor(8,3) = 14
        inverse_neighbor(9,3) = 11
        inverse_neighbor(10,3) = 12
        !'_62/___', '_63/___', '_64/___', '_65/___'
        inverse_neighbor(11,3) = 9
        inverse_neighbor(12,3) = 10
        inverse_neighbor(13,3) = 7
        inverse_neighbor(14,3) = 8
        !'_23/___', '_25/___'
        inverse_neighbor(15,3) = 18
        inverse_neighbor(16,3) = 17
        !'_43/___', '_45/___'
        inverse_neighbor(17,3) = 16
        inverse_neighbor(18,3) = 15
        !'123/___', '134/___', '145/___', '152/___'
        inverse_neighbor(19,3) = 25
        inverse_neighbor(20,3) = 26
        inverse_neighbor(21,3) = 23
        inverse_neighbor(22,3) = 24
        !'623/___', '634/___', '645/___', '652/___'
        inverse_neighbor(23,3) = 21
        inverse_neighbor(24,3) = 22
        inverse_neighbor(25,3) = 19
        inverse_neighbor(26,3) = 20
        !'__1/123', '__1/134', '__1/145', '__1/152'
        inverse_neighbor(27,3) = 47
        inverse_neighbor(28,3) = 48
        inverse_neighbor(29,3) = 49
        inverse_neighbor(30,3) = 50
        !'__2/123', '__2/623', '__2/152', '__2/652'
        inverse_neighbor(31,3) = 39
        inverse_neighbor(32,3) = 40
        inverse_neighbor(33,3) = 41
        inverse_neighbor(34,3) = 42
        !'__3/123', '__3/623', '__3/134', '__3/634'
        inverse_neighbor(35,3) = 45
        inverse_neighbor(36,3) = 46
        inverse_neighbor(37,3) = 43
        inverse_neighbor(38,3) = 44
        !'__4/134', '__4/634', '__4/145', '__4/645'
        inverse_neighbor(39,3) = 31
        inverse_neighbor(40,3) = 32
        inverse_neighbor(41,3) = 33
        inverse_neighbor(42,3) = 34
        !'__5/145', '__5/645', '__5/152', '__5/652'
        inverse_neighbor(43,3) = 37
        inverse_neighbor(44,3) = 38
        inverse_neighbor(45,3) = 35
        inverse_neighbor(46,3) = 36
        !'__6/623', '__6/634', '__6/645', '__6/652'
        inverse_neighbor(47,3) = 27
        inverse_neighbor(48,3) = 28
        inverse_neighbor(49,3) = 29
        inverse_neighbor(50,3) = 30
        !'_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152'
        inverse_neighbor(51,3) = 63
        inverse_neighbor(52,3) = 64
        inverse_neighbor(53,3) = 66!65
        inverse_neighbor(54,3) = 65!66
        inverse_neighbor(55,3) = 59
        inverse_neighbor(56,3) = 60
        inverse_neighbor(57,3) = 62
        inverse_neighbor(58,3) = 61
        !'_62/623', '_62/652', '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652'
        inverse_neighbor(59,3) = 55
        inverse_neighbor(60,3) = 56
        inverse_neighbor(61,3) = 58
        inverse_neighbor(62,3) = 57
        inverse_neighbor(63,3) = 51
        inverse_neighbor(64,3) = 52
        inverse_neighbor(65,3) = 54!53
        inverse_neighbor(66,3) = 53!54
        !'_23/123', '_23/623', '_25/152', '_25/652'
        inverse_neighbor(67,3) = 73
        inverse_neighbor(68,3) = 74
        inverse_neighbor(69,3) = 71
        inverse_neighbor(70,3) = 72
        !'_43/134', '_43/634', '_45/145', '_45/645'
        inverse_neighbor(71,3) = 69
        inverse_neighbor(72,3) = 70
        inverse_neighbor(73,3) = 67
        inverse_neighbor(74,3) = 68

        ghost_nodes_module_ready = .true.
    endif

end subroutine

end module module_MPI
