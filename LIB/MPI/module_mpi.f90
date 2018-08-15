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

    ! We require two stages: first, we fill all ghost nodes which are simple copy,
    ! then in the second stage we can use interpolation and fill the remaining ones.
    ! In order not to send ALL data in both stages, we allocate one buffer for each stage.
    integer(kind=ik) :: Nstages = 2

    ! send/receive buffer, integer and real
    ! allocate in init substep not in synchronize subroutine, to avoid slow down when using
    ! large numbers of processes and blocks per process, when allocating on every call to the routine
    integer(kind=ik), allocatable :: int_send_buffer(:,:), int_recv_buffer(:,:)
    integer(kind=ik), allocatable :: int_recv_counter(:,:), int_send_counter(:,:)
    integer(kind=ik), allocatable :: recv_counter(:,:), send_counter(:,:)
    integer(kind=ik), allocatable :: real_pos(:,:)
    real(kind=rk), allocatable    :: new_send_buffer(:,:), new_recv_buffer(:,:)

    ! an array to count how many messages we send to the other mpiranks.
    integer(kind=ik), allocatable :: communication_counter(:,:), int_pos(:,:)

    ! internally, we flatten the ghost nodes layers to a line. this is stored in
    ! this buffer (NOTE max size is (blocksize)*(ghost nodes size + 1)*(number of datafields))
    real(kind=rk), allocatable    :: line_buffer(:)

    ! restricted/predicted data buffer
    real(kind=rk), allocatable    :: res_pre_data(:,:,:,:), tmp_block(:,:,:,:)

    ! it is faster to use named consts than strings, although strings are nicer to read
    integer(kind=ik), PARAMETER   :: EXCLUDE_REDUNDANT = 1, INCLUDE_REDUNDANT = 2, ONLY_REDUNDANT = 3
    integer(kind=ik), PARAMETER   :: SENDER = 1, RECVER = 2, RESPRE = 3

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
    PUBLIC :: synchronize_ghosts_generic_sequence, init_ghost_nodes, check_unique_origin

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


!! initialize ghost nodes module. allocate buffers and create data bounds array,
!! which we use to rapidly identify a ghost nodes layer
subroutine init_ghost_nodes( params )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in) :: params
    ! local variables
    integer(kind=ik) :: buffer_N_int, buffer_N, Bs, g, Neqn, number_blocks, rank
    integer(kind=ik) :: ineighbor, Nneighbor, leveldiff, idata_bounds_type, Ncpu
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
        Ncpu            = params%number_procs

        if (rank==0) write(*,'("---------------------------------------------------------")')
        if (rank==0) write(*,'("                     GHOST-INIT ")')
        if (rank==0) write(*,'("---------------------------------------------------------")')

        if (g>=(Bs+1)/2) then
            call abort(921151369, "Young skywalker, you failed at set g>=(Bs+1)/2 which implies &
            & that the ghost nodes layer can span beyond an entire finer block. Either decrease &
            & number_ghost_nodes or increase number_block_nodes.")
        endif

        ! synchronize buffer length
        ! assume: all blocks are used, all blocks have external neighbors,
        ! max neighbor number: 2D = 12, 3D = 56
        ! max neighborhood size, 2D: (Bs+g+1)*(g+1)
        ! max neighborhood size, 3D: (Bs+g+1)*(g+1)*(g+1)
        if ( params%threeD_case ) then
            !---3d---3d---
            ! space dimensions: used in the static arrays as index
            dim = 3
            ! per neighborhood relation, we send 5 integers as metadata in the int_buffer
            ! at most, we can have 56 neighbors ACTIVE per block
            buffer_N_int = number_blocks * 56 * 5
            ! how many possible neighbor relations are there?
            Nneighbor = 74

            allocate( tmp_block( Bs+2*g, Bs+2*g, Bs+2*g, Neqn) )
        else
            !---2d---2d---
            ! space dimensions: used in the static arrays as index
            dim = 2
            ! per neighborhood relation, we send 5 integers as metadata in the int_buffer
            ! at most, we can have 12 neighbors ACTIVE per block
            buffer_N_int = number_blocks * 12 * 5
            ! how many possible neighbor relations are there?
            Nneighbor = 16

            allocate( tmp_block( Bs+2*g, Bs+2*g, 1, Neqn) )
        end if

        ! size of ghost nodes buffer. Note this contains only the ghost nodes layer
        ! for all my blocks. previous versions allocated one of those per "friend"
        buffer_N = number_blocks * Neqn * ((Bs+2*g)**dim - Bs**dim)

        !-----------------------------------------------------------------------
        ! allocate auxiliary memory
        !-----------------------------------------------------------------------
        ! allocate synch buffer
        if (rank==0) then
            write(*,'("GHOSTS-INIT: Attempting to allocate the ghost-sync-buffer.")')

            write(*,'("GHOSTS-INIT: buffer_N_int=",i12," buffer_N=",i12," Nstages=",i1)') &
            buffer_N_int, buffer_N, Nstages

            write(*,'("GHOSTS-INIT: On each MPIRANK, Int  buffer:", f9.4, "GB")') &
                2.0*dble(buffer_N_int)*dble(Ncpu)*dble(Nstages)*8e-9

            write(*,'("GHOSTS-INIT: On each MPIRANK, Real buffer:", f9.4, "GB")') &
                2.0*dble(buffer_N)*dble(Nstages)*8e-9
            write(*,'("---------------- allocating now ----------------")')
        endif

        ! wait now so that if allocation fails, we get at least the above info
        call MPI_barrier( WABBIT_COMM, status(1))

        allocate( int_send_buffer( 1:buffer_N_int, 1:Nstages), stat=status(1) )
        allocate( int_recv_buffer( 1:buffer_N_int, 1:Nstages), stat=status(2) )
        allocate( new_send_buffer( 1:buffer_N, 1:Nstages), stat=status(3) )
        allocate( new_recv_buffer( 1:buffer_N, 1:Nstages), stat=status(4) )

        if (maxval(status) /= 0) call abort(999999, "Buffer allocation failed. Not enough memory?")

        if (rank==0) then

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "new_send_buffer", shape(new_send_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "new_recv_buffer", shape(new_recv_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "int_send_buffer", shape(int_send_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
             "int_recv_buffer", shape(int_recv_buffer)

            write(*,'("GHOSTS-INIT: on each mpirank, new_send_buffer size is",f9.4," GB ")') &
             product(real(shape(new_send_buffer)))*8e-9

             write(*,'("GHOSTS-INIT: on each mpirank, new_recv_buffer size is",f9.4," GB ")') &
              product(real(shape(new_recv_buffer)))*8e-9

            write(*,'("GHOSTS-INIT: on each mpirank, int_send_buffer size is",f9.4," GB ")') &
             product(real(shape(int_send_buffer)))*4e-9

             write(*,'("GHOSTS-INIT: on each mpirank, int_recv_buffer size is",f9.4," GB ")') &
              product(real(shape(int_recv_buffer)))*4e-9
        endif

        ! this is a list of communications with all other procs
        allocate( communication_counter(1:Ncpu, 1:Nstages) )
        allocate( int_pos(1:Ncpu, 1:Nstages) )
        allocate( real_pos(1:Ncpu, 1:Nstages) )
        allocate( line_buffer( Neqn*(Bs+2*g)**(dim) ) )

        allocate( recv_counter(0:Ncpu-1, 1:Nstages), send_counter(0:Ncpu-1, 1:Nstages) )
        allocate( int_recv_counter(0:Ncpu-1, 1:Nstages), int_send_counter(0:Ncpu-1, 1:Nstages) )

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
