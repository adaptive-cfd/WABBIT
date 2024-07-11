module module_MPI

    use mpi
    use module_params
    use module_timing
    use module_wavelets
    use module_treelib
    use module_helpers

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas
    use module_forestMetaData

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! everything is save by default
    SAVE

    ! send/receive buffer, integer and real
    ! allocate in init substep not in synchronize subroutine, to avoid slow down when using
    ! large numbers of processes and blocks per process, when allocating on every call to the routine
    integer(kind=ik), PARAMETER   :: S_META_FULL = 7  ! how many metadata entries we collect in the logik-loop
#ifdef DEV
    integer(kind=ik), PARAMETER   :: S_META_SEND = 5  ! how many metadata entries will be send, plus rank
#else
    integer(kind=ik), PARAMETER   :: S_META_SEND = 4  ! how many metadata entries will be send
#endif
    integer(kind=ik), allocatable :: meta_send_counter(:), meta_recv_counter(:)
    integer(kind=ik), allocatable :: data_recv_counter(:), data_send_counter(:)
    integer(kind=ik), allocatable :: meta_send_all(:)
    real(kind=rk), allocatable    :: rData_sendBuffer(:), rData_recvBuffer(:)
    !
    integer(kind=ik), allocatable :: send_request(:)
    integer(kind=ik), allocatable :: recv_request(:)

    ! arrays to iterate over messages we send to the other mpiranks.
    integer(kind=ik), allocatable :: int_pos(:)
    integer(kind=ik), allocatable :: real_pos(:)

    ! internally, we flatten the ghost nodes layers to a line. this is stored in
    ! this buffer (NOTE max size is (blocksize)*(ghost nodes size + 1)*(number of datafields))
    real(kind=rk), allocatable    :: line_buffer(:)

    ! restricted/predicted data buffer
    real(kind=rk), allocatable    :: res_pre_data(:,:,:,:), tmp_block(:,:,:,:)

    ! it is faster to use named consts than strings, although strings are nicer to read
    integer(kind=ik), PARAMETER   :: SENDER = 1, RECVER = 2, RESPRE = 3

    ! We frequently need to know the indices of a patch for ghost nodes or parts of interor. Thus we save them
    ! once in a large, module-global array (which is faster than computing it every time with tons
    ! of IF-THEN clauses).
    ! This arrays indices are:
    ! ijkPatches([start,end], [dir (x,y,z)], [neighborhood / relation], [level-diff], [sender/receiver/up-downsampled])
    ! The third index is 1:74 for ghost relations, 0 for whole block, -1:-8 for SC decimation
    integer(kind=ik), dimension(1:2, 1:3, -8:74, -1:1, 1:3) :: ijkPatches

    ! it is useful to keep a named constant for the dimensionality here (we use
    ! it to access e.g. two/three D arrays in inverse_neighbor)
    integer(kind=ik) :: dim

    ! we use this flag to call the allocation routine only once.
    logical :: ghost_nodes_module_ready = .false.

    ! As we always loop over hvy_id, only one block needs to be filtered at a time, so if we filter a
    ! new block we delete the old filtered values from the previous one.
    real(kind=rk), allocatable :: hvy_restricted(:, :, :, :)
    real(kind=rk), allocatable :: hvy_predicted(:, :, :, :)
    ! We need to save which block currently is filtered
    integer(kind=ik) :: restricted_hvy_ID
    integer(kind=ik) :: predicted_hvy_ID



    ! two shift parameters (asymmetric and symmetric ) used for selection of interpolation
    ! bounds on sender side. used to avoid one-sided interpolation if desired. They're
    ! global so we can save them easily to the ghosts_bounds.dat file
    integer(kind=ik) :: A, S

!---------------------------------------------------------------------------------------------
! public parts of this module

    PUBLIC :: sync_ghosts_tree, sync_ghosts_RHS_tree, sync_TMP_from_MF, sync_TMP_from_all, sync_SCWC_from_MC
    PUBLIC :: blocks_per_mpirank, synchronize_lgt_data, reset_ghost_nodes, init_ghost_nodes, move_mallat_patch_block, family_setup_patches, xfer_ensure_correct_buffer_size
    PUBLIC :: coarse_extension_modify_tree, coarse_extension_reconstruct_tree, xfer_block_data, prepare_update_family_metadata

    ! we set up a table that gives us directly the inverse neighbor relations.
    ! it is filled (once) in init_ghost_nodes
    PUBLIC :: inverse_neighbor
    integer(kind=ik), dimension(1:74,2:3) :: inverse_neighbor


contains

#include "synchronize_ghosts_generic.f90"
#include "blocks_per_mpirank.f90"
#include "synchronize_lgt_data.f90"
#include "reset_ghost_nodes.f90"
#include "calc_data_bounds.f90"
#include "restrict_predict_data.f90"
#include "reconstruction_step.f90"
#include "xfer_block_data.f90"

!! initialize ghost nodes module. allocate buffers and create data bounds array,
!! which we use to rapidly identify a ghost nodes layer
subroutine init_ghost_nodes( params )
    implicit none
    type (type_params), intent(in) :: params
    ! local variables
    integer(kind=ik) :: buffer_N_int, buffer_N, g, Neqn, number_blocks, rank
    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: ineighbor, Nneighbor, leveldiff, Ncpu
    integer(kind=ik) ::  j, rx0, rx1, ry0, ry1, rz0, rz1, sx0, sx1, sy0, sy1, sz0, sz1
    integer(kind=ik) :: i, k, status(1:4)
    integer(kind=ik) :: ijkrecv(2,3)
    integer(kind=ik) :: ijkbuffer(2,3)
    integer(kind=ik) :: ijksend(2,3)

    ! on second call, nothing happens
    if (.not. ghost_nodes_module_ready) then
        number_blocks   = params%number_blocks
        Bs              = params%Bs
        g               = params%g
        ! HACK: in the design phase, we thought to always sync neqn commponents
        ! but in some cases we end up synching more!
        Neqn            = N_MAX_COMPONENTS
        rank            = params%rank
        Ncpu            = params%number_procs

        if (rank==0) then
            write(*,'("---------------------------------------------------------")')
            write(*,'("             ╱▔▔▔▔▔▔╲ ╭━━━━━━━━━━╮")')
            write(*,'("            ▕ ╭━╮╭━╮ ▏┃GHOST-INIT┃")')
            write(*,'("            ▕ ┃╭╯╰╮┃ ▏╰┳━━━━━━━━━╯")')
            write(*,'("            ▕ ╰╯╭╮╰╯ ▏ ┃ ")')
            write(*,'("            ▕   ┃┃   ▏━╯ ")')
            write(*,'("            ▕   ╰╯   ▏" )')
            write(*,'("            ▕╱╲╱╲╱╲╱╲▏")')
            write(*,'("---------------------------------------------------------")')
            write(*,'("GHOST-INIT: We can synchronize at most N_MAX_COMPONENTS=",i2)') N_MAX_COMPONENTS
            write(*,'("GHOST-INIT: g= ",i2, " , g_RHS= ", i2)') params%g, params%g_rhs
        endif

        if ( params%dim==3 ) then
            if (g>=(Bs(1)+1)/2 .or. g>=(Bs(2)+1)/2 .or. g>=(Bs(3)+1)/2) then
              call abort(921151369, "Young skywalker, you failed at set g>=(Bs+1)/2 (in at least one direction) which implies &
              & that the ghost nodes layer can span beyond an entire finer block. Either decrease &
              & number_ghost_nodes or increase number_block_nodes.")
            endif
        else
            if (g>=(Bs(1)+1)/2 .or. g>=(Bs(2)+1)/2) then
              call abort(821151369, "Young skywalker, you failed at set g>=(Bs+1)/2 (in at least one direction) which implies &
              & that the ghost nodes layer can span beyond an entire finer block. Either decrease &
              & number_ghost_nodes or increase number_block_nodes.")
            endif
        endif

#ifdef DEV
        if (number_blocks<1) then
            call abort(2422021, "Ghost setup called but number_blocks undefnd - call init_ghost_nodes AFTER allocate_tree. Maybe you forgot --memory?")
        endif
#endif

        ! synchronize buffer length
        ! assume: all blocks are used, all blocks have external neighbors,
        ! max number of active neighbors: 2D = 12, 3D = 56
        dim = params%dim
        if ( dim == 3 ) then
            !---3d---3d---
            ! per neighborhood relation, we send 6 integers as metadata in the int_buffer
            ! at most, we can have 56 neighbors ACTIVE per block
            buffer_N_int = number_blocks * 56 * 6
            ! how many possible neighbor relations are there?
            Nneighbor = 74
        else
            !---2d---2d---
            ! per neighborhood relation, we send 6 integers as metadata in the int_buffer
            ! at most, we can have 12 neighbors ACTIVE per block
            buffer_N_int = number_blocks * 12 * 6
            ! how many possible neighbor relations are there?
            Nneighbor = 16
        end if

        ! size of ghost nodes buffer. Note this contains only the ghost nodes layer
        ! for all my blocks
        if ( dim == 3 ) then
            buffer_N = number_blocks * Neqn * ( (Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g) - (Bs(1)*Bs(2)*Bs(3)) )
        else
            ! 2D case
            buffer_N = number_blocks * Neqn * ( (Bs(1)+2*g)*(Bs(2)+2*g) - (Bs(1)*Bs(2)) )
        end if

        !-----------------------------------------------------------------------
        ! allocate auxiliary memory
        !-----------------------------------------------------------------------
        ! allocate synch buffer
        if (rank==0) then
            write(*,'("GHOST-INIT: Attempting to allocate the ghost-sync-buffer.")')

            write(*,'("GHOST-INIT: buffer_N_int=",i12," buffer_N=",i12)') &
            buffer_N_int, buffer_N

            write(*,'("GHOST-INIT: On each MPIRANK, Int  buffer:", f9.4, "GB")') &
                2.0*dble(buffer_N_int)*8e-9

            write(*,'("GHOST-INIT: On each MPIRANK, Real buffer:", f9.4, "GB")') &
                2.0*dble(buffer_N)*8e-9
            write(*,'("---------------- allocating now ----------------")')
        endif

        ! wait now so that if allocation fails, we get at least the above info
        call MPI_barrier( WABBIT_COMM, status(1))

        allocate( meta_send_all(1:buffer_N_int), stat=status(1) )
        allocate( rData_sendBuffer(1:buffer_N+buffer_N_int+Ncpu), stat=status(3) )
        allocate( rData_recvBuffer(1:buffer_N+buffer_N_int+Ncpu), stat=status(4) )

        if (maxval(status) /= 0) call abort(999999, "Buffer allocation failed. Not enough memory?")

        if (rank==0) then
            write(*,'("GHOST-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
            "rData_sendBuffer", shape(rData_sendBuffer)

            write(*,'("GHOST-INIT: on each mpirank, Allocated ",A25," SHAPE=",7(i9,1x))') &
            "rData_recvBuffer", shape(rData_recvBuffer)

            write(*,'("GHOST-INIT: on each mpirank, rData_sendBuffer size is",f9.4," GB ")') &
            product(real(shape(rData_sendBuffer)))*8e-9

            write(*,'("GHOST-INIT: on each mpirank, rData_recvBuffer size is",f9.4," GB ")') &
            product(real(shape(rData_recvBuffer)))*8e-9
        endif

        ! thanks to the mix of Fortran and MPI this is quite a nightmare with 0- or 1-based arrays
        allocate( send_request(1:2*Ncpu) )
        allocate( recv_request(1:2*Ncpu) )
        allocate( line_buffer( 1:Neqn*product(Bs(1:params%dim)+2*g) ) )
        allocate( int_pos(0:Ncpu-1), real_pos(0:Ncpu-1))
        allocate( data_recv_counter(0:Ncpu-1), data_send_counter(0:Ncpu-1) )
        allocate( meta_recv_counter(0:Ncpu-1), meta_send_counter(0:Ncpu-1) )

        ! wait now so that if allocation fails, we get at least the above info
        call MPI_barrier( WABBIT_COMM, status(1))

        !-----------------------------------------------------------------------
        ! set up constant arrays
        !-----------------------------------------------------------------------
        ! We frequently need to know the indices of a ghost nodes patch. Thus we save them
        ! once in a module-global array (which is faster than computing it every time with tons
        ! of IF-THEN clauses).
        ! This arrays indices are:
        ! ijkPatches([start,end], [dir], [ineighbor], [leveldiff], [isendrecv])
        call ghosts_setup_patches(params, gminus=params%g, gplus=params%g, output_to_file=.true.)

        ! set up table with inverse neighbor relations
        if (.true.) then ! for folding only
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
        endif

        ! this routine is not performance-critical
        call MPI_barrier( WABBIT_COMM, status(1))

        if (rank==0) then
            write(*,'("---------------------------------------------------------")')
            write(*,'("                     GHOST-INIT complete :)")')
            write(*,'("---------------------------------------------------------")')
        endif

        ghost_nodes_module_ready = .true.
    endif

end subroutine

subroutine xfer_ensure_correct_buffer_size(params, hvy_block)
    implicit none
    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)

    integer(kind=ik) :: i,j,k, Bs(1:3), g, nx, ny, nz, nc

    if (allocated(res_pre_data)) deallocate(res_pre_data)
    if (allocated(tmp_block)) deallocate(tmp_block)

    Bs = params%Bs
    g  = params%g
    nx = size(hvy_block,1)
    ny = size(hvy_block,2)
    nz = size(hvy_block,3)
    nc = size(hvy_block,4)

    ! now we know how large the patches are we'd like to store in the RESPRE buffer
    i = maxval( ijkPatches(2,1,:,:,RESPRE) )*2
    j = maxval( ijkPatches(2,2,:,:,RESPRE) )*2
    k = maxval( ijkPatches(2,3,:,:,RESPRE) )*2

    allocate( res_pre_data( i, j, k, nc) )

    if (params%dim == 3) then
        allocate( tmp_block( Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, nc) )
    else
        allocate( tmp_block( Bs(1)+2*g, Bs(2)+2*g, 1, nc) )
    endif

    ! initialize that no block is currently filtered
    restricted_hvy_ID = -1
    predicted_hvy_ID = -1

    if (allocated(hvy_restricted)) then
        if (nc > size(hvy_restricted, 4)) deallocate(hvy_restricted)
    endif
    if (.not. allocated(hvy_restricted)) allocate(hvy_restricted(1:nx, 1:ny, 1:nz, 1:nc) )

    if (allocated(hvy_predicted)) then
        if (nc > size(hvy_predicted, 4)) deallocate(hvy_predicted)
    endif
    if (.not. allocated(hvy_restricted)) then
        if (params%dim == 3) allocate(hvy_restricted(1:nx*2, 1:ny*2, 1:nz*2, 1:nc) )
        if (params%dim == 2) allocate(hvy_restricted(1:nx*2, 1:ny*2, 1, 1:nc) )
    endif


end subroutine


! setup the ghost patches (i.e. the indices of the sender/receiver parts
! for any neighborhood relation), where a smaller subset of the ghost nodes
! can be synced (the g you pass must be <= params%g)
subroutine ghosts_setup_patches(params, gminus, gplus, output_to_file)
    implicit none
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in) :: gminus, gplus
    logical, intent(in) :: output_to_file

    integer(kind=ik) :: N_neighbors, level_diff, neighborhood, i, j, k
    integer(kind=ik) :: ijkrecv(2,3)
    integer(kind=ik) :: ijkbuffer(2,3)
    integer(kind=ik) :: ijksend(2,3)

    ! full reset of all patch definitions. Note we set 1 not 0
    ijkPatches = 1

    N_neighbors = 74
    if (params%dim==2) then
        N_neighbors=16
    endif

    ! set all neighbour relations
    do neighborhood = 1, N_neighbors
        do level_diff = -1, 1
            !---------larger ghost nodes (used for refinement and coarsening and maybe RHS)------------
            ! The receiver bounds are hard-coded with a huge amount of tedious indices.
            call set_recv_bounds( params, ijkrecv, neighborhood, level_diff, gminus, gplus)
            ijkPatches(1:2, 1:3, neighborhood, level_diff, RECVER) = ijkrecv

            ! Luckily, knowing the receiver bounds, we can compute the sender bounds as well as
            ! the indices in the buffer arrays, where we temporarily store patches, if they need to be
            ! up or down sampled.
            call compute_sender_buffer_bounds(params, ijkrecv, ijksend, ijkbuffer, neighborhood, level_diff, &
                 gminus, gplus, output_to_file )
            ijkPatches(1:2, 1:3, neighborhood, level_diff, SENDER) = ijksend
            ijkPatches(1:2, 1:3, neighborhood, level_diff, RESPRE) = ijkbuffer

            ! apply a shift (if we sync only a subset of the layer)
            !
            ! g g g g g g g g g g g g        g g g g g g g g g g g g
            ! g g g g g g g g g g g g        g g g g g g g g g g g g
            ! g g g g g g g g g g g g        g g S S S S S S S S g g
            ! g g g u u u u u u g g g        g g S u u u u u u S g g
            ! g g g u u u u u u g g g        g g S u u u u u u S g g
            ! g g g u u u u u u g g g        g g S u u u u u u S g g
            ! g g g u u u u u u g g g        g g S u u u u u u S g g
            ! g g g u u u u u u g g g        g g S u u u u u u S g g
            ! g g g g g g g g g g g g        g g S S S S S S S S g g
            ! g g g g g g g g g g g g        g g g g g g g g g g g g
            ! g g g g g g g g g g g g        g g g g g g g g g g g g
            !
            ! g: ghost u: interior s: synced
            ! here, params%g=3 and g=1

            ijkPatches(1:2,1:dim, neighborhood, level_diff,RECVER) = ijkPatches(1:2,1:dim, neighborhood, level_diff,RECVER) + (params%g-gminus)
            ijkPatches(1:2,1:dim, neighborhood, level_diff,SENDER) = ijkPatches(1:2,1:dim, neighborhood, level_diff,SENDER) + (params%g-gminus)
        enddo
    enddo


#ifdef DEV
    ! this output can be plotted using the python script
    if ((params%rank==0) .and. output_to_file) then
        open(16,file='ghost_bounds.dat',status='replace')
        write(16,'(i3)') params%Bs(1)
        write(16,'(i3)') params%Bs(2)
        write(16,'(i3)') params%Bs(3)
        write(16,'(i3)') params%g
        write(16,'(i3)') S
        write(16,'(i3)') A
        do neighborhood = 1, N_neighbors
            do level_diff = -1, 1
                do i = 1, 2
                    do j = 1, 3
                        do k = 1, 3
                            write(16,'(i3)') ijkPatches(i, j, neighborhood, level_diff, k)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        close(16)
    endif
#endif
end subroutine



! setup the family patches (i.e. the indices of the sender/receiver parts
! for any neighborhood relation), this depends on the filter being used
subroutine family_setup_patches(params, output_to_file)
    implicit none
    type (type_params), intent(in) :: params
    logical, intent(in) :: output_to_file

    integer(kind=ik) :: N_family, level_diff, family, i, j, k, g
    integer(kind=ik) :: ijkrecv(2,3)
    integer(kind=ik) :: ijkbuffer(2,3)
    integer(kind=ik) :: ijksend(2,3)

    ! full reset of all patch definitions. Note we set 1 not 0
    ijkPatches = 1

    N_family = 8
    if (params%dim==2) then
        N_family=4
    endif
    g = params%g

    ! set full block relation
    ! The receiver bounds are hard-coded with a huge amount of tedious indices, they are identical to senders
    call set_recv_bounds( params, ijkrecv, 0, 0, g, g)
    ijkPatches(1:2,1:3, 0, 0, RECVER) = ijkrecv
    ijkPatches(1:2,1:3, 0, 0, SENDER) = ijkrecv

    ! set mother/daughter relations
    ! sender / receiver relations between the level differences are inverted
    do family = -1, -N_family, -1
        ! Fine sender or receiver: affects sc in mallat-ordering, connected to relation -1
        ! This takes different border into account
        call set_recv_bounds( params, ijkrecv, family, -1, g, g)
        ijkPatches(1:2,1:3, family,  1, SENDER) = ijkrecv
        ijkPatches(1:2,1:3, family, -1, RECVER) = ijkrecv

        ! Coarse sender or receiver: affects specific part of block
        call set_recv_bounds( params, ijkrecv, family, 1, g, g)
        ijkPatches(1:2,1:3, family, -1, SENDER) = ijkrecv
        ijkPatches(1:2,1:3, family,  1, RECVER) = ijkrecv
    enddo

#ifdef DEV
    ! this output can be plotted using the python script
    if ((params%rank==0) .and. output_to_file) then
        open(16,file='family_bounds.dat',status='replace')
        write(16,'(i3)') params%Bs(1)
        write(16,'(i3)') params%Bs(2)
        write(16,'(i3)') params%Bs(3)
        write(16,'(i3)') params%g
        write(16,'(i3)') S
        write(16,'(i3)') A
        do family = -1, -N_family, -1
            do level_diff = -1, 1
                do i = 1, 2
                    do j = 1, 3
                        do k = 1, 3
                            write(16,'(i3)') ijkPatches(i, j, family, level_diff, k)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        close(16)
    endif
#endif

end subroutine

end module module_MPI
