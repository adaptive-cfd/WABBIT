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

    ! We need to set the BC into the ghost patches after synching, for now it is done in the synch routine itself so we need access to physics moduel to do that
    ! JB Todo: Check if this is clever or we should do it otherwise
    use module_physics_metamodule

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! everything is save by default
    SAVE

    ! send/receive buffer, integer and real
    ! allocate in init substep not in synchronize subroutine, to avoid slow down when using
    ! large numbers of processes and blocks per process, when allocating on every call to the routine
    integer(kind=ik), PARAMETER   :: S_META_FULL = 7  ! how many metadata entries we collect in the logic-loop
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
    ! ijkPatches([start,end], [dir (x,y,z)], [neighborhood / relation (with lvl_diff)], [sender/receiver/up-downsampled])
    ! The third index is 1:56*3 for ghost relations for all three level-diffs, 0 for whole block, -1:-16 for SC decimation (with lvl_diff)
    integer(kind=ik), allocatable :: ijkPatches(:, :, :, :)  ! size=(1:2, 1:3, -16:56*3, 1:3)


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
    PUBLIC :: coarse_extension_modify, coarse_extension_reconstruct_tree, xfer_block_data, prepare_update_family_metadata
    PUBLIC :: block_has_valid_neighbor, block_is_root, block_is_leaf


contains

#include "synchronize_ghosts_generic.f90"
#include "blocks_per_mpirank.f90"
#include "synchronize_lgt_data.f90"
#include "reset_ghost_nodes.f90"
#include "calc_data_bounds.f90"
#include "restrict_predict_data.f90"
#include "reconstruction_step.f90"
#include "xfer_block_data.f90"
#include "block_relations.f90"

!! initialize ghost nodes module. allocate buffers and create data bounds array,
!! which we use to rapidly identify a ghost nodes layer
subroutine init_ghost_nodes( params )
    implicit none
    type (type_params), intent(in) :: params
    ! local variables
    integer(kind=ik) :: buffer_N_int, buffer_N, g, Neqn, number_blocks, rank
    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: ineighbor, leveldiff, Ncpu
    integer(kind=ik) ::  j, rx0, rx1, ry0, ry1, rz0, rz1, sx0, sx1, sy0, sy1, sz0, sz1
    integer(kind=ik) :: i, k, status(1:4)
    integer(kind=ik) :: ijkrecv(2,3)
    integer(kind=ik) :: ijkbuffer(2,3)
    integer(kind=ik) :: ijksend(2,3)
    real(kind=rk)    :: memory_this, memory_total

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
            write(*,'("╔",78("═"),"╗")')
            write(*,'("║", 19(" "), A, A41)')  "╱▔▔▔▔▔▔╲ ╭━━━━━━━━━━╮", "║"
            write(*,'("║", 18(" "), A, A41)') "▕ ╭━╮╭━╮ ▏┃GHOST-INIT┃", "║"
            write(*,'("║", 18(" "), A, A41)') "▕ ┃╭╯╰╮┃ ▏╰┳━━━━━━━━━╯", "║"
            write(*,'("║", 18(" "), A, A41)') "▕ ╰╯╭╮╰╯ ▏ ┃          ", "║"
            write(*,'("║", 18(" "), A, A41)') "▕   ┃┃   ▏━╯          ", "║"
            write(*,'("║", 18(" "), A, A41)') "▕   ╰╯   ▏            ", "║"
            write(*,'("║", 18(" "), A, A41)') "▕╱╲╱╲╱╲╱╲▏            ", "║"
            write(*,'("╚",78("═"),"╝")')
            write(*,'("GHOST-INIT: We can synchronize at most N_MAX_COMPONENTS=",i2)') N_MAX_COMPONENTS
            write(*,'("GHOST-INIT: g= ",i2, " , g_RHS= ", i2)') params%g, params%g_rhs
        endif

        if (g>=(Bs(1)+1)/2 .or. g>=(Bs(2)+1)/2 .or. (g>=(Bs(3)+1)/2 .and. params%dim==3)) then
            call abort(921151369, "Young skywalker, you failed at set g>=(Bs+1)/2 (in at least one direction) which implies &
            & that the ghost nodes layer can span beyond an entire finer block. Either decrease &
            & number_ghost_nodes or increase number_block_nodes.")
        endif

#ifdef DEV
        if (number_blocks<1) then
            call abort(2422021, "Ghost setup called but number_blocks undefnd - call init_ghost_nodes AFTER allocate_tree. Maybe you forgot --memory?")
        endif
#endif

        ! synchronize buffer length
        ! assume: all blocks are used, all blocks have external neighbors,
        ! max number of active neighbors per stage: 2D = 12+8+8, 3D = 56+24+24, being on lvl J+1,J,J-1
        ! per neighborhood relation, we save S_META_FULL and send S_META_SEND integers as metadata in the int_buffer
        if ( params%dim == 3 ) then
            !---3d---3d---
            buffer_N_int = number_blocks * (56+24+24)
        else
            !---2d---2d---
            buffer_N_int = number_blocks * (12+8+8)
        end if

        ! size of ghost nodes buffer. Note this contains only the ghost nodes layer
        ! for all my blocks
        if ( params%dim == 3 ) then
            buffer_N = number_blocks * Neqn * ( (Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g) - (Bs(1)*Bs(2)*Bs(3)) )
        else
            ! 2D case
            buffer_N = number_blocks * Neqn * ( (Bs(1)+2*g)*(Bs(2)+2*g) - (Bs(1)*Bs(2)) )
        end if
        memory_total = 0.0

        !-----------------------------------------------------------------------
        ! allocate auxiliary memory
        !-----------------------------------------------------------------------
        ! allocate synch buffer
        if (rank==0) then
            write(*,'("GHOST-INIT: Allocating ghost-sync-buffers.")')

            ! write(*,'("GHOST-INIT: buffer_N_int=",i12," buffer_N=",i12)') &
            ! buffer_N_int, buffer_N

            ! write(*,'("GHOST-INIT: On each MPIRANK, Int  buffer:", f9.4, "GB")') &
            !     2.0*dble(buffer_N_int)*8e-9

            ! write(*,'("GHOST-INIT: On each MPIRANK, Real buffer:", f9.4, "GB")') &
            !     2.0*dble(buffer_N)*8e-9
            ! write(*,'(20("─"), A20, 40("─"))') "   allocating now   "
        endif

        ! wait now so that if allocation fails, we get at least the above info
        call MPI_barrier( WABBIT_COMM, status(1))

        allocate( meta_send_all(1:buffer_N_int*S_META_FULL), stat=status(1) )
        allocate( rData_sendBuffer(1:buffer_N+buffer_N_int*S_META_SEND+Ncpu), stat=status(2) )
        allocate( rData_recvBuffer(1:buffer_N+buffer_N_int*S_META_SEND+Ncpu), stat=status(3) )
        memory_this = (product(real(shape(meta_send_all)))+product(real(shape(rData_sendBuffer)))+product(real(shape(rData_recvBuffer))))*8.0e-9
        memory_total = memory_total + memory_this

        if (maxval(status(1:3)) /= 0) call abort(999999, "Buffer allocation failed. Not enough memory?")

        if (rank==0) then
            write(*,'("GHOST-INIT: ALLOCATED ",A19," MEM=",f8.4," GB per rank, shape=",7(i9,1x))') &
            "meta_send_all", product(real(shape(meta_send_all)))*8e-9, shape(meta_send_all)

            write(*,'("GHOST-INIT: ALLOCATED ",A19," MEM=",f8.4," GB per rank, shape=",7(i9,1x))') &
            "rData_sendBuffer", product(real(shape(rData_sendBuffer)))*8e-9, shape(rData_sendBuffer)

            write(*,'("GHOST-INIT: ALLOCATED ",A19," MEM=",f8.4," GB per rank, shape=",7(i9,1x))') &
            "rData_recvBuffer", product(real(shape(rData_recvBuffer)))*8e-9, shape(rData_recvBuffer)            
        endif

        ! allocate ijkPatches holding the indices of sender and receiver patches
        ! doesnt scale so size is not counted
        allocate( ijkPatches(1:2, 1:3, -16:56*3, 1:3), stat=status(1) )

        ! thanks to the mix of Fortran and MPI this is quite a nightmare with 0- or 1-based arrays
        allocate( send_request(1:2*Ncpu) )
        allocate( recv_request(1:2*Ncpu) )
        allocate( line_buffer( 1:Neqn*product(Bs(1:params%dim)+2*g) ) )
        allocate( int_pos(0:Ncpu-1), real_pos(0:Ncpu-1))
        allocate( data_recv_counter(0:Ncpu-1), data_send_counter(0:Ncpu-1) )
        allocate( meta_recv_counter(0:Ncpu-1), meta_send_counter(0:Ncpu-1) )
        memory_this = (10.0*real(Ncpu)+product(real(shape(line_buffer))))*8.0e-9
        memory_total = memory_total + memory_this

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
        call family_setup_patches(params, output_to_file=.true.)

        ! this routine is not performance-critical
        call MPI_barrier( WABBIT_COMM, status(1))

        if (rank==0) then
            write(*,'("INIT: Measured local (on 1 cpu) memory (only ghosts!):   ",g13.5," GB per rank")') memory_total
            write(*,'("╔", 78("═"), "╗")')
            write(*,'("║", 20(" "), A, A39)') "GHOST-INIT complete :)", "║"
            write(*,'("╚", 78("═"), "╝")')
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
    i = maxval( ijkPatches(2,1,:,RESPRE) )*2
    j = maxval( ijkPatches(2,2,:,RESPRE) )*2
    k = maxval( ijkPatches(2,3,:,RESPRE) )*2

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

    integer(kind=ik) :: N_neighbors, lvl_diff, neighborhood, i, j, k
    integer(kind=ik) :: ijk_recv(2,3), ijk_send(2,3), ijk_buff(2,3)
    logical :: debug_to_file
    character(len=80) :: debug_file_name

    ! full reset of all patch definitions. Note we set 1 not 0
    ijkPatches = 1

    N_neighbors = 56*3  ! 56 for all leveldifferences

    debug_to_file = .false.
#ifdef DEV
    debug_to_file = output_to_file
#endif
    if (debug_to_file) then
        write(debug_file_name, '(A, A5, 3(A, i0), A)') &
            "ghost_bounds.dat"
        open(16,file=debug_file_name,status='replace')
        write(16,'(3(A, i0), 2(A))') "% dim=", params%dim, ", Bs=", params%Bs(1), ", g=", params%g, ", CDF=", params%wavelet
        write(16,'(A, A13, 8(A15))') "% ", "Neighborhood" , "lvl_diff", "send/recv/buff", "idx_1", "idx_2", "idy_1", "idy_2", "idz_1", "idz_2"
    endif

    ! set all neighbour relations
    do neighborhood = 1, N_neighbors
        if (neighborhood <= 56) then
            lvl_diff = 0
        elseif (neighborhood <= 2*56) then
            lvl_diff = +1
        else
            lvl_diff = -1
        endif

        !---------larger ghost nodes (used for refinement and coarsening and maybe RHS)------------
        ! The receiver and sender bounds are hard-coded with their directions.
        call set_recv_bounds( params, ijk_recv, neighborhood, lvl_diff, gminus, gplus)
        call set_send_bounds( params, ijk_send, ijk_buff, neighborhood, lvl_diff, gminus, gplus)


        ! ! debugging borders to terminal
        ! if (params%rank == 0 .and. any(lvl_diff == (/ 0 /))) then
        !     if (sum(ijkrecv) /= 6) then
        !         write(*, '("R R", i0 ," N ", i2, " lvl_diff ", i2, " patch ", 6(i2, 1x), "  - prod ", i3)') &
        !             params%rank, neighborhood, lvl_diff, ijkrecv, product(ijkrecv(2, 1:params%dim) - ijkrecv(1, 1:params%dim))
        !     endif
        !     if (sum(ijksend) /= 6) then
        !         write(*, '("S R", i0 ," N ", i2, " lvl_diff ", i2, " patch ", 6(i2, 1x), "  - prod ", i3)') &
        !             params%rank, neighborhood, lvl_diff, ijksend, product(ijksend(2, 1:params%dim) - ijksend(1, 1:params%dim))
        !     endif
        !     if (sum(ijkbuff) /= 6) then
        !         write(*, '("B R", i0 ," N ", i2, " lvl_diff ", i2, " patch ", 6(i2, 1x), "  - prod ", i3)') &
        !             params%rank, neighborhood, lvl_diff, ijkbuff, product(ijkbuff(2, 1:params%dim) - ijkbuff(1, 1:params%dim))
        !     endif
        ! endif

        if (debug_to_file) then
            ! debugging borders to file
            ! 1: neighborhood, 2: lvl_diff, 3: send/receive/respre, 4-9: Indices
            write(16,'(9(i15))') neighborhood, lvl_diff, SENDER, ijk_send(:, :)
            write(16,'(9(i15))') neighborhood, lvl_diff, RECVER, ijk_recv(:, :)
            write(16,'(9(i15))') neighborhood, lvl_diff, RESPRE, ijk_buff(:, :)
        endif

        ijkPatches(1:2, 1:3, neighborhood, RECVER) = ijk_recv
        ijkPatches(1:2, 1:3, neighborhood, SENDER) = ijk_send
        ijkPatches(1:2, 1:3, neighborhood, RESPRE) = ijk_buff
    enddo

    if (debug_to_file) then
        close(16)
    endif
end subroutine



! setup the family patches (i.e. the indices of the sender/receiver parts
! for any neighborhood relation), this depends on the filter being used
subroutine family_setup_patches(params, output_to_file)
    implicit none
    type (type_params), intent(in) :: params
    logical, intent(in) :: output_to_file

    integer(kind=ik) :: N_family, lvl_diff, family, i, j, k, g, i_f
    integer(kind=ik) :: ijk_buff(2,3)
    logical :: debug_to_file
    character(len=80) :: debug_file_name

    ! full reset of all patch definitions. Note we set 1 not 0
    ijkPatches = 1

    N_family = 8
    if (params%dim==2) then
        N_family=4
    endif
    g = params%g

    debug_to_file = .false.
#ifdef DEV
        debug_to_file = output_to_file
#endif
    if (debug_to_file) then
        write(debug_file_name, '(A, A5, 3(A, i0), A)') &
            "family_bounds.dat"
        open(16,file=debug_file_name,status='replace')
        write(16,'(3(A, i0))') "% dim=", params%dim, ", Bs=", params%Bs(1), ", g=", params%g
        write(16,'(A, A13, 8(A15))') "% ", "Relation" , "lvl_diff", "idx_1", "idx_2", "idy_1", "idy_2", "idz_1", "idz_2"
    endif

    ! set full block relation
    ! The receiver bounds are hard-coded with a huge amount of tedious indices, they are identical to senders
    call set_recv_bounds( params, ijk_buff, 0, 0, g, g)
    ijkPatches(1:2,1:3, 0, RECVER) = ijk_buff
    ijkPatches(1:2,1:3, 0, SENDER) = ijk_buff

    if (debug_to_file) then
        ! debugging borders to file
        ! 1: neighborhood, 2: lvl_diff, 3: send/receive/respre, 4-9: Indices
        write(16,'(9(i15))') 0, 0, ijk_buff(:, :)
    endif

    ! set mother/daughter relations
    ! sender / receiver relations between the level differences are inverted
    ! lvl_diff = -1:
    ! Daughter / Finer sender or receiver: affects sc in mallat-ordering, actually all are the same
    ! lvl_diff = +1:
    ! Mother / Coarser sender or receiver: affects specific part of block
    do lvl_diff = +1, -1, -2
        do family = -1, -N_family, -1
            ! for lvl_diff = -1 we have the indices -9 until -16
            i_f = family
            if (lvl_diff == -1) i_f = family -8

            call set_recv_bounds( params, ijk_buff, i_f, lvl_diff, g, g)
            ijkPatches(1:2,1:3, i_f, SENDER) = ijk_buff
            ijkPatches(1:2,1:3, i_f, RECVER) = ijk_buff

            ! if (sum(ijk_buff) /= 6) then
            !     write(*, '("B R", i0 ," N ", i2, " lvl_diff ", i2, " patch ", 6(i2, 1x), "  - prod ", i3)') &
            !         params%rank, i_f, lvl_diff, ijk_buff, product(ijk_buff(2, 1:params%dim) - ijk_buff(1, 1:params%dim))
            ! endif

            if (debug_to_file) then
                ! debugging borders to file
                ! 1: neighborhood, 2: lvl_diff, 3: send/receive/respre, 4-9: Indices
                write(16,'(9(i15))') i_f, lvl_diff, ijk_buff(:, :)
            endif
        enddo
    enddo

    if (debug_to_file) then
        close(16)
    endif

end subroutine

end module module_MPI
