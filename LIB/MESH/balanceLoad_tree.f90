!> \image html balancing.svg "Load balancing" width=400
! ********************************************************************************************
subroutine balanceLoad_tree( params, hvy_block, tree_ID)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)      :: params                     !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    !=============================================================================

    integer(kind=ik)  :: rank, proc_dist_id, proc_data_id, ierr, number_procs, &
                         k, N, l, com_i, com_N, heavy_id, sfc_id, neq, &
                         lgt_free_id, hvy_free_id, hilbertcode(params%Jmax)
    integer(kind=tsize) :: treecode, hilbertcode2
    ! block distribution lists
    integer(kind=ik), allocatable, save :: opt_dist_list(:), sfc_com_list(:,:), sfc_sorted_list(:,:)
    real(kind=rk) :: t0, t1

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! on only one proc, no balancing is required
    if (params%number_procs == 1) return

    ! start time
    t0 = MPI_Wtime()

    ! MPI_parameters
    rank = params%rank
    number_procs = params%number_procs

    ! allocate block to proc lists
    if (.not.allocated(opt_dist_list)) allocate( opt_dist_list(1:number_procs))
    ! allocate sfc com list, maximal number of communications is when every proc wants to send all of his blocks
    ! NOTE: it is not necessary or wise to reset this array (it is large!)
    if (.not.allocated(sfc_com_list)) allocate( sfc_com_list( 3, number_procs*params%number_blocks ) )

    ! allocate space filling curve list, number of elements is the number of active blocks
    ! and for each block, we store the space-filling-curve-index and the lgt ID
    if (.not.allocated(sfc_sorted_list)) allocate( sfc_sorted_list( 3, size(lgt_block,1)) )

    ! number of blocks
    N = params%number_blocks
    neq = params%n_eqn



    !---------------------------------------------------------------------------------
    ! 1st: define how many blocks each mpirank should have.
    !---------------------------------------------------------------------------------

    ! optimal distribution of blocks per mpirank. The simple division of "num_blocks" by "number_procs" actually
    ! yields a double (since it is not guaranteed that all mpiranks hold the exact same number of blocks)
    ! using the integer division, decimal places are cut
    opt_dist_list(:) = lgt_n(tree_ID) / number_procs

    ! as this does not necessarily work out, distribute remaining blocks on the first CPUs
    opt_dist_list(1:mod(lgt_n(tree_ID), number_procs)) = opt_dist_list(1:mod(lgt_n(tree_ID), number_procs)) + 1

    ! some error control -> did we loose blocks? should never happen.
#ifdef DEV
    if ( sum(opt_dist_list) /= lgt_n(tree_ID)) then
        call abort(1028,"ERROR: while reading from file, we seem to have gained/lost some blocks during distribution...")
    end if
#endif

    !---------------------------------------------------------------------------------
    ! 2nd: calculate space filling curve index for all blocks
    !---------------------------------------------------------------------------------
    t1 = MPI_wtime()
    select case (params%block_distribution)
    case("sfc_z")
        !-----------------------------------------------------------
        ! Z - curve - position is simply treecode
        !-----------------------------------------------------------
        do k = 1, lgt_n(tree_ID)
            sfc_sorted_list(:, k) = -1
            treecode = get_tc(lgt_block(lgt_active(k, tree_ID), IDX_TC_1 : IDX_TC_2))
            sfc_sorted_list(1, k) = lgt_active(k, tree_ID)
            call set_tc(sfc_sorted_list(2:3, k), treecode)
        end do
    case("sfc_hilbert")
        !-----------------------------------------------------------
        ! Hilbert curve - we need to translate the treecode
        !-----------------------------------------------------------
        do k = 1, lgt_n(tree_ID)
            sfc_sorted_list(:, k) = -1

            ! transfer treecode to hilbertcode
            treecode = get_tc(lgt_block(lgt_active(k, tree_ID), IDX_TC_1 : IDX_TC_2))
            if (params%dim == 3) then
                call treecode_to_hilbertcode_3D( treecode, hilbertcode2, dim=params%dim, &
                    level=lgt_block(lgt_active(k, tree_ID), IDX_MESH_LVL), max_level=params%Jmax)
            else
                call treecode_to_hilbertcode_2D( treecode, hilbertcode2, dim=params%dim, &
                    level=lgt_block(lgt_active(k, tree_ID), IDX_MESH_LVL), max_level=params%Jmax)
            endif
 
            sfc_sorted_list(1, k) = lgt_active(k, tree_ID)
            call set_tc(sfc_sorted_list(2:3, k), hilbertcode2)
        end do

    case default
        call abort(210309, "The block_dist is unkown: "//params%block_distribution)
    end select

    ! sort sfc_list according to the first dimension, thus the position on
    ! the space filling curve (this was a bug, fixed: Thomas, 13/03/2018)
    if (lgt_n(tree_ID) > 1) then
        call quicksort(sfc_sorted_list, 1, lgt_n(tree_ID), 3)
    end if
    call toc( "balanceLoad_tree (SFC+sort)", 91, MPI_wtime()-t1 )

    !---------------------------------------------------------------------------------
    ! 3rd: plan communication (fill list of blocks to transfer)
    !---------------------------------------------------------------------------------
    t1 = MPI_wtime()
    ! proc_dist_id: process responsible for current part of sfc
    ! proc_data_id: process who stores data of sfc element

    ! we start the loop on the root rank (0), then assign the first elements
    ! of the SFC, then to second rank, etc. (thus: proc_dist_id is a loop variable)
    proc_dist_id = 0

    ! communication counter. each communication (=send and receive) is stored
    ! in a long list
    com_i = 0

    ! prepare lists for transfering of blocks
    ! COLLECTIVE OPERATION
    do k = 1, lgt_n(tree_ID)
        ! if the current owner of the SFC is supposed to have zero blocks
        ! then it does not really own this part of the SFC. So we look for the
        ! first rank which is supposed to hold at least one block, and declare it as owner
        ! of this part. NOTE: as we try to minimize communication during send/recv in
        ! load balancing, it may well be that the list of active mpiranks (ie those
        ! which have nonzero number of blocks) is non contiguous, i.e.
        ! opt_dist_list = 1 1 1 0 0 0 0 1 0 1
        ! can happen.
        do while ( opt_dist_list(proc_dist_id+1) == 0 )
            proc_dist_id = proc_dist_id + 1
        end do

        ! find out on which mpirank lies the block that we're looking at
        call lgt2proc( proc_data_id, sfc_sorted_list(1, k), params%number_blocks )

        ! does this block lie on the right mpirank, i.e., the current part of the
        ! SFC? if so, nothing needs to be done. otherwise, the following if is active
        if ( proc_dist_id /= proc_data_id ) then
            ! as this block is one the wrong rank, it will be sent away from its
            ! current owner (proc_data_id) to the owner of this part of the
            ! SFC (proc_dist_id)

            ! save this send+receive operation in the list of planned communications
            ! column
            !    1     sender mpirank
            !    2     receiver mpirank
            !    3     block light data id
            com_i = com_i + 1
            sfc_com_list(1, com_i) = proc_data_id           ! sender mpirank
            sfc_com_list(2, com_i) = proc_dist_id           ! receiver mpirank
            sfc_com_list(3, com_i) = sfc_sorted_list(1,k)   ! light id of block
        end if

        ! The opt_dist_list defines how many blocks this rank should have, and
        ! we just treated one (which either already was on the mpirank or will be on
        ! it after communication), so remove one item from the opt_dist_list
        opt_dist_list( proc_dist_id+1 ) = opt_dist_list( proc_dist_id+1 ) - 1

        ! if there is no more blocks to be checked, increase mpirank counter by one
        if ( opt_dist_list( proc_dist_id+1 ) == 0 ) then
            proc_dist_id = proc_dist_id + 1
        end if
    end do

    !---------------------------------------------------------------------------------
    ! 4th: actual communication (send/recv)
    !---------------------------------------------------------------------------------
    call block_xfer( params, sfc_com_list(1:3, 1:com_i), com_i, hvy_block, msg="balanceLoad_tree" )

    call toc( "balanceLoad_tree (comm)", 92, MPI_wtime()-t1 )

    ! the block xfer changes the light data, and afterwards active lists are outdated.
    call updateMetadata_tree(params, tree_ID)

    ! timing
    call toc( "balanceLoad_tree (TOTAL)", 90, MPI_wtime()-t0 )
end subroutine balanceLoad_tree
