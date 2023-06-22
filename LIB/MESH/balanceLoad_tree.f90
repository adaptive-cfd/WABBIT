!> \image html balancing.svg "Load balancing" width=400
! ********************************************************************************************
subroutine balanceLoad_tree( params, hvy_block, tree_ID, predictable_dist)

    implicit none

    type (type_params), intent(in)      :: params                     !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    !> if true balance the load will always give the same block distribution
    !> for the same treestructure (default is False)
    logical, intent(in),optional        :: predictable_dist
    !=============================================================================

    integer(kind=ik)  :: rank, proc_dist_id, proc_data_id, ierr, number_procs, &
                         k, N, l, com_i, com_N, heavy_id, sfc_id, neq, &
                         lgt_free_id, hvy_free_id, hilbertcode(params%Jmax)
    ! block distribution lists
    integer(kind=ik), allocatable, save :: opt_dist_list(:), dist_list(:), friends(:,:), &
                     affinity(:), sfc_com_list(:,:), sfc_sorted_list(:,:)
    real(kind=rk) :: t0, t1
    logical       :: is_predictable

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! check if argument is present or not
    is_predictable=.False.
    if (present(predictable_dist)) then
      is_predictable=predictable_dist
    endif

    if (params%number_procs == 1) then
        ! on only one proc, no balancing is required
        return
    endif

    ! start time
    t0 = MPI_Wtime()

    ! MPI_parameters
    rank = params%rank
    number_procs = params%number_procs

    ! allocate block to proc lists
    if (.not.allocated(opt_dist_list)) allocate( opt_dist_list(1:number_procs))
    if (.not.allocated(dist_list)) allocate( dist_list(1:number_procs))
    ! allocate sfc com list, maximal number of communications is when every proc wants to send all of his blocks
    ! NOTE: it is not necessary or wise to reset this array (it is large!)
    if (.not.allocated(sfc_com_list)) allocate( sfc_com_list( number_procs*params%number_blocks, 3 ) )

    ! allocate space filling curve list, number of elements is the number of active blocks
    ! and for each block, we store the space-filling-curve-index and the lgt ID
    if (.not.allocated(sfc_sorted_list)) allocate( sfc_sorted_list( size(lgt_block,1), 2) )

    ! number of blocks
    N = params%number_blocks
    neq = params%n_eqn



    !---------------------------------------------------------------------------------
    ! First step: define how many blocks each mpirank should have.
    !---------------------------------------------------------------------------------
    if (is_predictable) then
       call set_desired_num_blocks_per_rank(params, opt_dist_list, tree_ID)
    else
       call set_desired_num_blocks_per_rank(params, dist_list, opt_dist_list, tree_ID)
    endif

    ! at this point, we know how many blocks a mpirank has: "dist_list(myrank+1)"
    ! and how many it should have, if equally distributed: "opt_dist_list(myrank+1)"

    ! if (rank==0) then
    !     call append_t_file( "balancing.t", (/ real(maxval(opt_dist_list-dist_list),kind=rk),&
    !     real(minval(opt_dist_list-dist_list),kind=rk),&
    !     real(sum(abs(opt_dist_list-dist_list)),kind=rk) /))
    ! endif

    !---------------------------------------------------------------------------------
    ! 1st: calculate space filling curve index for all blocks
    !---------------------------------------------------------------------------------
    t1 = MPI_wtime()
    select case (params%block_distribution)
    case("sfc_z")
        !-----------------------------------------------------------
        ! Z - curve
        !-----------------------------------------------------------
        if (params%dim == 3) then
            do k = 1, lgt_n(tree_ID)
                call treecode_to_sfc_id_3D( sfc_id, lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), params%Jmax )
                sfc_sorted_list(k, 1:2) = (/sfc_id, lgt_active(k, tree_ID)/)
            end do
        else
            do k = 1, lgt_n(tree_ID)
                call treecode_to_sfc_id_2D( sfc_id, lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), params%Jmax )
                sfc_sorted_list(k, 1:2) = (/sfc_id, lgt_active(k, tree_ID)/)
            end do
        endif
    case("sfc_hilbert")
        !-----------------------------------------------------------
        ! Hilbert curve
        !-----------------------------------------------------------
        if (params%dim == 3) then
            do k = 1, lgt_n(tree_ID)
                ! transfer treecode to hilbertcode
                call treecode_to_hilbertcode_3D( lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), hilbertcode, params%Jmax)
                ! calculate sfc position from hilbertcode
                call treecode_to_sfc_id_3D( sfc_id, hilbertcode, params%Jmax )
                sfc_sorted_list(k, 1:2) = (/sfc_id, lgt_active(k, tree_ID)/)
            end do
        else
            do k = 1, lgt_n(tree_ID)
                ! transfer treecode to hilbertcode
                call treecode_to_hilbertcode_2D( lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), hilbertcode, params%Jmax)
                ! calculate sfc position from hilbertcode
                call treecode_to_sfc_id_2D( sfc_id, hilbertcode, params%Jmax )
                sfc_sorted_list(k, 1:2) = (/sfc_id, lgt_active(k, tree_ID)/)
            end do
        endif

    case default
        call abort(210309, "The block_dist is unkown"//params%block_distribution)

    end select


    ! sort sfc_list according to the first dimension, thus the position on
    ! the space filling curve (this was a bug, fixed: Thomas, 13/03/2018)
    if (lgt_n(tree_ID) > 1) then
        call quicksort_ik(sfc_sorted_list, 1, lgt_n(tree_ID), 1, 2)
    end if
    call toc( "balanceLoad_tree (SFC+sort)", MPI_wtime()-t1 )

    !---------------------------------------------------------------------------------
    ! 2nd: plan communication (fill list of blocks to transfer)
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
        call lgt2proc( proc_data_id, sfc_sorted_list(k,2), params%number_blocks )

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
            sfc_com_list(com_i, 1) = proc_data_id           ! sender mpirank
            sfc_com_list(com_i, 2) = proc_dist_id           ! receiver mpirank
            sfc_com_list(com_i, 3) = sfc_sorted_list(k,2)   ! light id of block
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
    ! 3rd: actual communication (send/recv)
    !---------------------------------------------------------------------------------
    call block_xfer( params, sfc_com_list, com_i, hvy_block, msg="balanceLoad_tree" )
    call toc( "balanceLoad_tree (comm)", MPI_wtime()-t1 )

    ! the block xfer changes the light data, and afterwards active lists are outdated.
    ! NOTE: an idea would be to also xfer the neighboring information (to save the updateNeighbors_tree
    ! call) but that is tricky: the neighbor list contains light ID of the neighbors, and those
    ! also change with the xfer.
    call updateMetadata_tree(params, tree_ID)

    ! timing
    call toc( "balanceLoad_tree (TOTAL)", MPI_wtime()-t0 )
end subroutine balanceLoad_tree
