!> \brief Get current distribution of blocks among CPUs, compute optimal #blocks for each mpirank
!>        ---------------------------------------------
!! \details
!! The function returns a list "blocks_per_rank" of size mpisize, where the number of blocks for each mpirank
!! is contained. Note 1-based indexing. \n
!! In addition we return "blocks_per_rank_optimal" of size mpisize, where we set the number of blocks for each rank
!! such that the distribution is as homogeneous as possible and required block transfer is minimized
!--------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine set_desired_num_blocks_per_rank1(params, blocks_per_rank_optimal, tree_ID) ! predictable
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    integer(kind=ik), intent(out)       :: blocks_per_rank_optimal(0:)          !> block distribution lists. Note zero-based indexing. (as in MPI library)
    integer(kind=ik), intent(in)        :: tree_ID

    ! optimal distribution of blocks per mpirank. The simple division of "num_blocks" by "number_procs" actually
    ! yields a double (since it is not guaranteed that all mpiranks hold the exact same number of blocks)
    ! using the integer division, decimal places are cut
    blocks_per_rank_optimal(:) = lgt_n(tree_ID) / params%number_procs

    ! as this does not necessarily work out, distribute remaining blocks on the first CPUs
    if (mod(lgt_n(tree_ID), params%number_procs) > 0) then
        blocks_per_rank_optimal(0:mod(lgt_n(tree_ID), params%number_procs)-1) = &
        blocks_per_rank_optimal(0:mod(lgt_n(tree_ID), params%number_procs)-1) + 1
    end if

    ! some error control -> did we loose blocks? should never happen.
#ifdef DEV
    if ( sum(blocks_per_rank_optimal) /= lgt_n(tree_ID)) then
        call abort(1028,"ERROR: while reading from file, we seem to have gained/lost some blocks during distribution...")
    end if
#endif

end subroutine


subroutine set_desired_num_blocks_per_rank2(params, blocks_per_rank, blocks_per_rank_optimal, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    integer(kind=ik), intent(out)       :: blocks_per_rank(:), blocks_per_rank_optimal(:)   !> block distribution lists. Note 1-based indexing.
    integer(kind=ik), intent(in)        :: tree_ID

    integer                             :: num_blocks, proc_id, avg_blocks, number_procs, rank, excess_blocks
    integer(kind=ik)                    :: ierr                                 ! MPI error variable

    rank                    = params%rank
    number_procs            = params%number_procs
    blocks_per_rank         = 0
    blocks_per_rank_optimal = 0

    ! save my number of active blocks
    blocks_per_rank(rank+1) = hvy_n(tree_ID)
    ! count number of active blocks and current block distribution
    call MPI_Allreduce(MPI_IN_PLACE, blocks_per_rank, number_procs, MPI_INTEGER4, MPI_SUM, WABBIT_COMM, ierr)

    ! count global number of blocks on all mpiranks
    num_blocks = lgt_n(tree_ID)

    ! optimal distribution of blocks per mpirank. The simple division of "num_blocks" by "number_procs" actually
    ! yields a double (since it is not guaranteed that all mpiranks hold the exact same number of blocks)
    ! using the integer division, decimal places are cut
    avg_blocks = num_blocks / number_procs
    blocks_per_rank_optimal(:) = avg_blocks

    ! some blocks are missing due to the integer division
    excess_blocks = num_blocks - sum(blocks_per_rank_optimal)


    ! distribute remaining blocks (the excess blocks, if we have some)
    do while ( excess_blocks>0 )

        ! first we try to be clever and increase the counter of "desired" blocks for
        ! procs that already have more blocks than they should (by one)
        do proc_id = 1, number_procs
            ! check if this proc_id already has more blocks than it is supposed to
            ! and if so, we attribute it one of the excess blocks
            if ( blocks_per_rank(proc_id) > avg_blocks) then
                blocks_per_rank_optimal(proc_id) = blocks_per_rank_optimal(proc_id) + 1
                ! we got rid of one excess block
                excess_blocks = excess_blocks - 1
                ! no more blocks to distribute?
                if (excess_blocks==0) exit
            end if
        end do

        ! no more blocks to distribute?
        if (excess_blocks==0) exit

        ! second, it may be that this is not enough: there are still bocks to be
        ! distributed. so now we repeat the loop, but look for mpiranks that already
        ! have enough blocks and give them one more.
        do proc_id = 1, number_procs
            if ( blocks_per_rank(proc_id) == avg_blocks) then
                blocks_per_rank_optimal(proc_id) = blocks_per_rank_optimal(proc_id) + 1
                excess_blocks = excess_blocks - 1
                if (excess_blocks==0) exit
            end if
        end do

        ! no more blocks to distribute?
        if (excess_blocks==0) exit

        ! third, it may still not be enough. so now all procs that currently have
        ! less than the average get one more block each
        do proc_id = 1, number_procs
            if ( blocks_per_rank(proc_id) < avg_blocks) then
                blocks_per_rank_optimal(proc_id) = blocks_per_rank_optimal(proc_id) + 1
                excess_blocks = excess_blocks - 1
                if (excess_blocks==0) exit
            end if
        end do
    end do ! end of excess block distribution

#ifdef DEV
    if (rank==0) then
        ! error checking. the sum of newly distributed blocks must of course be
        ! the same as the number we had before distribution
        if (sum(blocks_per_rank_optimal)/=num_blocks .or. maxval(abs(blocks_per_rank_optimal-avg_blocks))>1) then
            write(*,*) "something went wrong - during balancing, we lost or gained some blocks", excess_blocks
            write(*,*) "or we have more than +-1 block difference among them"
            write(*,*) blocks_per_rank_optimal
            call abort(11191,"ERROR lost some blocks")
        end if
    end if
#endif

end subroutine
