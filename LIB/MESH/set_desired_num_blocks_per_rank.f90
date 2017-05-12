!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name set_desired_num_blocks_per_rank.f90
!> \version 0.4
!> \author engels, msr
!
!> \brief create block distribution lists
!
!> 
!! input:    - params, light data, lists of active blocks \n
!! output:   - distribution arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 28/11/16 - create
!
! ********************************************************************************************

subroutine set_desired_num_blocks_per_rank(params, dist_list, opt_dist_list, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> block distribution lists
    integer(kind=ik), intent(out)       :: dist_list(:), opt_dist_list(:)

    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! loop variables
    integer                             :: k, num_blocks, proc_id, avg_blocks, number_procs, rank, excess_blocks, ierr

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    dist_list = 0
    opt_dist_list = 0

!---------------------------------------------------------------------------------------------
! main body

    ! count number of active blocks and current block distribution
    do k = 1, lgt_n
        ! count block for proc corresponing to light data
        ! the array dist_list has the length of the number of mpi procs and holds the number
        ! of blocks used on this process. From the light data ID "k", we can compute the mpirank
        ! of the process holding the block "proc_id"
        call lgt_id_to_proc_rank( proc_id, lgt_active(k), params%number_blocks )
        dist_list( proc_id + 1 ) = dist_list( proc_id + 1 ) + 1
    end do

    ! count global number of blocks on all mpiranks
    num_blocks = sum(dist_list)

    ! optimal distribution of blocks per mpirank. The simple division of "num_blocks" by "number_procs" actually
    ! yields a double (since it is not guaranteed that all mpiranks hold the exact same number of blocks)
    ! using the integer division, decimal places are cut
    avg_blocks = num_blocks / number_procs
    opt_dist_list(:) = avg_blocks

    ! some blocks are missing due to the integer division
    excess_blocks = num_blocks - sum(opt_dist_list)

    ! distribute remaining blocks (the excess blocks, if we have some)
    do while ( excess_blocks>0 )

        ! first we try to be clever and increase the counter of "desired" blocks for
        ! procs that already have more blocks than they should (by one)
        do proc_id = 1, number_procs
            ! check if this proc_id has more blocks than it is supposed to and if so, we attribute it one of the excess blocks
            if ( dist_list(proc_id) > avg_blocks) then
                opt_dist_list(proc_id) = opt_dist_list(proc_id) + 1
                ! we got rid of one excess block
                excess_blocks = excess_blocks - 1
                ! no more blocks to distribute?
                if (excess_blocks==0) exit
            end if
        end do

        ! no more blocks to distribute?
        if (excess_blocks==0) exit

        ! second, it may be that this is not enough: there are still bocks to be
        ! distributed. so now we repeat the loop, but look for mpiranks that have
        ! enough blocks and give them one more.
        do proc_id = 1, number_procs
            if ( dist_list(proc_id) == avg_blocks) then
                opt_dist_list(proc_id) = opt_dist_list(proc_id) + 1
                excess_blocks = excess_blocks - 1
                if (excess_blocks==0) exit
            end if
        end do

        ! no more blocks to distribute?
        if (excess_blocks==0) exit

        ! third, it may still not be enough...so just pick some
        do proc_id = 1, number_procs
            if ( dist_list(proc_id) < avg_blocks) then
                opt_dist_list(proc_id) = opt_dist_list(proc_id) + 1
                excess_blocks = excess_blocks - 1
                if (excess_blocks==0) exit
            end if
        end do

    end do ! end of excess block distribution

    if (rank==0) then
        ! error checking. the sum of newly distributed blocks must of course be
        ! the same as the number we had before distribution
        if (sum(opt_dist_list)/=num_blocks .or. maxval(abs(opt_dist_list-avg_blocks))>1) then
            write(*,*) "something went wrong - during balancing, we lost or gained some blocks", excess_blocks
            write(*,*) "or we have more than +-1 block difference among them"
            write(*,*) opt_dist_list
            stop
        end if
    end if

end subroutine set_desired_num_blocks_per_rank
