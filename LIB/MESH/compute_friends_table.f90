!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name compute_friends_table.f90
!> \version 0.4
!> \author engels, msr
!
!> \brief compute friends table. 
!> \details This is a (mpisize) x (mpisize) array of integers, and it counts
!! how many neighbor relations the mpiranks have to each other mpirank. the matrix is
!! (most of the time) symmetric, exceptions from symmetry can occur in the case of finer/coarser
!! relations.
!! In the friends array, we thus just know who's my most important neighbor. The array can look like this: (mpisize=4)
!! |     |     |     |    |
!! |-----|-----|-----|----|
!! | -1  |  3  |  2  |  4 |
!! |  3  | -1  |  22 |  12|
!! |  2  |  22 |  -1 |  27|
!! |  4  | 12  |  27 |  -1|
!! Note we set (-1) on the diagonal (since an mpirank is likely its own best friend, they're egocentric in that regard)
!! \n
!! input:    - params, neighbor lists \n
!! output:   - friends table \n
!
!> = log ======================================================================================
!! \n
!! 28/11/16 - create
!
! ********************************************************************************************

subroutine compute_friends_table(params, hvy_neighbor, friends, hvy_active, hvy_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> friends table
    integer(kind=ik),intent(inout)      :: friends(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> friends array buffer
    integer(kind=ik), allocatable       :: friends_loc(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, l, proc_id

    ! MPI variables
    integer(kind=ik)                    :: number_procs, rank, ierr

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! We need the buffer for the friends array:
    allocate( friends_loc( 1:number_procs, 1:number_procs ))
    friends_loc = 0
    friends = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active heavy data
    do k = 1, hvy_n
        ! loop over all  possible neighbor directions
        do l = 1, 16

            ! neighbor exists
            if ( hvy_neighbor( hvy_active(k), l) /= -1 ) then

                ! proc rank from ligth id
                call lgt_id_to_proc_rank( proc_id, hvy_neighbor( hvy_active(k), l), params%number_blocks )

                ! note each rank fills only one line in this array, the others are all zero
                friends_loc(rank+1, proc_id+1) = friends_loc(rank+1, proc_id+1)+1

            end if

        end do
    end do

    ! overwrite my own.. (rank is zero-based)
    friends_loc(rank+1,rank+1) = -1

    ! the friends must be known on all mpiranks
    call MPI_Allreduce(friends_loc, friends, number_procs**2, MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD, ierr)

    deallocate(friends_loc)

end subroutine compute_friends_table
