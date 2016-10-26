! ********************************
! WABBIT
! --------------------------------
!
! broadcast all light block data
!
! name: broadcast_light_data.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine broadcast_light_data()

    use mpi
    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                    :: k, N, rank, ierr

    N = blocks_params%number_max_blocks

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! loop over all blocks (light data)
    do k = 1, N
        ! broadcast data to all other procs
        ! integer data
        call MPI_Bcast(blocks(k)%treecode, 10, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%level, 1, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%neighbor_treecode, 8*10, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%neighbor_id, 8, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%neighbor2_treecode, 4*10, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%neighbor2_id, 4, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%neighbor_number, 8, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%refinement, 1, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%proc_rank, 1, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%proc_data_id, 1, MPI_INTEGER4, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        ! character data
        call MPI_Bcast(blocks(k)%neighbor_dir, 2*8, MPI_CHARACTER, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(blocks(k)%neighbor2_dir, 2*4, MPI_CHARACTER, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)
        ! logical data
        call MPI_Bcast(blocks(k)%active, 1, MPI_LOGICAL, blocks(k)%proc_rank, MPI_COMM_WORLD, ierr)

    end do

end subroutine broadcast_light_data
