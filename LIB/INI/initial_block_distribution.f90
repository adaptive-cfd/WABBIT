! ********************************
! WABBIT
! --------------------------------
!
! distribute initial blocks on processes
!
! name: initial_block_distribution.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine initial_block_distribution(distribution)

    use mpi
    use module_params
    use module_blocks

    implicit none

    character(len=80), intent(in)   :: distribution

    integer(kind=ik)                :: light_id, heavy_id, rank, n_procs, ierr, proc_id

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_procs, ierr)

    select case(distribution)
        case("equal")

            proc_id = 0
            ! loop over all blocks
            do light_id = 1, blocks_params%number_max_blocks

                ! block is active
                if ( blocks(light_id).active == .true. ) then

                    ! find free heavy id
                    if ( rank == proc_id ) then
                        ! find free id
                        call get_heavy_free_block(heavy_id)
                        ! save light data id in heavy data on specific proc
                        blocks_data(heavy_id)%block_id = light_id
                    end if
                    ! sychronize free id
                    call MPI_Bcast(heavy_id, 1, MPI_INTEGER4, proc_id, MPI_COMM_WORLD, ierr)

                    ! save heavy-id and proc-id in light data
                    blocks(light_id).proc_data_id   = heavy_id
                    blocks(light_id).proc_rank      = proc_id
                    ! increase proc-id
                    proc_id                         = proc_id + 1

                else
                    ! block is not active, so reset id's to heavy data
                    blocks(light_id).proc_data_id   = -1
                    blocks(light_id).proc_rank      = -1

                end if
                ! reset proc counter
                if ( proc_id==n_procs ) proc_id = 0

            end do

    end select

end subroutine initial_block_distribution
