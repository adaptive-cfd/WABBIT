!> \brief debug lgt_block data
!!  proc 0 gather all data and compare the data to his own light data \n
!!
!!  is currently unused, very helpful routine for development of block structures, lgt/hvy data
!!
!! input:    - params, light data \n
!! output:   - status of lgt_block synchronzation \n
! ********************************************************************************************

subroutine check_lgt_block_synchronization( params )
    implicit none

    type (type_params), intent(in)      :: params
    ! local light data array
    integer(kind=ik)                    :: my_lgt_block( size(lgt_block,1) , size(lgt_block,2)), lgt_block_0( size(lgt_block,1) , size(lgt_block,2))
    integer(kind=ik)                    :: ierr                                 ! MPI error variable
    integer(kind=ik)                    :: k, l, lgt_start, a                   ! loop variables
    integer(kind=ik), allocatable, save :: lgt_all(:,:,:), lgt_all2(:,:,:)      ! lgt data


    if (.not. allocated(lgt_all)) Then
        allocate(lgt_all(1:params%number_procs, size(lgt_block,1), size(lgt_block,2)))
        allocate(lgt_all2(1:params%number_procs, size(lgt_block,1), size(lgt_block,2)))
    endif


    lgt_all = 0
    lgt_all(params%rank+1,:,:) = lgt_block

    ! gather light data
    call MPI_Allreduce(lgt_all, lgt_all2, params%number_procs*size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_SUM, WABBIT_COMM, ierr)

    ! loop over all blocks
    do k = 1, size(lgt_block, 1)
        do l = 1, size(lgt_block, 2)

            ! compare light data on all ranks
            if ( maxval( abs(lgt_all2(:,k,l) - lgt_all2(1,k,l) ) ) > 0 .and. lgt_all2(1,k,1)/=-1) then
                ! error case
                write(*,'(80("!"))')
                write(*, '("ERROR: lgt_block is not synchron, lgt_block id: ",i5," position: ", i5)') k, l
                write(*,'(80("!"))')
                do a =1, params%number_procs
                    write(*,'("rank=",i3,"  :  ",40(i2,1x))') a, lgt_all2(a,k,:)
                enddo
                write(*,'(80("!"))')
                call abort(12399)
            end if

        end do
    end do

end subroutine check_lgt_block_synchronization



subroutine write_lgt_data(params, file)
    implicit none
    type (type_params), intent(in)      :: params
    character(len=*), intent(in) :: file
    integer(kind=ik) :: k

    if (params%rank /= 0) return

    open(14, file=file, status='replace')
    do k = 1, size(lgt_block,1)
        write(14,'(2000(i2,1x))') lgt_block(k,:)
    enddo
    close(14)

end subroutine



subroutine read_lgt_data(params, file)
    implicit none
    type (type_params), intent(in)      :: params
    character(len=*), intent(in) :: file
    integer(kind=ik) :: num_lines, num_cols

    call read_intarray_from_ascii_file_mpi(file, lgt_block, 0)
end subroutine



! This function works, but better use write_neighborhood_info(hvy_neighbor, dim) and loop over all heavy data as it is more descriptive
subroutine write_neighbors(params, file, tree_ID)
    implicit none
    type (type_params), intent(in)      :: params
    integer(kind=ik) , allocatable, save :: tmp(:,:), tmp2(:,:)
    integer(kind=ik), intent(in) :: tree_ID
    character(len=*), intent(in) :: file
    integer(kind=ik) :: k, lgt_start, rank, lgt_id, N

    N = size(hvy_neighbor,2)

    if (.not. allocated(tmp))  allocate(tmp(1:params%number_procs*params%number_blocks,N))
    if (.not. allocated(tmp2))  allocate(tmp2(1:params%number_procs*params%number_blocks,N))

    tmp = 0

    do k = 1, hvy_n(tree_ID)
        call hvy2lgt( lgt_id, hvy_active(k,tree_ID), params%rank, params%number_blocks )
        tmp(lgt_id,:) = hvy_neighbor(hvy_active(k,tree_ID),:)
    enddo

    call MPI_Allreduce(tmp, tmp2, size(tmp), MPI_INTEGER4, MPI_SUM, WABBIT_COMM, k)

    if (params%rank /= 0) return

    open(14, file=file, status='replace')
    do k = 1, size(tmp2,1)
        write(14,'(178(i6,1x))') tmp2(k,:)
    enddo
    close(14)
end subroutine
