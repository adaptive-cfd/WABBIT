!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name check_lgt_block_synchronization.f90
!> \version 0.4
!> \author msr
!
!> \brief debug lgt_block data
!
!>
!!  proc 0 gather all data and compare the data to his own light data \n
!!
!!  is currently unused, very helpful routine for development of block structures, lgt/hvy data
!!
!! input:    - params, light data \n
!! output:   - status of lgt_block synchronzation \n
!!
!!
!! = log ======================================================================================
!! \n
!! 29/11/16 - create
! ********************************************************************************************

subroutine check_lgt_block_synchronization( params, lgt_block)
    !---------------------------------------------------------------------------------------------
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)     :: lgt_block(:, :)

    ! local light data array
    integer(kind=ik)                    :: my_lgt_block( size(lgt_block,1) , size(lgt_block,2)), lgt_block_0( size(lgt_block,1) , size(lgt_block,2))

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! loop variables
    integer(kind=ik)                    :: k, l, lgt_start, a
    ! lgt data
    integer(kind=ik), allocatable, save :: lgt_all(:,:,:), lgt_all2(:,:,:)

    !---------------------------------------------------------------------------------------------

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



subroutine write_lgt_data(params, lgt_block, file)
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)     :: lgt_block(:, :)
    character(len=*), intent(in) :: file
    integer(kind=ik) :: k

    if (params%rank /= 0) return

    open(14, file=file, status='replace')
    do k = 1, size(lgt_block,1)
        write(14,'(2000(i2,1x))') lgt_block(k,:)
    enddo
    close(14)

end subroutine



subroutine read_lgt_data(params, lgt_block, file)
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    character(len=*), intent(in) :: file
    integer(kind=ik) :: num_lines, num_cols

    call read_intarray_from_ascii_file_mpi(file, lgt_block, 0)
end subroutine



subroutine write_neighbors(params, hvy_active, hvy_n, hvy_neighbor, file)
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    integer(kind=ik) , allocatable, save :: tmp(:,:), tmp2(:,:)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)     :: hvy_n
    character(len=*), intent(in) :: file
    integer(kind=ik) :: k, lgt_start, rank, lgt_id

    if (.not. allocated(tmp))  allocate(tmp(1:params%number_procs*params%number_blocks,16))
    if (.not. allocated(tmp2))  allocate(tmp2(1:params%number_procs*params%number_blocks,16))

    tmp = 0

    do k = 1, hvy_n
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
        tmp(lgt_id,:) = hvy_neighbor(hvy_active(k),:)
    enddo

    call MPI_Allreduce(tmp, tmp2, size(tmp), MPI_INTEGER4, MPI_SUM, WABBIT_COMM, k)

    if (params%rank /= 0) return

    open(14, file=file, status='replace')
    do k = 1, size(tmp2,1)
        write(14,'(16(i6,1x))') tmp2(k,:)
    enddo
    close(14)
end subroutine
