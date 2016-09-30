! ********************************
! WABBIT
! --------------------------------
!
! refine every block to create the
! wavelet safety zone
!
! name: refine_everywhere.f90
! date: 30.09.2016
! author: engels, msr
! version: 0.2
!
! ********************************

subroutine refine_everywhere()

    use module_params
    use module_blocks
    use module_interpolation

    implicit none

!    real(kind=rk), dimension(:,:), allocatable      :: u1, u2, u3, u5, u6
!    real(kind=rk)                                   :: detail1, detail2
!
    integer(kind=ik)    :: k, N
!    integer(kind=ik)                                :: dF, k, N, block_num, Bs, g, allocate_error, max_status
!    integer(kind=ik), dimension(:,:), allocatable   :: dF_status
!
    N           = blocks_params%number_max_blocks
!    Bs = blocks_params%size_block
!    g  = blocks_params%number_ghost_nodes

    ! set status "refine" for all active blocks
    do k = 1, N
        if (blocks(k)%active) then
            blocks(k)%refinement = 1
        end if
    end do

    ! check if block has reached maximal level
    call respect_min_max_treelevel()

    ! interpolate the new mesh
    call interpolate_mesh()


!    ! allocate memory for array with refinement status for all fields
!    allocate( dF_status(N, blocks_params%number_data_fields), stat=allocate_error )
!
!    dF_status = 0
!
!    ! synchronize ghostnodes
!    call synchronize_ghosts()
!
!    ! clear old refinement status for all blocks
!    do k = 1, N
!        block_num = blocks_params%active_list(k)
!        blocks(block_num)%refinement = 0
!        blocks(block_num)%data_fields(:)%detail = 0.0_rk
!    end do
!
!    ! loop over all fields
!    do dF = 1, blocks_params%number_data_fields
!        ! loop over all active blocks
!        do k = 1, N
!            block_num = blocks_params%active_list(k)
!            dF_status(k, dF) = +1
!        end do
!    end do
!
!    ! set block refinement status
!    do k = 1, N
!
!        block_num = blocks_params%active_list(k)
!
!        max_status = -99
!        ! loop over all data fields
!        do dF = 1, blocks_params%number_data_fields
!            ! block refinement status is the maximal status
!            max_status = max( max_status, dF_status(k, dF) )
!        end do
!
!        blocks(block_num)%refinement = max_status
!
!    end do
!
!    ! check if block has reached maximal level
!    call respect_min_max_treelevel()
!
!    deallocate( dF_status, stat=allocate_error )

end subroutine refine_everywhere
