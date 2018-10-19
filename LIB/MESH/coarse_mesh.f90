! ********************************************************************************************
!> \brief Apply mesh coarsening: Merge tagged blocks into new, coarser blocks
!
!> \details
!> \author engels, msr, Pkrah
!! \date 17/08/18   - non blocking mpi communication for sending and receiving blocks during
!!                    gather (PKrah, commit d48299f4231040f619b2f2af5f56bf4f72994ff5  )
!! \date 22/05/2017 - Rewrite using modular subroutines, works for 2/3d cases
!! \date 08/11/16 - switch to v0.4, split old interpolate_mesh subroutine into two refine/coarsen
!!            subroutines
! ********************************************************************************************
subroutine coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_active, hvy_n, hvy_tmp )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)

    ! loop variables
    integer(kind=ik)                    :: k, maxtl, N, j
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: light_ids(1:8), mpirank_owners(1:8)
    integer(kind=ik), allocatable, save :: xfer_list(:,:)
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: data_rank, n_xfer, ierr

!---------------------------------------------------------------------------------------------
! variables initialization

    maxtl = params%max_treelevel
    ! number of blocks to merge, 4 or 8
    N = 2**params%dim
    ! at worst every block is on a different rank
    if (.not. allocated(xfer_list)) allocate(xfer_list(size(lgt_block,1),3))

    ! transfer counter
    n_xfer = 0

!---------------------------------------------------------------------------------------------
! main body

    !---------------------------------------------------------------------------
    ! first, prepare for xfer (gather information: which blocks are sent where)
    !---------------------------------------------------------------------------
    do k = 1, lgt_n
        ! Check if the block will be coarsened.
        !
        ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
        ! the active list. NOTE HERE: due to previous iterations, some light data id are already
        ! marked for xfer, (and thus given the -7 status) but they are still in active block list
        ! so: check again if this block is REALLY active (the active list is updated later)
        !
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
        ! ensure_gradedness removes the -1 flag if not all sister blocks share it

        if ( lgt_block(lgt_active(k), 1) >= 0 .and. lgt_block(lgt_active(k), maxtl+IDX_REFINE_STS) == -1) then
            ! find all sisters (including the block in question, so four or eight blocks)
            ! their light IDs are in "light_ids" and ordered by their last treecode-digit
            call find_sisters( params, lgt_active(k), light_ids(1:N), lgt_block, lgt_n, lgt_sortednumlist )

            ! figure out on which rank the sisters lie, xfer them if necessary
            do j = 1, N
                call lgt_id_to_proc_rank( mpirank_owners(j), light_ids(j), params%number_blocks )
            enddo

            ! The merging will be done on the mpirank which holds the most of the sister blocks
            data_rank = most_common_element( mpirank_owners(1:N) )

            do j = 1, N
                if (mpirank_owners(j) /= data_rank) then
                    ! MPI xfer required. Add the xfer to the list
                    n_xfer = n_xfer + 1
                    xfer_list(n_xfer, 1) = mpirank_owners(j)  ! send from this rank ..
                    xfer_list(n_xfer, 2) = data_rank          ! ... to this rank
                    xfer_list(n_xfer, 3) = light_ids(j)       ! transfer this block
                endif

                ! don't forget: mark all 4/8 sisters as treated here, in order not to trigger this
                ! loop again: we use the temporary status -7
                lgt_block(light_ids(j), maxtl+IDX_REFINE_STS) = -7
            enddo
        endif
    enddo

    ! actual xfer
    call block_xfer( params, xfer_list, n_xfer, lgt_block, hvy_block, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_tmp )

    ! the active lists are outdates after the transfer: we need to create
    ! them or find_sisters will not be able to do its job
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, &
    hvy_n, lgt_sortednumlist, .true. )

    ! actual merging
    do k = 1, lgt_n
        ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
        ! the active list. NOTE HERE: due to previous loops, some light data id are already
        ! coarsened, (and thus given the -1) but they are still in active block list
        ! so: check again if this block is REALLY active (the active list is updated later)
        !
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
        ! ensure_gradedness removes the -1 flag if not all sister blocks share it
        if ( lgt_block(lgt_active(k), 1) >= 0 .and. lgt_block(lgt_active(k), maxtl+IDX_REFINE_STS) == -7) then
            ! merge the four blocks into one new block. Merging is done in two steps,
            ! first for light data (which all CPUS do redundantly, so light data is kept synched)
            ! Then only the responsible rank will perform the heavy data merging.
            call find_sisters( params, lgt_active(k), light_ids(1:N), lgt_block, lgt_n, lgt_sortednumlist )
            call merge_blocks( params, hvy_block, lgt_block, light_ids(1:N) )
            ! note the newly merged block has status 0

            mesh_has_changed = .true.
        endif
    enddo


end subroutine coarse_mesh


logical function mesh_adapted()
    implicit none

    mesh_adapted=mesh_has_changed
end function
