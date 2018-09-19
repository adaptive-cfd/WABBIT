! ********************************************************************************************
!> \file
!> \brief Apply mesh coarsening: Merge tagged blocks into new, coarser blocks
!
!> \details
!> \author engels, msr, Pkrah
!! \date 17/08/18   - non blocking mpi communication for sending and receiving blocks during
!!                    gather (PKrah, commit d48299f4231040f619b2f2af5f56bf4f72994ff5  )
!! \date 22/05/2017 - Rewrite using modular subroutines, works for 2/3d cases
!! \date 08/11/16 - switch to v0.4, split old interpolate_mesh subroutine into two refine/coarsen
!!            subroutines
!! \todo it would be faster if the merging would allready start, after all sisters are recieved
!! on the gather rank
! ********************************************************************************************



subroutine coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )

    implicit none
    !--------------------------------------------------------------
    type (type_params), intent(in)      :: params                  !< user defined parameter structure
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)         !< light data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)!< heavy data array - block data
    integer(kind=ik), intent(inout)     :: lgt_active(:)           !< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n                   !< number of active blocks (light data)
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)  !< sorted list of numerical treecodes, used for block finding
    !--------------------------------------------------------------
    ! loop variables and bounds
    integer(kind=ik)                    :: k, maxtl, N_sisters, j
    ! rank of all sister blocks
    integer(kind=ik),save,allocatable   :: rank_sisters(:)
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: light_ids(1:8,lgt_n)
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: gather_rank,  n_req, ierr,n_merge
    ! non blocking send receives use communication requests
    integer(kind=ik),save,allocatable   :: request(:)

    mesh_has_changed=.false.
    maxtL = params%max_treelevel

    if (params%threeD_case) then
      N_sisters = 8
    else
      N_sisters = 4
    endif

    ! at worst every block is on a different rank
    if (.not. allocated(request))       allocate(request(N_sisters*size(lgt_block,1)))
    if (.not. allocated(rank_sisters))  allocate(rank_sisters(N_sisters))

    n_req   = 0_ik
    n_merge = 0_ik
    ! loop over all active light data
    do k = 1, lgt_n
      ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
      ! the active list. NOTE HERE: due to previous loops, some light data id are already
      ! coarsened, (and thus given the -1) but they are still in active block list
      ! so: check again if this block is REALLY active (the active list is updated later)
      !
      ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
      ! ensure_gradedness removes the -1 flag if not all sister blocks share it (NOTE: Before 06/2018, we
      ! used a second status: -2 to mark blocks that can *really* be coarsened)
      if ( lgt_block(lgt_active(k), 1) >= 0 .and. lgt_block(lgt_active(k), maxtl+idx_refine_sts) == -1) then
          ! count the number of merges
          n_merge=n_merge+1

          ! find all sisters (including the block in question, so four blocks)
          ! their light IDs are in "light_ids" and ordered by their last treecode-digit
          call find_sisters( params, lgt_active(k), light_ids(1:N_sisters,n_merge), lgt_block, lgt_n, lgt_sortednumlist )

          do j=1, N_sisters
        	   !  Check which CPU holds this block. The CPU will also hold the merged, new block
             call lgt_id_to_proc_rank( rank_sisters(j), lgt_active(k), params%number_blocks)
          end do
          ! find the rank which holds the most blocks to define the gather rank
          gather_rank = most_common_element(rank_sisters)
          ! gather all four sisters on the process "gather_rank". The light_ids are updated in the routine
          ! and they are still in the same order (0,1,2,3)-sister. It is just that they are now on one CPU
          call gather_blocks_on_proc( params, hvy_block, lgt_block, gather_rank, light_ids(1:N_sisters,n_merge), request, n_req )
        end if
      end do

      ! In gather_block_on_proc we have initiated communications with irecv and isend, to get the
      ! block information for sisters on different ranks.
      ! Before the sisters can be merged we have to make sure, that sending and receiving is done!
      if (n_req > 0_ik ) call MPI_Waitall( n_req, request, MPI_STATUSES_IGNORE, ierr)

      ! After we have gathered all n_merge*N sister blocks, we merge them into n_merge parent blocks
      do k = 1, n_merge
          ! merge the four blocks into one new block. Merging is done in two steps,
          ! first for light data (which all CPUS do redundantly, so light data is kept synched)
          ! Then only the responsible rank will perform the heavy data merging.
          call merge_blocks( params, hvy_block, lgt_block, light_ids(1:N_sisters,k) )
      end do

      if (n_merge> 0_ik) then
        ! the mesh changes if we merge blocks.
        ! we set this flag to true, so that other functions outside this module can have
        ! the information that the mesh has changed
        mesh_has_changed=.true.

       ! free all the light ids of blocks, which have been gathered by a non blocking recv/send.
       ! these blocks where marked with refine status 55 in gahter_blocks_on_proc
       do k = 1, lgt_n
           if ( lgt_block(lgt_active(k), maxtl + idx_refine_sts) == 55) then
            lgt_block(lgt_active(k),:)=-1
           endif
        end do
     endif

end subroutine coarse_mesh



!> this function returns true if the mesh has changed during the last iteration
!> in the time loop
logical function mesh_adapted()
  implicit none

  mesh_adapted=mesh_has_changed
end function
