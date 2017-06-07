!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name coarse_mesh.f90
!> \version 0.5
!> \author engels, msr
!
!> \brief Apply mesh coarsening: Merge tagged blocks into new, coarser blocks
!
!!
!!
!! = log ======================================================================================
!! \n
!! 22/05/2017 Rewrite using modular subroutines, works for 2/3d cases
!! 08/11/16 - switch to v0.4, split old interpolate_mesh subroutine into two refine/coarsen
!            subroutines
! ********************************************************************************************

subroutine coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

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

    ! loop variables
    integer(kind=ik)                    :: k, maxtl
    ! list of block ids, proc ranks
    integer(kind=ik), allocatable       :: light_ids(:)
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: data_rank


!---------------------------------------------------------------------------------------------
! variables initialization

    maxtL = params%max_treelevel

    if (params%threeD_case) then
      allocate( light_ids(1:8) )
    else
      allocate( light_ids(1:4) )
    endif

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active light data
    do k = 1, lgt_n
      ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
      ! the active list. NOTE HERE: due to previous loops, some light data id are already
      ! coarsened, (and thus given the -1) but they are still in active block list
      ! so: check again if this block is REALLY active (the active list is updated later)
      !
      ! SECOND condition: block wants to coarsen, definetly, i.e. it has the status -2. Note -1
      ! is just a temporary state and four blocks, all -1, are given the -2 state
      ! is ensure_completeness finds all four. So that means, by extension, we previously
      ! searched all 4 sister blocks, and we found them all.
      if ( lgt_block(lgt_active(k), 1) >= 0 .and. lgt_block(lgt_active(k), maxtl+2) == -2) then
          ! Check which CPU holds this block. The CPU will also hold the merged, new block
          call lgt_id_to_proc_rank( data_rank, lgt_active(k), params%number_blocks )

          ! find all sisters (including the block in question, so four blocks)
          ! their light IDs are in "light_ids" and ordered by their last treecode-digit
          call find_sisters( params, lgt_active(k), light_ids, lgt_block, lgt_n, lgt_sortednumlist )

          ! gather all four sisters on the process "datarank". The light_ids are updated in the routine
          ! and they are still in the same order (0,1,2,3)-sister. It is just that they are now on one CPU
          call gather_blocks_on_proc( params, hvy_block, lgt_block, data_rank, light_ids )

          ! merge the four blocks into one new block. Merging is done in two steps,
          ! first for light data (which all CPUS do redundantly, so light data is kept synched)
          ! Then only the responsible rank will perform the heavy data merging.
          call merge_blocks( params, hvy_block, lgt_block, light_ids )
        endif
    end do

    deallocate( light_ids )

end subroutine coarse_mesh
