!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name merge_blocks.f90
!> \version 0.5
!> \author engels
!
!> \brief Merge 4 or 8 blocks into a new one (restriction operator).
!
!> \details
!! This routine merges either 4 or 8 blocks into one new block with half the resolution
!! It is supposed that all blocks are on the same mpirank, thus gather_blocks_on_proc has
!! to be called prior to merging.
!! \n
!! Note ghost nodes are not copied (thus synching is required afterwards)
!! \n
!! Note we keep the light data synchronized among CPUS, so that after moving, all CPU are up-to-date with their light data.
!! However, the active lists are outdated after this routine.
! ********************************************************************************************
subroutine merge_blocks( params, hvy_block, lgt_block, lgt_blocks_to_merge )
  implicit none

  !> user defined parameter structure
  type (type_params), intent(in)      :: params
  !> light data array
  integer(kind=ik), intent(inout)     :: lgt_block(:, :)
  !> heavy data array - block data
  real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
  !> list of blocks to merge, can contain 4 or 8 blocks
  integer(kind=ik), intent(inout)     :: lgt_blocks_to_merge(:)

  ! number of blocks to be merged, can be 4 or 8
  integer(kind=ik) :: N_merge
  ! what CPU is responsible for merging:
  integer(kind=ik) :: data_rank(8)
  ! list of block ids, proc ranks
  integer(kind=ik) :: heavy_ids(8)

  integer(kind=ik) :: i1, i2, im, i, g, level, lgt_merge_id, maxtL, hvy_merge_id
  integer(kind=ik), dimension(3) ::  bound1, bound2, boundm, Bs

!---------------------------------------------------------------------------------------------
! variables initialization

  ! number of blocks to be merged
  N_merge = size(lgt_blocks_to_merge,1)
  Bs = params%Bs
  g  = params%n_ghosts
  maxtL = params%max_treelevel
  ! level of merged block
  level = lgt_block( lgt_blocks_to_merge(1), maxtL + idx_mesh_lvl )
!  write(*,*) "level= ",level
  if ( N_merge /= 4 .and. N_merge /= 8) then
    call abort('You try to merge neither n=4 or 8 blocks...this cannot work.')
  endif

  ! Check which CPU holds the blocks. The CPU will also hold the merged, new block
  do i = 1, N_merge
    call lgt_id_to_proc_rank( data_rank(i), lgt_blocks_to_merge(i), params%number_blocks )
  enddo

  ! Check if all blocks lie on the same rank
  if ( maxval(data_rank(1:N_merge)-data_rank(1)) /= 0 ) then
    call abort("You try to merge blocks on different ranks, but you must call gather_ranks before.")
  endif

  do i = 1, size(lgt_blocks_to_merge)
      if (lgt_block(lgt_blocks_to_merge(i),level) /= i-1) then
          call abort(647483," You try to merge blocks which do not belong together")
      endif
      do i1 = 1, level-1
          if (lgt_block(lgt_blocks_to_merge(i),i1) /= lgt_block(lgt_blocks_to_merge(1),i1)) then
              call abort(647483," You try to merge blocks which do not belong together")
          endif
      enddo
  enddo


!---------------------------------------------------------------------------------------------
! main body

  ! merge the specified blocks into one new block. Merging is done in two steps,
  ! first for light data (which all CPUS do redundantly, so light data is kept synched)
  ! Then only the responsible rank will perform the heavy data merging.

  ! a) light data (collective operation)
  ! fetch a free light ID for the merged blocks
  call get_free_local_light_id( params, data_rank(1), lgt_block, lgt_merge_id)
  ! create light data entry for the new block
  lgt_block( lgt_merge_id, : ) = -1
  lgt_block( lgt_merge_id, 1:level-1 ) = lgt_block( lgt_blocks_to_merge(1), 1:level-1 )
  lgt_block( lgt_merge_id, maxtl+ idx_mesh_lvl ) = level-1
  lgt_block( lgt_merge_id, maxtl+ idx_refine_sts ) = 0


  ! b) heavy data merging (individual operation)
  if ( data_rank(1) == params%rank) then
      ! loop over the sisters and store the corresponding mpirank and heavy index
      ! in a list (still, we are dealing with all four sisters)
      do i = 1, N_merge
        call lgt_id_to_hvy_id( heavy_ids(i), lgt_blocks_to_merge(i), data_rank(1), params%number_blocks )
      enddo

      ! get heavy id of merge block
      call lgt_id_to_hvy_id( hvy_merge_id, lgt_merge_id, data_rank(1), params%number_blocks )

      do i = 1,3
        bound1(i) = g+1
        bound2(i) = Bs(i)+g
        boundm(i) = (Bs(i)-1)/2 + 1 + g
      enddo

      if (N_merge == 4) then
        ! ************ 2D case ***********************
        ! sister 0
        hvy_block(bound1(1):boundm(1), bound1(2):boundm(2), :, :, hvy_merge_id) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, :,:, heavy_ids(1) )
        ! sister 1
        hvy_block(bound1(1):boundm(1), boundm(2):bound2(2), :, :, hvy_merge_id) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, :,:, heavy_ids(2) )
        ! sister 2
        hvy_block(boundm(1):bound2(1), bound1(2):boundm(2), :, :, hvy_merge_id) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, :,:, heavy_ids(3) )
        ! sister 3
        hvy_block(boundm(1):bound2(1), boundm(2):bound2(2), :, :, hvy_merge_id) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, :,:, heavy_ids(4) )

      elseif (N_merge == 8) then
        ! ************ 3D case ***********************
        ! sister 0
        hvy_block(bound1(1):boundm(1), bound1(2):boundm(2), bound1(3):boundm(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(1) )
        ! sister 1
        hvy_block(bound1(1):boundm(1), boundm(2):bound2(2), bound1(3):boundm(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(2) )
        ! sister 2
        hvy_block(boundm(1):bound2(1), bound1(2):boundm(2), bound1(3):boundm(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(3) )
        ! sister 3
        hvy_block(boundm(1):bound2(1), boundm(2):bound2(2), bound1(3):boundm(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(4) )
        ! sister 4
        hvy_block(bound1(1):boundm(1), bound1(2):boundm(2), boundm(3):bound2(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(5) )
        ! sister 5
        hvy_block(bound1(1):boundm(1), boundm(2):bound2(2), boundm(3):bound2(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(6) )
        ! sister 6
        hvy_block(boundm(1):bound2(1), bound1(2):boundm(2), boundm(3):bound2(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(7) )
        ! sister 7
        hvy_block(boundm(1):bound2(1), boundm(2):bound2(2), boundm(3):bound2(3), :, hvy_merge_id ) = hvy_block( g+1:Bs(1)+g:2, g+1:Bs(2)+g:2, g+1:Bs(3)+g:2, :, heavy_ids(8) )

      endif

      ! delete old heavy data (recall that they are all on the same CPU)
      do i = 1, N_merge
          hvy_block( :,:,:,:,heavy_ids(i) ) = 5.0e5_rk
      enddo
  end if


  ! merging is complete now, remove the original blocks from light data:
  do i = 1, N_merge
    lgt_block( lgt_blocks_to_merge(i), : ) = -1
  enddo
end subroutine merge_blocks
