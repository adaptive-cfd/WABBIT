!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name gather_blocks_on_proc.f90
!> \version 0.5
!> \author engels
!
!> \brief Gather a specified list of blocks on one specified CPU. Used for e.g. for block merging.
!
!> \details
!! This routine gathers all block specified in the list "lgt_blocks_to_gather" on one CPU, namely the CPU
!! "gather_rank". It does not hurt if some or all block are already on the CPU, nothing is done in that case.
!! Ghost nodes are transferred as well, works for 2D/3D data of any size and for any number of blocks to be transferred
!! \n
!! Note we keep the light data synchronized among CPUS, so that after moving, all CPU are up-to-date with their light data.
!! However, the active lists are outdated after this routine.
!!
!! Uses blocking communication.
! ********************************************************************************************
subroutine gather_blocks_on_proc( params, hvy_block, lgt_block, gather_rank, lgt_blocks_to_gather )
  implicit none

  !> user defined parameter structure
  type (type_params), intent(in)      :: params
  !> light data array
  integer(kind=ik), intent(inout)     :: lgt_block(:, :)
  !> heavy data array - block data
  real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
  !> list of blocks (lightIDs) to be gathered. on output, the new IDs of
  !! these blocks are returned.
  integer(kind=ik), intent(inout)     :: lgt_blocks_to_gather(:)
  !> on which MPIRANK will be gather the desired blocks?
  integer(kind=ik), intent(in)        :: gather_rank

  integer(kind=ik) :: i, hvy_id, owner_rank, myrank, npoints
  integer(kind=ik) :: lgt_free_id, tag, hvy_free_id, ierr

  myrank = params%rank

  ! look at all blocks in the gather list
  do i = 1, size(lgt_blocks_to_gather)
      ! which mpirank owns the block?
      call lgt_id_to_proc_rank( owner_rank, lgt_blocks_to_gather(i), params%number_blocks )

      ! if gather_rank == owner_rank, nothing to be done, as block is already on
      ! right processor
      if ( gather_rank /= owner_rank) then
          ! use a unique tag for message exchange
          tag = lgt_blocks_to_gather(i)

          ! get a free place to gather the block ("lgt_free_id") on the gather_rank
          ! for light data synching, it is good to have all procs do that
          call get_free_local_light_id( params, gather_rank, lgt_block, lgt_free_id)

          ! Am I the target rank who receives all the data?
          if ( myrank == gather_rank) then
              !------------------------
              ! RECV CASE
              !------------------------
              ! get hvy id where to store the data
              call lgt_id_to_hvy_id( hvy_free_id, lgt_free_id, myrank, params%number_blocks )

              npoints = size(hvy_block,1)*size(hvy_block,2)*size(hvy_block,3)*size(hvy_block,4)
              call MPI_recv( hvy_block(:,:,:,:,hvy_free_id), npoints, MPI_REAL8, owner_rank, tag, MPI_COMM_WORLD, ierr)

          elseif ( myrank == owner_rank) then
              ! Am I the owner of this block, so will I have to send data?
              !------------------------
              ! SEND CASE
              !------------------------
              ! what heavy ID (on this proc) does the block have?
              call lgt_id_to_hvy_id( hvy_id, lgt_blocks_to_gather(i), owner_rank, params%number_blocks )

              npoints = size(hvy_block,1)*size(hvy_block,2)*size(hvy_block,3)*size(hvy_block,4)
              call MPI_send( hvy_block(:,:,:,:,hvy_id), npoints, MPI_REAL8, gather_rank, tag, MPI_COMM_WORLD, ierr)

              ! delete old heavy data
              hvy_block(:,:,:,:,hvy_id) = 4.0e11_rk
          endif

          ! even if I am not concerned with sending or recv, the light data changes, assuming the copy went through.
          ! if it did not, code hangs anyways. so here assume it worked and on all CPU just copy the light data
          lgt_block( lgt_free_id, : ) = lgt_block( lgt_blocks_to_gather(i), : )
          ! we must also delete the original block (also done on ALL CPUS)
          lgt_block( lgt_blocks_to_gather(i), : ) = -1
          ! we return the new index of the blocks we moved
          lgt_blocks_to_gather(i) = lgt_free_id
      endif
  enddo

! NOTE; even though light data is synched, active lists etc have to be created afterwards!
end subroutine gather_blocks_on_proc
