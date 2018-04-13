!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_equidistant_base_mesh.f90
!> \version 0.5
!> \author engels
!
!> \brief This routine tries to create an equidistant mesh on the specified level Jmin, so all blocks on
!! this level are set to active and their treecodes are stored.
!! Since the grid changes, the neighbor relations and active-lists are updated as well.
!
!> \note In almost all cases, you'll want to call reset_grid before call this routine, as any existing
!! blocks are not deleted here. Not calling reset_grid is possible but will cause strange behavior and
!! probably is unintended. This routine warns if active blocks are found. \n
!!
!!
!! = log ======================================================================================
!! \n
!! 03 Apr 2017 - create
!
! ********************************************************************************************

subroutine create_equidistant_base_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, Jmin, verbosity )

  !---------------------------------------------------------------------------------------------
  ! variables

  implicit none

  !> user defined parameter structure
  type (type_params), intent(in)      :: params
  !> light data array
  integer(kind=ik), intent(inout)     :: lgt_block(:, :)
  !> heavy data array - block data
  real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
  !> heavy data array - neighbor data
  integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
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
  !> what level to initialize?
  integer(kind=ik), intent(in)        :: Jmin

  !> write output
  logical, intent(in)                 :: verbosity

  ! MPI error variable
  integer(kind=ik)                    :: ierr

  ! loop control variables in space
  integer(kind=ik)                    :: ix, iy, iz, icpu, nx, ny, nz

  integer(kind=ik)                    :: num_blocks, number_procs
  integer(kind=ik)                    :: d

  ! on my section, which is the first and last light id?
  integer(kind=ik)                    :: lgt_id_first, lgt_id_last, lgt_id, heavy_id
  integer(kind=ik),allocatable        :: blocks_per_rank_list(:)
  integer(kind=ik),allocatable        :: treecode(:), my_lgt_block(:,:)

  !-----------------------------------------------------------------------------
  ! initialization
  !-----------------------------------------------------------------------------

  ! data demensionality
  if ( params%threeD_case) then
    d = 3
    nz = 2**Jmin
  else
    d = 2
    nz = 1
  endif
  ! therefore, the total number of blocks (on all cpus) is
  num_blocks = (2**Jmin)**d
  ! shortcut for number of cpu
  number_procs = params%number_procs
  ! reset light data list (note setting -1 is required for get_free_light_id to work!)
  lgt_block = -1
  ! number of blocks in x,y direction
  nx = 2**Jmin
  ny = 2**Jmin

  if ( (params%rank == 0) .and. verbosity ) then
    write(*,'(A)') "INIT: initializing an equidistant grid..."
    write(*,'(" Jmin=",i2," Nblocks=",i6," (on all cpus)")') Jmin, num_blocks
    write(*,'(" on this level, we have (",i3," x ",i3," x ",i3,") Blocks")') nx, ny, nz
    ! check if there is already some data on the grid (see routine header description)
    if ( maxval(lgt_block(:,1))>=0 ) then
      write(*,'(A)') "ERROR: CREATE-EQUIDISTANT MESH is called with NON_EMPTY DATA!!!!!"
    endif
  endif

  !-----------------------------------------------------------------------------
  ! Nblocks per CPU
  !-----------------------------------------------------------------------------
  ! determine how many blocks each processor has to hold so that the total number
  ! is equal to the desired one on the coarsest grid. Note in many situations, the
  ! initial number of blocks on the coarsest level is quite low, so many cpu may end up
  ! not holding any active blocks, but the code will modify this later, once new blocks
  ! arise.
  ! this list contains (on each mpirank) the number of blocks for each mpirank. note
  ! zero indexing as required by MPI
  allocate( blocks_per_rank_list( 0:number_procs-1 ) )
  ! set list to the average value
  blocks_per_rank_list = (num_blocks - mod(num_blocks,number_procs)) / number_procs
  ! as this does not necessairly work out, distribute remaining blocks on the first CPUs
  if (mod(num_blocks, number_procs) > 0) then
    blocks_per_rank_list(0:mod(num_blocks, number_procs)-1) = blocks_per_rank_list(1:mod(num_blocks, number_procs)) + 1
  end if
  ! some error control -> did we loose blocks? should never happen.
  if ( sum(blocks_per_rank_list) /= num_blocks) then
    call abort("ERROR: on the coarsest grid, we seem to have gained/lost some blocks during distribution...")
  end if


  !-----------------------------------------------------------------------------
  ! Generate and distribute the blocks
  !-----------------------------------------------------------------------------
  ! allocate treecode
  allocate( treecode( params%max_treelevel ) )

  ! loop over blocks in x,y,z directions (in the 2d case, 3rd loop degenerates)
  do ix = 1, nx
    do iy = 1, ny
      do iz = 1, nz
        ! for each of the points (ix,iy,iz), find an mpirank to hold it.
        do icpu = 0, number_procs -1
          ! can the current cpu "icpu" still accept more blocks?
          if ( blocks_per_rank_list(icpu) > 0 ) then
            ! yes, it can.
            ! remove one block from the list
            blocks_per_rank_list(icpu) = blocks_per_rank_list(icpu) - 1

            ! create the new block on that cpu.
            if (params%rank == icpu) then
              ! for this new block, compute the treecode
              if (params%threeD_case) then
                call encoding_3D( treecode, ix, iy, iz, num_blocks, params%max_treelevel )
              else
                call encoding_2D( treecode, ix, iy, 2**Jmin, 2**Jmin, params%max_treelevel )
              endif

              ! on my section of the global light data list, which is the first and last light id I hold?
              call hvy_id_to_lgt_id( lgt_id_first, 1, params%rank, params%number_blocks )
              call hvy_id_to_lgt_id( lgt_id_last, params%number_blocks, params%rank, params%number_blocks )

              ! find and set free heavy data id, note: look for free id in light data, but use
              ! search routine only on MY LOCAL light data -> so, returned id works directly on heavy data
              call get_free_light_id( heavy_id, lgt_block( lgt_id_first:lgt_id_last, 1 ), params%number_blocks )

              ! look for global light id of this block
              call hvy_id_to_lgt_id( lgt_id, heavy_id, params%rank, params%number_blocks )

              ! save treecode in global light id list (NOTE: we need to sync that as only one proc did it..)
              lgt_block( lgt_id, 1:params%max_treelevel ) = treecode
              lgt_block( lgt_id, params%max_treelevel+1 ) = Jmin
              lgt_block( lgt_id, params%max_treelevel+2 ) = 0

              ! reset block data to zero
              hvy_block(:,:,:,:,heavy_id) = real(params%rank, kind=rk)
            end if
            ! as this block is now given to one cpu, we can leave the loop over
            ! cpus and take care of the next one.
            exit
          else
            ! it cannot -> choose another cpu
            cycle
          end if
        end do
      end do
    end do
  end do


  ! synchronize light data. This is necessary as all CPUs above created their blocks locally.
  ! As they all pass the same do loops, the counter array blocks_per_rank_list does not have to
  ! be synced. However, the light data has to.
  allocate ( my_lgt_block(1:size(lgt_block,1),1:size(lgt_block,2)) )
  my_lgt_block = lgt_block
  lgt_block = 0
  call MPI_Allreduce(my_lgt_block, lgt_block, size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_MAX, params%WABBIT_COMM, ierr)
  deallocate( my_lgt_block )

  ! as the grid has changed (we created it here), we now update the heavy and light
  ! active lists, as well as neighbor relations.
  ! update list of sorted nunmerical treecodes, used for finding blocks
  call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
  ! update neighbor relations
  call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

end subroutine create_equidistant_base_mesh
