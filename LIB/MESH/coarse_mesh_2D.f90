!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name coarse_mesh_2D.f90
!> \version 0.4
!> \author msr
!
!> \brief coarse the mesh: \n
!! every proc work on light data array \n
!
!> 
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4, split old interpolate_mesh subroutine into two refine/coarsen
!            subroutines
! ********************************************************************************************

subroutine coarse_mesh_2D( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )

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
    integer(kind=ik)                    :: k, i, N, lgt_merge_id, hvy_merge_id
    ! grid parameter
    integer(kind=ik)                    :: Bs, g, i1, im, i2
    ! treecode varaible
    integer(kind=ik)                    :: me_treecode(params%max_treelevel)
    ! max treecode level
    integer(kind=ik)                    :: maxtL
    ! mesh level
    integer(kind=ik)                    :: level
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: light_ids(4), proc_rank(4), heavy_ids(4)
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: data_rank
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    sub_t0 = MPI_Wtime()
    N = params%number_blocks
    maxtL = params%max_treelevel
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

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
              call lgt_id_to_proc_rank( data_rank, lgt_active(k), N )

              ! get treecodes for current block and sister-blocks
              me_treecode = lgt_block( lgt_active(k), 1:maxtL)
              ! current block level
              level = lgt_block( lgt_active(k), maxtL+1 )

              ! find all sisters (including the block in question, so four blocks)
              ! their light IDs are in "light_ids" and ordered by their last treecode-digit
              call find_sisters( params, lgt_active(k), light_ids, lgt_block, lgt_n, lgt_sortednumlist )

              ! gather all four sisters on the process "datarank". The light_ids are updated in the routine
              ! and they are still in the same order (0,1,2,3)-sister. It is just that they are now on one CPU
              call gather_blocks_on_proc( params, hvy_block, lgt_block, data_rank, light_ids )

              ! merge the four blocks into one new block. Merging is done in two steps,
              ! first for light data (which all CPUS do redundantly, so light data is kept synched)
              ! Then only the responsible rank will perform the heavy data merging.

              ! a) light data
              ! fetch a free light ID for the merged blocks
              call get_free_local_light_id( params, data_rank, lgt_block, lgt_merge_id)
              ! create light data entry for the new block
              lgt_block( lgt_merge_id, : ) = -1
              lgt_block( lgt_merge_id, 1:level-1 ) = me_treecode(1:level-1)
              lgt_block( lgt_merge_id, maxtl+1 ) = level-1
              lgt_block( lgt_merge_id, maxtl+2 ) = 0

              ! b) heavy data merging
              if (data_rank == params%rank) then
                  ! loop over the sisters and store the corresponding mpirank and heavy index
                  ! in a list (still, we are dealing with all four sisters)
                  do i = 1, 4
                    call lgt_id_to_proc_rank( proc_rank(i), light_ids(i), N )
                    call lgt_id_to_hvy_id( heavy_ids(i), light_ids(i), proc_rank(i), N )
                  enddo

                  ! get heavy id of merge block
                  call lgt_id_to_hvy_id( hvy_merge_id, lgt_merge_id, data_rank, N )

                  i1 = g+1
                  i2 = Bs+g
                  im = (Bs-1)/2 + 1 + g
                  ! sister 0
                  hvy_block(i1:im, i1:im, :, :, hvy_merge_id) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, :,:, heavy_ids(1) )
                  ! sister 1
                  hvy_block(i1:im, im:i2, :, :, hvy_merge_id) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, :,:, heavy_ids(2) )
                  ! sister 2
                  hvy_block(im:i2, i1:im, :, :, hvy_merge_id) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, :,:, heavy_ids(3) )
                  ! sister 3
                  hvy_block(im:i2, im:i2, :, :, hvy_merge_id) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, :,:, heavy_ids(4) )

                  ! delete old heavy data (recall that they are all on the same CPU)
                  do i = 1, 4
                      hvy_block( :,:,:,:,heavy_ids(i) ) = 5.0e5_rk
                  enddo
              endif

              ! merging is complete now, remove the original blocks from light data:
              do i = 1, 4
                lgt_block( light_ids(i), : ) = -1
                lgt_block( light_ids(i), maxtl+1 ) = -1
                lgt_block( light_ids(i), maxtl+2 ) = 0
              enddo
            endif
    end do


    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "coarse_mesh" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "coarse_mesh"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine coarse_mesh_2D
