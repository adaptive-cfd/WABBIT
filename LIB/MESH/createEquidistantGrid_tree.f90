!> \brief This routine tries to create an equidistant mesh on the specified level Jmin, so all blocks on
!! this level are set to active and their treecodes are stored.
!! Since the grid changes, the neighbor relations and active-lists are updated as well.
! ********************************************************************************************

subroutine createEquidistantGrid_tree( params, hvy_block, Jmin, verbosity, tree_ID )
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params                         !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    integer(kind=ik), intent(in)        :: Jmin                           !> what level to initialize?
    logical, intent(in)                 :: verbosity                      !> write output
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik)                    :: ierr                           ! MPI error variable
    integer(kind=ik)                    :: ix, iy, iz, icpu, nx, ny, nz   ! loop control variables in space
    integer(kind=ik)                    :: num_blocks, number_procs, k
    integer(kind=ik)                    :: d

    ! on my section, which is the first and last light id?
    integer(kind=ik)                    :: lgt_id_first, lgt_id_last, lgt_id, heavy_id
    integer(kind=ik),allocatable        :: blocks_per_rank_list(:)
    integer(kind=ik),allocatable        :: treecode(:)

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    !-----------------------------------------------------------------------------
    ! initialization
    !-----------------------------------------------------------------------------

    ! data dimensionality
    d = params%dim
    if (params%dim == 3) then
        nz = 2**Jmin
    else
        nz = 1
    endif

    ! therefore, the total number of blocks (on all cpus) is
    num_blocks = (2**Jmin)**d
    ! shortcut for number of cpu
    number_procs = params%number_procs

    call reset_tree(params, verbosity, tree_ID )

    ! number of blocks in x,y direction
    nx = 2**Jmin
    ny = 2**Jmin

    if ( (params%rank == 0) .and. verbosity ) then
        write(*,'(A)') "EQUI: initializing an equidistant grid..."
        write(*,'("EQUI: Jmin=",i2," Nblocks=",i6," (on all cpus)")') Jmin, num_blocks
        write(*,'("EQUI: On this level, we have (",i3," x ",i3," x ",i3,") Blocks")') nx, ny, nz
        write(*,'("EQUI: tree_ID=",i2)') tree_ID

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
    allocate( treecode( params%Jmax ) )
    ! loop over blocks in x,y,z directions (in the 2d case, 3rd loop degenerates)
    ! NOTE: This ordering is necessary for POSTPROCESSING flusi to wabbit!
    do ix = nx, 1, -1
        do iy = ny,1,-1
            do iz = nz,1,-1
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

                            call encoding(treecode, (/ix,iy,iz/), d, num_blocks, Jmin)

                            ! on my section of the global light data list, which is the first and last light id I hold?
                            call hvy2lgt( lgt_id_first, 1, params%rank, params%number_blocks )
                            call hvy2lgt( lgt_id_last, params%number_blocks, params%rank, params%number_blocks )

                            ! get a free block on this rank
                            call get_free_local_light_id( params, icpu, lgt_id)
                            ! and the corresponding heavy id
                            call lgt2hvy( heavy_id, lgt_id, icpu, params%number_blocks )

                            ! save treecode in global light id list (NOTE: we need to sync that as only one proc did it..)

                            lgt_block( lgt_id, : ) = -1
                            lgt_block( lgt_id, 1:Jmin ) = treecode(1:Jmin)
                            lgt_block( lgt_id, params%Jmax+IDX_MESH_LVL ) = Jmin
                            lgt_block( lgt_id, params%Jmax+IDX_REFINE_STS ) = 0
                            lgt_block( lgt_id, params%Jmax+IDX_TREE_ID ) = tree_ID
                        end if
                        ! as this block is now given to one cpu, we can leave the loop over
                        ! cpus and take care of the next one.
                        exit
                    else
                        ! it cannot -> choose another cpu
                        cycle
                    end if
                end do ! loop over cpu

            end do
        end do
    end do


    ! synchronize light data. This is necessary as all CPUs above created their blocks locally.
    ! As they all pass the same do loops, the counter array blocks_per_rank_list does not have to
    ! be synced. However, the light data has to.
    call synchronize_lgt_data( params, refinement_status_only=.false. )

    ! as the grid has changed (we created it here), we now update the heavy and light
    ! active lists, as well as neighbor relations.
    call updateMetadata_tree(params, tree_ID)

    ! the grid we created above is not at all balanced - on the contrary, it fills up CPU1, then CPU2
    ! correct this by balancing:
    call balanceLoad_tree(params, hvy_block, tree_ID)
end subroutine createEquidistantGrid_tree
