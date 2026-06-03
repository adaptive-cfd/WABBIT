!> \image html balancing.svg "Load balancing" width=400
! ********************************************************************************************
subroutine balanceLoad_tree( params, hvy_block, tree_ID, balanceMode, balance_ref, balance_level, &
                             balance_name, time, hvy_tmp, Jmin_set, full_tree_grid )
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)      :: params                     !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    
    !> Balance mode: controls the balancing strategy
    !! "standard"           - distribute all blocks evenly along SFC (default)
    !! "refine_predictive"  - blocks with +1 status count as 2^dim blocks (for pre-refinement balancing)
    !! "selective"          - only balance blocks matching filter, others stay on current rank
    character(len=*), intent(in), optional :: balanceMode
    
    !> Refinement status filter: balance only blocks with these ref values (array input, e.g., (/ +1, 0 /))
    integer(kind=ik), intent(in), optional :: balance_ref(:)
    
    !> Level filter: balance only blocks on these levels (array input, e.g., (/ Jmax /) or (/ Jmax-1, Jmax /))
    integer(kind=ik), intent(in), optional :: balance_level(:)
    
    !> balanceLoad is always the same, however we want to differentiate between them, so we can assign a name
    character(len=*), intent(in), optional :: balance_name
    !> current simulation time, used for development output only
    real(kind=rk), intent(in), optional :: time
    !> optional second heavy data array - transferred with same pattern as hvy_block
    real(kind=rk), intent(inout), optional :: hvy_tmp(:, :, :, :, :)
    !> Jmin value override for metadata update after balancing
    integer(kind=ik), intent(in), optional :: Jmin_set
    !> full tree grid flag, needs to be set to true if balancing is done on a full tree grid instead of a leaf-only grid
    logical, intent(in), optional :: full_tree_grid
    !=============================================================================

    integer(kind=ik)  :: rank, mpirank_shouldBe, mpirank_currently, ierr, number_procs, &
        k, N, l, com_i, com_N, sfc_id, lgt_free_id, hvy_free_id, hilbertcode(params%Jmax), hvy_id, lgt_ID, Nblocks_toDistribute, &
        num_blocks_count(3), i_var, balance_id, sfc_i, Jmin_set_apply
    integer(kind=tsize) :: treecode, hilbertcode2
    character(len=cshort) :: mode
    ! block distribution lists
    integer(kind=ik), allocatable, save :: blocksPerRank_balanced(:), sfc_com_list(:,:), sfc_sorted_list(:,:)
    integer(kind=ik), allocatable :: unassigned_per_rank(:)
    real(kind=rk) :: t0, t1
    character(len=cshort) :: format_string, string_prepare
    character(len=4) :: string_kind(3)
    logical :: full_tree_grid_apply


    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! on only one proc, no balancing is required
    if (params%number_procs == 1) return

    ! start time
    t0 = MPI_Wtime()

    ! MPI_parameters
    rank = params%rank
    number_procs = params%number_procs

    Jmin_set_apply = params%Jmin
    if (present(Jmin_set)) Jmin_set_apply = Jmin_set
    full_tree_grid_apply = .false.
    if (present(full_tree_grid)) full_tree_grid_apply = full_tree_grid

    !---------------------------------------------------------------------------
    ! Set balance mode
    !---------------------------------------------------------------------------
    mode = "standard"  ! default
    if (present(balanceMode)) mode = balanceMode

    !---------------------------------------------------------------------------
    ! Translate balanceMode to balance_id to avoid string comparisons in loops
    ! (following pattern from prepare_ghost_synch_metadata)
    !---------------------------------------------------------------------------
    ! First digit: base mode
    select case(trim(mode))
    case("standard")
        balance_id = 1
        ! Standard balancing: distribute all blocks evenly along SFC
        
    case("refine_predictive")
        balance_id = 2
        ! Predictive refinement: blocks with +1 status count as 2^dim blocks
        
    case("selective")
        balance_id = 3
        ! Selective: only balance blocks matching filter, others stay put
        
    case default
        call abort(260112, "balanceLoad_tree: Unknown balance mode: " // trim(mode))
    end select

    ! Second digit: filter type
    if (present(balance_ref) .and. present(balance_level)) then
        ! Both ref and level filter active
        balance_id = balance_id + 10*3
    elseif (present(balance_ref)) then
        ! Only ref filter active
        balance_id = balance_id + 10*1
    elseif (present(balance_level)) then
        ! Only level filter active
        balance_id = balance_id + 10*2
    endif

    ! Validation: selective mode requires filters
    if (balance_id == 3) then
        call abort(260112, "balanceLoad_tree: selective mode requires balance_ref or balance_level!")
    endif

    ! development counters, just to see how many blocks are send/received/kept
    num_blocks_count(1:3) = 0

    ! allocate block to proc lists
    if (.not.allocated(blocksPerRank_balanced)) allocate( blocksPerRank_balanced(1:number_procs))
    ! allocate sfc com list, maximal number of communications is when every proc wants to send all of his blocks
    ! NOTE: it is not necessary or wise to reset this array (it is large!)
    if (.not.allocated(sfc_com_list)) allocate( sfc_com_list( 3, number_procs*params%number_blocks ) )

    ! allocate space filling curve list, number of elements is the number of active blocks
    ! and for each block, we store the space-filling-curve-index and the lgt ID
    ! NOTE: it is not necessary or wise to reset this array (it is large!)
    if (.not.allocated(sfc_sorted_list)) allocate( sfc_sorted_list( 3, size(lgt_block,1)) )

    ! number of blocks
    N = params%number_blocks

    !---------------------------------------------------------------------------------
    ! Calculate how many blocks to distribute among mpiranks
    !---------------------------------------------------------------------------------
    select case(mod(balance_id, 10))  ! First digit = base mode
    case(1)  ! standard
        ! Standard balancing: distribute all blocks that are on the grid right now
        Nblocks_toDistribute = lgt_n(tree_ID)
        
    case(2)  ! refine_predictive
        ! Predictive refinement: we distribute blocks such that the distribution 
        ! AFTER refinement is balanced. Blocks with status +1 will be refined, 
        ! resulting in 2**dim blocks, all others remain unchanged.
        Nblocks_toDistribute = 0
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)
            if ( lgt_block(lgt_id , IDX_REFINE_STS) == +1) then
                ! block will be refined, resulting in 2**D blocks
                Nblocks_toDistribute = Nblocks_toDistribute + 2**params%dim
            else
                ! block will not be refined, resulting in 1 block
                Nblocks_toDistribute = Nblocks_toDistribute + 1
            endif        
        enddo
        
    case(3)  ! selective
        ! Selective balancing: count assigned blocks globally and unassigned blocks per rank
        ! Use lgt2proc to determine which rank owns each unassigned block
        Nblocks_toDistribute = 0
        allocate(unassigned_per_rank(0:number_procs-1))
        unassigned_per_rank = 0
        
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)
            if (block_matches_filter(lgt_id, balance_id, balance_ref, balance_level)) then
                ! Assigned block: will be balanced via SFC
                Nblocks_toDistribute = Nblocks_toDistribute + 1
            else
                ! Unassigned block: count per rank using lgt2proc
                call lgt2proc(mpirank_currently, lgt_id, params%number_blocks)
                unassigned_per_rank(mpirank_currently) = unassigned_per_rank(mpirank_currently) + 1
            endif
        enddo
    end select

    !---------------------------------------------------------------------------------
    ! 1st: define how many blocks each mpirank should have.
    !---------------------------------------------------------------------------------

    ! Calculate balanced distribution: distribute evenly across all ranks
    ! In selective mode: only assigned blocks, overflow handling in plan_communication_selective
    ! In standard/refine_predictive: all blocks
    blocksPerRank_balanced(:) = Nblocks_toDistribute / number_procs

    ! as this does not necessarily work out, distribute remaining blocks on the first CPUs
    blocksPerRank_balanced(1:mod(Nblocks_toDistribute, number_procs)) = blocksPerRank_balanced(1:mod(Nblocks_toDistribute, number_procs)) + 1

    ! some error control -> did we loose blocks? should never happen.
#ifdef DEV
    if ( sum(blocksPerRank_balanced) /= Nblocks_toDistribute) then
        call abort(1028,"ERROR: we seem to have gained/lost some blocks during distribution...")
    end if
#endif

    !---------------------------------------------------------------------------------
    ! 2nd: calculate space filling curve index for blocks to be balanced
    !---------------------------------------------------------------------------------
    t1 = MPI_wtime()
    sfc_i = 0  ! counter for blocks added to SFC list    
    select case (params%block_distribution)
    case("sfc_z")
        !-----------------------------------------------------------
        ! Z - curve - position is simply treecode
        !-----------------------------------------------------------
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)
            
            ! In selective mode: skip blocks that don't match filter
            if (mod(balance_id, 10) == 3) then  ! selective mode
                if (.not. block_matches_filter(lgt_id, balance_id, balance_ref, balance_level)) cycle
            endif
            
            sfc_i = sfc_i + 1
            sfc_sorted_list(:, sfc_i) = -1
            sfc_sorted_list(1, sfc_i) = lgt_id
            treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
            call set_tc(sfc_sorted_list(2:3, sfc_i), treecode)
        end do
    case("sfc_hilbert")
        !-----------------------------------------------------------
        ! Hilbert curve - we need to translate the treecode
        !-----------------------------------------------------------
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)
            
            ! In selective mode: skip blocks that don't match filter
            if (mod(balance_id, 10) == 3) then  ! selective mode
                if (.not. block_matches_filter(lgt_id, balance_id, balance_ref, balance_level)) cycle
            endif
            
            ! transfer treecode to hilbertcode
            treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
            if (params%dim == 3) then
                call treecode_to_hilbertcode_3D( treecode, hilbertcode2, dim=params%dim, &
                    level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)
            else
                call treecode_to_hilbertcode_2D( treecode, hilbertcode2, dim=params%dim, &
                    level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)
            endif
 
            sfc_i = sfc_i + 1
            sfc_sorted_list(:, sfc_i) = -1
            sfc_sorted_list(1, sfc_i) = lgt_id
            call set_tc(sfc_sorted_list(2:3, sfc_i), hilbertcode2)
        end do

    case default
        call abort(210309, "The block_dist is unkown: "//params%block_distribution)
    end select

    ! sort sfc_list according to the first dimension, thus the position on
    ! the space filling curve (this was a bug, fixed: Thomas, 13/03/2018)
    ! For non-selective balancing, sfc_i == lgt_n(tree_ID)
    if (sfc_i > 1) then
        call quicksort(sfc_sorted_list, 1, sfc_i, 3)
    end if
    call toc( "balanceLoad_tree (SFC+sort)", 91, MPI_wtime()-t1 )

    !---------------------------------------------------------------------------------
    ! 3rd: plan communication (fill list of blocks to transfer)
    !---------------------------------------------------------------------------------
    t1 = MPI_wtime()
    
    ! Dispatch to appropriate communication planning routine based on balance mode
    select case(mod(balance_id, 10))  ! First digit = base mode
    case(1, 2)  ! standard or refine_predictive
        call plan_communication_standard(params, sfc_sorted_list, sfc_i, blocksPerRank_balanced, &
                                         balance_id, sfc_com_list, com_i, num_blocks_count)
    case(3)  ! selective
        call plan_communication_selective(params, tree_ID, sfc_sorted_list, sfc_i, blocksPerRank_balanced, &
                                          unassigned_per_rank, balance_id, balance_ref, balance_level, &
                                          sfc_com_list, com_i, num_blocks_count)
        ! Clean up
        deallocate(unassigned_per_rank)
    end select
    
    call toc( "balanceLoad_tree (comm)", 92, MPI_wtime()-t1 )

    !---------------------------------------------------------------------------------
    ! 4th: actual communication (send/recv)
    !---------------------------------------------------------------------------------
    call block_xfer( params, sfc_com_list(1:3, 1:com_i), com_i, hvy_block, msg="balanceLoad_tree", hvy_tmp=hvy_tmp )

    ! the block xfer changes the light data, and afterwards active lists are outdated.
    call updateMetadata_tree(params, tree_ID, search_overlapping=full_tree_grid_apply, Jmin_set=Jmin_set_apply)

    ! development output, gather information on rank 0 and print to file
    ! reuse blocksPerRank_balanced for this
    if (params%debug_balanceLoad) then
        string_kind = (/ "send", "keep", "recv" /)
        
        do i_var = 1, 3
            ! gather number of blocks send/received/kept on rank 0
            call MPI_GATHER(num_blocks_count(i_var), 1, MPI_INTEGER, blocksPerRank_balanced, 1, MPI_INTEGER, 0, WABBIT_COMM, ierr)
            if (params%rank == 0) then
                if (present(balance_name)) then
                    string_prepare = trim(adjustl(balance_name))//","
                else
                    string_prepare = "balanceLoad,"
                endif
                string_prepare = trim(adjustl(string_prepare))//trim(adjustl(string_kind(i_var)))//","
                if (present(time)) then
                    write(string_prepare,'(A,es12.5,",")') trim(adjustl(string_prepare)), time
                else
                    string_prepare = trim(adjustl(string_prepare))//"-1.0E+00,"  ! set negative time, just to have csv with the same length in every row
                endif
                ! Single IO operation with dynamic format
                open(unit=99, file='debug_balanceLoad.csv', status="unknown", position="append")
                if (params%number_procs == 1) then
                    write(99, '(A, i0)') trim(adjustl(string_prepare)), blocksPerRank_balanced(1)
                else
                    write(format_string, '("(A, i0,",i0,"("","",i0))")') params%number_procs - 1
                    write(99, format_string) trim(adjustl(string_prepare)), blocksPerRank_balanced(1), blocksPerRank_balanced(2:params%number_procs)
                endif
                close(99)
            endif
        enddo
    endif

    ! timing
    call toc( "balanceLoad_tree (TOTAL)", 90, MPI_wtime()-t0 )

contains

    !> \brief Plan communication for selective balancing mode
    !> \details Three-phase algorithm for selective load balancing:
    !!
    !! OVERVIEW:
    !! Selective mode balances only blocks matching filter criteria (refinement status/level)
    !! via SFC, while non-matching blocks stay on their current ranks. This can create
    !! overflow: a rank may receive too many assigned blocks + its unassigned blocks exceed capacity.
    !!
    !! PHASE 1: SFC Distribution of Assigned Blocks
    !! - Walk through SFC-sorted list of assigned blocks (filtered by balance_ref/balance_level)
    !! - Assign blocks to ranks sequentially along SFC for spatial locality
    !! - Track how many assigned blocks each rank will receive in blocks_assigned_per_rank
    !! - This may create overflow on some ranks (assigned + unassigned > capacity)
    !!
    !! PHASE 2: Compute Overflow and Free Space
    !! - overflow_or_free(rank) = capacity - (assigned + unassigned)
    !! - Positive values: free space (can receive overflow blocks)
    !! - Negative values: overflow (must evacuate unassigned blocks)
    !! - All ranks compute identical arrays (collective loops over global metadata)
    !!
    !! PHASE 3: Proportional Redistribution of Overflow
    !! - Evacuate unassigned blocks from overflowing ranks
    !! - Distribute proportionally to ranks with free space using efficient crawling:
    !!   * ratio = total_free / total_overflow
    !!   * Each overflow block i gets target position int(i * ratio) in virtual free space
    !!   * Crawl through ranks to find which rank owns that position
    !!   * Crawling state never resets (monotonic target positions) → O(n) complexity
    !! - Processing in lgt_active order ensures deterministic distribution
    !!
    !! NOTE: All loops are COLLECTIVE operations using global metadata (lgt_block, lgt_active)
    !! so all ranks compute identical block assignments without MPI synchronization
    subroutine plan_communication_selective(params, tree_ID, sfc_sorted_list, N_sfc_blocks, &
                                           blocksPerRank_balanced, unassigned_per_rank, &
                                           balance_id, balance_ref, balance_level, &
                                           sfc_com_list, com_N, num_blocks_count)
        implicit none
        type (type_params), intent(in)      :: params
        integer(kind=ik), intent(in)        :: tree_ID
        integer(kind=ik), intent(in)        :: sfc_sorted_list(:,:)
        integer(kind=ik), intent(in)        :: N_sfc_blocks        !> number of assigned blocks in SFC list
        integer(kind=ik), intent(inout)     :: blocksPerRank_balanced(:)
        integer(kind=ik), intent(in)        :: unassigned_per_rank(0:) !> unassigned blocks per rank
        integer(kind=ik), intent(in)        :: balance_id
        integer(kind=ik), intent(in), optional :: balance_ref(:), balance_level(:)
        integer(kind=ik), intent(inout)     :: sfc_com_list(:,:)
        integer(kind=ik), intent(out)       :: com_N               !> number of communications planned
        integer(kind=ik), intent(inout)     :: num_blocks_count(3) !> debug counters
        
        ! Local variables
        integer(kind=ik) :: k, lgt_id, mpirank, mpisize, ierr, rank_i, mpirank_currently
        integer(kind=ik) :: mpirank_shouldBe, rank_should_have, total_overflow, total_free, i_overflow, i_free, i_proc
        integer(kind=ik) :: overflow_block_counter, target_free_index, free_index_cumulative, target_rank
        integer(kind=ik), allocatable :: blocks_assigned_per_rank(:)   ! how many assigned blocks each rank should get
        integer(kind=ik), allocatable :: overflow_or_free(:)            ! overflow (negative) or free space (positive) per rank
        real(kind=rk) :: ratio  ! proportional distribution factor
        
        mpirank = params%rank
        mpisize = params%number_procs
        
        ! Allocate arrays for MPI communication
        allocate(blocks_assigned_per_rank(0:mpisize-1))
        allocate(overflow_or_free(0:mpisize-1))
        
        blocks_assigned_per_rank = 0
        overflow_or_free = 0  ! Initialize explicitly for safety
        com_N = 0
        
        !-------------------------------------------------------------------------
        ! PHASE 1: Distribute assigned blocks via SFC
        !-------------------------------------------------------------------------
        ! Walk through the SFC-sorted list of assigned blocks and assign them
        ! to ranks sequentially, consuming one slot from blocksPerRank_balanced
        
        mpirank_shouldBe = 0
        do k = 1, N_sfc_blocks
            lgt_id = sfc_sorted_list(1, k)  ! Column 1 = lgt_id, columns 2-3 = treecode
            ! Determine current owner of this block
            call lgt2proc(mpirank, lgt_id, params%number_blocks)
            
            ! Find next rank with available slots
            do while (blocksPerRank_balanced(mpirank_shouldBe+1) <= 0)
                mpirank_shouldBe = mpirank_shouldBe + 1
            end do
            
            ! Assign this block to mpirank_shouldBe
            blocks_assigned_per_rank(mpirank_shouldBe) = blocks_assigned_per_rank(mpirank_shouldBe) + 1
            blocksPerRank_balanced(mpirank_shouldBe+1) = blocksPerRank_balanced(mpirank_shouldBe+1) - 1
            
            ! Plan communication if block needs to move
            if (mpirank /= mpirank_shouldBe) then
                com_N = com_N + 1
                ! store: sender, receiver, lgt_id
                sfc_com_list(1, com_N) = mpirank
                sfc_com_list(2, com_N) = mpirank_shouldBe
                sfc_com_list(3, com_N) = lgt_id
                
                if (params%rank == mpirank) num_blocks_count(1) = num_blocks_count(1) + 1  ! sent
                if (params%rank == mpirank_shouldBe) num_blocks_count(3) = num_blocks_count(3) + 1  ! received
            else
                if (params%rank == mpirank) num_blocks_count(2) = num_blocks_count(2) + 1  ! kept
            end if
        end do
        
        !-------------------------------------------------------------------------
        ! PHASE 2: Compute overflow and free space globally
        !-------------------------------------------------------------------------
        ! After Phase 1, we know:
        ! - blocks_assigned_per_rank(rank): how many assigned blocks each rank will have
        !   (computed identically on all ranks via collective loop over global sfc_sorted_list)
        ! - unassigned_per_rank(rank): how many unassigned blocks each rank has
        !   (global, from earlier MPI_ALLREDUCE in parent routine)
        ! Optimally, all unassigned blocks would stay on their current rank. However, some ranks
        ! may have received too many assigned blocks in Phase 1, resulting in overflow.
        ! We need to move those overflow blocks to ranks with free space.
        ! overflow_or_free = capacity - total = params%number_blocks - (assigned + unassigned)
        !   > 0: free space (can receive blocks)
        !   < 0: overflow (needs to send blocks away)
        
        ! Compute overflow/free space for all ranks (no MPI needed - all ranks have identical data)
        overflow_or_free(0:mpisize-1) = params%number_blocks - (blocks_assigned_per_rank(0:mpisize-1) + unassigned_per_rank(0:mpisize-1))

        ! total overflow is sum of negative entries in overflow_or_free
        total_overflow = 0
        total_free = 0
        do rank_i = 0, mpisize-1
            if (overflow_or_free(rank_i) < 0) total_overflow = total_overflow + abs(overflow_or_free(rank_i))
            if (overflow_or_free(rank_i) > 0) total_free = total_free + overflow_or_free(rank_i)
        end do

        ! Check that total overflow is lower or equal to total free space
        if (total_overflow > total_free) then
            call abort(260112, "balanceLoad_tree: selective mode overflow exceeds free space! This is really weird")
        end if
        
        !-------------------------------------------------------------------------
        ! PHASE 3: Reassign overflow blocks from ranks with overflow to ranks with free space
        !-------------------------------------------------------------------------
        ! Algorithm: Proportional distribution with efficient single-pass crawling
        ! 
        ! 1. Ratio = total_free / total_overflow determines how many free slots per overflow block
        ! 2. Each overflow block i gets target position: int(i * ratio) in virtual free space array
        ! 3. Crawl through ranks to find which rank owns that target position
        ! 4. IMPORTANT: Crawling state (target_rank, free_index_cumulative) is NEVER reset between blocks
        !    - This works because target positions increase monotonically (0, ratio, 2*ratio, ...)
        !    - Results in O(n) complexity instead of O(n²)
        ! 5. Blocks are processed in lgt_active order, making distribution deterministic across all ranks
        
        if (total_overflow == 0) then
            ! No overflow, all blocks can stay on their current ranks
            ! Just count kept blocks for debug output
            do k = 1, lgt_n(tree_ID)
                lgt_id = lgt_active(k, tree_ID)
                if (.not. block_matches_filter(lgt_id, balance_id, balance_ref, balance_level)) then
                    call lgt2proc(mpirank_currently, lgt_id, params%number_blocks)
                    if (params%rank == mpirank_currently) num_blocks_count(2) = num_blocks_count(2) + 1  ! kept
                endif
            enddo
        else
            ! Overflow exists, need to redistribute
            ! Initialize crawling state (these variables are NEVER reset during the loop)
            ratio = real(total_free, kind=rk) / real(total_overflow, kind=rk)  ! proportionality factor
            overflow_block_counter = 0         ! how many overflow blocks we've placed so far (0-indexed)
            free_index_cumulative = 0          ! accumulated free space as we crawl through ranks
            target_rank = 0                    ! current rank we're examining (increments only, never resets)
            i_overflow = 0                     ! overflow blocks remaining on current rank (initialized to 0)
            i_proc = -1                        ! track current rank in block loop (initialized to invalid rank)
            
            ! Loop over all blocks in lgt_active order to find unassigned blocks on overflowing ranks
            do k = 1, lgt_n(tree_ID)
                lgt_id = lgt_active(k, tree_ID)
                call lgt2proc(mpirank_currently, lgt_id, params%number_blocks)
                
                ! Track when we encounter a new rank
                if (mpirank_currently /= i_proc) then
                    i_proc = mpirank_currently
                    i_overflow = max(0, -overflow_or_free(i_proc))  ! how many overflow blocks this rank needs to send
                end if
                
                ! Check if this is an unassigned block
                if (.not. block_matches_filter(lgt_id, balance_id, balance_ref, balance_level)) then
                    if (i_overflow > 0) then
                        ! This unassigned block needs to be evacuated from this overflowing rank
                        
                        ! Calculate target position in virtual free space array
                        target_free_index = int(real(overflow_block_counter, kind=rk) * ratio)
                        
                        ! Crawl forward through ranks to find which rank owns this position
                        do while (target_rank < mpisize)
                            if (overflow_or_free(target_rank) > 0) then
                                ! This rank has free space
                                if (free_index_cumulative + overflow_or_free(target_rank) > target_free_index) then
                                    ! Found it! This block goes to target_rank
                                    i_free = target_rank
                                    exit
                                else
                                    ! Not yet, accumulate and continue
                                    free_index_cumulative = free_index_cumulative + overflow_or_free(target_rank)
                                    target_rank = target_rank + 1
                                endif
                            else
                                ! No free space on this rank, skip it
                                target_rank = target_rank + 1
                            endif
                        end do
                        
                        ! Safety check
                        if (target_rank >= mpisize) then
                            call abort(260114, "balanceLoad_tree: selective mode crawling exceeded rank bounds!")
                        endif
                        
                        ! Plan communication: send from mpirank_currently to i_free
                        com_N = com_N + 1
                        sfc_com_list(1, com_N) = mpirank_currently
                        sfc_com_list(2, com_N) = i_free
                        sfc_com_list(3, com_N) = lgt_id
                        
                        ! Update counters
                        overflow_block_counter = overflow_block_counter + 1  ! next overflow block
                        i_overflow = i_overflow - 1                          ! one less to evacuate from this rank
                        
                        ! Debug output: track send/receive counts for this rank only
                        if (params%rank == mpirank_currently) num_blocks_count(1) = num_blocks_count(1) + 1  ! sent from this rank
                        if (params%rank == i_free) num_blocks_count(3) = num_blocks_count(3) + 1  ! received by this rank
                    else
                        ! This unassigned block can stay on this rank (no overflow)
                        if (params%rank == mpirank_currently) num_blocks_count(2) = num_blocks_count(2) + 1  ! kept
                    end if
                endif
            enddo
        endif
        
        deallocate(blocks_assigned_per_rank, overflow_or_free)
        
    end subroutine plan_communication_selective


    !> \brief Plan communication for standard/refine_predictive balancing modes
    !> \details Walk through SFC and assign contiguous segments to ranks sequentially
    subroutine plan_communication_standard(params, sfc_sorted_list, N_sfc_blocks, blocksPerRank_balanced, &
                                          balance_id, sfc_com_list, com_N, num_blocks_count)
        implicit none
        type (type_params), intent(in)      :: params
        integer(kind=ik), intent(in)        :: sfc_sorted_list(:,:)
        integer(kind=ik), intent(in)        :: N_sfc_blocks        !> number of blocks in SFC list
        integer(kind=ik), intent(inout)     :: blocksPerRank_balanced(:)
        integer(kind=ik), intent(in)        :: balance_id
        integer(kind=ik), intent(inout)     :: sfc_com_list(:,:)
        integer(kind=ik), intent(out)       :: com_N               !> number of communications planned
        integer(kind=ik), intent(inout)     :: num_blocks_count(3) !> debug counters
        
        integer(kind=ik) :: k, mpirank_shouldBe, mpirank_currently, lgt_id
        
        ! mpirank_shouldBe: rank responsible for current part of SFC
        ! mpirank_currently: rank currently holding the block
        
        ! we start the loop on the root rank (0), then assign the first elements
        ! of the SFC, then to second rank, etc. (thus: mpirank_shouldBe is a loop variable)
        mpirank_shouldBe = 0
        
        ! communication counter. each communication (=send and receive) is stored in a long list
        com_N = 0
        
        ! prepare lists for transfering of blocks
        ! ATTENTION: this loop is unusual because it does not loop over the active lists, but the sfc_sorted_list
        ! hence you find no lgt_active(k, tree_ID) here.
        ! COLLECTIVE OPERATION
        do k = 1, N_sfc_blocks
            ! Determine which mpirank should own the current portion of the SFC.
            !
            ! if the current owner of the SFC is supposed to have zero blocks
            ! then it does not really own this part of the SFC. So we look for the
            ! first rank which is supposed to hold at least one block, and declare it as owner
            ! of this part. NOTE: as we try to minimize communication during send/recv in
            ! load balancing, it may well be that the list of active mpiranks (ie those
            ! which have nonzero number of blocks) is non contiguous, i.e.
            ! blocksPerRank_balanced = 1 1 1 0 0 0 0 1 0 1
            ! can happen.
            do while ( blocksPerRank_balanced(mpirank_shouldBe+1) <= 0 )
                ! why <= 0 ? Because in refine_predictive mode, we may end up with a slight imbalance
                ! as we distribute chunks of 2**D or 1 future block - the desired number of blocks is not automatically
                ! divisible by this number
                mpirank_shouldBe = mpirank_shouldBe + 1
            end do
            
#ifdef DEV
            ! Should never happen...
            if (mpirank_shouldBe>params%number_procs-1) then
                call abort(262981,"Very unusual, we still have blocks leftover but no more mpiranks for distribution")
            endif
#endif
            
            ! find out on which mpirank lies the block that we're looking at
            call lgt2proc( mpirank_currently, sfc_sorted_list(1, k), params%number_blocks )
            
            ! does this block lie on the right mpirank, i.e., the current part of the
            ! SFC? if so, nothing needs to be done. otherwise, the following if is active
            if ( mpirank_shouldBe /= mpirank_currently ) then
                ! as this block is one the wrong rank, it will be sent away from its
                ! current owner (mpirank_currently) to the owner of this part of the
                ! SFC (mpirank_shouldBe)
                
                ! save this send+receive operation in the list of planned communications
                ! column
                !    1     sender mpirank
                !    2     receiver mpirank
                !    3     block light data id
                com_N = com_N + 1
                sfc_com_list(1, com_N) = mpirank_currently      ! sender mpirank
                sfc_com_list(2, com_N) = mpirank_shouldBe       ! receiver mpirank
                sfc_com_list(3, com_N) = sfc_sorted_list(1,k)   ! light id of block
            end if
            
            if (params%debug_balanceLoad) then
                ! Development checks - count how many blocks are send, received or kept
                if ( params%rank == mpirank_currently .and. mpirank_currently /= mpirank_shouldBe ) then
                    ! block is send away from this mpirank
                    num_blocks_count(1) = num_blocks_count(1) + 1
                elseif ( params%rank == mpirank_currently .and. mpirank_currently == mpirank_shouldBe ) then
                    ! block is kept on this mpirank
                    num_blocks_count(2) = num_blocks_count(2) + 1
                elseif ( params%rank == mpirank_shouldBe .and. mpirank_currently /= mpirank_shouldBe ) then
                    ! block is received by this mpirank
                    num_blocks_count(3) = num_blocks_count(3) + 1
                endif
            endif
            
            ! The blocksPerRank_balanced defines how many blocks this rank should have, and
            ! we just treated one (which either already was on the mpirank or will be on
            ! it after communication), so remove one item from the blocksPerRank_balanced
            lgt_id = sfc_sorted_list(1,k)
            if ((lgt_block(lgt_id , IDX_REFINE_STS) == +1).and.(balance_id == 2)) then
                ! block will be refined, resulting in 2**D blocks (refine_predictive mode)
                blocksPerRank_balanced( mpirank_shouldBe+1 ) = blocksPerRank_balanced( mpirank_shouldBe+1 ) - 2**params%dim
            else
                ! either standard balancing or block will not be refined, resulting in 1 block in both cases
                blocksPerRank_balanced( mpirank_shouldBe+1 ) = blocksPerRank_balanced( mpirank_shouldBe+1 ) - 1
            endif
            
        end do
        
    end subroutine plan_communication_standard

    !> \brief Helper function to check if a block matches the filter criteria
    !> \details Checks refinement status and/or level filters based on balance_id encoding
    !> Following pattern from GS_iteration_ref for array input handling
    logical function block_matches_filter(lgt_id, balance_id, balance_ref, balance_level)
        implicit none
        integer(kind=ik), intent(in) :: lgt_id, balance_id
        integer(kind=ik), intent(in), optional :: balance_ref(:), balance_level(:)
        
        integer(kind=ik) :: ref, level
        
        ! Default: no filter, all blocks match
        block_matches_filter = .true.
        
        ! Check filter type (second digit of balance_id)
        select case(balance_id / 10)
        case(0)
            ! No filter - all blocks match
            block_matches_filter = .true.
            
        case(1)
            ! Refinement status filter only
            ref = lgt_block(lgt_id, IDX_REFINE_STS)
            ! Use "any" to check if ref matches any value in the array (like GS_iteration_ref uses "all")
            block_matches_filter = any(ref == balance_ref)
            
        case(2)
            ! Level filter only
            level = lgt_block(lgt_id, IDX_MESH_LVL)
            ! Use "any" to check if level matches any value in the array
            block_matches_filter = any(level == balance_level)
            
        case(3)
            ! Both filters (AND logic)
            ref = lgt_block(lgt_id, IDX_REFINE_STS)
            level = lgt_block(lgt_id, IDX_MESH_LVL)
            block_matches_filter = any(ref == balance_ref) .and. any(level == balance_level)
        end select
    end function block_matches_filter

end subroutine balanceLoad_tree
