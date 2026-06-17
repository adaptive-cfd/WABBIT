subroutine post_meta_analyse(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)      :: file, operator
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, tc_length, hvy_id
    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik)       :: blocks_with_all_sisters, blocks_with_incomplete_sisters
    integer(kind=ik), allocatable :: blocks_per_level(:)
    integer(kind=ik)       :: level, Jmax

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1

    real(kind=rk), dimension(3)        :: domain

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_command_argument(1, operator)
    call get_command_argument(2, file)

    ! does the user need help?
    if (file=='--help' .or. file=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: meta analysis tools"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Read in a data field (2D or 3D) and performs various"
            write(*,*) " meta-analyses on the grid structure."
            write(*,*) " "
            write(*,*) " Available operators:"
            write(*,*) " "
            write(*,*) " --analyse-sisters"
            write(*,*) "   Checks if all sister blocks exist for every block."
            write(*,*) "   This is done by checking hvy_family for each block."
            write(*,*) "   If any sister is missing (entry -1), the block gets"
            write(*,*) "   status 1. If all sisters exist, the block gets status -1."
            write(*,*) " "
            write(*,*) " --analyse-levels"
            write(*,*) "   Counts the number of blocks on each refinement level"
            write(*,*) "   and prints a histogram of the block distribution."
            write(*,*) " "
            write(*,*) " Usage: "
            write(*,*) " ./wabbit-post --analyse-sisters file.h5"
            write(*,*) " ./wabbit-post --analyse-levels file.h5"
        end if
        return
    endif

    call check_file_exists(trim(file))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    ! unused so just fill any value
    params%order_discretization = "FD_2nd_central"
    params%order_predictor = "multiresolution_2nd"
    params%wavelet = "CDF20"
    params%g = 2_ik
    params%Jmax = tc_length+1
    params%n_eqn = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs

    allocate(params%butcher_tableau(1,1))
    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )

    ! init wavelet
    call setup_wavelet(params)

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

    ! read mesh and field
    call readHDF5vct_tree((/file/), params, hvy_block, tree_ID, synchronize_ghosts=.false.)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted numerical treecodes, used for finding blocks
    call updateMetadata_tree( params, tree_ID )

    ! ------------------------------------------------
    ! Perform the requested analysis
    ! ------------------------------------------------
    if (operator == "--analyse-sisters") then
        ! ------------------------------------------------
        ! Check for sister completeness
        ! ------------------------------------------------
        ! Loop over all heavy blocks and check if all sisters exist
        ! by examining hvy_family(hvy_id, 2:1+2**params%dim)
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            
            ! Check if any sister is missing (value -1 in hvy_family means non-existent)
            if ( any(hvy_family(hvy_id, 2:1+2**params%dim) == -1) ) then
                ! At least one sister is missing
                lgt_block(lgt_id, IDX_REFINE_STS) = 1
            else
                ! All sisters exist
                lgt_block(lgt_id, IDX_REFINE_STS) = -1
            endif
        enddo

        ! Synchronize refinement status across all MPI ranks
        call synchronize_lgt_data( params, refinement_status_only=.true. )

        ! Count blocks with all sisters (status -1) and without (status 1 or other)
        blocks_with_all_sisters = 0
        blocks_with_incomplete_sisters = 0

        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)
            if (lgt_block(lgt_id, IDX_REFINE_STS) == -1) then
                blocks_with_all_sisters = blocks_with_all_sisters + 1
            else
                blocks_with_incomplete_sisters = blocks_with_incomplete_sisters + 1
            endif
        enddo

        ! Output results
        if (params%rank == 0) then
            write(*,*) "============================================"
            write(*,*) " Sister Completeness Check Results"
            write(*,*) "============================================"
            write(*,'(A,I0)') " Total blocks: ", blocks_with_all_sisters + blocks_with_incomplete_sisters
            write(*,'(A,I0)') " Blocks with all sisters: ", blocks_with_all_sisters
            write(*,'(A,I0)') " Blocks without all sisters: ", blocks_with_incomplete_sisters
            write(*,'(A,F7.2,A)') " Percentage with all sisters: ", &
                100.0_rk * dble(blocks_with_all_sisters) / dble(blocks_with_all_sisters + blocks_with_incomplete_sisters), "%"
            write(*,*) "============================================"
        endif

    elseif (operator == "--analyse-levels") then
        ! ------------------------------------------------
        ! Analyse block distribution across levels
        ! ------------------------------------------------
        Jmax = params%Jmax
        allocate(blocks_per_level(0:Jmax))
        blocks_per_level = 0

        ! Count blocks on each level
        do k = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k, tree_ID)
            level = lgt_block(lgt_id, IDX_MESH_LVL)
            blocks_per_level(level) = blocks_per_level(level) + 1
        enddo

        ! Output results
        if (params%rank == 0) then
            write(*,*) "============================================"
            write(*,*) " Block Distribution by Refinement Level"
            write(*,*) "============================================"
            write(*,'(A,I0)') " Total blocks: ", lgt_n(tree_ID)
            write(*,'(A,I0)') " Maximum level (Jmax): ", Jmax
            write(*,*) "--------------------------------------------"
            write(*,'(A)') " Level  |  Number of Blocks  |  Percentage"
            write(*,*) "--------------------------------------------"
            do level = 0, Jmax
                if (blocks_per_level(level) > 0) then
                    write(*,'(I5,A,I12,A,F10.2,A)') level, "   |   ", blocks_per_level(level), &
                        "     |  ", 100.0_rk * dble(blocks_per_level(level)) / dble(lgt_n(tree_ID)), "%"
                endif
            enddo
            write(*,*) "============================================"
        endif

        deallocate(blocks_per_level)

    else
        call abort(202102161, "Unknown operator: "//trim(operator)//". Use --help for available options.")
    endif

end subroutine post_meta_analyse
