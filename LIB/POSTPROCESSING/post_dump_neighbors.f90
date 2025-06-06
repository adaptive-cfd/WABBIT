
subroutine post_dump_neighbors(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)      :: file, operator
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, tc_length
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id

    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_command_argument(1, operator)
    call get_command_argument(2, file)

    ! does the user need help?
    if (file=='--help' .or. file=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: dump neighbors"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Read in a data field (2D or 3D) and computes the neighbor relations"
            write(*,*) " (updateNeighbors_tree) for all blocks. then, dumps hvy_neighbor to "
            write(*,*) " ascii file. Useful for development"
            write(*,*) " Best used in serial monoproc: then, hvy_id == lgt_id"
        end if
        return
    endif

    call check_file_exists(trim(file))

    if (params%number_procs /= 1) then
        call abort( 202101281, "use on monoproc")
    endif


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
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call updateMetadata_tree( params, tree_ID )

    ! write neighbor_blocks file
    if (params%rank == 0) then
        open(16,file='hvy_neighbors.dat',status='replace')
        write(16,'(A)') "% lgt_id + 168*neighbors"
        do k=1,hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)

            write(16,'(169(i0, 1x))') lgt_id, hvy_neighbor(k, :)
        enddo
        close(16)
    endif
end subroutine post_dump_neighbors
