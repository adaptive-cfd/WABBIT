! extract 2D slice from 3D data (currently only x=const, i.e. y-z plane)
!
! Algorithm:
!
! load data
! for all blocks:
!       check if the block contains the x_quey position
!       if yes:
!           interpolate the data along the x-axis to match x_query
!           save this interpolated slice

subroutine post_extract_slice(params)
    use module_globals
    use module_params
    use module_mesh
    use module_helpers
    use mpi
    use module_forestMetaData

    implicit none

    character(len=cshort)                       :: fname, fname_out
    type (type_params), intent(inout)       :: params

    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
    integer(kind=ik)                        :: tree_ID=1, hvy_id

    integer(kind=ik)                        :: lgt_id, k, iteration, tc_length, dim, level
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), dimension(3)             :: x0, dx
    real(kind=rk), dimension(3)             :: domain
    real(kind=rk), dimension(4)             :: xi
    real(kind=rk)                           :: time, cut_pos
    real(kind=rk), allocatable              :: hvy_block_2Dslice(:, :, :, :, :)
    integer(kind=tsize)                     :: tc_3D, tc_2D
    integer(kind=ik), allocatable           :: lgt_block_2Dslice(:,:)
    integer(kind=ik), allocatable           :: lgt_active_2Dslice(:,:), hvy_active_2Dslice(:,:)

    integer(kind=ik)  :: Nblocks, Nblocks_total, Nblocks_max, g, cut_dim
    real(kind=rk)    :: cut_pos_normalized
    integer(kind=ik) :: icut, ixyz(1:3),mpicode, i_hvy, i_lgt
    logical          :: help1, help2, help3

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_cmd_arg( "--help", help1, default=.false. )
    call get_cmd_arg( "--h", help2, default=.false. )
    call get_cmd_arg( "-h", help3, default=.false. )

    ! does the user need help?
    if (help1 .or. help2 .or. help3 .and. params%rank==0) then
        if (params%rank==0) then
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " Wabbit extract slice: extract a 2D slice from a 3D dataset"
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " ./wabbit-post --extract-slice --file=FILE.h5 --output=SLICE.h5 --dim=1 --pos=0.5 --wavelet=CDF20"
            write(*,'(A)') ""
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " --file                - input file (3D data)"
            write(*,'(A)') " --output=slice_00.h5  - output file (2D slice)"
            write(*,'(A)') " --dim                 - dimension along which the slice is extracted (1:x, 2:y, 3:z)"
            write(*,'(A)') " --pos                 - position of the slice along the specified dimension"
            write(*,'(A)') " --wavelet=CDF20       - wavelet used for synching the data"
            write(*,'(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    call get_cmd_arg("--file", fname, "none")
    call get_cmd_arg("--output", fname_out, "slice_00.h5")
    call get_cmd_arg("--dim", cut_dim, -1)
    call get_cmd_arg("--pos", cut_pos, -1.0_rk)
    if (cut_dim<1 .or. cut_dim>3) call abort(260707, "--dim must be 1 (x), 2 (y) or 3 (z)")
    if (cut_pos<0.0_rk) call abort(260707, "--pos must be >= 0.0")

    call check_file_exists( fname )

    ! get some parameters from the file
    call read_attributes(fname, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%dim = dim
    params%Bs = Bs
    params%Jmax = tc_length
    params%n_eqn = 1
    ! for interpolation of the slice (in the x-direction) we need at most 2 ghost nodes.
    g = 2
    params%g = g
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )
    params%block_distribution = "sfc_hilbert"
    params%order_predictor = "multiresolution_4th"
    allocate(params%symmetry_vector_component(1:3))
    params%symmetry_vector_component(1:3) = "0"
    params%wavelet = "CDF20"
    call setup_wavelet(params)

    call allocate_forest(params, hvy_block)

    ! read data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID)

    call updateMetadata_tree(params, tree_ID)

    call sync_ghosts_tree( params, hvy_block, tree_ID )

    Nblocks = 0
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if ((x0(cut_dim) <= cut_pos) .and. (cut_pos < x0(cut_dim)+real(Bs(cut_dim)-1,kind=rk)*dx(cut_dim))) then
            Nblocks = Nblocks +1
        endif
    end do
    ! we need to allreaduce Nblocks inplace to get how many blocks will really be sliced in total, and how many at maximum on one rank
    call MPI_Allreduce(Nblocks, Nblocks_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpicode)
    call MPI_Allreduce(Nblocks, Nblocks_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpicode)
    if (params%rank==0) write(*,'(A,I0,A,I0,A,I0,A)') "out of the total ", lgt_n(tree_ID)," blocks, ", Nblocks_total, " are concerned by slicing, with at max ", Nblocks_max, " on one rank."

    allocate( hvy_block_2Dslice(1:Bs(2)+2*g, 1:Bs(3)+2*g, 1, 1, 1:Nblocks_max) )
    allocate( lgt_block_2Dslice(1:Nblocks_max*params%number_procs, 1:EXTRA_LGT_FIELDS))
    allocate( lgt_active_2Dslice(1:Nblocks_max*params%number_procs, 1), hvy_active_2Dslice(1:Nblocks_max*params%number_procs, 1))

    hvy_block_2Dslice = 0.0_rk
    lgt_block_2Dslice = -1

    i_hvy = 1  ! this is the index where we will write the next block in the 2D slice
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        tc_3D = get_tc(lgt_block(lgt_id, IDX_TC_1:IDX_TC_2))
        level = lgt_block(lgt_id, IDX_MESH_LVL)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if ((x0(cut_dim) <= cut_pos) .and. (cut_pos < x0(cut_dim)+real(Bs(cut_dim)-1,kind=rk)*dx(cut_dim))) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! floor index: in our notation, this is point 1 (the second point)
            ! where: | is location of interpolation
            ! -----o-----o-|---o-----o----- (points)
            !
            ! -----1-----2-----3-----4------ (indiced in fortran, one-based indexing [the last argument to lagrange_polynomial])
            !
            ! -----0-----1-----2-----3----- (x_normalized)
            !            ^
            !           ix
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            icut = floor( (cut_pos-x0(cut_dim))/dx(cut_dim) ) + g+1

            ! coordinate in the local coordinate system (above), where the point IX corresonds
            ! to 1.
            cut_pos_normalized = 1 + ( cut_pos - x0(cut_dim)+real(icut-g-1, kind=rk)*dx(cut_dim) )

            if ((icut<g+1) .or. (icut>Bs(cut_dim)+g)) then
                write(*,*) icut
                write(*,*) x0, ":", dx, ":", cut_pos
                call abort(716717, "well")
            endif

            xi = (/0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk/)

            ! interpolate the slice along the x-direction
            if (cut_dim==1) then
                hvy_block_2Dslice(:, :, 1, 1, i_hvy) = &
                    +lagrange_polynomial(cut_pos_normalized, xi, 1)*hvy_block( icut-1, :, :, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 2)*hvy_block( icut  , :, :, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 3)*hvy_block( icut+1, :, :, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 4)*hvy_block( icut+2, :, :, 1, hvy_id)
            elseif (cut_dim==2) then
                hvy_block_2Dslice(:, :, 1, 1, i_hvy) = &
                    +lagrange_polynomial(cut_pos_normalized, xi, 1)*hvy_block( :, icut-1, :, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 2)*hvy_block( :, icut  , :, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 3)*hvy_block( :, icut+1, :, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 4)*hvy_block( :, icut+2, :, 1, hvy_id)
            elseif (cut_dim==3) then
                hvy_block_2Dslice(:, :, 1, 1, i_hvy) = &
                    +lagrange_polynomial(cut_pos_normalized, xi, 1)*hvy_block( :, :, icut-1, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 2)*hvy_block( :, :, icut  , 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 3)*hvy_block( :, :, icut+1, 1, hvy_id) &
                    +lagrange_polynomial(cut_pos_normalized, xi, 4)*hvy_block( :, :, icut+2, 1, hvy_id)
            endif

            ! from the treecode, we compute the ix,iy,iz of the 2D slice
            call decoding_b(ixyz, tc_3D, dim=params%dim, level=level, max_level=tc_length)

            ! note 2D data is (x,y) but our slice is (y,z) hence the oddity (iy,iz,1)
            ! NOTE: in paraview I saw that I had to invert iy,iz to iz,iy although I do not
            ! copletely understand why.
            if (cut_dim==1) then
                call encoding_b((/ixyz(2),ixyz(3)/), tc_2D, dim=2, level=level, max_level=tc_length)
            elseif (cut_dim==2) then
                call encoding_b((/ixyz(1),ixyz(3)/), tc_2D, dim=2, level=level, max_level=tc_length)
            elseif (cut_dim==3) then
                call encoding_b((/ixyz(1),ixyz(2)/), tc_2D, dim=2, level=level, max_level=tc_length)
            endif

            ! for writing in the 2D slice, we need to compute the new lgt_id
            call hvy2lgt(i_lgt, i_hvy, params%rank, Nblocks_max)

            ! copy computed treecode and som other information to lgt_block for the 2D slice
            lgt_block_2Dslice(i_lgt, :)   = -1
            ! lgt_block_2Dslice(i_lgt, 1:level) = treecode(1:level)
            lgt_block_2Dslice(i_lgt, IDX_MESH_LVL)   = level
            lgt_block_2Dslice(i_lgt, IDX_TREE_ID)    = 1
            lgt_block_2Dslice(i_lgt, IDX_REFINE_STS) = 0
            call set_tc(lgt_block_2Dslice( i_lgt, IDX_TC_1:IDX_TC_2), tc_2D)
            i_hvy = i_hvy + 1
        endif
    end do

    ! save 2D slice to disk, change code to now expect 2D data
    params%dim = 2
    ! ToDo: Julius - tbh I am not sure if this is correct for data where Bs(1) != Bs(2) != Bs(3)
    if (cut_dim==1) then
        params%Bs = (/Bs(2), Bs(3), 1/)
        params%domain_size = (/domain(2), domain(3), 0.0_rk/)
    elseif (cut_dim==2) then
        params%Bs = (/Bs(1), Bs(3), 1/)
        params%domain_size = (/domain(1), domain(3), 0.0_rk/)
    elseif (cut_dim==3) then
        params%Bs = (/Bs(1), Bs(2), 1/)
        params%domain_size = (/domain(1), domain(2), 0.0_rk/)
    endif

    ! -------------------
    ! reinitialize all arrays with our 2D data
    deallocate(lgt_block, lgt_active, hvy_active)
    allocate(lgt_block(size(lgt_block_2Dslice,1), size(lgt_block_2Dslice,2)))
    allocate(lgt_active(size(lgt_active_2Dslice,1), size(lgt_active_2Dslice,2)))
    allocate(hvy_active(size(hvy_active_2Dslice,1), size(hvy_active_2Dslice,2)))

    ! we only need to set lgt_block and the amount of blocks, the rest will be updated later on
    lgt_block  = lgt_block_2Dslice
    params%number_blocks = Nblocks_max
    ! -------------------
    ! now correct meta information and everything
    call synchronize_lgt_data( params, refinement_status_only=.false. )  ! this syncs the lgt_block between the processors
    call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)  ! this updates the active lists and neighbor/family relations

    call saveHDF5_tree( fname_out, time, iteration, 1, params, hvy_block_2Dslice, tree_ID)

end subroutine
