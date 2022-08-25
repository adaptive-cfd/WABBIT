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
    use module_precision
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

    integer(hsize_t), dimension(4)          :: size_field
    integer(hid_t)                          :: file_id
    integer(kind=ik)                        :: lgt_id, k, nz, iteration, tc_length, dim, level
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), dimension(3)             :: x0, dx
    real(kind=rk), dimension(3)             :: domain
    real(kind=rk), dimension(4)             :: xi
    real(kind=rk)                           :: time, x_query
    real(kind=rk), allocatable              :: hvy_block_2Dslice(:, :, :, :, :)
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik), allocatable           :: tree(:), sum_tree(:), blocks_per_rank(:), treecode(:)
    integer(kind=ik), allocatable           :: lgt_block_2Dslice(:,:)
    integer(kind=ik), allocatable           :: lgt_active_2Dslice(:,:), hvy_active_2Dslice(:,:)

    integer(kind=ik)  :: Nblocks, g
    real(kind=rk)    :: x,y,z, x_query_normalized
    real(kind=rk)    :: maxi,mini,squari,meani,qi
    real(kind=rk)    :: maxl,minl,squarl,meanl,ql
    integer(kind=ik) :: ix,iy,iz,mpicode, ioerr, i

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_cmd_arg("--file", fname, "not given")
    call get_cmd_arg("--output", fname_out, "slice_00.h5")
    call get_cmd_arg("--x", x_query, 0.0_rk)

    call check_file_exists( fname )

    if (params%number_procs > 1) then
        ! it shouldn't be difficult, buut I did not do it. problem: active lists for 2D slice
        ! are no longer simple vectors 1:N
        call abort(21060301, "This routine cannot be run in parallel right now.")
    endif

    ! get some parameters from the file
    call read_attributes(fname, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%dim = dim
    params%Bs = Bs
    params%max_treelevel = tc_length
    params%n_eqn = 1
    ! for interpolation of the slice (in the x-direction) we need at most 2 ghost nodes.
    g = 2
    params%n_ghosts = g
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%number_blocks = lgt_n(tree_ID)
    params%block_distribution = "sfc_hilbert"
    params%order_predictor = "multiresolution_4th"
    allocate(params%symmetry_vector_component(1:3))
    params%symmetry_vector_component(1:3) = "0"
    ! hack to save on memory
    N_MAX_COMPONENTS = 1

    allocate(treecode(1:tc_length))

    call allocate_forest(params, hvy_block)

    ! read data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID)

    call updateMetadata_tree(params, tree_ID)

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    Nblocks = 0
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if ((x0(1) <= x_query) .and. (x_query < x0(1)+real(Bs(1)-1,kind=rk)*dx(1))) then
            Nblocks = Nblocks +1
        endif
    end do

    write(*,*) "out of the total ", lgt_n(tree_ID)," blocks, ", Nblocks, " are concerned by slicing"

    allocate( hvy_block_2Dslice(1:Bs(2)+2*g, 1:Bs(3)+2*g, 1, 1, 1:Nblocks) )
    allocate( lgt_block_2Dslice(1:Nblocks, 1:tc_length+EXTRA_LGT_FIELDS))
    allocate( lgt_active_2Dslice(1:Nblocks, 1), hvy_active_2Dslice(1:Nblocks, 1))

    hvy_block_2Dslice = 0.0_rk
    lgt_block_2Dslice = -1

    i = 1
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if ((x0(1) <= x_query) .and. (x_query < x0(1)+real(Bs(1)-1,kind=rk)*dx(1))) then
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
            ix = floor( (x_query-x0(1))/dx(1) ) + g+1

            ! coordinate in the local coordinate system (above), where the point IX corresonds
            ! to 1.
            x_query_normalized = 1 + ( x_query - x0(1)+real(ix-g-1, kind=rk)*dx(1) )

            if ((ix<g+1) .or. (ix>Bs(1)+g)) then
                write(*,*) ix
                write(*,*) x0, ":", dx, ":", x_query
                call abort(716717, "well")
            endif

            xi = (/0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk/)

            ! interpolate the slice along the x-direction
            hvy_block_2Dslice(:, :, 1, 1, i) = &
            +lagrange_polynomial(x_query_normalized, xi, 1)*hvy_block( ix-1, :, :, 1, hvy_id) &
            +lagrange_polynomial(x_query_normalized, xi, 2)*hvy_block( ix  , :, :, 1, hvy_id) &
            +lagrange_polynomial(x_query_normalized, xi, 3)*hvy_block( ix+1, :, :, 1, hvy_id) &
            +lagrange_polynomial(x_query_normalized, xi, 4)*hvy_block( ix+2, :, :, 1, hvy_id)

            ! this is just the inverse of what get_block_spacing_origin2 computes
            ! we use it to create the 2D treecode
            iy = 1 + nint( x0(2) / (dx(2)*real(Bs(2)-1, kind=rk)) )
            iz = 1 + nint( x0(3) / (dx(3)*real(Bs(3)-1, kind=rk)) )

            level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)

            ! note 2D data is (x,y) but our slice is (y,z) hence the oddity (iy,iz,1)
            ! NOTE: in paraview I saw that I had to invert iy,iz to iz,iy although I do not
            ! copletely understand why.
            call encoding_revised(treecode, (/iz,iy/), 2, level)

            ! copy computed treecode and som other information to lgt_block for the 2D slice
            lgt_block_2Dslice(i, 1:level) = treecode(1:level)
            lgt_block_2Dslice(i, tc_length+IDX_MESH_LVL)   = level
            lgt_block_2Dslice(i, tc_length+IDX_TREE_ID)    = 1
            lgt_block_2Dslice(i, tc_length+IDX_REFINE_STS) = 0

            i = i + 1
        endif
    end do

    ! save 2D slice to disk, change code to now expect 2D data
    params%dim = 2
    params%Bs = (/Bs(2), Bs(3), 1/)
    params%domain_size = (/domain(2), domain(3), 0.0_rk/)

    do i = 1, Nblocks
        lgt_active_2Dslice(i, tree_ID) = i
        hvy_active_2Dslice(i, tree_ID) = i
    enddo


    ! -------------------
    deallocate(lgt_block, lgt_active, hvy_active)
    allocate(lgt_block(size(lgt_block_2Dslice,1), size(lgt_block_2Dslice,2)))
    allocate(lgt_active(size(lgt_active_2Dslice,1), size(lgt_active_2Dslice,2)))
    allocate(hvy_active(size(hvy_active_2Dslice,1), size(hvy_active_2Dslice,2)))

    lgt_block  = lgt_block_2Dslice
    lgt_active = lgt_active_2Dslice
    hvy_active = hvy_active_2Dslice
    hvy_n(tree_ID) = Nblocks
    lgt_n(tree_ID) = Nblocks
    ! -------------------

    call saveHDF5_tree( fname_out, time, iteration, 1, params, hvy_block_2Dslice, tree_ID)

end subroutine
