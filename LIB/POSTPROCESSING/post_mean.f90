subroutine post_mean(params)
    use module_precision
    use module_params
    use module_mesh
    use mpi
    use module_forestMetaData

    implicit none
    character(len=cshort)                   :: fname, fname_out                 !> name of the file
    type (type_params), intent(inout)       :: params                           !> parameter struct

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id

    integer(hsize_t), dimension(4)          :: size_field
    integer(hid_t)                          :: file_id
    integer(kind=ik)                        :: lgt_id, k, nz, iteration, dim, g
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), dimension(3)             :: x0, dx
    real(kind=rk), dimension(3)             :: domain
    real(kind=rk)                           :: time
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik), allocatable           :: tree(:), sum_tree(:), blocks_per_rank(:)

    real(kind=rk)    :: x,y,z
    real(kind=rk)    :: maxi,mini,squari,meani,qi,inti
    real(kind=rk)    :: maxl,minl,squarl,meanl,ql
    integer(kind=ik) :: ix,iy,iz,mpicode, ioerr, rank, i, tc_length

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    !-----------------------------------------------------------------------------------------------------
    rank = params%rank
    !-----------------------------------------------------------------------------------------------------
    call get_command_argument(2,fname)
    if (fname =='--help' .or. fname=='--h') then
        if (rank==0) then
            write(*,*) "WABBIT postprocessing: compute mean value of field. We also output the volume integral"
            write(*,*) "(Lx*Ly*Lz*mean(field) to the file results.txt"
            write(*,*) "mpi_command -n number_procs ./wabbit-post --mean filename.h5 result.txt"
        end if
        return
    endif

    if (rank==0) write (*,*) "Computing spatial mean of file: "//trim(adjustl(fname))
    call check_file_exists( fname )

    ! add some parameters from the file
    call read_attributes(fname, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%Bs = Bs
    params%n_eqn = 1
    params%n_ghosts = 2_ik
    params%order_predictor = "multiresolution_2nd"
    params%max_treelevel = tc_length
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%number_blocks = ceiling( real(lgt_n(tree_ID)) / real(params%number_procs) )

    call allocate_forest(params, hvy_block)

    ! read input data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call updateMetadata_tree( params, tree_ID )

    ! compute an additional quantity that depends also on the position
    ! (the others are translation invariant)
    if (params%dim == 3) then
        nz = Bs(3)
    else
        nz = 1
    end if

    meanl = 0.0_rk

    g = params%n_ghosts
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if (params%dim == 3) then
            meanl = meanl + sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 1, hvy_id))*dx(1)*dx(2)*dx(3)
        else
            meanl = meanl + sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, 1, hvy_id))*dx(1)*dx(2)
        endif
    end do

    call MPI_REDUCE(meanl,meani,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,WABBIT_COMM,mpicode)

    if (params%dim == 3) then
        meani = meani / (params%domain_size(1)*params%domain_size(2)*params%domain_size(3))
        inti = meani*product(domain)
    else
        meani = meani / (params%domain_size(1)*params%domain_size(2))
        inti = meani*product(domain(1:2))
    endif

    if (rank == 0) then
        write(*,*) "Computed mean value is: ", meani
        write(*,*) "Computed integral value is: ", inti

        ! write volume integral to disk
        call get_command_argument(3,fname_out)
        open(14,file=fname_out, status='replace')
        write(14,*) inti
        close(14)
    endif
end subroutine post_mean
