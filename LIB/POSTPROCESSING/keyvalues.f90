!> \brief loads the specified *.h5 file and creates a *.key file that contains
!! min / max / mean / L2 norm of the field data. This is used for testing
!! so that we don't need to store entire fields but rather the *.key only
!*********************************************************************************************

subroutine keyvalues(fname, params)

    use module_IO
    use module_precision
    use module_params
    use module_mesh
    use module_forestMetaData
    use mpi

    implicit none
    character(len=*), intent(in)            :: fname                            !> name of the file
    type (type_params), intent(inout)       :: params                           !> parameter struct
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                        :: tree_ID=1, hvy_id

    integer(hsize_t), dimension(4)          :: size_field
    integer(hid_t)                          :: file_id
    integer(kind=ik)                        :: lgt_id, k, nz, iteration, tc_length, dim
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), dimension(3)             :: x0, dx
    real(kind=rk), dimension(3)             :: domain
    real(kind=rk)                           :: time
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik), allocatable           :: tree(:), sum_tree(:), blocks_per_rank(:)
    integer(kind=ik)                        :: sum_curve(2)
    character(len=12)                       :: curves(2)
    real(kind=rk)                           :: x,y,z
    real(kind=rk)                           :: maxi,mini,squari,meani,qi
    real(kind=rk)                           :: maxl,minl,squarl,meanl,ql
    integer(kind=ik)                        :: ix,iy,iz,mpicode, ioerr, rank, i

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    rank = params%rank
    curves(1) = 'sfc_hilbert'
    curves(2) = 'sfc_z'

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    !-----------------------------------------------------------------------------------------------------
    if (fname=='--help' .or. fname=='--h') then
        if (rank==0) then
            write(*,*) "WABBIT postprocessing routine to load a specified *.h5 file and create a *.key file that contains  min / max / mean / L2 norm of the field data."
            write(*,*) "mpi_command -n number_procs ./wabbit-post --keyvalues filename.h5"
        end if
        return
    endif

    call check_file_exists( fname )

    if (rank==0) write (*,'("analyzing file ",a20," for keyvalues")') trim(adjustl(fname))

    ! get some parameters from the file
    call read_attributes(fname, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%dim = dim
    params%Bs = Bs
    params%max_treelevel = tc_length
    params%n_eqn = 1
    params%n_ghosts = 0
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    ! make sure there is enough memory allocated
    params%number_blocks = ceiling( real(lgt_n(tree_ID)) / real(params%number_procs) )
    params%order_predictor = "multiresolution_2nd"

    call allocate_forest(params, hvy_block)

    ! the work array needs to be allocated as balance load requires a buffer.
    ! it can However be smaller than what is allocated in allocate_grid.
    allocate( hvy_tmp( size(hvy_block,1), size(hvy_block,2), size(hvy_block,3), size(hvy_block,4), size(hvy_block,5)) )

    ! read data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID, verbosity=.true.)

    call updateMetadata_tree(params, tree_ID)

    ! compute an additional quantity that depends also on the position
    ! (the others are translation invariant)
    if (params%dim==2) Bs(3)=1

    allocate(tree(1:params%max_treelevel))
    allocate(sum_tree(1:params%max_treelevel))
    allocate(blocks_per_rank(1:params%number_procs))

    ! test all existing space filling curves
    do i = 1,size(curves)
        tree = 0_ik
        params%block_distribution=trim(curves(i))

        call balanceLoad_tree(params, hvy_block, tree_ID)


        call MPI_ALLGATHER(hvy_n,1,MPI_INTEGER,blocks_per_rank,1,MPI_INTEGER, &
        WABBIT_COMM,mpicode)

        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)

            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)

            tree = tree + (sum(blocks_per_rank(1:rank))+k)*lgt_block(lgt_id,1:params%max_treelevel)
        end do

        call MPI_REDUCE(tree,sum_tree, params%max_treelevel, MPI_INTEGER, &
        MPI_SUM, 0, WABBIT_COMM, mpicode)
        sum_curve(i) = sum(sum_tree)
    end do


    maxl = 0.0_rk
    minl = 99e99_rk
    squarl = 0.0_rk
    meanl = 0.0_rk
    ql = 0.0_rk

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        do iz = 1, Bs(3)
            do iy = 1, Bs(2)
                do ix = 1, Bs(1)
                    x = dble(ix-1)*dx(1)/dble(Bs(1)) + x0(1)
                    y = dble(iy-1)*dx(2)/dble(Bs(2)) + x0(2)
                    z = dble(iz-1)*dx(3)/dble(Bs(3)) + x0(3)

                    ql = ql + hvy_block(ix,iy,iz,1,hvy_id)*(x+y+z)
                end do
            end do
        end do

        maxl = max(maxl,maxval(hvy_block(:,:,:,:,hvy_id)))
        minl = min(minl,minval(hvy_block(:,:,:,:,hvy_id)))
        squarl = squarl + sum(hvy_block(:,:,:,:,hvy_id)**2)
        meanl  = meanl +sum(hvy_block(:,:,:,:,hvy_id))
    end do

    call MPI_REDUCE(ql,qi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(maxl,maxi,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(minl,mini,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(squarl,squari,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(meanl,meani,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,WABBIT_COMM,mpicode)

    qi = qi / lgt_n(tree_ID)
    squari = squari / lgt_n(tree_ID)
    meani = meani / lgt_n(tree_ID)

    if (rank == 0) then
        open  (59, file=fname(1:index(fname,'.'))//'key', &
        status = 'replace', action='write', iostat=ioerr)
        write (59,'(6(es15.8,1x), 2(i10,1x))') time, maxi, mini, meani, squari, qi , &
        sum_curve(1), sum_curve(2)
        write (*,'(A)') "Result:"
        write (* ,'(6(A15,1x),2(A12,1x))') "time","maxval","minval","meanval","sumsquares", &
        "Q-integral", curves(1), curves(2)
        write (* ,'(6(es15.8,1x),2(i12,1x))') time, maxi, mini, meani, squari, qi, &
        sum_curve(1), sum_curve(2)
        write (*,'(A)') "These values can be used to compare two HDF5 files"
        close (59)
    endif

    deallocate( hvy_tmp )

end subroutine keyvalues
