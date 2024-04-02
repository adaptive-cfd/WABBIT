!> \brief loads the specified *.h5 file and creates a *.key file that contains
!! min / max / mean / L2 norm of the field data. This is used for testing
!! so that we don't need to store entire fields but rather the *.key only
!*********************************************************************************************

subroutine keyvalues(fname, params)

    use module_globals
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
    integer(kind=ik), allocatable           :: tree(:), tree_temp(:), sum_tree(:), blocks_per_rank(:)
    integer(kind=ik)                        :: sum_curve(2)
    character(len=12)                       :: curves(2)
    real(kind=rk)                           :: x,y,z
    real(kind=rk)                           :: val_max, val_min, val_squares, val_mean, val_Q, val_grid
    integer(kind=ik)                        :: ix,iy,iz,mpicode, ioerr, rank, i, g

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
    params%Jmax = tc_length
    params%n_eqn = 1
    params%g = 2
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    ! make sure there is enough memory allocated
    params%number_blocks = 2*ceiling( real(lgt_n(tree_ID)) / real(params%number_procs) )
    ! actually unused, but initialization needs to be done
    params%wavelet = "CDF20"
    call setup_wavelet(params)
    g = params%g

    call allocate_forest(params, hvy_block)

    ! the work array needs to be allocated as balance load requires a buffer.
    ! it can However be smaller than what is allocated in allocate_grid.
    allocate( hvy_tmp( size(hvy_block,1), size(hvy_block,2), size(hvy_block,3), size(hvy_block,4), size(hvy_block,5)) )

    ! read data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID, verbosity=.true., synchronize_ghosts=.false.)

    call updateMetadata_tree(params, tree_ID)

    allocate(tree(1:params%Jmax))
    allocate(tree_temp(1:params%Jmax))
    allocate(sum_tree(1:params%Jmax))
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
            
            ! Julius: CHANGE_LGT_BLOCK
            ! Probably this can be changed to a implementation without arrays but then I have to take
            ! a further look at it to stay consistent with old results
            tree_temp(:) = -1_ik
            call tcb2array( get_tc(lgt_block(lgt_id,IDX_TC_1 : IDX_TC_2)), &
                tree_temp, params%dim, lgt_block(lgt_id,IDX_MESH_LVL), params%Jmax)
            tree = tree + (sum(blocks_per_rank(1:rank))+k)*tree_temp(1:params%Jmax)
        end do

        call MPI_REDUCE(tree,sum_tree, params%Jmax, MPI_INTEGER, &
        MPI_SUM, 0, WABBIT_COMM, mpicode)
        sum_curve(i) = sum(sum_tree)
    end do


    val_max     = 0.0_rk
    val_min     = 99e99_rk
    val_squares = 0.0_rk
    val_mean    = 0.0_rk
    val_q       = 0.0_rk
    val_grid    = 0.0_rk

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if (params%dim == 3) then
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        x = dble(ix-1)*dx(1) + x0(1)
                        y = dble(iy-1)*dx(2) + x0(2)
                        z = dble(iz-1)*dx(3) + x0(3)

                        val_squares = val_squares + product(dx(1:3))*hvy_block(ix,iy,iz,1,hvy_id)**2
                        val_mean    = val_mean + product(dx(1:3))*hvy_block(ix,iy,iz,1,hvy_id)
                        val_Q       = val_Q + product(dx(1:3))*hvy_block(ix,iy,iz,1,hvy_id)*(x+y+z)
                        val_max     = max(val_max, hvy_block(ix,iy,iz,1,hvy_id))
                        val_min     = min(val_min, hvy_block(ix,iy,iz,1,hvy_id))
                    end do
                end do
            end do
        else
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    x = dble(ix-1)*dx(1) + x0(1)
                    y = dble(iy-1)*dx(2) + x0(2)

                    val_squares = val_squares + product(dx(1:2))*hvy_block(ix,iy,1,1,hvy_id)**2
                    val_mean    = val_mean + product(dx(1:2))*hvy_block(ix,iy,1,1,hvy_id)
                    val_Q       = val_Q + product(dx(1:2))*hvy_block(ix,iy,1,1,hvy_id)*(x+y)
                    val_max     = max(val_max, hvy_block(ix,iy,1,1,hvy_id))
                    val_min     = min(val_min, hvy_block(ix,iy,1,1,hvy_id))
                end do
            end do
        endif

        val_grid = val_grid + sum(x0)
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE,val_Q,1,MPI_DOUBLE_PRECISION,MPI_SUM,WABBIT_COMM,mpicode)
    call MPI_ALLREDUCE(MPI_IN_PLACE,val_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,WABBIT_COMM,mpicode)
    call MPI_ALLREDUCE(MPI_IN_PLACE,val_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,WABBIT_COMM,mpicode)
    call MPI_ALLREDUCE(MPI_IN_PLACE,val_squares,1,MPI_DOUBLE_PRECISION,MPI_SUM,WABBIT_COMM,mpicode)
    call MPI_ALLREDUCE(MPI_IN_PLACE,val_mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,WABBIT_COMM,mpicode)
    call MPI_ALLREDUCE(MPI_IN_PLACE,val_grid,1,MPI_DOUBLE_PRECISION,MPI_SUM,WABBIT_COMM,mpicode)

    val_mean = val_mean / product(params%domain_size(1:params%dim))

    if (rank == 0) then
        open  (59, file=fname(1:index(fname,'.'))//'key', &
        status = 'replace', action='write', iostat=ioerr)
        write (59,'(7(es15.8,1x), 2(i10,1x))') time, val_max, val_min, val_mean, val_squares, val_Q, val_grid , &
        sum_curve(1), sum_curve(2)
        write (*,'(A)') "Result:"
        write (* ,'(7(A15,1x),2(A12,1x))') "time","maxval","minval","meanval","sumsquares", &
        "Q-integral","grid-integral", curves(1), curves(2)
        write (* ,'(7(es15.8,1x),2(i12,1x))') time, val_max, val_min, val_mean, val_squares, val_Q, val_grid, &
        sum_curve(1), sum_curve(2)
        write (*,'(A)') "These values can be used to compare two HDF5 files"
        close (59)
    endif

    deallocate( hvy_tmp )

end subroutine keyvalues
