!> \file
! WABBIT
!> \name keyvalues.f90
!> \version 0.5
!> \author sm, engels
!
!> \brief loads the specified *.h5 file and creates a *.key file that contains
!! min / max / mean / L2 norm of the field data. This is used for testing
!! so that we don't need to store entire fields but rather the *.key only
!! \version 10/1/18 - create commit b2719e1aa2339f4f1f83fb29bd2e4e5e81d05a2a
!! \version 4/4/18 - add space filling curve test commit 47b4c1e24eabd112c950928e77e134046fc05d9a
!*********************************************************************************************

subroutine keyvalues(fname, params)

    use module_IO
    use module_precision
    use module_params
    use module_mesh
    use mpi

    implicit none
    !> name of the file
    character(len=*), intent(in)            :: fname
    !> parameter struct
    type (type_params), intent(inout)       :: params
    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    real(kind=rk), allocatable              :: hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable           :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(hsize_t), dimension(4)          :: size_field
    integer(hid_t)                          :: file_id
    integer(kind=ik)                        :: lgt_id, k, Bs, nz, iteration, lgt_n, hvy_n, tc_length, dim
    real(kind=rk), dimension(3)             :: x0, dx
    real(kind=rk), dimension(3)             :: domain
    real(kind=rk)                           :: time
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik), allocatable           :: tree(:), sum_tree(:), blocks_per_rank(:)

    integer(kind=ik)  :: sum_curve(2)
    character(len=12) :: curves(2)
    real(kind=rk)    :: x,y,z
    real(kind=rk)    :: maxi,mini,squari,meani,qi
    real(kind=rk)    :: maxl,minl,squarl,meanl,ql
    integer(kind=ik) :: ix,iy,iz,mpicode, ioerr, rank, i
    !-----------------------------------------------------------------------------------------------------
    rank = params%rank
    curves(1) = 'sfc_hilbert'
    curves(2) = 'sfc_z'
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
    call read_attributes(fname, lgt_n, time, iteration, domain, Bs, tc_length, dim)
    if (dim==3) then
        params%threeD_case = .true.
    else
        params%threeD_case = .false.
    end if
    params%Bs = Bs
    params%max_treelevel = tc_length
    params%n_eqn = 1
    params%nr_ghosts = 0
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    ! make sure there is enough memory allocated
    params%number_blocks = (dim**2) * (lgt_n/params%number_procs)

    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, .false.)

    ! the work array needs to be allocated as balance load requires a buffer.
    ! it can However be smaller than what is allocated in allocate_grid.
    allocate( hvy_work( size(hvy_block,1), size(hvy_block,2), size(hvy_block,3), size(hvy_block,4), size(hvy_block,5)  ) )

    call read_mesh(fname, params, lgt_n, hvy_n, lgt_block)
    call read_field(fname, 1, params, hvy_block, hvy_n )

    call create_active_and_sorted_lists( params, lgt_block, &
    lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
    ! compute an additional quantity that depends also on the position
    ! (the others are translation invariant)
    Bs = params%Bs
    if (params%threeD_case) then
        nz = Bs
    else
        nz = 1
    end if

    allocate(tree(1:params%max_treelevel))
    allocate(sum_tree(1:params%max_treelevel))
    allocate(blocks_per_rank(1:params%number_procs))

    ! test all existing space filling curves
    do i=1,size(curves)
        tree = 0_ik
        params%block_distribution=trim(curves(i))
        call balance_load(params, lgt_block, hvy_block, hvy_neighbor, &
        lgt_active, lgt_n, hvy_active, hvy_n, hvy_work)
        call create_active_and_sorted_lists( params, lgt_block, &
        lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
        lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        call MPI_ALLGATHER(hvy_n,1,MPI_INTEGER,blocks_per_rank,1,MPI_INTEGER, &
        WABBIT_COMM,mpicode)
        do k=1,hvy_n
            call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
            tree = tree + (sum(blocks_per_rank(1:rank))+k)*lgt_block(lgt_id,1:params%max_treelevel)
        end do
        call MPI_REDUCE(tree,sum_tree, params%max_treelevel, MPI_INTEGER, &
        MPI_SUM,0,WABBIT_COMM,mpicode)
        sum_curve(i) = sum(sum_tree)
    end do


    maxl = 0.0_rk
    minl = 99e99_rk
    squarl = 0.0_rk
    meanl = 0.0_rk
    ql = 0.0_rk

    do k = 1,hvy_n
        call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        do iz = 1, nz
            do iy = 1, Bs
                do ix = 1, Bs
                    x = dble(ix-1)*dx(1)/dble(Bs) + x0(1)
                    y = dble(iy-1)*dx(2)/dble(Bs) + x0(2)
                    z = dble(iz-1)*dx(3)/dble(nz) + x0(3)

                    ql = ql + hvy_block(ix,iy,iz,1,hvy_active(k))*(x+y+z)
                end do
            end do
        end do
        maxl = max(maxl,maxval(hvy_block(:,:,:,:,hvy_active(k))))
        minl = min(minl,minval(hvy_block(:,:,:,:,hvy_active(k))))
        squarl = squarl + sum(hvy_block(:,:,:,:,hvy_active(k))**2)
        meanl  = meanl +sum(hvy_block(:,:,:,:,hvy_active(k)))
    end do

    call MPI_REDUCE(ql,qi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(maxl,maxi,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(minl,mini,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(squarl,squari,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,WABBIT_COMM,mpicode)
    call MPI_REDUCE(meanl,meani,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,WABBIT_COMM,mpicode)

    qi = qi / lgt_n
    squari = squari / lgt_n
    meani = meani / lgt_n

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

    deallocate( hvy_work )

end subroutine keyvalues
