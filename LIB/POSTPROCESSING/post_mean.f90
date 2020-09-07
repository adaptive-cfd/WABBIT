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
!*********************************************************************************************

subroutine post_mean(params)
    use module_IO
    use module_precision
    use module_params
    use module_mesh
    use mpi

    implicit none
    !> name of the file
    character(len=80)            :: fname, fname_out
    !> parameter struct
    type (type_params), intent(inout)       :: params
    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(hsize_t), dimension(4)          :: size_field
    integer(hid_t)                          :: file_id
    integer(kind=ik)                        :: lgt_id, k, nz, iteration, lgt_n, hvy_n, dim, g
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
    call read_attributes(fname, lgt_n, time, iteration, domain, Bs, tc_length, dim)

    params%Bs = Bs
    params%n_eqn = 1
    params%n_ghosts = 2_ik
    params%order_predictor = "multiresolution_2nd"
    params%max_treelevel = tc_length
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%number_blocks = 2_ik*lgt_n/params%number_procs

    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
    hvy_active, lgt_sortednumlist)

    call read_mesh(fname, params, lgt_n, hvy_n, lgt_block)
    call read_field(fname, 1, params, hvy_block, hvy_n )

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)
    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    ! compute an additional quantity that depends also on the position
    ! (the others are translation invariant)
    if (params%dim == 3) then
        nz = Bs(3)
    else
        nz = 1
    end if

    meanl = 0.0_rk

    g = params%n_ghosts
    do k = 1,hvy_n
        call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        if (params%dim == 3) then
            meanl = meanl + sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, 1, hvy_active(k)))*dx(1)*dx(2)*dx(3)
        else
            meanl = meanl + sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, 1, hvy_active(k)))*dx(1)*dx(2)
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
