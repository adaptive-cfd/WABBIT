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

subroutine keyvalues(fname, params, help)

    use module_IO
    use module_precision
    use module_params
    use module_initialization, only: allocate_grid
    use module_mesh
    use mpi

    implicit none
    !> name of the file
    character(len=*), intent(in)            :: fname
    !> parameter struct
    type (type_params), intent(inout)       :: params
    !> help flag
    logical, intent(in)                     :: help
    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(kind=ik), allocatable           :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable              :: real_send_buffer(:,:), real_receive_buffer(:,:)
    integer(hsize_t), dimension(4)          :: size_field    
    integer(hid_t)                          :: file_id
    integer(kind=ik)                        :: lgt_id, k, Bs, nz, iteration, lgt_n, hvy_n
    real(kind=rk), dimension(3)             :: x0, dx
    real(kind=rk), dimension(3)             :: domain
    real(kind=rk)                           :: time
    integer(hsize_t), dimension(2)          :: dims_treecode

    real(kind=rk)    :: x,y,z
    real(kind=rk)    :: maxi,mini,squari,meani,qi
    real(kind=rk)    :: maxl,minl,squarl,meanl,ql
    integer(kind=ik) :: ix,iy,iz,mpicode, ioerr
!-----------------------------------------------------------------------------------------------------

    if (help) then
        if (params%rank==0) then
            write(*,*) "WABBIT postprocessing routine to load a specified *.h5 file and create a *.key file that contains  min / max / mean / L2 norm of the field data."
            write(*,*) "mpi_command -n number_procs ./wabbit-post 2[3]D --keyvalues filename.h5"
        end if
    else
        call check_file_exists( fname )

        if (params%rank==0) write (*,'("analyzing file ",a20," for keyvalues")') trim(adjustl(fname))

        ! ! get some parameters from the file
        call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)
        if ( params%threeD_case ) then
            call get_size_datafield(4, file_id, "blocks", size_field)
        else
            call get_size_datafield(3, file_id, "blocks", size_field(1:3))
        end if
        call get_size_datafield(2, file_id, "block_treecode", dims_treecode)
        params%max_treelevel = int(dims_treecode(1), kind=ik)
        call close_file_hdf5(file_id)
        params%number_block_nodes = int(size_field(1),kind=ik)
        params%number_data_fields = 1
        params%number_ghost_nodes = 0
        call get_attributes(fname, lgt_n, time, iteration, domain)
        params%Lx = domain(1)
        params%Ly = domain(2)
        if (params%threeD_case) params%Lz = domain(3)
        params%number_blocks = lgt_n/params%number_procs + mod(lgt_n,params%number_procs)
        call allocate_grid( params, lgt_block, hvy_block, hvy_work,&
            hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist,&
            int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
        call read_mesh(fname, params, lgt_n, hvy_n, lgt_block)
        call read_field(fname, 1, params, hvy_block, hvy_n )
        call create_active_and_sorted_lists( params, lgt_block, &
            lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        ! compute an additional quantity that depends also on the position
        ! (the others are translation invariant)
        Bs = params%number_block_nodes
        if (params%threeD_case) then 
            nz = Bs
        else
            nz = 1
        end if
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

        call MPI_REDUCE(ql,qi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpicode)
        call MPI_REDUCE(maxl,maxi,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpicode)
        call MPI_REDUCE(minl,mini,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,mpicode)
        call MPI_REDUCE(squarl,squari,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpicode)
        call MPI_REDUCE(meanl,meani,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpicode)

        qi = qi / lgt_n
        squari = squari / lgt_n
        meani = meani / lgt_n

        if (params%rank == 0) then
            open  (59, file=fname(1:index(fname,'.'))//'key', &
                status = 'replace', action='write', iostat=ioerr)
            write (59,'(6(es15.8,1x))') time, maxi, mini, meani, squari, qi
            write (*,'(A)') "Result:"
            write (* ,'(6(A15,1x))') "time","maxval","minval","meanval","sumsquares","Q-integral"
            write (* ,'(6(es15.8,1x))') time, maxi, mini, meani, squari, qi
            write (*,'(A)') "These values can be used to compare two HDF5 files"
            close (59)
        endif
    end if
end subroutine keyvalues
