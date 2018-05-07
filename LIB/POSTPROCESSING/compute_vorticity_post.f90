!> \file
! WABBIT
!> \name compute_vorticity_post.f90
!> \version 0.5
!> \author sm
!
!> \brief postprocessing routine for subsequent vorticity calculation from datafields ux, uy (, uz) saved in .h5 files
! = log ======================================================================================
!
!> \version 02/02/18 - create commit 13cb3d25ab12e20cb38e5b87b9a1e27a8fe387e8
!-----------------------------------------------------------------------------------------------------

subroutine compute_vorticity_post(help, params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_initialization, only: allocate_grid, allocate_com_arrays
    use module_mpi
    use module_operators

    implicit none

    !> help flag
    logical, intent(in)                :: help
    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: file_ux, file_uy, file_uz
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, lgt_n, hvy_n, max_neighbors
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=1), allocatable          :: hvy_synch(:, :, :, :)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    integer(kind=ik), allocatable      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable         :: real_send_buffer(:,:), real_receive_buffer(:,:)
    integer(hsize_t), dimension(4)     :: size_field
    character(len=80)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik), allocatable      :: com_matrix(:,:,:)
    integer(kind=ik), allocatable      :: com_lists(:, :, :, :)
    integer(hsize_t), dimension(2)     :: dims_treecode

!-----------------------------------------------------------------------------------------------------

    if (help) then
        if (params%rank==0) then
            write(*,*) "wabbit postprocessing routine for subsequent vorticity calculation"
            write(*,*) "mpi_command -n number_procs ./wabbit-post 2D --vorticity source_ux.h5 source_uy.h5 derivative-order(2 or 4)"
            write(*,*) "mpi_command -n number_procs ./wabbit-post 3D --vorticity source_ux.h5 source_uy.h5 source_uz.h5 derivative-order(2 or 4)"
        end if
    else
        ! get values from command line (filename and level for interpolation)
        call get_command_argument(3, file_ux)
        call check_file_exists(trim(file_ux))
        call get_command_argument(4, file_uy)
        call check_file_exists(trim(file_uy))
        if (params%threeD_case) then
            call get_command_argument(5, file_uz)
            call check_file_exists(trim(file_uz))
            call get_command_argument(6, order)
        else
            call get_command_argument(5, order)
        end if
        if (order == "4") then
            params%order_discretization = "FD_4th_central_optimized"
            params%number_ghost_nodes = 4_ik
        elseif (order == "2") then
            params%order_discretization = "FD_2nd_central"
            params%number_ghost_nodes = 2_ik
        else
            call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")
        end if
        params%order_predictor = "multiresolution_4th"

        ! get some parameters from one of the files (they should be the same in all of them)
        call open_file_hdf5( trim(adjustl(file_ux)), file_id, .false.)
        if ( params%threeD_case ) then
            call get_size_datafield(4, file_id, "blocks", size_field)
            params%number_data_fields  = 3
            max_neighbors = 74
        else
            call get_size_datafield(3, file_id, "blocks", size_field(1:3))
            params%number_data_fields  = 2
            max_neighbors = 12
        end if
        call get_size_datafield(2, file_id, "block_treecode", dims_treecode)
        params%max_treelevel = int(dims_treecode(1), kind=ik)
        call close_file_hdf5(file_id)
        call get_attributes(file_ux, lgt_n, time, iteration, domain)
        ! only lgt_n/number_procs blocks necessary (since we do not want to refine)
        params%number_blocks = lgt_n/params%number_procs
        if (params%rank==0) params%number_blocks = params%number_blocks + &
            mod(lgt_n, params%number_procs)
        params%Lx = domain(1)
        params%Ly = domain(2)
        if (params%threeD_case) params%Lz = domain(3)
        params%number_block_nodes = int(size_field(1),kind=ik)
        params%mpi_data_exchange = "Non_blocking_Isend_Irecv"

        ! allocate data
        call allocate_grid( params, lgt_block, hvy_block, hvy_work, hvy_synch, &
            hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist,&
            int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
        ! allocate communication arrays
        call allocate_com_arrays(params, com_lists, com_matrix)
        ! read mesh and field
        call read_mesh(file_ux, params, lgt_n, hvy_n, lgt_block)
        call read_field(file_ux, 1, params, hvy_block, hvy_n)
        call read_field(file_uy, 2, params, hvy_block, hvy_n)
        if (params%threeD_case) call read_field(file_uz, 3, params, hvy_block, hvy_n)
        ! create lists of active blocks (light and heavy data)
        ! update list of sorted nunmerical treecodes, used for finding blocks
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
            lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
            lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
            hvy_active, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:), &
            com_matrix, .true., int_send_buffer, int_receive_buffer, &
            real_send_buffer, real_receive_buffer )

        ! calculate vorticity from velocities
        do k=1,hvy_n
            call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
            if (params%threeD_case) then
                call compute_vorticity(hvy_block(:,:,:,1,hvy_active(k)), &
                    hvy_block(:,:,:,2,hvy_active(k)), hvy_block(:,:,:,3,hvy_active(k)),&
                    dx, params%number_block_nodes, params%number_ghost_nodes,&
                    params%order_discretization, hvy_work(:,:,:,1:3,hvy_active(k)))
            else
                call compute_vorticity(hvy_block(:,:,:,1,hvy_active(k)), &
                    hvy_work(:,:,:,2,hvy_active(k)), hvy_work(:,:,:,3,hvy_active(k)),&
                    dx, params%number_block_nodes, params%number_ghost_nodes, &
                    params%order_discretization, hvy_work(:,:,:,:,hvy_active(k)))
            end if
        end do
        write( fname,'(a, "_", i12.12, ".h5")') 'vorx', nint(time * 1.0e6_rk)
        call write_field(fname, time, iteration, 1, params, lgt_block,&
            hvy_work, lgt_active, lgt_n, hvy_n)
        if (params%threeD_case) then
            write( fname,'(a, "_", i12.12, ".h5")') 'vory', nint(time * 1.0e6_rk)
            call write_field(fname, time, iteration, 2, params, lgt_block,&
                hvy_work, lgt_active, lgt_n, hvy_n)
            write( fname,'(a, "_", i12.12, ".h5")') 'vorz', nint(time * 1.0e6_rk)
            call write_field(fname, time, iteration, 3, params, lgt_block, &
                hvy_work, lgt_active, lgt_n, hvy_n)
        end if
    end if
end subroutine compute_vorticity_post
