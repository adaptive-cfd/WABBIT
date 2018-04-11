!> \file
! WABBIT
!> \name sparse_to_dense.f90
!> \version 0.5
!> \author sm
!
!> \brief postprocessing routine for interpolation of a given field to the desired level
!
! = log ======================================================================================
!> \date  31/01/18 - create hashcode: commit 13cb3d25ab12e20cb38e5b87b9a1e27a8fe387e8
!-----------------------------------------------------------------------------------------------------

subroutine sparse_to_dense(help, params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_initialization, only: allocate_grid, allocate_com_arrays
    use module_mpi

    implicit none

    !> help flag
    logical, intent(in)                :: help
    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: file_in
    character(len=80)      :: file_out
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration

    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(kind=ik), allocatable           :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable              :: real_send_buffer(:,:), real_receive_buffer(:,:)
    integer(hsize_t), dimension(4)          :: size_field
    integer(kind=ik)                        :: hvy_n, lgt_n, max_neighbors, level, k
    integer(hid_t)                          :: file_id
    character(len=2)                        :: level_in, order
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik), allocatable           :: com_matrix(:,:,:)
    integer(kind=ik), allocatable           :: com_lists(:, :, :, :)
    integer(kind=ik)                        :: treecode_size
!-----------------------------------------------------------------------------------------------------

    if (help .eqv. .true.) then
        if (params%rank==0) then
            write(*,*) "postprocessing subroutine to refine/coarse mesh to a uniform grid (up and downsampling ensured). command line:"
            write(*,*) "mpi_command -n number_procs ./wabbit-post 2D --sparse-to-dense source.h5 target.h5 target_treelevel order-predictor(2 or 4)"
        end if
    else
        ! get values from command line (filename and level for interpolation)
        call get_command_argument(3, file_in)
        call check_file_exists(trim(file_in))
        call get_command_argument(4, file_out)
        call get_command_argument(5, level_in)
        read(level_in,*) level
        call get_command_argument(6, order)
        if (order == "4") then
            params%order_predictor = "multiresolution_4th"
            params%number_ghost_nodes = 4_ik
        elseif (order == "2") then
            params%order_predictor = "multiresolution_2nd"
            params%number_ghost_nodes = 2_ik
        else
            call abort(392,"ERROR: chosen predictor order invalid or not (yet) implemented. choose between 4 (multiresolution_4th) and 2 (multiresolution_2nd)")
        end if

        ! get some parameters from file
        call open_file_hdf5( trim(adjustl(file_in)), file_id, .false.)
        if ( params%threeD_case ) then
            call get_size_datafield(4, file_id, "blocks", size_field)
        else
            call get_size_datafield(3, file_id, "blocks", size_field(1:3))
        end if
        call get_size_datafield(2, file_id, "block_treecode", dims_treecode)
        call close_file_hdf5(file_id)
        call get_attributes(file_in, lgt_n, time, iteration, domain)
        ! set max_treelevel for allocation of hvy_block
        if (dims_treecode(1)<=level) then
            params%max_treelevel = level
        else
            params%max_treelevel = int(dims_treecode(1), kind=ik)
        end if
        params%min_treelevel = level
        params%Lx = domain(1)
        params%Ly = domain(2)
        if (params%threeD_case) then
            params%number_blocks = (8_ik**params%max_treelevel)/params%number_procs&
                + mod(8_ik**params%max_treelevel,params%number_procs)
            params%Lz = domain(3)
            max_neighbors = 74
        else
            params%number_blocks = (4_ik**params%max_treelevel)/params%number_procs&
                + mod(4_ik**params%max_treelevel,params%number_procs)
            max_neighbors = 12
        end if
        params%number_block_nodes = int(size_field(1),kind=ik)
        params%number_data_fields  = 1
        params%mpi_data_exchange = "Non_blocking_Isend_Irecv"
        ! allocate data
        call allocate_grid(params, lgt_block, hvy_block, &
            hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, .true., hvy_work, &
            int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer)
        ! allocate communication arrays
        call allocate_com_arrays(params, com_lists, com_matrix)
        ! read field
        call read_mesh(file_in, params, lgt_n, hvy_n, lgt_block)
        call read_field(file_in, 1, params, hvy_block, hvy_n)
        ! create lists of active blocks (light and heavy data)
        ! update list of sorted nunmerical treecodes, used for finding blocks
        call create_active_and_sorted_lists( params, lgt_block, &
            lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, &
            lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        ! balance the load
        params%block_distribution="sfc_hilbert"
        call balance_load(params, lgt_block, hvy_block,&
            hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n)

        ! create lists of active blocks (light and heavy data) after load balancing (have changed)
        call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
            lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active,&
            lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor,&
            hvy_active, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:), &
            com_matrix, .true., int_send_buffer, int_receive_buffer, &
            real_send_buffer, real_receive_buffer )
        ! refine/coarse to attain desired level, respectively
        !coarsen
        do while (max_active_level( lgt_block, lgt_active, lgt_n )>level)
            ! set refinement status to -1 (coarsen) everywhere
            lgt_block(:, params%max_treelevel +2) = -1
            ! check where coarsening is actually needed
            call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
            ! this might not be necessary since we start from an admissible grid
            call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )
            call ensure_completeness( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )
            call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )
            call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, &
                hvy_active, hvy_n, lgt_sortednumlist, .true. )
            ! update neighbor relations
            call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n,&
                lgt_sortednumlist, hvy_active, hvy_n )
        end do
        ! refine
        do while (min_active_level( lgt_block, lgt_active, lgt_n )<level)
            ! check where refinement is actually needed
            do k=1, lgt_n
                if (treecode_size(lgt_block(lgt_active(k),:), params%max_treelevel) < level)&
                    lgt_block(lgt_active(k), params%max_treelevel +2) = 1
            end do
            call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )
            if ( params%threeD_case ) then
                ! 3D:
                call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
            else
                ! 2D:
                call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:),&
                    hvy_active, hvy_n )
            end if
            call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
                lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
            ! update neighbor relations
            call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
                lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
            call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor,&
                hvy_active, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:),&
                com_matrix, .true., int_send_buffer, int_receive_buffer, &
                real_send_buffer, real_receive_buffer )
        end do

        call balance_load( params, lgt_block, hvy_block, &
            hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )
        call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
            lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        call write_field(file_out, time, iteration, 1, params, lgt_block, &
            hvy_block, lgt_active, lgt_n, hvy_n)

        if (params%rank==0 ) then
            write(*,'("Wrote data of input-file",a20," now on uniform grid (level",i3, ") to file",a20)'), &
                trim(file_in), level, trim(file_out)
             write(*,'("Minlevel:", i3," Maxlevel:" i3, " (should be identical now)")'), &
                 min_active_level( lgt_block, lgt_active, lgt_n ),&
                 max_active_level( lgt_block, lgt_active, lgt_n )
        end if
    end if
end subroutine sparse_to_dense
