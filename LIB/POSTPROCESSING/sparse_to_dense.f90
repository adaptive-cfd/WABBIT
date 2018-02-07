!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_attributes.f90
!> \version 0.5
!> \author sm
!
!> \brief postprocessing routine for interpolation of a given field to the desired level
!
!>
!! input: \n
!!           - help flag, parameter array \n
!! = log ======================================================================================
!! \n
!! 31/01/18 - create

subroutine sparse_to_dense(help, params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_initialization
    use module_mpi

    implicit none

    logical, intent(in)                :: help
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

    if (help .eqv. .true.) then
        if (params%rank==0) then
            write(*,*) "postprocessing subroutine to refine/coarse mesh to a uniform grid (up and downsampling ensured). command line:"
            write(*,*) "./wabbit-post 2D --sparse-to-dense source.h5 target.h5 target_treelevel order-predictor(2 or 4)"
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
            call abort("chosen predictor order invalid or not (yet) implemented. choose between 4 (multiresolution_4th) and 2 (multiresolution_2nd)")
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
            params%number_blocks = 8_ik**params%max_treelevel
            params%Lz = domain(3)
            max_neighbors = 74
        else
            params%number_blocks = 4_ik**params%max_treelevel
            max_neighbors = 12
        end if
        params%number_block_nodes = int(size_field(1),kind=ik)
        params%number_data_fields  = 1
        params%mpi_data_exchange = "Non_blocking_Isend_Irecv"
        ! allocate data
        call allocate_grid( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
        ! allocate communication arrays
        call allocate_com_arrays(params, com_lists, com_matrix)
        ! read field
        call read_mesh_and_attributes(file_in, params, lgt_n, hvy_n, lgt_block, time, iteration)
        call read_field(file_in, 1, params, hvy_block, hvy_n)
        ! create lists of active blocks (light and heavy data)
        ! update list of sorted nunmerical treecodes, used for finding blocks
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

        ! ! balance the load
        ! params%block_dist="sfc_hilbert"
        ! if (params%threeD_case) then
        !     call balance_load_3D(params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n)
        ! else
        !     call balance_load_2D(params, lgt_block, hvy_block(:,:,1,:,:), hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n)
        ! end if

        ! ! create lists of active blocks (light and heavy data)
        ! ! update list of sorted nunmerical treecodes, used for finding blocks
        ! call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        ! ! update neighbor relations
        ! call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        ! refine/coarse to attain desired level, respectively
        !coarsen 
        do while (max_active_level( lgt_block, lgt_active, lgt_n )>level)
            lgt_block(:, params%max_treelevel +2) = -1
            call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
            call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )
            call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        end do
        ! refine
        do while (min_active_level( lgt_block, lgt_active, lgt_n )<level)
            ! update neighbor relations
            call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
            call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:), com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
            do k=1, lgt_n
                if (treecode_size(lgt_block(lgt_active(k),:), params%max_treelevel) < level) lgt_block(lgt_active(k), params%max_treelevel +2) = 1
            end do
            if ( params%threeD_case ) then
                ! 3D:
                call refinement_execute_3D( params, lgt_block, hvy_block(:,:,:,:,:), hvy_active, hvy_n )
            else
                ! 2D:
                call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_active, hvy_n )
            end if
            call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        end do
        call write_field(file_out, time, iteration, 1, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n)

        if (params%rank==0 ) write(*,'("Wrote data of input-file",a20," now on uniform grid (level",i3, ") to file",a20)'), trim(file_in), level, trim(file_out)
    end if
end subroutine sparse_to_dense
