!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name compute_vorticity_post.f90
!> \version 0.5
!> \author sm
!
!> \brief postprocessing routine for subsequent computation of the vorticity
!
!>
!! input: \n
!!           - help flag, parameter array \n
!! = log ======================================================================================
!! \n
!! 02/02/18 - create

subroutine compute_vorticity_post(help, params)
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
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration
    character(len=2)       :: order

    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(kind=ik), allocatable           :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable              :: real_send_buffer(:,:), real_receive_buffer(:,:)
    integer(hsize_t), dimension(4)          :: size_field

    if (help .eqv. .true.) then
        if (params%rank==0) then
            write(*,*) ""
            write(*,*) "./wabbit-post 2D --compute-vorticity source.h5 derivative-order(2 or 4)"
        end if
    else
        ! get values from command line (filename and level for interpolation)
        call get_command_argument(3, file_ux)
        call check_file_exists(trim(file_ux))
        call get_command_argument(4, file_uy)
        call check_file_exists(trim(file_uy))
        if (params%threeD_case) then
            call abort(888,"3d not yet implemented")
        else
            call get_command_argument(4, order)
        end if
        if (order .eqv. "4") then
            params%order_discretization = "FD_4th_central_optimized"
            params%number_ghost_nodes = 4_ik
        elseif (order .eqv. "2") then
            params%order_discretization = "FD_2nd_central"
            params%number_ghost_nodes = 2_ik
        else
            call abort("chosen discretization order invalid or not yet implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")
        end if

        ! get some parameters from one of the files
        call open_file_hdf5( trim(adjustl(file_ux)), file_id, .false.)
        if ( params%threeD_case ) then
            call get_size_datafield(4, file_id, "blocks", size_field)
            params%number_data_fields  = 3
        else
            call get_size_datafield(3, file_id, "blocks", size_field(1:3))
            params%number_data_fields  = 2
        end if
        call close_file_hdf5(file_id)
        call get_attributes(file_ux, lgt_n, time, iteration, domain)
        params%Lx = domain(1)
        params%Ly = domain(2)
        if (params%threeD_case) then
            params%number_blocks = 8_ik**params%max_treelevel
            params%Lz = domain(3)
        else
            params%number_blocks = 4_ik**params%max_treelevel
        end if
        ! set some parameters
        params%number_block_nodes = int(size_field(1),kind=ik)
        ! allocate data
        call allocate_grid( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
        ! read field
        call read_mesh_and_attributes(file_in, params, lgt_n, hvy_n, lgt_block, time, iteration)
        call read_field(file_ux, 1, params, hvy_block, hvy_n)
        call read_field(file_uy, 2, params, hvy_block, hvy_n)

        call compute_vorticity(hvy_block(:,:,1,:,, v, w, dx, params%number_block_nodes, params%number_ghost_nodes, params%order_discretization, hvy_work(:,:,1,...))
    end if
    call write_field(fname, time, iteration, 1, params, lgt_block, hvy_work(:,:,:,:,:), lgt_active, lgt_n, hvy_n)
end subroutine compute_vorticity_post
