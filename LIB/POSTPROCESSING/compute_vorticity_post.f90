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
    use module_mpi
    use module_operators

    implicit none

    !> help flag
    logical, intent(in)                :: help
    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: file_ux, file_uy, file_uz
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, lgt_n, hvy_n, Bs, tc_length, dim
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=80)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    !-----------------------------------------------------------------------------------------------------

    if (help) then
        if (params%rank==0) then
            write(*,*) "wabbit postprocessing routine for subsequent vorticity calculation"
            write(*,*) "mpi_command -n number_procs ./wabbit-post 2D --vorticity source_ux.h5 source_uy.h5 derivative-order(2 or 4)"
            write(*,*) "mpi_command -n number_procs ./wabbit-post 3D --vorticity source_ux.h5 source_uy.h5 source_uz.h5 derivative-order(2 or 4)"
        end if
        return
    endif

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

    ! decide which order for discretization and predictor is used. Note predictor
    ! is used in ghost nodes sync'ing
    if (order == "4") then
        params%order_discretization = "FD_4th_central_optimized"
        params%order_predictor = "multiresolution_4th"
        params%number_ghost_nodes = 4_ik

    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%order_predictor = "multiresolution_2nd"
        params%number_ghost_nodes = 2_ik

    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")

    end if

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_ux, lgt_n, time, iteration, domain, Bs, tc_length, dim)
    params%max_treelevel = tc_length
    params%number_data_fields = dim
    params%Lx = domain(1)
    params%Ly = domain(2)
    params%Lz = domain(3)
    params%number_block_nodes = Bs
    allocate(params%butcher_tableau(1,1))
    params%non_uniform_mesh_correction = .true. ! This is an important switch for the OLD ghost nodes.
    ! only (4* , for safety) lgt_n/number_procs blocks necessary (since we do not want to refine)
    !> \todo change that for 3d case
    params%number_blocks = 4_ik*lgt_n/params%number_procs
    if (params%rank==0) params%number_blocks = params%number_blocks + &
    mod(lgt_n, params%number_procs)

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist, .true., hvy_work)

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

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

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
            hvy_block(:,:,:,2,hvy_active(k)), hvy_work(:,:,:,3,hvy_active(k)),&
            dx, params%number_block_nodes, params%number_ghost_nodes, &
            params%order_discretization, hvy_work(:,:,:,:,hvy_active(k)))
        end if
    end do
    write( fname,'(a, "_", i12.12, ".h5")') 'vorx', nint(time * 1.0e6_rk)
    call write_field(fname, time, iteration, 1, params, lgt_block,&
    hvy_work, lgt_active, lgt_n, hvy_n, hvy_active )
    if (params%threeD_case) then
        write( fname,'(a, "_", i12.12, ".h5")') 'vory', nint(time * 1.0e6_rk)
        call write_field(fname, time, iteration, 2, params, lgt_block,&
        hvy_work, lgt_active, lgt_n, hvy_n,  hvy_active)
        write( fname,'(a, "_", i12.12, ".h5")') 'vorz', nint(time * 1.0e6_rk)
        call write_field(fname, time, iteration, 3, params, lgt_block, &
        hvy_work, lgt_active, lgt_n, hvy_n, hvy_active)
    end if
end subroutine compute_vorticity_post
