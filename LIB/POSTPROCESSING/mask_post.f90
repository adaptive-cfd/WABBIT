subroutine mask_post(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators
    use module_acm, only : create_mask_3D

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: fname_ini, fname_mask, fname_grid
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, lgt_n, hvy_n, tc_length, g
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=80)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    real(kind=rk), allocatable :: us(:,:,:,:)
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)
    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h') then
        if (params%rank==0) then
            write(*,*) "wabbit postprocessing routine for subsequent vorticity calculation"
            write(*,*) " ./wabbit-post --mask PARAMS.ini grid_000.h5 mask_000.h5"
        end if
        return
    endif

    call get_command_argument(3, fname_grid)
    call get_command_argument(4, fname_mask)

    write(*,*) fname_ini, fname_grid, fname_mask

    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, fname_ini )
    ! have the pysics module read their own parameters
    call init_physics_modules( params, fname_ini )
    ! initalize debugging ( this is mainly time measurements )
    call allocate_init_timing( params )


    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname_grid, lgt_n, time, iteration, domain, Bs, tc_length, params%dim)

    ! only (4* , for safety) lgt_n/number_procs blocks necessary (since we do not want to refine)
    !> \todo change that for 3d case
    params%number_blocks = 4*lgt_n/params%number_procs
    Bs = params%Bs
    g = params%n_ghosts

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist, hvy_work=hvy_work)

    ! read mesh and field
    call read_mesh(fname_grid, params, lgt_n, hvy_n, lgt_block)
    call read_field(fname_grid, 1, params, hvy_block, hvy_n)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    if (.not. allocated(us)) allocate(us(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, 1:3))

    ! calculate vorticity from velocities
    do k = 1, hvy_n
        call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        call create_mask_3D( time, x0, dx, params%Bs, params%n_ghosts, &
        hvy_work(:,:,:,1, hvy_active(k)), us )
    end do

    call write_field(fname_mask, time, iteration, 1, params, lgt_block,&
    hvy_work, lgt_active, lgt_n, hvy_n, hvy_active )

end subroutine mask_post
