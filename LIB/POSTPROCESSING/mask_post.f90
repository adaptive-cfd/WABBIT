subroutine mask_post(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators
    use module_acm, only : create_mask_3D_ACM

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)      :: fname_ini, fname_mask, fname_grid
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, tc_length, g
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), allocatable      :: hvy_n(:), lgt_n(:)
    integer(kind=ik)                   :: tree_ID=1, hvy_id

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    real(kind=rk), allocatable :: us(:,:,:,:)
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

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
    call read_attributes(fname_grid, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, &
    params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! only (4* , for safety) lgt_n/number_procs blocks necessary (since we do not want to refine)
    !> \todo change that for 3d case
    params%number_blocks = 4*lgt_n(tree_ID)/params%number_procs
    Bs = params%Bs
    g = params%n_ghosts

    ! allocate data
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_work=hvy_work, hvy_n=hvy_n, lgt_n=lgt_n)

    ! read mesh and field
    call read_mesh(fname_grid, params, lgt_n, hvy_n, lgt_block, tree_ID)
    call read_field(fname_grid, 1, params, hvy_block, hvy_n, tree_ID)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call createActiveSortedLists_tree( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID)

    ! update neighbor relations
    call updateNeighbors_tree( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID)

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    if (.not. allocated(us)) allocate(us(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, 1:3))

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        call create_mask_3D_ACM( time, x0, dx, params%Bs, params%n_ghosts, &
        hvy_work(:,:,:,1, hvy_id), us )
    end do

    call saveHDF5_tree(fname_mask, time, iteration, 1, params, lgt_block,&
    hvy_work, lgt_active, lgt_n, hvy_n, hvy_active, tree_ID )

end subroutine mask_post
