
subroutine post_filtertest(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_bridge_interface, only : initialize_communicator
    use module_time_step, only: filter_wrapper

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)      :: file_ux, infile
    real(kind=rk)          :: time,x ,y
    integer, allocatable :: lgt_n(:)
    integer, allocatable :: hvy_n(:)
    integer(kind=ik) :: iteration, k, lgt_id, tc_length, tree_N, i, l, ix, iy
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
    character(len=cshort)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    call get_command_argument(2, infile)
    call get_command_argument(3, file_ux)

    call check_file_exists(file_ux)

    call ini_file_to_params( params, infile )
    call initialize_communicator(params)
    call init_physics_modules( params, infile, params%N_mask_components )

    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp, hvy_mask, hvy_n, lgt_n)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call reset_forest(params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    call read_attributes(file_ux, lgt_n(tree_ID_flow), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    call read_mesh(file_ux, params, lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), lgt_block)

    do k = 1, params%n_eqn
        call read_field(file_ux, k, params, hvy_block, hvy_n(tree_ID_flow) )
    enddo

    call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active(:,tree_ID_flow), &
    lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), 1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! to get a mesh where coarsening is possible evrywhere refine one time
call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), &
lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), "everywhere", tree_ID=tree_ID_flow )

call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k = 1, hvy_n(tree_ID_flow)
    !!!!!!!
    ! call random_data(hvy_block(:,:,:,:,hvy_active(k,tree_ID_flow)))
    !!!!!!!!!

    call hvy2lgt( lgt_id, hvy_active(k,tree_ID_flow), params%rank, params%number_blocks )
    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

    do iy = 1, params%Bs(2)+2*params%n_ghosts
        y = real(iy-(params%n_ghosts+1), kind=rk) * dx(2) + x0(2)
        do ix = 1, params%Bs(1)+2*params%n_ghosts
            x = real(ix-(params%n_ghosts+1), kind=rk) * dx(1) + x0(1)

            ! use cos functions because theyre symmetric (symmetry BC)
            hvy_block(ix, iy, :, :, hvy_active(k,tree_ID_flow)) = cos(14._rk*x/params%domain_size(1) * 2.0_rk*pi) * cos(14._rk*y/params%domain_size(2) * 2.0_rk*pi)
        enddo
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

    call write_field("input_000020000000.h5", time, iteration, 1, params, lgt_block, hvy_block, lgt_active(:,tree_ID_flow),&
    lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), hvy_active(:,tree_ID_flow) )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

    call filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), hvy_neighbor)

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

    call write_field("filter_000020000000.h5", time, iteration, 1, params, lgt_block, hvy_block, lgt_active(:,tree_ID_flow),&
    lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), hvy_active(:,tree_ID_flow) )

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    call reset_forest(params, lgt_block, lgt_active, lgt_n,hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    call read_attributes("input_000020000000.h5", lgt_n(tree_ID_flow), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    call read_mesh("input_000020000000.h5", params, lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), lgt_block)

    do k = 1, params%n_eqn
        call read_field("input_000020000000.h5", k, params, hvy_block, hvy_n(tree_ID_flow) )
    enddo

    call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active(:,tree_ID_flow), &
    lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), 1)

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

    call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_ID_flow), &
    lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), &
    1, "everywhere", hvy_tmp=hvy_tmp, ignore_maxlevel=.true. )

    write(*,*) "adapted to=", lgt_n(tree_ID_flow)

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

    call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), &
    lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), "everywhere", tree_ID=tree_ID_flow )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

    write(*,*) "refined to=", lgt_n(tree_ID_flow)

    call write_field("coarsenrefine_000020000000.h5", time, iteration, 1, params, lgt_block, hvy_block, lgt_active(:,tree_ID_flow),&
    lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), hvy_active(:,tree_ID_flow) )


end subroutine
