!> \brief postprocessing routine for subsequent vorticity calculation from datafields ux, uy (, uz) saved in .h5 files
!-----------------------------------------------------------------------------------------------------

subroutine postGhostSyncTest(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_operators

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)                  :: file_in, file_out, file_uz, operator
    real(kind=rk)                      :: time, x, y
    integer(kind=ik)                   :: iteration, k, lgt_id, lgt_n, hvy_n, tc_length, g, ix, iy, hvy_id
    integer(kind=ik), dimension(3)     :: Bs
    character(len=2)                   :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=cshort)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: nwork

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, operator)
    call get_command_argument(2, file_in)
    call get_command_argument(3, file_out)


    ! does the user need help?
    if (file_in=='--help' .or. file_in=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " ./wabbit-post --ghost-test ux_000.h5 outfile_00.h5"
            write(*,*) " Note modified redundant points"
            write(*,*) "-----------------------------------------------------------"

        end if
        return
    endif
!
    call check_file_exists(trim(file_in))



    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_in, lgt_n, time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    NiterationsGhosts = 30000
    params%order_discretization = "FD_4th_central_optimized"
    params%order_predictor = "multiresolution_4th"
    params%n_ghosts = 6_ik

    params%iter_ghosts = .true.
    params%ghost_nodes_redundant_point_coarseWins = .false.
    params%wavelet_transform_type="biorthogonal"
    params%wavelet="CDF4,4"

    params%max_treelevel = tc_length
    params%n_eqn = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))

    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "x"
    params%symmetry_vector_component(2) = "y"
    if (params%dim==3) then
        params%symmetry_vector_component(3) = "z"
    endif

    Bs = params%Bs
    g  = params%n_ghosts

    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n)/real(params%number_procs) )

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist)

    ! read mesh
    call read_mesh(file_in, params, lgt_n, hvy_n, lgt_block)

    ! read actual data
    call read_field(file_in, 1, params, hvy_block, hvy_n)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    ! loop over my active heavy data
    do k = 1, hvy_n
        hvy_id = hvy_active(k)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! incl ghost nodes
        do iy = 1, Bs(2)+2*g
            do ix = 1, Bs(1)+2*g
                x = dble(ix-(g+1)) * dx(1) + x0(1) - 0.75_rk*0.5_rk
                y = dble(iy-(g+1)) * dx(2) + x0(2) - 0.5_rk
                hvy_block(ix,iy,:,:,hvy_id) = exp( -( (x)**2 + (y)**2 ) / 1.0e-4_rk )
            end do
        end do
    enddo

    call write_field("in_00.h5", time, iteration, 1, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_active )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    call write_field("out_00.h5", time, iteration, 1, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_active )
end subroutine postGhostSyncTest
