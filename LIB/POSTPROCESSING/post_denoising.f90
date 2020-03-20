subroutine post_denoising(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_initialization
    use module_helpers
    use module_forest

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: indicator="threshold-state-vector", file_in, file_out
    real(kind=rk)          :: time, eps=-1.0_rk, Z, x0(1:3), dx(1:3)
    integer(kind=ik), allocatable           :: lgt_block(:, :), hvy_n(:), lgt_n(:)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:,:)
    integer(kind=ik)                        :: max_neighbors, level, k, tc_length, hvy_id, lgt_id, iter
    integer(hid_t)                          :: file_id
    character(len=80)                       :: order
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: treecode_size, number_dense_blocks, i, l, dim, iteration
    integer(kind=ik)                        :: lgt_n_tmp, tree_n, g, mpierr, Bs(1:3)
    !-----------------------------------------------------------------------------------------------------

    call get_command_argument(2, file_in)
    if (file_in == '--help' .or. file_in == '--h') then
        if ( params%rank==0 ) then
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "                DENOISING ‹(•¿•)›"
            write(*,*) "--------------------------------------------------------------"
            write(*,*) " --eps-norm=L2"
            write(*,*) " --wavelet=CDF44"
            write(*,*) " --input=data_00000.h5 "
            write(*,*) " --output=data_00000.h5 "
            write(*,*) "-------------------------------------------------------------"
            write(*,*)
        end if
        return
    end if


    !----------------------------------
    ! read parameters
    !----------------------------------
    call get_cmd_arg_str( "--eps-norm", params%eps_norm, default="L2" )
    call get_cmd_arg_str( "--input", file_in, default="" )
    call get_cmd_arg_str( "--output", file_out, default="" )
    call get_cmd_arg_str( "--indicator", indicator, default="threshold-state-vector" )
    call get_cmd_arg_str( "--wavelet", order, default="CDF44" )
    call get_cmd_arg_dbl( "--eps", params%eps, default=1.0e-2_rk )


    ! Check parameters for correct inputs:
    if (order == "CDF20") then
        params%harten_multiresolution = .true.
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik

    elseif (order == "CDF40") then
        params%harten_multiresolution = .true.
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik

    elseif (order == "CDF44") then
        params%harten_multiresolution = .false.
        params%wavelet_transform_type = 'biorthogonal'
        params%order_predictor = "multiresolution_4th"
        params%wavelet='CDF4,4'
        params%n_ghosts = 6_ik
    else
        call abort(20030202, "The --order parameter is not correctly set [CDF40, CDF20, CDF44]")

    end if

    g = params%n_ghosts

    params%coarsening_indicator = indicator
    params%forest_size = 3
    params%n_eqn = 1
    allocate(params%field_names(params%n_eqn))
    allocate(params%threshold_state_vector_component(params%n_eqn))
    params%eps_normalized = .false.
    params%threshold_state_vector_component = .true.
    params%block_distribution = "sfc_hilbert"
    allocate(params%butcher_tableau(1,1))



    !-------------------------------------------
    ! check and find common params in all h5-files
    !-------------------------------------------
    call read_attributes(file_in, lgt_n_tmp, time, iteration, params%domain_size, params%Bs, params%max_treelevel, params%dim)
    Bs = params%Bs

    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (times 3 because we need total, coheren, incohrent)
    params%number_blocks = ceiling(  3.0*real(lgt_n_tmp)/real(params%number_procs) )


!*****************
    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    hvy_neighbor = -1
    lgt_n = 0 ! reset number of active light blocks
    hvy_n = 0
    tree_n = 0 ! reset number of trees in forest

    call read_field2tree(params, (/file_in/), 1, 1, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    ! synching is required for the adaptation step
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,1), hvy_n(1) )

    ! initial guess for incohrent is total
    call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_dest=2, tree_id_source=1)



    ! call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,2), lgt_n(2), &
    !     lgt_sortednumlist(:,:,2), hvy_active(:,2), hvy_n(2), 2, indicator, hvy_tmp)


    do iter = 1, 10
        !=======================================================================
        ! from current incoh part, compute threshold
        Z = 0.0_rk
        do k = 1, hvy_n(2)
            hvy_id = hvy_active(k, 2)
            call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
            Z = Z + sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, 1, hvy_id)**2 )*dx(1)*dx(2)
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, Z, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        Z = Z / (2.0_rk*product( params%domain_size(1:params%dim) ))
        params%eps =  sqrt( (4.0_rk/3.0_rk) * Z * log( (Bs(1)-1)**params%dim * 2.0_rk**( dble(params%dim*params%max_treelevel))) )
        if (params%rank==0) write(*,*) "eps=", params%eps
        !=======================================================================

        call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_dest=2, tree_id_source=1)

        call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,2), lgt_n(2), &
        lgt_sortednumlist(:,:,2), hvy_active(:,2), hvy_n(2), 2, indicator, hvy_tmp)

        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

        call substract_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, 2, 1)
        ! now tree=2 is the incoherent part

        ! call multiply_tree_with_scalar(params, hvy_block, hvy_active, hvy_n, 2, -1.0_rk)
    enddo

    call write_tree_field( "y_000000000000.h5", params, lgt_block, lgt_active, hvy_block, &
            lgt_n, hvy_n, hvy_active, 1, 2, time, iteration )

    !
    ! call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    ! lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID_flow, params%coarsening_indicator, hvy_block ) !???? hvy_tmp


! fuer CVS brauchen wir dann noch \epsilon was von Z (der enstrophy) und der feinsten aufloesung abhaengt. fuer L^2 normalisierte wavelets ist
! der threshold:
!
! \epsilon = \sqrt{2/3 \sigma^2 \ln N}
!
! wobei \sigma^2 die varianz (= 2 Z) der incoh. vorticity ist.
! typischerweise erhaelt man diese mit 1-3 iterationen.
! als ersten schritt koennen wir einfach Z der totalen stroemung nehmen.
! N ist die maximale aufloesung, typicherweise 2^{d J}.



    ! call write_field(file_out, time, iteration, 1, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_active )

end subroutine post_denoising
