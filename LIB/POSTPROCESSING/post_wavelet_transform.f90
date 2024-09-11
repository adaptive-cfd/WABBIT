
subroutine post_wavelet_transform(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: fname_in, fname_out, operator
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgt_ID, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_ID, ix, iy, iz
    logical                            :: full_tree, split_components

    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id, pos_underscore
    real(kind=rk), dimension(3)        :: domain

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, operator)
    call get_command_argument(2, fname_in)

    ! does the user need help?
    if (operator=='--help' .or. operator=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: compute wavelet decomposition or wavelet reconstruction"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " compute wavelet transformation for an existing field"
            write(*,*) ""
            write(*,*) " ./wabbit-post --wavelet-decompose ux_000.h5 --wavelet=CDF40"
            write(*,*) " ./wabbit-post --wavelet-reconstruct ux_000.h5 --wavelet=CDF40"
            write(*,*) ""
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --wavelet=CDF44"
            write(*,*) " --full_tree=0         - Save or read in full tree format or only leaf-layer"
            write(*,*) " --split_components=0  - For decomposition: Save individual files for each wavelet component"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif

    call get_cmd_arg( "--full_tree", full_tree, default=.false. )
    call get_cmd_arg( "--split-components", split_components, default=.false.)
    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )

    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    call setup_wavelet(params)

    params%forest_size = 1
    params%n_eqn = 1

    call check_file_exists(trim(fname_in))


    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname_in, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    params%Jmax = tc_length
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))

    Bs = params%Bs
    g  = params%g


    ! for full wavelet operations we need 8/7 for 3D or 4/3 for 2D as blocks
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk))

    if (.not. split_components) then
        params%n_eqn = 1
    else
        params%n_eqn = 2**params%dim
    endif

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

    ! read input data
    call readHDF5vct_tree( (/ fname_in /), params, hvy_block, tree_ID)

    call sync_ghosts_tree( params, hvy_block, tree_ID)

    ! decompose or reconstruct
    select case(operator)
    case("--wavelet-decompose")
        call wavelet_decompose_full_tree(params, hvy_block, tree_ID, hvy_tmp)

        if (.not. full_tree) then
            call prune_fulltree2leafs(params, tree_ID)
        endif

        if (split_components) then
            ! split components so that each 2x2 block of points has either SC or a WC and we save them individually, easier for visualization
            do k = 1, hvy_n(tree_ID)
                hvy_ID = hvy_active(k, tree_ID)

                ! loop over all 2x2 cluster of points
                do ix = 1, size(hvy_block, 1), 2
                    do iy = 1, size(hvy_block, 2), 2
                        do iz = 1, size(hvy_block, 3), 2
                            if (params%dim == 2) then
                                ! set components WX, WY, WXY
                                hvy_block(ix:ix+1, iy:iy+1, iz, 2, hvy_id) = hvy_block(ix+1, iy  , iz, 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz, 3, hvy_id) = hvy_block(ix  , iy+1, iz, 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz, 4, hvy_id) = hvy_block(ix+1, iy+1, iz, 1, hvy_id)
                                ! set SC
                                hvy_block(ix:ix+1, iy:iy+1, iz, 1, hvy_id) = hvy_block(ix  , iy  , iz, 1, hvy_id)
                            else
                                ! set components WX, WY, WXY, WZ, WXZ, WYZ, WXYZ
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 2, hvy_id) = hvy_block(ix+1, iy  , iz  , 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 3, hvy_id) = hvy_block(ix  , iy+1, iz  , 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 4, hvy_id) = hvy_block(ix+1, iy+1, iz  , 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 5, hvy_id) = hvy_block(ix  , iy  , iz+1, 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 6, hvy_id) = hvy_block(ix+1, iy  , iz+1, 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 7, hvy_id) = hvy_block(ix  , iy+1, iz+1, 1, hvy_id)
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 8, hvy_id) = hvy_block(ix+1, iy+1, iz+1, 1, hvy_id)
                                ! set SC
                                hvy_block(ix:ix+1, iy:iy+1, iz:iz+1, 1, hvy_id) = hvy_block(ix  , iy  , iz  , 1, hvy_id)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif

    case("--wavelet-reconstruct")
        ! we need the saved values from hvy_tmp so currently we cannot use this at all
        call abort(240909, "Reconstruction currently is not implemented as I don't know what to do with the saved values in hvy_tmp")

        ! this currently only works if we read in the full tree
        if (.not. full_tree) then
            call abort(240909, "Reconstruction currently only works when providing the values in full tree format")
        endif

        call wavelet_reconstruct_full_tree(params, hvy_block, hvy_tmp, tree_ID)
    end select

    ! append ending before last underscore or if it doesnt exist before last .
    pos_underscore = index(fname_in, "_", back=.true.)
    if (pos_underscore < 0) pos_underscore = index(fname_in, '_', back=.true.)

    ! save files - either one for all or one for each component for visiblity
    if (.not. split_components) then
        write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WD", fname_in(pos_underscore:LEN_TRIM(fname_in))
        call saveHDF5_tree(fname_out, time, iteration, 1, params, hvy_block, tree_ID )
    else
        ! save SC, WX, WY, WXY
        write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-SC", fname_in(pos_underscore:LEN_TRIM(fname_in))
        call saveHDF5_tree(fname_out, time, iteration, 1, params, hvy_block, tree_ID )
        write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WX", fname_in(pos_underscore:LEN_TRIM(fname_in))
        call saveHDF5_tree(fname_out, time, iteration, 2, params, hvy_block, tree_ID )
        write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WY", fname_in(pos_underscore:LEN_TRIM(fname_in))
        call saveHDF5_tree(fname_out, time, iteration, 3, params, hvy_block, tree_ID )
        write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WXY", fname_in(pos_underscore:LEN_TRIM(fname_in))
        call saveHDF5_tree(fname_out, time, iteration, 4, params, hvy_block, tree_ID )
        if (params%dim == 3) then
            ! save XZ, WXZ, WYZ, WXYZ
            write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WZ", fname_in(pos_underscore:LEN_TRIM(fname_in))
            call saveHDF5_tree(fname_out, time, iteration, 5, params, hvy_block, tree_ID )
            write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WXZ", fname_in(pos_underscore:LEN_TRIM(fname_in))
            call saveHDF5_tree(fname_out, time, iteration, 6, params, hvy_block, tree_ID )
            write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WYZ", fname_in(pos_underscore:LEN_TRIM(fname_in))
            call saveHDF5_tree(fname_out, time, iteration, 7, params, hvy_block, tree_ID )
            write(fname_out, '(A, A, A)') fname_in(1:pos_underscore-1), "-WXYZ", fname_in(pos_underscore:LEN_TRIM(fname_in))
            call saveHDF5_tree(fname_out, time, iteration, 8, params, hvy_block, tree_ID )
        endif
    endif


end subroutine
