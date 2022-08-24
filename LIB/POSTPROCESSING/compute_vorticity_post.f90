!> \brief postprocessing routine for subsequent vorticity calculation from datafields ux, uy (, uz) saved in .h5 files
!-----------------------------------------------------------------------------------------------------

subroutine compute_vorticity_post(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_ux, file_uy, file_uz, operator
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgt_id, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs
    character(len=2)                   :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: nwork

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
    call get_command_argument(2, file_ux)
    ! does the user need help?
    if (file_ux=='--help' .or. file_ux=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: vorticity / divergence / Q-criterion"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Computes either quantity from velocity files. Output is stored"
            write(*,*) " in predefined files."
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --vorticity"
            write(*,*) "./wabbit-post --vorticity source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER]"
            write(*,*) " Computes 3 (3D) or 1 (2D) vorticity component, saves in "
            write(*,*) " vorx_*.h5 [vory_*.h5] [vorz_*.h5]"
            write(*,*) " order = 2 or 4"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --vor-abs"
            write(*,*) "./wabbit-post --vor-abs source_ux.h5 source_uy.h5 source_uz.h5 [ORDER]"
            write(*,*) " Computes vorticity magnitude of 3D velocity field, saves in "
            write(*,*) " vorabs_*.h5"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --divergence"
            write(*,*) "./wabbit-post --divergence source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER]"
            write(*,*) " Computes divergence of 2D/3D velocity field, saves in "
            write(*,*) " divu_*.h5"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --Q"
            write(*,*) "./wabbit-post --Q source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER]"
            write(*,*) " Computes Q-criterion of 2D/3D velocity field, saves in "
            write(*,*) " Qcrit_*.h5"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(file_ux))

    call get_command_argument(3, file_uy)
    call check_file_exists(trim(file_uy))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_ux, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    if (params%dim == 3) then
        call get_command_argument(4, file_uz)
        call check_file_exists(trim(file_uz))
        call get_command_argument(5, order)
    else
        call get_command_argument(4, order)
    end if

    ! decide which order for discretization and predictor is used. Note predictor
    ! is used in ghost nodes sync'ing
    if (order == "4") then
        params%order_discretization = "FD_4th_central_optimized"
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik

    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik

    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")

    end if

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
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )

    nwork = 1
    if (operator == "--vorticity") then
        nwork = 3
    endif

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=nwork)

    ! read input data
    if (params%dim == 3) then
        call readHDF5vct_tree( (/file_ux, file_uy/), params, hvy_block, tree_ID)
    else
        call readHDF5vct_tree( (/file_ux, file_uy, file_uz/), params, hvy_block, tree_ID)
    end if

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call createActiveSortedLists_tree( params, tree_ID)

    ! update neighbor relations
    call updateNeighbors_tree( params, tree_ID )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    ! calculate vorticity from velocities
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if (operator == "--vorticity") then
            if (params%dim == 3) then
                call compute_vorticity( hvy_block(:,:,:,1,hvy_id), &
                hvy_block(:,:,:,2,hvy_id), hvy_block(:,:,:,3,hvy_id),&
                dx, Bs, g,&
                params%order_discretization, hvy_tmp(:,:,:,1:3,hvy_id))
            else
                call compute_vorticity(hvy_block(:,:,:,1,hvy_id), &
                hvy_block(:,:,:,2,hvy_id), hvy_block(:,:,:,1,hvy_id),&
                dx, Bs, g, &
                params%order_discretization, hvy_tmp(:,:,:,:,hvy_id))
            end if

        elseif (operator=="--vor-abs") then
            call compute_vorticity_abs(hvy_block(:,:,:,1,hvy_id), &
            hvy_block(:,:,:,2,hvy_id), hvy_block(:,:,:,3,hvy_id),&
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvy_id))

        elseif (operator == "--divergence") then
            call divergence( hvy_block(:,:,:,1,hvy_id), &
            hvy_block(:,:,:,2,hvy_id), &
            hvy_block(:,:,:,3,hvy_id),&
            dx, Bs, g, &
            params%order_discretization, hvy_tmp(:,:,:,1,hvy_id))

        elseif (operator == "--Q") then
            call compute_Qcriterion( hvy_block(:,:,:,1,hvy_id), &
            hvy_block(:,:,:,2,hvy_id), &
            hvy_block(:,:,:,3,hvy_id),&
            dx, Bs, g, &
            params%order_discretization, hvy_tmp(:,:,:,1,hvy_id))

        else
            call abort(1812011, "operator is neither --vorticity --vor-abs --divergence --Q")

        endif
    end do

    if (operator == "--vorticity") then
        write( fname,'(a, "_", i12.12, ".h5")') 'vorx', nint(time * 1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

        if (params%dim == 3) then
            write( fname,'(a, "_", i12.12, ".h5")') 'vory', nint(time * 1.0e6_rk)
            call saveHDF5_tree(fname, time, iteration, 2, params, hvy_tmp, tree_ID)
            write( fname,'(a, "_", i12.12, ".h5")') 'vorz', nint(time * 1.0e6_rk)
            call saveHDF5_tree(fname, time, iteration, 3, params, hvy_tmp, tree_ID)
        end if

    elseif (operator == "--vor-abs") then
        write( fname,'(a, "_", i12.12, ".h5")') 'vorabs', nint(time * 1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

    elseif (operator=="--divergence") then
        write( fname,'(a, "_", i12.12, ".h5")') 'divu', nint(time * 1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

    elseif (operator=="--Q") then
        write( fname,'(a, "_", i12.12, ".h5")') 'Qcrit', nint(time * 1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )
    endif
end subroutine compute_vorticity_post
