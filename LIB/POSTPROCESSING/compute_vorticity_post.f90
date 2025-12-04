!> \brief postprocessing routine for subsequent vorticity calculation from datafields ux, uy (, uz) saved in .h5 files
!-----------------------------------------------------------------------------------------------------

subroutine compute_vorticity_post(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_ux, file_uy, file_uz, operator
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs
    character(len=cshort)                   :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID

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
    if (file_ux=='--help' .or. file_ux=='--h' .or. file_ux=='-h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit postprocessing: vorticity / divergence / Q-criterion"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Computes either quantity from velocity files. Output is stored"
            write(*, '(A)') " in predefined files."
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " --vorticity"
            write(*, '(A)') "./wabbit-post --vorticity source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER]"
            write(*, '(A)') " Computes 3 (3D) or 1 (2D) vorticity component, saves in "
            write(*, '(A)') " vorx_*.h5 [vory_*.h5] [vorz_*.h5]"
            write(*, '(A)') " order = 2 or 4"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " --vor-abs"
            write(*, '(A)') "./wabbit-post --vor-abs source_ux.h5 source_uy.h5 source_uz.h5 [ORDER]"
            write(*, '(A)') " Computes vorticity magnitude of 3D velocity field, saves in "
            write(*, '(A)') " vorabs_*.h5"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " --divergence"
            write(*, '(A)') "./wabbit-post --divergence source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER]"
            write(*, '(A)') " Computes divergence of 2D/3D velocity field, saves in "
            write(*, '(A)') " divu_*.h5"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " --Q"
            write(*, '(A)') "./wabbit-post --Q source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER]"
            write(*, '(A)') " Computes Q-criterion of 2D/3D velocity field, saves in "
            write(*, '(A)') " Qcrit_*.h5"
            write(*, '(A)') "-----------------------------------------------------------"
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
    if (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%wavelet = "CDF20"  ! wavelet used for adaptive syncing
        params%g = 1_ik  ! 1 from stencil, 1 from wavelet
    elseif (order == "3_comp_12") then
        params%order_discretization = "FD_3rd_comp_1_2"
        params%wavelet = "CDF40"  ! wavelet used for adaptive syncing
        params%g = 3_ik  ! 2 from stencil, 3 from wavelet
    elseif (order == "4") then
        params%order_discretization = "FD_4th_central"
        params%wavelet = "CDF40"  ! wavelet used for adaptive syncing
        params%g = 3_ik  ! 2 from stencil, 3 from wavelet
    elseif (order == "4op") then
        params%order_discretization = "FD_4th_central_optimized"
        params%wavelet = "CDF40"  ! wavelet used for adaptive syncing
        params%g = 3_ik  ! 3 from stencil, 3 from wavelet
    elseif (order == "4_comp_13") then
        params%order_discretization = "FD_4th_comp_1_3"
        params%wavelet = "CDF40"  ! wavelet used for adaptive syncing
        params%g = 3_ik  ! 3 from stencil, 3 from wavelet
        if (operator == "--dissipation") params%g = 4_ik  ! dissipation uses L=DG with different stencil width
    elseif (order == "4_comp_04") then
        params%order_discretization = "FD_4th_comp_0_4"
        params%wavelet = "CDF40"  ! wavelet used for adaptive syncing
        params%g = 4_ik  ! 4 from stencil, 3 from wavelet
    elseif (order == "5_comp_23") then
        params%order_discretization = "FD_5th_comp_2_3"
        params%wavelet = "CDF60"  ! wavelet used for adaptive syncing
        params%g = 5_ik  ! 3 from stencil, 5 from wavelet
    elseif (order == "6") then
        params%order_discretization = "FD_6th_central"
        params%wavelet = "CDF60"  ! wavelet used for adaptive syncing
        params%g = 5_ik  ! 3 from stencil, 5 from wavelet
    elseif (order == "6_comp_24") then
        params%order_discretization = "FD_6th_comp_2_4"
        params%wavelet = "CDF60"  ! wavelet used for adaptive syncing
        params%g = 5_ik  ! 4 from stencil, 5 from wavelet
        if (operator == "--dissipation") params%g = 6_ik  ! dissipation uses L=DG with different stencil width
    elseif (order == "6_comp_15") then
        params%order_discretization = "FD_6th_comp_1_5"
        params%wavelet = "CDF60"  ! wavelet used for adaptive syncing
        params%g = 5_ik  ! 5 from stencil, 5 from wavelet
        if (operator == "--dissipation") params%g = 6_ik  ! dissipation uses L=DG with different stencil width
    elseif (order == "6_comp_06") then
        params%order_discretization = "FD_6th_comp_0_6"
        params%wavelet = "CDF60"  ! wavelet used for adaptive syncing
        params%g = 6_ik  ! 6 from stencil, 5 from wavelet
    elseif (order == "7_comp_34") then
        params%order_discretization = "FD_7th_comp_3_4"
        params%wavelet = "CDF80"  ! wavelet used for adaptive syncing
        params%g = 7_ik  ! 4 from stencil, 7 from wavelet
    elseif (order == "8") then
        params%order_discretization = "FD_8th_central"
        params%wavelet = "CDF80"  ! wavelet used for adaptive syncing
        params%g = 7_ik  ! 4 from stencil, 7 from wavelet
    elseif (order == "8_comp_35") then
        params%order_discretization = "FD_8th_comp_3_5"
        params%wavelet = "CDF80"  ! wavelet used for adaptive syncing
        params%g = 7_ik  ! 5 from stencil, 7 from wavelet
    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 2 (FD_2nd_central), 4 (FD_4th_central), 6 (FD_6th_central), 4op (FD_4th_central_optimized), 3_comp_12 (FD_3rd_cop_1_2), 4_comp_13 (FD_4th_comp_1_3), 4_comp_04 (FD_4th_comp_0_4), 5_comp_23 (FD_5th_comp_2_3), 6_comp_24 (FD_6th_comp_2_4), 6_comp_15 (FD_6th_comp_1_5), 6_comp_06 (FD_6th_comp_0_6), 7_comp_34 (FD_7th_comp_3_4), 8 (FD_8th_central), 8_comp_35 (FD_8th_comp_3_5)")
    end if

    params%Jmax = tc_length
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
    g  = params%g


    call setup_wavelet(params)

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
        call readHDF5vct_tree( (/file_ux, file_uy, file_uz/), params, hvy_block, tree_ID)
    else
        call readHDF5vct_tree( (/file_ux, file_uy/), params, hvy_block, tree_ID)
    end if

    call sync_ghosts_tree( params, hvy_block, tree_ID)

    ! calculate vorticity from velocities
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)

        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgtID, x0, dx )

        if (operator == "--vorticity") then
            call compute_vorticity(hvy_block(:,:,:,1:params%dim,hvyID), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,:,hvyID))

        elseif (operator=="--vor-abs") then
            call compute_vorticity_abs(hvy_block(:,:,:,1:3,hvyID), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        elseif (operator=="--helicity") then
            call compute_helicity(hvy_block(:,:,:,1:3,hvyID), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1:3,hvyID))

        elseif (operator=="--hel-abs") then
            call compute_helicity_abs(hvy_block(:,:,:,1:3,hvyID), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        elseif (operator == "--divergence") then
            call compute_divergence( hvy_block(:,:,:,1:params%dim,hvyID), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        ! elseif (operator == "--laplace") then
        !     call compute_divergence( hvy_block(:,:,:,1,hvyID), &
        !     hvy_block(:,:,:,2,hvyID), &
        !     hvy_block(:,:,:,3,hvyID),&
        !     dx, Bs, g, &
        !     params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        elseif (operator == "--Q") then
            call compute_Qcriterion( hvy_block(:,:,:,1:params%dim,hvyID), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        elseif (operator == "--copy") then
            hvy_tmp(:,:,:,1,hvyID) = hvy_block(:,:,:,1,hvyID)
        
        elseif (operator == "--dissipation") then
            call compute_dissipation( hvy_block(:,:,:,1:params%dim,hvyID), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        else
            call abort(1812011, "operator is neither --vorticity --vor-abs --divergence --Q")

        endif
    end do

    if (operator == "--vorticity") then
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'vorx', int(time), nint((time-int(time))*1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

        if (params%dim == 3) then
            write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'vory', int(time), nint((time-int(time))*1.0e6_rk)
            call saveHDF5_tree(fname, time, iteration, 2, params, hvy_tmp, tree_ID)
            write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'vorz', int(time), nint((time-int(time))*1.0e6_rk)
            call saveHDF5_tree(fname, time, iteration, 3, params, hvy_tmp, tree_ID)
        end if

    elseif (operator == "--vor-abs") then
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'vorabs', int(time), nint((time-int(time))*1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )
    
    elseif (operator=="--helicity") then
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'helx', int(time), nint((time-int(time))*1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'hely', int(time), nint((time-int(time))*1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 2, params, hvy_tmp, tree_ID )
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'helz', int(time), nint((time-int(time))*1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 3, params, hvy_tmp, tree_ID )
    
    elseif (operator=="--hel-abs") then
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'helabs', int(time), nint((time-int(time))*1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

    elseif (operator=="--divergence") then
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'divu', int(time), nint((time-int(time))*1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

    elseif (operator=="--Q") then
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'Qcrit', int(time), nint((time-int(time))*1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )
    
    elseif (operator=="--dissipation") then
        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'dissipation', int(time), nint((time-int(time))*1.0e6_rk)

        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

    elseif (operator=="--copy") then

        write( fname,'(a, "_", i6.6, i6.6, ".h5")') 'copyUx', int(time), nint((time-int(time))*1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )
    endif
end subroutine compute_vorticity_post
