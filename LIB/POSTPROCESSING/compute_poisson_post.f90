!> \brief postprocessing routine for Poisson-based operations: pressure from velocity, velocity from vorticity
!-----------------------------------------------------------------------------------------------------

subroutine compute_poisson_post(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_poisson
    use module_fft
    use module_nspp

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_in1, file_in2, file_in3, operator
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs
    character(len=cshort)              :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_mask(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID

    character(len=cshort)              :: fname, fname_mask
    logical                            :: read_mask
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: n_hvy_tmp, mpicode

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and operator)
    call get_command_argument(1, operator)
    call get_command_argument(2, file_in1)
    
    ! does the user need help?
    if (file_in1=='--help' .or. file_in1=='--h' .or. file_in1=='-h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit postprocessing: Poisson operations"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Computes derived quantities using Poisson solvers."
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " --press-from-vel"
            write(*, '(A)') "./wabbit-post --press-from-vel source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER]"
            write(*, '(A)') " Computes pressure from velocity field, saves in "
            write(*, '(A)') " pressure_*.h5"
            write(*, '(A)') " Note: Requires physics='NSPP' setup for RHS computation"
            write(*, '(A)') " Optional: --mask=mask_file.h5 to include mask/penalization"
            write(*, '(A)') " Optional: --wavelet=CDF44 to override default wavelet"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " --vel-from-vor"
            write(*, '(A)') "./wabbit-post --vel-from-vor source_vorx.h5 [source_vory.h5] [source_vorz.h5] [ORDER]"
            write(*, '(A)') " Reconstructs velocity from vorticity, saves in "
            write(*, '(A)') " ux_*.h5 [uy_*.h5] [uz_*.h5]"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " ORDER: 2, 3_comp_12, 4, 4op, 4_comp_13, 4_comp_04, 5_comp_23,"
            write(*, '(A)') "        6, 6_comp_24, 6_comp_15, 6_comp_06, 7_comp_34, 8, 8_comp_35"
            write(*, '(A)') " (default: 4)"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(file_in1))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_in1, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! Now read additional files based on operator and dimension
    if (operator == "--vel-from-vor" .and. params%dim == 2) then
        ! 2D vorticity: only one file needed, next argument is ORDER
        call get_command_argument(3, order)
    else
        ! All other cases: need file_in2
        call get_command_argument(3, file_in2)
        call check_file_exists(trim(file_in2))
        
        if (params%dim == 3) then
            call get_command_argument(4, file_in3)
            call check_file_exists(trim(file_in3))
            call get_command_argument(5, order)
        else
            call get_command_argument(4, order)
        end if
    endif

    call get_cmd_arg( "--wavelet", params%wavelet, default="---" )
    call get_cmd_arg( "--mask", fname_mask, default="---" )
    
    ! Check if mask file should be read
    read_mask = (fname_mask /= "---")
    if (read_mask) then
        call check_file_exists(trim(fname_mask))
        if (params%rank == 0) then
            write(*,'(A)') "Reading mask from: "//trim(fname_mask)
        endif
    endif

    ! decide which order for discretization and predictor is used
    if (order == "" .or. order == "4") then
        params%order_discretization = "FD_4th_central"
        params%poisson_order = "FD_4th_central"
        if (params%wavelet == "---") params%wavelet = "CDF40"
    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%poisson_order = "FD_2nd_central"
        if (params%wavelet == "---") params%wavelet = "CDF20"
    elseif (order == "3_comp_12") then
        params%order_discretization = "FD_3rd_comp_1_2"
        params%poisson_order = "FD_3rd_comp_1_2"
        if (params%wavelet == "---") params%wavelet = "CDF40"
    elseif (order == "4op") then
        params%order_discretization = "FD_4th_central_optimized"
        params%poisson_order = "FD_4th_central_optimized"
        if (params%wavelet == "---") params%wavelet = "CDF40"
    elseif (order == "4_comp_13") then
        params%order_discretization = "FD_4th_comp_1_3"
        params%poisson_order = "FD_4th_comp_1_3"
        if (params%wavelet == "---") params%wavelet = "CDF40"
        params%g = 3_ik
    elseif (order == "4_comp_04") then
        params%order_discretization = "FD_4th_comp_0_4"
        params%poisson_order = "FD_4th_comp_0_4"
        if (params%wavelet == "---") params%wavelet = "CDF40"
    elseif (order == "5_comp_14") then
        params%order_discretization = "FD_5th_comp_1_4"
        params%poisson_order = "FD_5th_comp_1_4"
        if (params%wavelet == "---") params%wavelet = "CDF60"
    elseif (order == "5_comp_23") then
        params%order_discretization = "FD_5th_comp_2_3"
        params%poisson_order = "FD_5th_comp_2_3"
        if (params%wavelet == "---") params%wavelet = "CDF60"
    elseif (order == "6") then
        params%order_discretization = "FD_6th_central"
        params%poisson_order = "FD_6th_central"
        if (params%wavelet == "---") params%wavelet = "CDF60"
    elseif (order == "6_comp_24") then
        params%order_discretization = "FD_6th_comp_2_4"
        params%poisson_order = "FD_6th_comp_2_4"
        if (params%wavelet == "---") params%wavelet = "CDF60"
    elseif (order == "6_comp_15") then
        params%order_discretization = "FD_6th_comp_1_5"
        params%poisson_order = "FD_6th_comp_1_5"
        if (params%wavelet == "---") params%wavelet = "CDF60"
    elseif (order == "6_comp_06") then
        params%order_discretization = "FD_6th_comp_0_6"
        params%poisson_order = "FD_6th_comp_0_6"
        if (params%wavelet == "---") params%wavelet = "CDF60"
    elseif (order == "7_comp_34") then
        params%order_discretization = "FD_7th_comp_3_4"
        params%poisson_order = "FD_7th_comp_3_4"
        if (params%wavelet == "---") params%wavelet = "CDF80"
    elseif (order == "8") then
        params%order_discretization = "FD_8th_central"
        params%poisson_order = "FD_8th_central"
        if (params%wavelet == "---") params%wavelet = "CDF80"
    elseif (order == "8_comp_35") then
        params%order_discretization = "FD_8th_comp_3_5"
        params%poisson_order = "FD_8th_comp_3_5"
        if (params%wavelet == "---") params%wavelet = "CDF80"
    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 2, 3_comp_12, 4, 4op, 4_comp_13, 4_comp_04, 5_comp_23, 6, 6_comp_24, 6_comp_15, 6_comp_06, 7_comp_34, 8, 8_comp_35")
    end if

    params%Jmax = tc_length
    params%n_eqn = params%dim+1
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

    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1
    params%N_mask_components = 6

    ! we need to initialize all poisson parameters
    params%poisson_Jmin = 0_ik
    params%poisson_cycle_end_criteria = "tolerance"
    ! params%poisson_cycle_it = 6_ik
    params%poisson_cycle_tol_rel = 1e-8_rk
    params%poisson_cycle_tol_abs = 1.0e-20_rk
    params%poisson_cycle_max_it = 50_ik
    params%poisson_GS_it = 10_ik
    params%poisson_Sync_it = 2_ik
    params%poisson_coarsest = "FFT"
    params%poisson_balanceLoad = 0_ik
    params%FFT_accuracy = "spectral"  ! FD or spectral

    ! for full wavelet operations we need 8/7 for 3D or 4/3 for 2D as blocks, there can be an imbalance thats why we add another factor and for low block amounts we just add a few
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk) * 1.2_rk)+4

    call setup_wavelet(params, params%g)
    call setup_Laplacian_stencils(params, params%g)
    call fft_initialize(params)

    Bs = params%Bs
    g  = params%g

    ! Determine required array sizes based on operator
    if (operator == "--press-from-vel") then
        ! pressure computation needs physics RHS
        params%physics_type = "NSPP"
        params%n_eqn_rhs = params%dim  ! velocity components
        n_hvy_tmp = 2*params%dim  ! sufficient for pressure computation
        
        ! Setup minimal NSPP parameters manually (no INI file in postprocessing)
        call get_cmd_arg_dbl( "--nu", params_nspp%nu, default=0.0e-2_rk )
        call get_cmd_arg_dbl( "--C_eta", params_nspp%C_eta, default=1.0e-3_rk )
        
        ! Initialize essential NSPP module parameters
        params_nspp%dim = params%dim
        params_nspp%discretization = params%order_discretization
        params_nspp%compute_flow = .true.
        params_nspp%n_ghosts = g
        params_nspp%penalization = (fname_mask /= "---")  ! True if mask file provided
        params_nspp%C_eta_temp = params_nspp%C_eta
        params_nspp%domain_size = params%domain_size
        params_nspp%periodic_BC = params%periodic_BC
        params_nspp%symmetry_BC = params%symmetry_BC
        params_nspp%skew_symmetry = .true.

        ! this is quite a hack
        params_nspp%initialized = .true.
        
        ! Get MPI info
        call MPI_COMM_SIZE (WABBIT_COMM, params_nspp%mpisize, mpicode)
        call MPI_COMM_RANK (WABBIT_COMM, params_nspp%mpirank, mpicode)
        
        if (params%rank == 0) then
            write(*,'(A)') "NSPP parameters for pressure computation:"
            write(*,'(A,es12.4)') "  nu = ", params_nspp%nu
            write(*,'(A,L1)')     "  penalization = ", params_nspp%penalization
            if (params_nspp%penalization) write(*,'(A,es12.4)') "  C_eta = ", params_nspp%C_eta
        endif
        
    elseif (operator == "--vel-from-vor") then
        ! velocity from vorticity needs sufficient tmp space
        if (params%dim == 3) then
            n_hvy_tmp = 9  ! Option A: solve all 3 components efficiently
        else
            n_hvy_tmp = 3  ! 2D only needs stream function solve
        endif
        
    else
        call abort(300126, "Unknown operator. Use --press-from-vel or --vel-from-vor")
    endif

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, &
                         neqn_hvy_tmp=n_hvy_tmp, hvy_mask=hvy_mask)

    ! Initialize mask to zero if not reading from file
    if (.not. read_mask) then
        hvy_mask = 0.0_rk
        ! color to 1
        hvy_mask(:,:,:,5,:) = 1.0_rk
    endif

    ! read input data
    if (strings_are_similar(operator, "--press-from-vel")) then
        ! Read velocity components
        if (params%dim == 3) then
            call readHDF5vct_tree( (/file_in1, file_in2, file_in3/), params, hvy_block, tree_ID)
        else
            call readHDF5vct_tree( (/file_in1, file_in2/), params, hvy_block, tree_ID)
        end if
        
        ! Read mask if provided
        if (read_mask) then
            ! read in only actual mask for now, in theory we should be able to read in solid velocity as well
            call readHDF5vct_tree( (/fname_mask/), params, hvy_mask(:,:,:,1:1,:), tree_ID)
        endif
        
    elseif (strings_are_similar(operator, "--vel-from-vor")) then
        ! Read vorticity components
        if (params%dim == 3) then
            call readHDF5vct_tree( (/file_in1, file_in2, file_in3/), params, hvy_block, tree_ID)
        else
            ! 2D: only one vorticity component
            call readHDF5vct_tree( (/file_in1/), params, hvy_block, tree_ID)
        end if
    endif

    ! Perform the computation
    if (strings_are_similar(operator, "--press-from-vel")) then
        !-----------------------------------------------------------------------
        ! Compute pressure from velocity
        !-----------------------------------------------------------------------
        call pressure_from_velocity(params, time, hvy_block, hvy_tmp, hvy_mask, tree_ID)
        
        ! Save pressure field
        write( fname,'(a, "_", a, ".h5")') 'p', &
               trim(adjustl(timestr(time)))
        call saveHDF5_tree(fname, time, iteration, params%dim+1, params, hvy_block, tree_ID)
        
    elseif (strings_are_similar(operator, "--vel-from-vor")) then
        !-----------------------------------------------------------------------
        ! Compute velocity from vorticity
        !-----------------------------------------------------------------------
        call velocity_from_vorticity(params, hvy_block, hvy_tmp, tree_ID)
        
        ! Save velocity components
        write( fname,'(a, "_", a, ".h5")') 'ux', &
               trim(adjustl(timestr(time)))
        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_block, tree_ID)
        
        write( fname,'(a, "_", a, ".h5")') 'uy', &
               trim(adjustl(timestr(time)))
        call saveHDF5_tree(fname, time, iteration, 2, params, hvy_block, tree_ID)
        
        if (params%dim == 3) then
            write( fname,'(a, "_", a, ".h5")') 'uz', &
                   trim(adjustl(timestr(time)))
            call saveHDF5_tree(fname, time, iteration, 3, params, hvy_block, tree_ID)
        endif
    endif

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_mask=hvy_mask)
    
    if (params%rank == 0) then
        write(*,'(A)') "-----------------------------------------------------------"
        write(*,'(A,A,A)') " Successfully completed ", trim(operator), " operation"
        write(*,'(A)') "-----------------------------------------------------------"
    endif

end subroutine compute_poisson_post
