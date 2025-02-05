!> \brief postprocessing routine for derivative calculation from a datafield saved in a .h5 file
!-----------------------------------------------------------------------------------------------------

subroutine post_derivative(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file, operator
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs
    character(len=2)                   :: order, der_dim, der_order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID, der_dim_i, der_order_i

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
    call get_command_argument(2, file)
    ! does the user need help?
    if (file=='--help' .or. file=='--h' .or. file=='-h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit postprocessing: Derivative"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Computes a derivative quantity from a file. Output is stored"
            write(*, '(A)') " in predefined files."
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') "./wabbit-post --derivative source.h5 [DER-DIM] [DER-ORDER] [CON-ORDER]"
            write(*, '(A)') " Computes first order derivative in one direction."
            write(*, '(A)') "    der_dim = 1, 2 or 3     for derivative in x-, y- or z-direction"
            write(*, '(A)') "    der_order = 1 or 2      for first or second order derivative"
            write(*, '(A)') "    con_order = 2, 4 or 6   as convergence order of the stencil"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(file))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    call get_command_argument(3, der_dim)
    call get_command_argument(4, der_order)
    call get_command_argument(5, order)


    ! decide which order for discretization and predictor is used. Note predictor
    ! is used in ghost nodes sync'ing
    if (order == "4") then
        params%order_discretization = "FD_4th_central"
        params%g = 4_ik
        params%wavelet="CDF40"
    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%g = 2_ik
        params%wavelet="CDF20"
    elseif (order == "6") then
        params%order_discretization = "FD_6th_central"
        params%g = 7_ik
        params%wavelet="CDF60"
    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 6 (FD_6th_central), 4 (FD_4th_central) and 2 (FD_2nd_central)")
    end if

    ! translate der_order and der_dim and check if everything is alright
    if (der_dim == "1") then
        der_dim_i = 1
    elseif (der_dim == "2") then
        der_dim_i = 2
    elseif (der_dim == "3" .and. params%dim==3) then
        der_dim_i = 3
    else
        call abort(250124, "Cannot compute derivative in this direction, choose between 1 for x-, 2 for y- and 3 for z-direction (3 only for 3D)")
    endif
    if (der_order == "1") then
        der_order_i = 1
    elseif (der_order == "2") then
        der_order_i = 2
    else
        call abort(250124, "Cannot compute this derivative, choose between 1 for first order and 2 for second order derivative")
    endif

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

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=nwork)

    ! read input data
    call readHDF5vct_tree( (/file/), params, hvy_block, tree_ID)

    call sync_ghosts_tree( params, hvy_block, tree_ID)

    ! calculate vorticity from velocities
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)

        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgtID, x0, dx )

            call compute_derivative( hvy_block(:,:,:,1,hvyID), dx, Bs, g, der_dim_i, der_order_i, params%order_discretization, hvy_tmp(:,:,:,1,hvyID))            
    end do

    if (der_order_i == 1 .and. der_dim_i == 1) write( fname,'(a, "_", i12.12, ".h5")') 'dx', nint(time * 1.0e6_rk)
    if (der_order_i == 1 .and. der_dim_i == 2) write( fname,'(a, "_", i12.12, ".h5")') 'dy', nint(time * 1.0e6_rk)
    if (der_order_i == 1 .and. der_dim_i == 3) write( fname,'(a, "_", i12.12, ".h5")') 'dz', nint(time * 1.0e6_rk)
    if (der_order_i == 2 .and. der_dim_i == 1) write( fname,'(a, "_", i12.12, ".h5")') 'dxdx', nint(time * 1.0e6_rk)
    if (der_order_i == 2 .and. der_dim_i == 2) write( fname,'(a, "_", i12.12, ".h5")') 'dydy', nint(time * 1.0e6_rk)
    if (der_order_i == 2 .and. der_dim_i == 3) write( fname,'(a, "_", i12.12, ".h5")') 'dzdz', nint(time * 1.0e6_rk)

    call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

end subroutine post_derivative