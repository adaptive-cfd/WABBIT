subroutine post_compression_unit_test(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: operator
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: k, lgtID, tc_length, g, j
    integer(kind=ik), dimension(3)     :: Bs

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID

    character(len=cshort)              :: fname
    real(kind=rk) :: norm_L2, norm_Linfty, error_L2, error_Linfty
    real(kind=rk), dimension(3)        :: dx, x0
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: ierr, n_eps, Nb
    real(kind=rk) :: eps(1:51)

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, operator)
    ! does the user need help?
    if (operator=='--help' .or. operator=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: wavelet compression test"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " For every value of C_eps:"
            write(*,*) " "
            write(*,*) " 1) create the reference data (a bunch of Gauss blobs) on an equidistant grid"
            write(*,*) " 2) compress the grid using the thresholding rule"
            write(*,*) " 3) re-refine the mesh to its original resolution"
            write(*,*) " 4) compute the error by comparing to the analytical solution"
            write(*,*) " "
            write(*,*) " Notes:"
            write(*,*) " This routine is somehow redundant with the python compression"
            write(*,*) " test, but it is much faster because we do not involve I/O."
            write(*,*) " Works in 2D and 3D"
            write(*,*) " "
            write(*,*) " Parameters:"
            write(*,*) " --memory"            
            write(*,*) " --wavelet"
            write(*,*) " --eps-normalized"
            write(*,*) " --eps-norm"
            write(*,*) " --Bs"
            write(*,*) " --Jmax"
            write(*,*) " --dim"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif

    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )
    call get_cmd_arg( "--eps-normalized", params%eps_normalized, default=.true. )
    call get_cmd_arg( "--eps-norm", params%eps_norm, default="Linfty" )
    call get_cmd_arg( "--Bs", params%Bs(1), default=18 )
    call get_cmd_arg( "--Bs", params%Bs(2), default=18 )
    call get_cmd_arg( "--Bs", params%Bs(3), default=18 )
    call get_cmd_arg( "--Jmax", params%Jmax, default=4 )
    call get_cmd_arg( "--dim", params%dim, default=2 )

    ! fixed parameters
    params%domain_size = 2.0_rk
    params%n_eqn = 1
    params%forest_size = 1
    params%physics_type = "postprocessing"
    params%coarsening_indicator = "threshold-state-vector"
    params%force_maxlevel_dealiasing = .false.
    params%Jmin = 1

    call setup_wavelet(params, params%g)
    
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    allocate(params%threshold_state_vector_component(1:params%n_eqn))

    params%threshold_state_vector_component = .true.
    
    Bs = params%Bs
    g  = params%g
    
    
    eps(1:51) = (/1.00000000e-10, 1.58489319e-10, 2.51188643e-10, 3.98107171e-10, &
    6.30957344e-10, 1.00000000e-09, 1.58489319e-09, 2.51188643e-09, &
    3.98107171e-09, 6.30957344e-09, 1.00000000e-08, 1.58489319e-08, &
    2.51188643e-08, 3.98107171e-08, 6.30957344e-08, 1.00000000e-07, &
    1.58489319e-07, 2.51188643e-07, 3.98107171e-07, 6.30957344e-07, &
    1.00000000e-06, 1.58489319e-06, 2.51188643e-06, 3.98107171e-06, &
    6.30957344e-06, 1.00000000e-05, 1.58489319e-05, 2.51188643e-05, &
    3.98107171e-05, 6.30957344e-05, 1.00000000e-04, 1.58489319e-04, &
    2.51188643e-04, 3.98107171e-04, 6.30957344e-04, 1.00000000e-03, &
    1.58489319e-03, 2.51188643e-03, 3.98107171e-03, 6.30957344e-03, &
    1.00000000e-02, 1.58489319e-02, 2.51188643e-02, 3.98107171e-02, &
    6.30957344e-02, 1.00000000e-01, 1.58489319e-01, 2.51188643e-01, &
    3.98107171e-01, 6.30957344e-01, 1.00000000e+00/)
    n_eps = size(eps, dim=1)

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=params%n_eqn)
    
    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    if (params%rank==0) then
        open(14, file="error.csv", status='replace')
        write (14,'((A,";",A,";",A,";",A))') "eps", "error_L2", "error_Linfty", "NB_compressed"
        close(14)

        open(14, file="header.txt", status='replace')
        write (14,'("Bs=",i3," Jmax=",i2," dim=",i1," eps-norm=",A," eps-normalized=",L)') Bs(1), params%Jmax, params%dim, &
        adjustl(trim(params%eps_norm)), params%eps_normalized
        close(14)
    endif

    !--------------------------------------------------------------------------------------------
    ! actual routine follows
    !--------------------------------------------------------------------------------------------


    do j = 1, n_eps
        params%eps = eps(j)
        
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        ! (1) create the data on the finest, equidistant grid
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        call createEquidistantGrid_tree( params, hvy_block, params%Jmax, .true., tree_ID )

        do k = 1, hvy_n(tree_ID)
            hvyID = hvy_active(k,tree_ID)    
            call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
            call get_block_spacing_origin( params, lgtID, x0, dx )
    
            call set_block_testing_data(params, hvy_block(:,:,:,1,hvyID), x0, dx)    
        end do

        if (j==1) then
            call saveHDF5_tree("input_0000.h5", 0.0_rk, 1_ik, 1, params, hvy_block, tree_ID)
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        ! (2) compress it with the wavelet and threshold paras%eps
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID))
        call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp)
        
        Nb = lgt_n(tree_ID)
        
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        ! (3) refine if back to the original resolution
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID))
        call refineToEquidistant_tree(params, hvy_block, hvy_tmp, tree_ID, params%Jmax)

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        ! (4) compute the error (knowing the analytical solution - this is very fast)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        error_L2     = 0.0_rk
        error_Linfty = 0.0_rk
        norm_L2      = 0.0_rk
        norm_Linfty  = 0.0_rk

        do k = 1, hvy_n(tree_ID)
            hvyID = hvy_active(k,tree_ID)    
            call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
            call get_block_spacing_origin( params, lgtID, x0, dx )
    
            ! reference data in hvy_tmp
            call set_block_testing_data(params, hvy_tmp(:,:,:,1,hvyID), x0, dx)    

            if (params%dim == 3) then
                error_L2     = error_L2 + sum( (hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g,1,hvyID)-hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g,1,hvyID))**2 )
                error_Linfty = max( error_Linfty, maxval(abs(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g,1,hvyID)-hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g,1,hvyID))) )

                norm_L2     = norm_L2 + sum( (hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g,1,hvyID))**2 )
                norm_Linfty = max( norm_Linfty, maxval(abs(hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g,1,hvyID))) )
            else
                error_L2     = error_L2 + sum( (hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1,1,hvyID)-hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, 1,1,hvyID))**2 )
                error_Linfty = max( error_Linfty, maxval(abs(hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1,1,hvyID)-hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, 1,1,hvyID))) )

                norm_L2     = norm_L2 + sum( (hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, 1,1,hvyID))**2 )
                norm_Linfty = max( norm_Linfty, maxval(abs(hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, 1,1,hvyID))) )
            endif
        end do

        call MPI_Allreduce(MPI_IN_PLACE, error_L2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, norm_L2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        error_L2 = sqrt(error_L2) / sqrt(norm_L2)

        call MPI_Allreduce(MPI_IN_PLACE, error_Linfty, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, norm_Linfty, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, ierr)
        error_Linfty = error_Linfty / norm_Linfty

        if (params%rank==0) then
            write(*,*) "Nb=", Nb
            write(*,'("eps=",es12.4," err_L2=",es12.4," err_Linfty=",es12.4,1x,A)') params%eps, error_L2, error_Linfty, params%wavelet

            open(14,file='error.csv',status='unknown',position='append')
            write(14,'(es15.8,";",es15.8,";",es15.8,";",i9)') params%eps, error_L2, error_Linfty, Nb
            close(14)
        endif
    enddo

end subroutine 