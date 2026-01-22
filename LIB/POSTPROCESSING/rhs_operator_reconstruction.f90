subroutine rhs_operator_reconstruction(params)
    use mpi
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_acm
    use module_time_step
    use module_ini_files_parser_mpi
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: file, mode, OPERATOR
    real(kind=rk) :: time, x, y, dx_fine, u_dx, u_dxdx, dx_inv, val, x2, y2, nu
    integer(kind=ik) :: iteration, k, tc_length, iblock, ix, iy, &
    g, iz, level,tree_ID_tmp,tree_ID_rhs_u,tree_ID_rhs_u_ei
    integer(kind=ik) :: ixx, iyy, ix2, iy2, nx_fine, ixx2,iyy2, n_nonzero
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :)
    integer :: hvy_id, lgt_id, fsize, j, tree_ID_u, tree_ID_ei
    character(len=cshort)              :: fname, fname_ini
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    character(len=1) :: dir
    logical, parameter :: verbosity=.False.
    integer(kind=ik) :: Nlgtn,Nhvyn
    integer(kind=ik) , save, allocatable :: lgt_active_ref(:,:), lgt_block_ref(:,:)
    integer(kind=ik) , save :: lgt_n_ref(2)=0_ik
    type(inifile) :: ini_file
    real(kind=rk) :: dt,CFL_number
    logical :: error_OOM

    ! following are to suppress warnings
    character(len=80) :: def_wavelet, def_operator
    logical :: def_adapt

    if (params%number_procs>1) call abort(2205121, "OperatorReconstruction is a serial routine...")

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, mode)

    ! does the user need help?
    if (mode=='--help' .or. mode=='--h' .or. mode=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --OP-rhs u_001.h5 inifile.ini --memory=10GB"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " This function computes the linearised RHS operator A (Ngrid x Ngrid)"
            write(*,*) " arround the given statevector (here u_001)!"
            write(*,*) " Mathematically speaking:"
            write(*,*) " A = [r_1, r_2, ..., r_Ngrid], where r_i = rhs(u_001 + e_i) - rhs(u_001)"
            write(*,*) " and e_i is the i-th standard basis vector of the corresponding grid"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    ! read ini parameter for physics rhs
    call get_command_argument( 3, fname_ini )
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, fname_ini )
    call init_physics_modules( params, fname_ini, params%N_mask_components )

    call get_command_argument(2, file)
    call check_file_exists(file)

    ! get some parameters from the grid file
    call read_attributes(file, Nlgtn, time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    params%Jmax = tc_length + 2
    params%forest_size = 10
    fsize = params%forest_size
    params%domain_size = domain
    params%Bs = Bs
    ! call read_ini_file_mpi(ini_file, fname_ini, .true.)
    ! call read_param_mpi(ini_file, 'ConvectionDiffusion', 'nu', nu )
    !---------------------------------------------------------------------------
    ! Adjustable PARAMETERS
    !---------------------------------------------------------------------------
    !OPERATOR = "RHS" ! decide between RHS or Identity
    OPERATOR = "EVOLVE" ! decide between RHS or Identity
    ! it gives warning when using same values as in_out and default so lets transcribe first
    def_wavelet = params%wavelet
    def_operator = OPERATOR
    def_adapt = params%adapt_tree
    call get_cmd_arg( "--wavelet", params%wavelet, default= def_wavelet )    ! if not set we take the parameter from the ini file
    call get_cmd_arg( "--operator", OPERATOR, default= def_operator )               ! if not set we take the parameter from the ini file
    call get_cmd_arg( "--adapt", params%adapt_tree, default= def_adapt )! if not set we take the parameter from the ini file
    CFL_number = params%CFL
    !params%adapt_tree = .False.
    !---------------------------------------------------------------------------

    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, hvy_block, hvy_work,hvy_mask=hvy_mask, hvy_tmp=hvy_tmp)
    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes(params)

    call reset_forest(params)
    lgt_n = 0 ! reset number of active light blocks
    hvy_n = 0
    tree_n = 0 ! reset number of trees in forest




    open(17, file=trim(adjustl(file))//'.info.txt', status='replace')
    if (params%adapt_tree) then
        write(17,'(A,1x,A,1x,A,1x,A,1x,A,1x,A," g=",i1," Bs=",i2," nu=",e12.4," CFL=",e12.4)') trim(params%order_discretization), &
        trim(params%order_predictor), " ",trim(params%wavelet)," ",trim(OPERATOR), &
        params%g, params%Bs(1), nu, &
        CFL_number
    else
        write(17,'(A,1x,A,1x,A," g=",i1," Bs=",i2," nu=",e12.4," CFL=",e12.4)') trim(params%order_discretization), &
        " static grid ",trim(OPERATOR), params%g, params%Bs(1), nu, CFL_number
    endif
    close(17)

    !---------------------------------------------------------------------------

    if ((params%order_discretization == "FD_4th_central_optimized").and.(params%g<4)) then
        call abort(33,"not enough g")
    endif


    g = params%g
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    tree_ID_u = 1
    call readHDF5vct_tree((/file/), params, hvy_block, tree_ID_u, verbosity=.true.) ! read both grid and data in this tree

    ! this hack ensures, on mono_CPU, that later on, refine+coarsening always ends up in the same order in hvy_actve
    if (params%adapt_tree) then
        call sync_ghosts_tree(params, hvy_block, tree_ID_u)
        call refine_tree( params, hvy_block, "everywhere", tree_ID=tree_ID_u, error_OOM=error_OOM )
        call sync_ghosts_tree(params, hvy_block, tree_ID_u)
        call adapt_tree( time, params, hvy_block, tree_ID_u, "everywhere", hvy_tmp )
    endif

    !---------------------------------------------------------------------------
    tree_ID_rhs_u = 2
    call readHDF5vct_tree( (/file/), params, hvy_block, tree_ID_rhs_u)
    call multiply_tree_with_scalar(params, hvy_block, tree_ID_rhs_u, 0.0_rk, verbosity) ! effectively, we just read the grid in this tree
    !---------------------------------------------------------------------------
    tree_ID_ei = 3
    call readHDF5vct_tree( (/file/), params, hvy_block, tree_ID_ei) ! effectively, we just read the grid in this tree
    call multiply_tree_with_scalar(params, hvy_block, tree_ID_ei, 0.0_rk, verbosity)
    !---------------------------------------------------------------------------
    tree_ID_rhs_u_ei = 4
    call readHDF5vct_tree( (/file/), params, hvy_block, tree_ID_rhs_u_ei) ! effectively, we just read the grid in this tree
    call multiply_tree_with_scalar(params, hvy_block, tree_ID_rhs_u_ei, 0.0_rk, verbosity)

    dx_fine = (2.0_rk**-maxActiveLevel_tree(tree_ID_u))*domain(2)/real((Bs(2)), kind=rk)
    nx_fine = nint(domain(2)/dx_fine)

    write(*,*) "nx_fine=", nx_fine
    write(*,*) "n_ghost=", g
    write(*,*) "nblocks=", lgt_n, "bs=", Bs, "npoints (op. matrix size!)=", lgt_n*bs(1)*bs(2)

    !---------------------------------------------------------------------------
    ! save the grid (for plotting in python)
    !---------------------------------------------------------------------------
    open(19, file=trim(adjustl(file))//'.operator_grid_points.txt', status='replace')
    do iblock = 1, hvy_n(tree_ID_u)
        call hvy2lgt(lgt_id, hvy_active(iblock,tree_ID_u), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )
        level = lgt_block(lgt_id, IDX_MESH_LVL)

        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                ixx = nint(x/dx_fine)+1
                iyy = nint(y/dx_fine)+1

                write(19,*) ixx, iyy, x, y, level
            enddo
        enddo
    enddo
    close(19)

    !---------------------------------------------------------------------------
    ! compute operator matrix
    !---------------------------------------------------------------------------
    open(17, file=trim(adjustl(file))//'.operator_matrix.txt', status='replace')

    ! we now calculate 1/h[rhs(u + h* e_i) - rhs(u)]1/h[rhs(u + h* e_i) - rhs(u)]

    !---------------------------------------------------------------------------
    ! store_reference_mesh
    !---------------------------------------------------------------------------
    call copy_tree(params, hvy_block, tree_ID_rhs_u, tree_ID_u)

    if ( params%adapt_tree ) then
        tree_ID_tmp = tree_ID_rhs_u

        call store_ref_meshes( lgt_block_ref, lgt_active_ref, lgt_n_ref, tree_ID_u, tree_ID_ei)
        Nhvyn = hvy_n(tree_ID_u)
        !---------------------------------------------------------------------------
        ! refine grid once
        !
        ! The total operator we investigate here is
        !
        ! R * E * C
        !
        ! and this refinement step is the first "R"
        !---------------------------------------------------------------------------

        call sync_ghosts_tree( params, hvy_block, tree_ID_tmp )
        call refine_tree( params, hvy_block, "everywhere", tree_ID=tree_ID_tmp, error_OOM=error_OOM )

        call sync_ghosts_tree(params, hvy_block, tree_ID_tmp)
    endif

    ! -----------------------------------------------
    ! Operator "E" from R * E * C
    ! -----------------------------------------------
    hvy_work(:,:,:,:,:,1) = 0
    if (OPERATOR == "RHS") then
        call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, tree_ID_rhs_u)
    else if( OPERATOR == "EVOLVE") then
        iteration = 0
        time = 0.0_rk
        dt = 1
        call timeStep_tree(time, dt, iteration, params, hvy_block, hvy_work, hvy_mask, hvy_tmp, tree_ID_rhs_u)
        do k = 1, hvy_n(tree_ID_rhs_u)
            hvy_id = hvy_active(k,tree_ID_rhs_u)
            hvy_work(:,:,:,:,hvy_id,1) = hvy_block(:,:,:,:,hvy_id)
        enddo
    else !! this operator is used for the identity
        do k = 1, hvy_n(tree_ID_rhs_u)
            hvy_id = hvy_active(k,tree_ID_rhs_u)
            hvy_work(:,:,:,:,hvy_id,1) = hvy_block(:,:,:,:,hvy_id)
        enddo
    endif
    !!!! only 2D!!!!
    iz = 1

    do iblock = 1, hvy_n(tree_ID_u)
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                !write(*,*) "--------------------point---------------------------"
                !---------------------------------------------------------------
                ! reset entire grid to zeros (do not care about performance, just reset all)
                call copy_tree(params, hvy_block, tree_ID_ei, tree_ID_u)
                !---------------------------------------------------------------
                ! set this one point we're looking at to 1
                call hvy2lgt(lgt_id, hvy_active(iblock,tree_ID_ei), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, x0, dx )

                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
                if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk

                ! save its indices on the fine grid
                ixx = nint(x/dx_fine)+1
                iyy = nint(y/dx_fine)+1


                ! set the one
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                hvy_block(ix, iy, iz, 1, hvy_active(iblock,tree_ID_ei)) = hvy_block(ix, iy, iz, 1, hvy_active(iblock,tree_ID_ei)) + 1.0_rk
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if ( params%adapt_tree ) then
                    call sync_ghosts_tree( params, hvy_block, tree_ID_ei )
                    call refine_tree( params, hvy_block, "everywhere", tree_ID=tree_ID_ei, error_OOM=error_OOM )
                endif

                call updateMetadata_tree(params, tree_ID_ei)
                !---------------------------------------------------------------
                ! synchronize ghosts (important! if e.g. coarseWins is actiyve and you happen to set the redundant value of a refined block, its overwritten to be zero again)
                ! Note: this also applies to coarse block bordering on a coarse block, if its ID is lower.
                ! In fact, each point is then computed only once. Note: if you set the point on a high lgt_id, then
                ! it will be "sync'ed down" to lower light IDs, so you can find the point more than once
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call sync_ghosts_tree(params, hvy_block, tree_ID_ei)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Now compute rhs(u+he_i) of the derivative 1/h[rhs(u + h* e_i) - rhs(u)]
                ! call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
                !     lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor, tree_ID_ei)

                if (OPERATOR == "RHS") then
                    call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, tree_ID_ei)
                else if( OPERATOR == "EVOLVE") then
                    iteration = 0
                    time = 0
                    call timeStep_tree(time, dt, iteration, params, hvy_block, hvy_work, hvy_mask, hvy_tmp, tree_ID_ei)
                    do k = 1, hvy_n(tree_ID_ei)! call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
                        hvy_work(:,:,:,:,hvy_active(k,tree_ID_ei),1) = hvy_block(:,:,:,:,hvy_active(k,tree_ID_ei))
                    enddo
                else !! this operator is used for the identity
                    do k = 1, hvy_n(tree_ID_ei)! call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
                        hvy_work(:,:,:,:,hvy_active(k,tree_ID_ei),1) = hvy_block(:,:,:,:,hvy_active(k,tree_ID_ei))
                    enddo
                endif
                ! This second sync step also synchronizes the derivative we computed previously
                ! Note that on a coarse/fine interface, wabbit computes two values for the derivative
                ! on the coarse and fine level. This synchronizing step lets us keep only either of those,
                ! depending on fineWins or coarseWins
                call sync_ghosts_tree(params, hvy_work(:,:,:,:,:,1), tree_ID_ei)
                call sync_ghosts_tree(params, hvy_work(:,:,:,:,:,1), tree_ID_rhs_u)

                call delete_tree(params, tree_ID_rhs_u_ei)
                call substract_two_trees(params, hvy_work(:,:,:,:,:,1), hvy_tmp, tree_ID1=tree_ID_ei, tree_ID2=tree_ID_rhs_u, dest_tree_ID=tree_ID_rhs_u_ei)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !call sync_ghosts(params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID_rhs), hvy_n(tree_ID_rhs))
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !write(*,*) "--------------------point-3s--------------------------"

                if ( params%adapt_tree ) then
                    call updateMetadata_tree(params, tree_ID_rhs_u_ei)

                    call sync_ghosts_tree(params, hvy_work(:,:,:,:,:,1), tree_ID_rhs_u_ei)
                    ! call adapt_tree( time, params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, lgt_active(:,tree_ID_rhs_u_ei), &
                    !  lgt_n(tree_ID_rhs_u_ei), lgt_sortednumlist(:,:,tree_ID_rhs_u_ei), hvy_active(:,tree_ID_rhs_u_ei), hvy_n(tree_ID_rhs_u_ei), tree_ID_rhs_u_ei, "everywhere", hvy_tmp, external_loop=.true. )

                    call coarse_tree_2_reference_mesh(params, lgt_block_ref, lgt_active_ref(:,1), lgt_n_ref(1), &
                    hvy_work(:,:,:,:,:,1), hvy_tmp, tree_ID_rhs_u_ei, verbosity)

                    call updateMetadata_tree(params, tree_ID_rhs_u_ei)

                    call sync_ghosts_tree(params, hvy_work(:,:,:,:,:,1), tree_ID_rhs_u_ei)
                endif
                ! save operator line to text file.
                ! note: unfortunately, we use the index on the finest level, i.e., we temporarily
                ! create a matrix N_max**2 by N_max**2, where N_max is Npoints on the finest level.
                ! Many points do not exist; they are on coarse levels. however, this is a problem
                ! for the python script, because it first reads the entire matrix, then removes zero cols/rows.
                do k = 1, hvy_n(tree_ID_rhs_u_ei)
                    call hvy2lgt(lgt_id, hvy_active(k,tree_ID_rhs_u_ei), params%rank, params%number_blocks)
                    call get_block_spacing_origin( params, lgt_id, x0, dx )

                    do iy2 = g+1, Bs(2)+g
                        do ix2 = g+1, Bs(1)+g
                            x = dble(ix2-(g+1)) * dx(1) + x0(1)
                            y = dble(iy2-(g+1)) * dx(2) + x0(2)

                            if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
                            if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk

                            ixx2 = nint(x/dx_fine)+1
                            iyy2 = nint(y/dx_fine)+1

                            val = hvy_work(ix2, iy2, iz, 1, hvy_active(k,tree_ID_rhs_u_ei),1) ! u_dx

                            if (abs(val) > 1.0e-13) then
                                ! this point is a nonzero value
                                write(17,'(i6,1x,i6,1x,es15.8)') ixx+(iyy-1)*nx_fine, ixx2+(iyy2-1)*nx_fine, val
                                write(*,*) "python col=", ixx2+(iyy2-1)*nx_fine -1 , "row=", ixx+(iyy-1)*nx_fine -1, val, "xy=",x,y, "block", k
                            endif
                            hvy_work(ix2, iy2, iz, 1, hvy_active(k,tree_ID_rhs_u_ei),1) = 0
                        end do
                    end do
                end do



            enddo
        enddo
    enddo
    close(17)

end subroutine
