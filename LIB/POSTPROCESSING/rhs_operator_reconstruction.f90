subroutine rhs_operator_reconstruction(params)
    use mpi
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_acm
    use module_time_step
    use module_ini_files_parser_mpi

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: file, mode, OPERATOR
    real(kind=rk) :: time, x, y, dx_fine, u_dx, u_dxdx, dx_inv, val, x2, y2, nu
    integer(kind=ik) :: iteration, k, tc_length, tree_N, iblock, ix, iy, &
    g, iz, a1, b1, a2, b2, level,tree_id_tmp,tree_id_rhs_u,tree_id_rhs_u_ei
    integer(kind=ik) :: ixx, iyy, ix2, iy2, nx_fine, ixx2,iyy2, n_nonzero
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), allocatable      :: lgt_n(:), hvy_n(:)
    integer :: hvy_id, lgt_id, fsize, j, tree_id_u, tree_id_ei
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


    params%max_treelevel = tc_length + 2
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
    call get_cmd_arg( "--wavelet", params%wavelet, default= params%wavelet )    ! if not set we take the parameter from the ini file
    call get_cmd_arg( "--operator", OPERATOR, default= OPERATOR )               ! if not set we take the parameter from the ini file
    call get_cmd_arg( "--adapt", params%adapt_mesh, default= params%adapt_mesh )! if not set we take the parameter from the ini file
    CFL_number = params%CFL
    !params%adapt_mesh = .False.
    params%iter_ghosts = .false.
    !---------------------------------------------------------------------------

    if (params%wavelet=="CDF44") then
        params%wavelet_transform_type = "biorthogonal"
        params%n_ghosts = 6_ik
    else
      params%wavelet_transform_type = "harten-multiresolution"
      params%n_ghosts = 4_ik
    endif

    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        hvy_active, lgt_sortednumlist, hvy_work,hvy_mask=hvy_mask, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)
    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    call reset_forest(params, lgt_block, lgt_active, lgt_n,hvy_active, hvy_n, &
    lgt_sortednumlist,tree_n)
    lgt_n = 0 ! reset number of active light blocks
    hvy_n = 0
    tree_n = 0 ! reset number of trees in forest




    open(17, file=trim(adjustl(file))//'.info.txt', status='replace')
    if (params%adapt_mesh) then
      write(17,'(A,1x,A,1x,A,1x,A,1x,A,1x,A," g=",i1," Bs=",i2, " coarseWins=",L1," nu=",e12.4," CFL=",e12.4)') trim(params%order_discretization), &
      trim(params%order_predictor), " ",trim(params%wavelet)," ",trim(OPERATOR), &
      params%n_ghosts, params%Bs(1), params%ghost_nodes_redundant_point_coarseWins, nu, &
      CFL_number
    else
      write(17,'(A,1x,A,1x,A," g=",i1," Bs=",i2," nu=",e12.4," CFL=",e12.4)') trim(params%order_discretization), &
      " static grid ",trim(OPERATOR), params%n_ghosts, params%Bs(1), nu, CFL_number
    endif
    close(17)

    !---------------------------------------------------------------------------

    if ((params%order_discretization == "FD_4th_central_optimized").and.(params%n_ghosts<4)) then
        call abort(33,"not enough g")
    endif


    g = params%n_ghosts
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    tree_id_u = 1
    call read_field2tree(params, (/file/), 1, tree_id_u, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    call create_active_and_sorted_lists(params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_id_u)
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_u),&
    lgt_n(tree_id_u), lgt_sortednumlist(:,:,tree_id_u), hvy_active(:,tree_id_u), hvy_n(tree_id_u) )

    ! this hack ensures, on mono_CPU, that later on, refine+coarsening always ends up in the same order in hvy_actve
    if (params%adapt_mesh) then
      call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_u), hvy_n(tree_id_u))
      call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_id_u), lgt_n(tree_id_u), &
      lgt_sortednumlist(:,:,tree_id_u), hvy_active(:,tree_id_u), hvy_n(tree_id_u), "everywhere", tree_ID=tree_id_u )

        call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_u), hvy_n(tree_id_u))
        call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_id_u), &
        lgt_n(tree_id_u), lgt_sortednumlist(:,:,tree_id_u), hvy_active(:,tree_id_u), hvy_n(tree_id_u), tree_id_u, "everywhere", hvy_tmp, external_loop=.true. )

    endif

    !---------------------------------------------------------------------------
    tree_id_rhs_u = 2
    call read_field2tree(params, (/file/), 1, tree_id_rhs_u, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    call create_active_and_sorted_lists(params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_id_rhs_u)
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_rhs_u),&
    lgt_n(tree_id_rhs_u), lgt_sortednumlist(:,:,tree_id_rhs_u), hvy_active(:,tree_id_rhs_u), hvy_n(tree_id_rhs_u) )
    call multiply_tree_with_scalar(params, hvy_block, hvy_active, hvy_n, tree_id_rhs_u, &
        0.0_rk, verbosity)
    !---------------------------------------------------------------------------
    tree_id_ei = 3
    call read_field2tree(params, (/file/), 1, tree_id_ei, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    call create_active_and_sorted_lists(params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_id_ei)
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_ei),&
    lgt_n(tree_id_ei), lgt_sortednumlist(:,:,tree_id_ei), hvy_active(:,tree_id_ei), hvy_n(tree_id_ei) )
    call multiply_tree_with_scalar(params, hvy_block, hvy_active, hvy_n, tree_id_ei, &
        0.0_rk, verbosity)
    !---------------------------------------------------------------------------
    tree_id_rhs_u_ei = 4
    call read_field2tree(params, (/file/), 1, tree_id_rhs_u_ei, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    call create_active_and_sorted_lists(params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_id_rhs_u_ei)
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_rhs_u_ei),&
    lgt_n(tree_id_rhs_u_ei), lgt_sortednumlist(:,:,tree_id_rhs_u_ei), hvy_active(:,tree_id_rhs_u_ei), hvy_n(tree_id_rhs_u_ei) )
    call multiply_tree_with_scalar(params, hvy_block, hvy_active, hvy_n, tree_id_rhs_u_ei, &
        0.0_rk, verbosity)

    dx_fine = (2.0_rk**-max_active_level(lgt_block, lgt_active(:,tree_id_u), lgt_n(tree_id_u)))*domain(2)/real((Bs(2)-1), kind=rk)
    nx_fine = nint(domain(2)/dx_fine)
    write(*,*) "nx_fine=", nx_fine
    write(*,*) "n_ghost=", g
    write(*,*) "nblocks=", lgt_n, "bs=", Bs, "npoints (op. matrix size!)=", lgt_n*bs(1)*bs(2)

    !---------------------------------------------------------------------------
    ! save the grid (for plotting in python)
    !---------------------------------------------------------------------------
    open(19, file=trim(adjustl(file))//'.operator_grid_points.txt', status='replace')
    do iblock = 1, hvy_n(tree_id_u)
        call hvy2lgt(lgt_id, hvy_active(iblock,tree_id_u), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)
        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
            do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
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
      ! store_rewrite(*,*) "hvyn old/new",Nhvyn, hvy_n(tree_id_rhs)ference_mesh
      !---------------------------------------------------------------------------
    call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
          hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_rhs_u, tree_id_u)

    if ( params%adapt_mesh ) then
      tree_id_tmp = tree_id_rhs_u

      call store_ref_meshes(lgt_block,     lgt_active,     lgt_n,  &
                lgt_block_ref, lgt_active_ref, lgt_n_ref, tree_id_u, tree_id_ei)
      Nhvyn = hvy_n(tree_id_u)
      !---------------------------------------------------------------------------
      ! refine grid ones
      !---------------------------------------------------------------------------
      call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_tmp), hvy_n(tree_id_tmp) )

      call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_id_tmp), lgt_n(tree_id_tmp), &
          lgt_sortednumlist(:,:,tree_id_tmp), hvy_active(:,tree_id_tmp), hvy_n(tree_id_tmp), "everywhere", tree_ID=tree_id_tmp )
      call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_tmp), hvy_n(tree_id_tmp))
    endif

    hvy_work(:,:,:,:,:,1) = 0
    if (OPERATOR == "RHS") then
      call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
         lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor, tree_id_rhs_u)
    else if( OPERATOR == "EVOLVE") then
      iteration = 0
      time = 0.0_rk
      dt = 1
      call time_stepper(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, hvy_mask, hvy_tmp, &
          hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist, tree_id_rhs_u)
      do k = 1, hvy_n(tree_id_rhs_u)
                hvy_id = hvy_active(k,tree_id_rhs_u)
                hvy_work(:,:,:,:,hvy_id,1) = hvy_block(:,:,:,:,hvy_id)
      enddo
    else !! this operator is used for the identity
      do k = 1, hvy_n(tree_id_rhs_u)
            hvy_id = hvy_active(k,tree_id_rhs_u)
            hvy_work(:,:,:,:,hvy_id,1) = hvy_block(:,:,:,:,hvy_id)
      enddo
    endif
    !!!! only 2D!!!!
    iz = 1

    do iblock = 1, hvy_n(tree_id_u)
        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
            do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                !write(*,*) "--------------------point---------------------------"
                !---------------------------------------------------------------
                ! reset entire grid to zeros (do not care about performance, just reset all)
                call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                    hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_ei, tree_id_u)
                !---------------------------------------------------------------
                ! set this one point we're looking at to 1
                call hvy2lgt(lgt_id, hvy_active(iblock,tree_id_ei), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
                if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk

                ! save its indices on the fine grid
                ixx = nint(x/dx_fine)+1
                iyy = nint(y/dx_fine)+1


                ! set the one
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                hvy_block(ix, iy, iz, 1, hvy_active(iblock,tree_id_ei)) = hvy_block(ix, iy, iz, 1, hvy_active(iblock,tree_id_ei)) + 1.0_rk
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if ( params%adapt_mesh ) then
                  call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_ei), hvy_n(tree_id_ei) )

                  call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_id_ei), lgt_n(tree_id_ei), &
                      lgt_sortednumlist(:,:,tree_id_ei), hvy_active(:,tree_id_ei), hvy_n(tree_id_ei), "everywhere", tree_ID=tree_id_ei )

                endif

                call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_ei), lgt_n(tree_id_ei), &
                          lgt_sortednumlist(:,:,tree_id_ei), hvy_active(:,tree_id_ei), hvy_n(tree_id_ei), tree_id_ei)
                !---------------------------------------------------------------
                ! synchronize ghosts (important! if e.g. coarseWins is actiyve and you happen to set the redundant value of a refined block, its overwritten to be zero again)
                ! Note: this also applies to coarse block bordering on a coarse block, if its ID is lower.
                ! In fact, each point is then computed only once. Note: if you set the point on a high lgt_id, then
                ! it will be "sync'ed down" to lower light IDs, so you can find the point more than once
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id_ei), hvy_n(tree_id_ei))
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Now compute rhs(u+he_i) of the derivative 1/h[rhs(u + h* e_i) - rhs(u)]
                ! call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
                !     lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor, tree_id_ei)

                if (OPERATOR == "RHS") then
                  call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
                     lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor, tree_id_ei)
                else if( OPERATOR == "EVOLVE") then
                  iteration = 0
                  time = 0
                  call time_stepper(time, dt, iteration, params, lgt_block, hvy_block, hvy_work, hvy_mask, hvy_tmp, &
                     hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist, tree_id_ei)
                     do k = 1, hvy_n(tree_id_ei)! call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
                       hvy_work(:,:,:,:,hvy_active(k,tree_id_ei),1) = hvy_block(:,:,:,:,hvy_active(k,tree_id_ei))
                     enddo
                else !! this operator is used for the identity
                  do k = 1, hvy_n(tree_id_ei)! call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, lgt_block, &
                    hvy_work(:,:,:,:,hvy_active(k,tree_id_ei),1) = hvy_block(:,:,:,:,hvy_active(k,tree_id_ei))
                  enddo
                endif
                ! This second sync step also synchronizes the derivative we computed previously
                ! Note that on a coarse/fine interface, wabbit computes two values for the derivative
                ! on the coarse and fine level. This synchronizing step lets us keep only either of those,
                ! depending on fineWins or coarseWins
                call sync_ghosts(params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_id_ei), hvy_n(tree_id_ei))
                call sync_ghosts(params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_id_rhs_u), hvy_n(tree_id_rhs_u))
                call delete_tree(params, lgt_block, lgt_active, lgt_n, hvy_active,hvy_n, tree_id_rhs_u_ei)
                call substract_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                    hvy_work(:,:,:,:,:,1), hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1=tree_id_ei, tree_id2=tree_id_rhs_u, dest_tree_id=tree_id_rhs_u_ei)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !call sync_ghosts(params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_id_rhs), hvy_n(tree_id_rhs))
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !write(*,*) "--------------------point-3s--------------------------"

                if ( params%adapt_mesh ) then
                  call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_rhs_u_ei), lgt_n(tree_id_rhs_u_ei), &
                            lgt_sortednumlist(:,:,tree_id_rhs_u_ei), hvy_active(:,tree_id_rhs_u_ei), hvy_n(tree_id_rhs_u_ei), tree_id_rhs_u_ei)
                  call sync_ghosts(params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_id_rhs_u_ei), hvy_n(tree_id_rhs_u_ei))
                  ! call adapt_mesh( time, params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, lgt_active(:,tree_id_rhs_u_ei), &
                  !  lgt_n(tree_id_rhs_u_ei), lgt_sortednumlist(:,:,tree_id_rhs_u_ei), hvy_active(:,tree_id_rhs_u_ei), hvy_n(tree_id_rhs_u_ei), tree_id_rhs_u_ei, "everywhere", hvy_tmp, external_loop=.true. )
                  call coarse_tree_2_reference_mesh(params, tree_n, &
                        lgt_block, lgt_active(:,tree_id_rhs_u_ei), lgt_n(tree_id_rhs_u_ei), lgt_sortednumlist(:,:,tree_id_rhs_u_ei), &
                        lgt_block_ref, lgt_active_ref(:,1),lgt_n_ref(1), &
                        hvy_work(:,:,:,:,:,1), hvy_active(:,tree_id_rhs_u_ei), hvy_n(tree_id_rhs_u_ei), hvy_tmp, hvy_neighbor, tree_id_rhs_u_ei, verbosity)

                  call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active(:,tree_id_rhs_u_ei), lgt_n(tree_id_rhs_u_ei), &
                                lgt_sortednumlist(:,:,tree_id_rhs_u_ei), hvy_active(:,tree_id_rhs_u_ei), hvy_n(tree_id_rhs_u_ei), tree_id_rhs_u_ei)

                  call sync_ghosts(params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_id_rhs_u_ei), hvy_n(tree_id_rhs_u_ei))
                endif
                ! save operator line to text file.
                ! note: unfortunately, we use the index on the finest level, i.e., we temporarily
                ! create a matrix N_max**2 by N_max**2, where N_max is Npoints on the finest level.
                ! Many points do not exist; they are on coarse levels. however, this is a problem
                ! for the python script, because it first reads the entie matrix, then removes zero cols/rows.
                do k = 1, hvy_n(tree_id_rhs_u_ei)
                    call hvy2lgt(lgt_id, hvy_active(k,tree_id_rhs_u_ei), params%rank, params%number_blocks)
                    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                    do iy2 = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                        do ix2 = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                            x = dble(ix2-(g+1)) * dx(1) + x0(1)
                            y = dble(iy2-(g+1)) * dx(2) + x0(2)

                            if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
                            if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk

                            ixx2 = nint(x/dx_fine)+1
                            iyy2 = nint(y/dx_fine)+1

                            val = hvy_work(ix2, iy2, iz, 1, hvy_active(k,tree_id_rhs_u_ei),1) ! u_dx

                            if (abs(val) > 1.0e-13) then
                                ! this point is a nonzero value
                                write(17,'(i6,1x,i6,1x,es15.8)') ixx+(iyy-1)*nx_fine, ixx2+(iyy2-1)*nx_fine, val
                                write(*,*) "python col=", ixx2+(iyy2-1)*nx_fine -1 , "row=", ixx+(iyy-1)*nx_fine -1, val, "xy=",x,y, "block", k
                            endif
                            hvy_work(ix2, iy2, iz, 1, hvy_active(k,tree_id_rhs_u_ei),1) = 0
                        end do
                    end do
                end do



            enddo
        enddo
    enddo
    close(17)

end subroutine
