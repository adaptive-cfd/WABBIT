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
    character(len=2)                   :: order

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
        params%g = 4_ik

    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%order_predictor = "multiresolution_2nd"
        params%g = 2_ik

    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")

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

    params%wavelet="CDF40"
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

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID))

    ! calculate vorticity from velocities
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)

        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgtID, x0, dx )

        if (operator == "--vorticity") then
            if (params%dim == 3) then
                call compute_vorticity( hvy_block(:,:,:,1,hvyID), &
                hvy_block(:,:,:,2,hvyID), hvy_block(:,:,:,3,hvyID),&
                dx, Bs, g,&
                params%order_discretization, hvy_tmp(:,:,:,1:3,hvyID))
            else
                call compute_vorticity(hvy_block(:,:,:,1,hvyID), &
                hvy_block(:,:,:,2,hvyID), hvy_block(:,:,:,1,hvyID),&
                dx, Bs, g, &
                params%order_discretization, hvy_tmp(:,:,:,:,hvyID))
            end if

        elseif (operator=="--vor-abs") then
            call compute_vorticity_abs(hvy_block(:,:,:,1,hvyID), &
            hvy_block(:,:,:,2,hvyID), hvy_block(:,:,:,3,hvyID),&
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        elseif (operator == "--divergence") then
            call divergence( hvy_block(:,:,:,1,hvyID), &
            hvy_block(:,:,:,2,hvyID), &
            hvy_block(:,:,:,3,hvyID),&
            dx, Bs, g, &
            params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        elseif (operator == "--Q") then
            call compute_Qcriterion( hvy_block(:,:,:,1,hvyID), &
            hvy_block(:,:,:,2,hvyID), &
            hvy_block(:,:,:,3,hvyID),&
            dx, Bs, g, &
            params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

        elseif (operator == "--copy") then
            hvy_tmp(:,:,:,1,hvyID) = hvy_block(:,:,:,1,hvyID)

        else
            call abort(1812011, "operator is neither --vorticity --vor-abs --divergence --Q")

        endif
    end do

    ! unique grid has to sync before saving to HDF5 as we store the 1st ghost node as well for visualization
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID))

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

    elseif (operator=="--copy") then

        write( fname,'(a, "_", i12.12, ".h5")') 'copyUx', nint(time * 1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )
    endif
end subroutine compute_vorticity_post


subroutine wavelet_test(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: fname1, fname_out, fname2
    real(kind=rk)                      :: time, x, y
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g, ix,iy,iz
    integer(kind=ik), dimension(3)     :: Bs
    character(len=2)                   :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID, dim

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
    call get_command_argument(2, fname1)
    call get_command_argument(3, fname2)

    ! does the user need help?
    if (fname1=='--help' .or. fname1=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: wavelet-test"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wavelet test on equidistant grid"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(fname1))
    call check_file_exists(trim(fname2))



    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname1, lgt_n(tree_ID), time, iteration, params%domain_size, &
    params%Bs, params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    params%wavelet = "CDF44"
    params%g = 8_ik

    params%n_eqn = 2
    allocate(params%butcher_tableau(1,1))

    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "0"

    Bs = params%Bs
    g  = params%g

    call setup_wavelet(params)

    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=params%n_eqn, hvy_work=hvy_work, nrhs_slots1=1)

    ! read input data
    call readHDF5vct_tree( (/fname1, fname2/), params, hvy_block, tree_ID)

    ! call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    ! hvy_active(:,tree_ID), hvy_n(tree_ID) )


    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgtID, x0, dx )

        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1) - 0.40625
                y = dble(iy-(g+1)) * dx(2) + x0(2) - 0.40625

              if (x<-params%domain_size(1)/2.0) x = x + params%domain_size(1)
              if (x>params%domain_size(1)/2.0) x = x - params%domain_size(1)
              if (y<-params%domain_size(2)/2.0) y = y + params%domain_size(2)
              if (y>params%domain_size(2)/2.0) y = y - params%domain_size(2)

                hvy_block(ix,iy,:,:,hvyID) = exp( -( (x)**2 + (y)**2 ) / 0.001 ) !+ 1.0e-2*rand_nbr()
                hvy_block(ix,iy,:,:,hvyID) = rand_nbr()
            enddo
        enddo

        ! call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    end do

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active(:,tree_ID), hvy_n(tree_ID) )

    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        hvy_tmp(:,:,:,:,hvyID) = hvy_block(:,:,:,:,hvyID)
    end do

    do k=1,50
    call substitution_step( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID), inputDataSynced=.false. )
enddo

    !------------------------------------------------------------------------
    ! do k = 1, hvy_n(tree_ID)
    !     hvyID = hvy_active(k,tree_ID)
    !     call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    ! end do
    !
    ! call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    ! hvy_active(:,tree_ID), hvy_n(tree_ID) )
    !
    ! do k = 1, hvy_n(tree_ID)
    !     hvyID = hvy_active(k,tree_ID)
    !     call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvyID))
    ! end do
    !------------------------------------------------------------------------

    call saveHDF5_tree("out1_111.h5", time, iteration, 1, params, hvy_block, tree_ID )
    call saveHDF5_tree("inp1_111.h5", time, iteration, 1, params, hvy_tmp, tree_ID )

    ! call saveHDF5_tree("out2_111.h5", time, iteration, 2, params, hvy_block, tree_ID )
    ! call saveHDF5_tree("inp2_111.h5", time, iteration, 2, params, hvy_tmp, tree_ID )

    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        hvy_tmp(:,:,:,:,hvyID) = abs( hvy_block(:,:,:,:,hvyID) - hvy_tmp(:,:,:,:,hvyID) )
    end do

    call saveHDF5_tree("delta1_111.h5", time, iteration, 1, params, hvy_tmp, tree_ID )
    ! call saveHDF5_tree("delta2_111.h5", time, iteration, 2, params, hvy_tmp, tree_ID )

    ! call substitution_step( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )
    ! call saveHDF5_tree("twice1_111.h5", time, iteration, 1, params, hvy_block, tree_ID )
end subroutine





subroutine wavelet_test_coarsening(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: fname1, fname_out, fname2
    real(kind=rk)                      :: time, x, y
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g, ix,iy,iz
    integer(kind=ik), dimension(3)     :: Bs
    character(len=2)                   :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID, dim

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: nwork, ierr, iter , kk, Jmax, nx,ny,nz,&
    nc,neighborhood, Nwcl, Nwcr, lgtID_neighbor, Nscl, Nscr
    logical :: coarsen, block1, block2, block3
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: sc, wcx, wcy, wcxy, tmp_reconst

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! does the user need help?
    if (fname1=='--help' .or. fname1=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: wavelet-test"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wavelet test on equidistant grid"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif


    ! decide which order for discretization and predictor is used. Note predictor
    ! is used in ghost nodes sync'ing
    params%wavelet = "CDF22"
    params%g = 3_ik

    params%wavelet = "CDF44"
    params%g = 7_ik

    ! params%wavelet = "CDF42"
    ! params%g = 5_ik


    params%Bs = 32
    params%Jmax = 4
    Jmax = params%Jmax
    params%dim = 2
    params%domain_size = (/1.0_rk, 1.0_rk, 1.0_rk/)
    params%n_eqn = 2
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "0"
    params%number_blocks = 3* 2**(params%dim*params%Jmax)
    time = 0.0_rk
    Bs = params%Bs
    g  = params%g

    call setup_wavelet(params)

    !----------------------------------------------------------------------------
    ! create some data on fine equidistant grid
    !----------------------------------------------------------------------------
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=params%n_eqn, hvy_work=hvy_work, nrhs_slots1=1)

    ! call createEquidistantGrid_tree( params, params%Jmax, .true., tree_ID )
    !
    ! ! create just some data...
    ! do k = 1, hvy_n(tree_ID)
    !     hvyID = hvy_active(k,tree_ID)
    !     call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
    !     call get_block_spacing_origin( params, lgtID, x0, dx )
    !
    !     do iy = g+1, Bs(2)+g
    !         do ix = g+1, Bs(1)+g
    !             x = dble(ix-(g+1)) * dx(1) + x0(1) - 0.5_rk
    !             y = dble(iy-(g+1)) * dx(2) + x0(2) - 0.5_rk
    !
    !             if (x<-params%domain_size(1)/2.0) x = x + params%domain_size(1)
    !             if (x>params%domain_size(1)/2.0) x = x - params%domain_size(1)
    !             if (y<-params%domain_size(2)/2.0) y = y + params%domain_size(2)
    !             if (y>params%domain_size(2)/2.0) y = y - params%domain_size(2)
    !
    !             ! hvy_block(ix,iy,:,:,hvyID) = exp( -( (x)**2 + (y)**2 ) / 0.01 ) + 1.0e-1*rand_nbr()
    !             hvy_block(ix,iy,:,:,hvyID) = rand_nbr()
    !         enddo
    !     enddo
    ! end do
    !
    ! call saveHDF5_tree("coarsening_331.h5", time, iteration, 1, params, hvy_block, tree_ID )




    select case(params%wavelet)
    case ("CDF44")
        Nscl = g+5 ! dictated by support of h_tilde (HD) filter for SC
        Nscr = g+5
        Nwcl = Nscl+3 ! chosen such that g_tilde (GD) not not see the copied SC
        Nwcr = Nscr+5
    case ("CDF42")
        Nscl = g+3 ! dictated by support of h_tilde (HD) filter for SC
        Nscr = g+3
        Nwcl = Nscl+3 ! chosen such that g_tilde (GD) not not see the copied SC
        Nwcr = Nscr+5
    case ("CDF22")
        Nscl = g+1 ! dictated by support of h_tilde (HD) filter for SC
        Nscr = g+1
        Nwcl = Nscl + 0 ! chosen such that g_tilde (GD) not not see the copied SC
        Nwcr = Nscr + 2
    case default
        call abort(6252, "not yet implemented waveler rsub")
    end select


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Case A: the EXACT coarse extension
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    call readHDF5vct_tree( (/"coarsening_331.h5"/), params, hvy_block, tree_ID)

! ****''''*** here the loop begins
do iter= 1, 1

    ! compute FWT
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    end do

    ! call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    ! flag for coarsening
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        ! keep blocks 3001,
        !             0331,
        !             2121
        coarsen = .true.
        block1 = ( (lgt_block(lgtID,1)==3).and.(lgt_block(lgtID,2)==0).and.(lgt_block(lgtID,3)==0).and.(lgt_block(lgtID,4)==1) )
        block2 = ( (lgt_block(lgtID,1)==0).and.(lgt_block(lgtID,2)==3).and.(lgt_block(lgtID,3)==3).and.(lgt_block(lgtID,4)==1) )
        block3 = ( (lgt_block(lgtID,1)==2).and.(lgt_block(lgtID,2)==1).and.(lgt_block(lgtID,3)==2).and.(lgt_block(lgtID,4)==1) )
        coarsen = (.not. (block1 .or. block2 .or. block3) )

        if (coarsen) then
            lgt_block( lgtID, Jmax+ IDX_REFINE_STS ) = -1
        else
            lgt_block( lgtID, Jmax+ IDX_REFINE_STS ) = 0
        endif
    end do

    call ensureGradedness_tree( params, tree_ID )


    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)

    if (.not. allocated(sc  )) allocate(  sc(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcx )) allocate( wcx(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcy )) allocate( wcy(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcxy)) allocate(wcxy(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )




    ! removal of WC
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)

        ! on blocks that are to be coarsened ...
        if (lgt_block( lgtID, Jmax+ IDX_REFINE_STS ) == -1) then
            ! block will coarsen - delete all WC on this block
            hvy_block( (g+2):(Bs(1)+g):2, (g+1):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+1):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+2):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
        endif

        ! ... as well as blocks the have neighbors which will be coarsened
        ! (proper "coarse extension")
        if (lgt_block( lgtID, Jmax+IDX_REFINE_STS ) == 0) then

            !....................
            ! unpack/inflated Mallat ordering
            sc   = 0.0_rk
            wcx  = 0.0_rk
            wcy  = 0.0_rk
            wcxy = 0.0_rk
            ! copy from Spaghetti to inflated Mallat ordering
            if (modulo(g, 2) == 0) then
                ! even g
                sc(   1:nx:2, 1:ny:2, :, :) = hvy_block(1:nx:2, 1:ny:2, :, 1:nc, hvyID)
                wcx(  1:nx:2, 1:ny:2, :, :) = hvy_block(2:nx:2, 1:ny:2, :, 1:nc, hvyID)
                wcy(  1:nx:2, 1:ny:2, :, :) = hvy_block(1:nx:2, 2:ny:2, :, 1:nc, hvyID)
                wcxy( 1:nx:2, 1:ny:2, :, :) = hvy_block(2:nx:2, 2:ny:2, :, 1:nc, hvyID)
            else
                ! odd g
                sc(   2:nx-1:2, 2:ny-1:2, :, :) = hvy_block(2:nx-1:2, 2:ny-1:2, :, 1:nc, hvyID)
                wcx(  2:nx-1:2, 2:ny-1:2, :, :) = hvy_block(3:nx:2, 2:ny-1:2  , :, 1:nc, hvyID)
                wcy(  2:nx-1:2, 2:ny-1:2, :, :) = hvy_block(2:nx-1:2, 3:ny:2  , :, 1:nc, hvyID)
                wcxy( 2:nx-1:2, 2:ny-1:2, :, :) = hvy_block(3:nx:2, 3:ny:2    , :, 1:nc, hvyID)
            endif
            !....................

            do neighborhood = 1, 8
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                    ! this neighbor will coarsen
                    if (lgt_block(lgtID_neighbor, Jmax+IDX_REFINE_STS)==-1) then
                        select case(neighborhood)
                        case(1)
                            wcx(1:Nwcl, :, :, 1:nc) = 0.0_rk
                            wcy(1:Nwcl, :, :, 1:nc) = 0.0_rk
                            wcxy(1:Nwcl, :, :, 1:nc) = 0.0_rk
                        case(2)
                            wcx(:, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                            wcy(:, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                            wcxy(:, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                        case(3)
                            wcx(nx-Nwcr:nx, :, :, 1:nc) = 0.0_rk
                            wcy(nx-Nwcr:nx, :, :, 1:nc) = 0.0_rk
                            wcxy(nx-Nwcr:nx, :, :, 1:nc) = 0.0_rk
                        case(4)
                            wcx(:, 1:Nwcl, :, 1:nc) = 0.0_rk
                            wcy(:, 1:Nwcl, :, 1:nc) = 0.0_rk
                            wcxy(:, 1:Nwcl, :, 1:nc) = 0.0_rk
                        case(5)
                            wcx(1:Nwcl, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                            wcy(1:Nwcl, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                            wcxy(1:Nwcl, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                        case(6)
                            wcx(1:Nwcl, 1:Nwcl, :, 1:nc) = 0.0_rk
                            wcy(1:Nwcl, 1:Nwcl, :, 1:nc) = 0.0_rk
                            wcxy(1:Nwcl, 1:Nwcl, :, 1:nc) = 0.0_rk
                        case(7)
                            wcx(nx-Nwcr:ny, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                            wcy(nx-Nwcr:ny, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                            wcxy(nx-Nwcr:ny, ny-Nwcr:ny, :, 1:nc) = 0.0_rk
                        case(8)
                            wcx(nx-Nwcr:ny, 1:Nwcl, :, 1:nc) = 0.0_rk
                            wcy(nx-Nwcr:ny, 1:Nwcl, :, 1:nc) = 0.0_rk
                            wcxy(nx-Nwcr:ny, 1:Nwcl, :, 1:nc) = 0.0_rk
                        end select
                    endif
                endif
            enddo

            ! repack to Spaghetti-ordering
            hvy_block(:,:, :, 1:nc, hvyID) = 0.0_rk
            if (modulo(g, 2) == 0) then
                hvy_block(1:nx:2, 1:ny:2, :, 1:nc, hvyID) =   sc(1:nx:2, 1:ny:2, :, 1:nc)
                hvy_block(2:nx:2, 1:ny:2, :, 1:nc, hvyID) =  wcx(1:nx:2, 1:ny:2, :, 1:nc)
                hvy_block(1:nx:2, 2:ny:2, :, 1:nc, hvyID) =  wcy(1:nx:2, 1:ny:2, :, 1:nc)
                hvy_block(2:nx:2, 2:ny:2, :, 1:nc, hvyID) = wcxy(1:nx:2, 1:ny:2, :, 1:nc)
            else
                hvy_block(2:nx-1:2, 2:ny-1:2, :, 1:nc, hvyID) =   sc(2:nx-1:2, 2:ny-1:2, :, 1:nc)
                hvy_block(3:nx:2, 2:ny-1:2  , :, 1:nc, hvyID) =  wcx(2:nx-1:2, 2:ny-1:2, :, 1:nc)
                hvy_block(2:nx-1:2, 3:ny:2  , :, 1:nc, hvyID) =  wcy(2:nx-1:2, 2:ny-1:2, :, 1:nc)
                hvy_block(3:nx:2, 3:ny:2    , :, 1:nc, hvyID) = wcxy(2:nx-1:2, 2:ny-1:2, :, 1:nc)
            endif

        endif
    end do


    ! wc are modified - we have to sync
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )


    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvyID))
    end do


    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    call executeCoarsening_tree( params, hvy_block, tree_ID, .true. )
    call updateMetadata_tree(params, tree_ID)

    call saveHDF5_tree("coarsening_332.h5", time, iteration, 1, params, hvy_block, tree_ID )

    call substitution_step( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID), inputDataSynced=.true. )
    call saveHDF5_tree("coarsening_333.h5", time, iteration, 1, params, hvy_block, tree_ID )

enddo
! end of loop here


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Case B: the actual coarse extension substituion
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


call readHDF5vct_tree( (/"coarsening_331.h5"/), params, hvy_block, tree_ID)

! ****''''*** here the loop begins
do iter= 1, 1

    ! compute FWT
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    end do

    ! call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    ! flag for coarsening
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        ! keep blocks 3001,
        !             0331,
        !             2121
        coarsen = .true.
        block1 = ( (lgt_block(lgtID,1)==3).and.(lgt_block(lgtID,2)==0).and.(lgt_block(lgtID,3)==0).and.(lgt_block(lgtID,4)==1) )
        block2 = ( (lgt_block(lgtID,1)==0).and.(lgt_block(lgtID,2)==3).and.(lgt_block(lgtID,3)==3).and.(lgt_block(lgtID,4)==1) )
        block3 = ( (lgt_block(lgtID,1)==2).and.(lgt_block(lgtID,2)==1).and.(lgt_block(lgtID,3)==2).and.(lgt_block(lgtID,4)==1) )
        coarsen = (.not. (block1 .or. block2 .or. block3) )

        if (coarsen) then
            lgt_block( lgtID, Jmax+ IDX_REFINE_STS ) = -1
        else
            lgt_block( lgtID, Jmax+ IDX_REFINE_STS ) = 0
        endif
    end do

    call ensureGradedness_tree( params, tree_ID )


    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)

    if (.not. allocated(sc  )) allocate(  sc(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcx )) allocate( wcx(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcy )) allocate( wcy(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcxy)) allocate(wcxy(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )




    ! removal of WC
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)

        ! on blocks that are to be coarsened ...
        if (lgt_block( lgtID, Jmax+ IDX_REFINE_STS ) == -1) then
            ! block will coarsen - delete all WC on this block
            hvy_block( (g+2):(Bs(1)+g):2, (g+1):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+1):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+2):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
        endif
    end do


    ! wc are modified - we have to sync
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )


    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvyID))
    end do


    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    call executeCoarsening_tree( params, hvy_block, tree_ID, .true. )
    call updateMetadata_tree(params, tree_ID)

    ! do kk = 1, 50
        call substitution_step( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID), inputDataSynced=.true. )
    ! enddo
    call saveHDF5_tree("coarsening_334.h5", time, iteration, 1, params, hvy_block, tree_ID )

enddo
! end of loop here



end subroutine
