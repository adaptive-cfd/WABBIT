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
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central) and 2 (FD_2nd_central)")

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

        ! elseif (operator == "--laplace") then
        !     call divergence( hvy_block(:,:,:,1,hvyID), &
        !     hvy_block(:,:,:,2,hvyID), &
        !     hvy_block(:,:,:,3,hvyID),&
        !     dx, Bs, g, &
        !     params%order_discretization, hvy_tmp(:,:,:,1,hvyID))

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
    use module_unit_test

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: fname1, fname_out, fname2
    real(kind=rk)                      :: time, x, y
    integer(kind=ik)                   :: iteration=-99, lgtID, tc_length, g, ix,iy,iz, ic, k
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
    integer(kind=tsize)                :: treecode, treecode_check
    real(kind=rk), allocatable :: wc(:,:,:,:,:)
    logical :: coarsen

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

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

    !!!!!!!!!!!!!!!!!!!!!!
    ! yes works in 2D or 3D
    params%dim = 2
    !!!!!!!!!!!!!!!!!!!!!
    params%wavelet = "CDF46"
    params%Bs = 32
    params%n_eqn = 1
    params%domain_size=1.0_rk
    params%Jmax=3
    params%Jmin=1
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "0"

    call setup_wavelet(params, params%g)

    Bs = params%Bs
    g  = params%g

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=params%n_eqn, hvy_work=hvy_work, nrhs_slots1=1)

    call unitTest_waveletDecomposition( params, hvy_block, hvy_work, hvy_tmp, tree_ID=1)
    call reset_tree(params, .true., tree_ID)

    call createEquidistantGrid_tree( params, hvy_block, 3, .true., tree_ID )

    ! flag for coarsening
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        ! keep a single block ( which will be 8 blocks due to completeness)
        ! block has ID 033 = 0*1 + 3*2^n_dim + 3*(2^n_dim)^2
        treecode_check = 3*2**params%dim + 3* 2**(2*params%dim)
        treecode = get_tc(lgt_block(lgtID, IDX_TC_1 : IDX_TC_2))
        coarsen = .not. (treecode == treecode_check)
        ! coarsen = .not. ( (lgt_block(lgtID,1)==0).and.(lgt_block(lgtID,2)==3).and.(lgt_block(lgtID,3)==3) )

        if (coarsen) then
            lgt_block( lgtID, IDX_REFINE_STS ) = -1
        else
            lgt_block( lgtID, IDX_REFINE_STS ) = 0
        endif
    end do

    call ensureGradedness_tree( params, tree_ID )
    call executeCoarsening_tree( params, hvy_block, tree_ID )
    call updateMetadata_tree(params, tree_ID) ! because we do not call adapt_mesh here

    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        if (params%dim == 3) then
            do ic = 1, params%n_eqn
                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            hvy_block(ix,iy,iz,ic,hvyID) = rand_nbr()
                        end do
                    end do
                end do
            end do
        else
            do ic = 1, params%n_eqn
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        hvy_block(ix,iy,:,ic,hvyID) = rand_nbr()
                    end do
                end do
            end do
        endif
    end do


    call saveHDF5_tree("data_0001.h5", 1.0_rk, iteration, 1, params, hvy_block, tree_ID )

    call coarseExtensionUpdate_tree( params, hvy_block, hvy_tmp, tree_ID )

    call saveHDF5_tree("data_0002.h5", 2.0_rk, iteration, 1, params, hvy_block, tree_ID )

    do k = 1, 20
        call coarseExtensionUpdate_tree( params, hvy_block, hvy_tmp, tree_ID )
    enddo

    call saveHDF5_tree("data_0003.h5", 3.0_rk, iteration, 1, params, hvy_block, tree_ID )


end subroutine





!     do k = 1, hvy_n(tree_ID)
!         hvyID = hvy_active(k,tree_ID)
!         call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
!
!         ! unpack/inflated Mallat ordering
!         call spaghetti2inflatedMallat_block(params, hvy_block(:,:,:,:,hvyID), wc)
!
!         ! on blocks that are to be coarsened ...
!         if (lgt_block( lgtID, IDX_REFINE_STS ) == -1) then
!             ! block will coarsen - delete all WC on this block
!             wc(:,:,:,:,2:8) = 0.0_rk
!         endif
!
!         ! ... as well as blocks the have neighbors which will be coarsened
!         ! (proper "coarse extension")
!         if (lgt_block( lgtID, IDX_REFINE_STS ) == 0) then
!
!             do neighborhood = 1, 8
!                 ! neighbor exists ?
!                 if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
!                     ! neighbor light data id
!                     lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
!                     ! this neighbor will coarsen
!                     if (lgt_block(lgtID_neighbor, IDX_REFINE_STS)==-1) then
!         ! call coarseExtensionManipulateWC_block(params, wc, neighborhood)
! !! here, unusual neighbors 1-8 are used, and we consequently cannot use the above routine
!                         ! select case(neighborhood)
!                         ! case(1)
!                         !     wc(1:params%Nwcl, :, :, 1:nc, 2) = 0.0_rk
!                         !     wc(1:params%Nwcl, :, :, 1:nc, 3) = 0.0_rk
!                         !     wc(1:params%Nwcl, :, :, 1:nc, 4) = 0.0_rk
!                         ! case(2)
!                         !     wc(:, ny-params%Nwcr:ny, :, 1:nc, 2) = 0.0_rk
!                         !     wc(:, ny-params%Nwcr:ny, :, 1:nc, 3) = 0.0_rk
!                         !     wc(:, ny-params%Nwcr:ny, :, 1:nc, 4) = 0.0_rk
!                         ! case(3)
!                         !     wc(nx-params%Nwcr:nx, :, :, 1:nc, 2) = 0.0_rk
!                         !     wc(nx-params%Nwcr:nx, :, :, 1:nc, 3) = 0.0_rk
!                         !     wc(nx-params%Nwcr:nx, :, :, 1:nc, 4) = 0.0_rk
!                         ! case(4)
!                         !     wc(:, 1:params%Nwcl, :, 1:nc, 2) = 0.0_rk
!                         !     wc(:, 1:params%Nwcl, :, 1:nc, 3) = 0.0_rk
!                         !     wc(:, 1:params%Nwcl, :, 1:nc, 4) = 0.0_rk
!                         ! case(5)
!                         !     wc(1:params%Nwcl, ny-params%Nwcr:ny, :, 1:nc, 2) = 0.0_rk
!                         !     wc(1:params%Nwcl, ny-params%Nwcr:ny, :, 1:nc, 3) = 0.0_rk
!                         !     wc(1:params%Nwcl, ny-params%Nwcr:ny, :, 1:nc, 4) = 0.0_rk
!                         ! case(6)
!                         !     wc(1:params%Nwcl, 1:params%Nwcl, :, 1:nc, 2) = 0.0_rk
!                         !     wc(1:params%Nwcl, 1:params%Nwcl, :, 1:nc, 3) = 0.0_rk
!                         !     wc(1:params%Nwcl, 1:params%Nwcl, :, 1:nc, 4) = 0.0_rk
!                         ! case(7)
!                         !     wc(nx-params%Nwcr:ny, ny-params%Nwcr:ny, :, 1:nc, 2) = 0.0_rk
!                         !     wc(nx-params%Nwcr:ny, ny-params%Nwcr:ny, :, 1:nc, 3) = 0.0_rk
!                         !     wc(nx-params%Nwcr:ny, ny-params%Nwcr:ny, :, 1:nc, 4) = 0.0_rk
!                         ! case(8)
!                         !     wc(nx-params%Nwcr:ny, 1:params%Nwcl, :, 1:nc, 2) = 0.0_rk
!                         !     wc(nx-params%Nwcr:ny, 1:params%Nwcl, :, 1:nc, 3) = 0.0_rk
!                         !     wc(nx-params%Nwcr:ny, 1:params%Nwcl, :, 1:nc, 4) = 0.0_rk
!                         ! end select
!                     endif
!                 endif
!             enddo
!
!             ! repack to Spaghetti-ordering
!             call inflatedMallat2spaghetti_block(params, wc, hvy_block(:,:, :,:,hvyID))
!         endif
!     end do






! subroutine wavelet_test(params)
!     use module_globals
!     use module_mesh
!     use module_params
!     use module_mpi
!     use module_operators
!     use module_forestMetaData
!
!     implicit none
!
!     !> parameter struct
!     type (type_params), intent(inout)  :: params
!     character(len=cshort)              :: fname1, fname_out, fname2
!     real(kind=rk)                      :: time, x, y
!     integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g, ix,iy,iz
!     integer(kind=ik), dimension(3)     :: Bs
!     character(len=2)                   :: order
!
!     real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
!     real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
!     real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
!     integer(kind=ik)                   :: tree_ID=1, hvyID, dim
!
!     character(len=cshort)              :: fname
!     real(kind=rk), dimension(3)        :: dx, x0
!     integer(hid_t)                     :: file_id
!     real(kind=rk), dimension(3)        :: domain
!     integer(kind=ik)                   :: nwork
!
!     ! this routine works only on one tree
!     allocate( hvy_n(1), lgt_n(1) )
!
!     ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
!     ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
!     ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
!     ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
!     ! to the last subroutine.)  -Thomas
!
!
!     !-----------------------------------------------------------------------------------------------------
!     ! get values from command line (filename and level for interpolation)
!     call get_command_argument(2, fname1)
!     call get_command_argument(3, fname2)
!
!     ! does the user need help?
!     if (fname1=='--help' .or. fname1=='--h') then
!         if (params%rank==0) then
!             write(*,*) "-----------------------------------------------------------"
!             write(*,*) " Wabbit postprocessing: wavelet-test"
!             write(*,*) "-----------------------------------------------------------"
!             write(*,*) " Wavelet test on equidistant grid"
!             write(*,*) "-----------------------------------------------------------"
!         end if
!         return
!     endif
!
!     call check_file_exists(trim(fname1))
!     call check_file_exists(trim(fname2))
!
!     ! get some parameters from one of the files (they should be the same in all of them)
!     call read_attributes(fname1, lgt_n(tree_ID), time, iteration, params%domain_size, &
!     params%Bs, params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
!
!
!     params%wavelet = "CDF44"
!     params%g = 8_ik
!
!     params%n_eqn = 2
!     allocate(params%butcher_tableau(1,1))
!
!     allocate(params%symmetry_vector_component(1:params%n_eqn))
!     params%symmetry_vector_component(1) = "0"
!
!     Bs = params%Bs
!     g  = params%g
!
!     call setup_wavelet(params)
!
!     ! no refinement is made in this postprocessing tool; we therefore allocate about
!     ! the number of blocks in the file (and not much more than that)
!     params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )
!
!     ! allocate data
!     call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=params%n_eqn, hvy_work=hvy_work, nrhs_slots1=1)
!
!     ! read input data
!     call readHDF5vct_tree( (/fname1, fname2/), params, hvy_block, tree_ID)
!
!     ! call sync_ghosts_tree( params, hvy_block, tree_ID )
!
!
!     do k = 1, hvy_n(tree_ID)
!         hvyID = hvy_active(k,tree_ID)
!         call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
!         call get_block_spacing_origin( params, lgtID, x0, dx )
!
!         do iy = g+1, Bs(2)+g
!             do ix = g+1, Bs(1)+g
!                 x = dble(ix-(g+1)) * dx(1) + x0(1) - 0.40625
!                 y = dble(iy-(g+1)) * dx(2) + x0(2) - 0.40625
!
!               if (x<-params%domain_size(1)/2.0) x = x + params%domain_size(1)
!               if (x>params%domain_size(1)/2.0) x = x - params%domain_size(1)
!               if (y<-params%domain_size(2)/2.0) y = y + params%domain_size(2)
!               if (y>params%domain_size(2)/2.0) y = y - params%domain_size(2)
!
!                 hvy_block(ix,iy,:,:,hvyID) = exp( -( (x)**2 + (y)**2 ) / 0.001 ) !+ 1.0e-2*rand_nbr()
!                 hvy_block(ix,iy,:,:,hvyID) = rand_nbr()
!             enddo
!         enddo
!
!         ! call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
!     end do
!
!     call sync_ghosts_tree( params, hvy_block, tree_ID )
!
!     do k = 1, hvy_n(tree_ID)
!         hvyID = hvy_active(k,tree_ID)
!         hvy_tmp(:,:,:,:,hvyID) = hvy_block(:,:,:,:,hvyID)
!     end do
!
!     do k=1,50
!     call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID), inputDataSynced=.false. )
! enddo
!
!     !------------------------------------------------------------------------
!     ! do k = 1, hvy_n(tree_ID)
!     !     hvyID = hvy_active(k,tree_ID)
!     !     call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
!     ! end do
!     !
!     ! call sync_ghosts_tree( params, hvy_block, tree_ID )
!     !
!     ! do k = 1, hvy_n(tree_ID)
!     !     hvyID = hvy_active(k,tree_ID)
!     !     call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvyID))
!     ! end do
!     !------------------------------------------------------------------------
!
!     call saveHDF5_tree("out1_111.h5", time, iteration, 1, params, hvy_block, tree_ID )
!     call saveHDF5_tree("inp1_111.h5", time, iteration, 1, params, hvy_tmp, tree_ID )
!
!     ! call saveHDF5_tree("out2_111.h5", time, iteration, 2, params, hvy_block, tree_ID )
!     ! call saveHDF5_tree("inp2_111.h5", time, iteration, 2, params, hvy_tmp, tree_ID )
!
!     do k = 1, hvy_n(tree_ID)
!         hvyID = hvy_active(k,tree_ID)
!         hvy_tmp(:,:,:,:,hvyID) = abs( hvy_block(:,:,:,:,hvyID) - hvy_tmp(:,:,:,:,hvyID) )
!     end do
!
!     call saveHDF5_tree("delta1_111.h5", time, iteration, 1, params, hvy_tmp, tree_ID )
!     ! call saveHDF5_tree("delta2_111.h5", time, iteration, 2, params, hvy_tmp, tree_ID )
!
!     ! call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )
!     ! call saveHDF5_tree("twice1_111.h5", time, iteration, 1, params, hvy_block, tree_ID )
! end subroutine





subroutine wavelet_test_coarsening(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_unit_test

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
    nc, neighborhood, lgtID_neighbor, ii
    integer(kind=tsize)                :: treecode, treecode1, treecode2, treecode3
    logical :: coarsen, block1, block2, block3
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: tmp_reconst
    real(kind=rk), allocatable, dimension(:,:,:,:,:), save :: wc

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
    params%Jmax = 3
    Jmax = params%Jmax
    params%dim = 3
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

    call unitTest_waveletDecomposition(params, hvy_block, hvy_work, hvy_tmp, tree_ID)

    ! call createEquidistantGrid_tree( params, 4, .true., tree_ID )
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



    call readHDF5vct_tree( (/"coarsening_331.h5"/), params, hvy_block, tree_ID)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Case A: the EXACT coarse extension
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ****''''*** here the loop begins
do iter= 1, 1

    ! compute FWT
    call sync_ghosts_tree( params, hvy_block, tree_ID )

    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    end do

    ! flag for coarsening
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        ! keep blocks 3001 = 3*1                                  + 1* (2**n_dim)**3
        !             0331 =       3* 2**n_dim + 3* (2**n_dim)**2 + 1* (2**n_dim)**3
        !             2121 = 2*1 + 1* 2**n_dim + 2* (2**n_dim)**2 + 1* (2**n_dim)**3
        coarsen = .true.
        treecode = get_tc(lgt_block(lgtID, IDX_TC_1 : IDX_TC_2))
        treecode1 = 3                                           + 1* 2**(3*params%dim)
        treecode2 =     3* 2**params%dim + 3* 2**(2*params%dim) + 1* 2**(3*params%dim)
        treecode3 = 2 + 1* 2**params%dim + 2* 2**(2*params%dim) + 1* 2**(3*params%dim)
        ! block1 = ( (lgt_block(lgtID,1)==3).and.(lgt_block(lgtID,2)==0).and.(lgt_block(lgtID,3)==0).and.(lgt_block(lgtID,4)==1) )
        ! block2 = ( (lgt_block(lgtID,1)==0).and.(lgt_block(lgtID,2)==3).and.(lgt_block(lgtID,3)==3).and.(lgt_block(lgtID,4)==1) )
        ! block3 = ( (lgt_block(lgtID,1)==2).and.(lgt_block(lgtID,2)==1).and.(lgt_block(lgtID,3)==2).and.(lgt_block(lgtID,4)==1) )
        coarsen = (.not. ((treecode == treecode1) .or. (treecode == treecode2) .or. (treecode == treecode3)) )

        if (coarsen) then
            lgt_block( lgtID, IDX_REFINE_STS ) = -1
        else
            lgt_block( lgtID, IDX_REFINE_STS ) = 0
        endif
    end do

    call ensureGradedness_tree( params, tree_ID )


    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)

    if (.not. allocated(wc  )) allocate( wc(1:nx, 1:ny, 1:nz, 1:nc, 1:8) )
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )




    ! removal of WC
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)

        ! on blocks that are to be coarsened ...
        if (lgt_block( lgtID, IDX_REFINE_STS ) == -1) then
            ! block will coarsen - delete all WC on this block
            hvy_block( (g+2):(Bs(1)+g):2, (g+1):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+1):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+2):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
        endif

        ! ... as well as blocks the have neighbors which will be coarsened
        ! (proper "coarse extension")
        if (lgt_block( lgtID, IDX_REFINE_STS ) == 0) then
            ! unpack/inflated Mallat ordering
            call spaghetti2inflatedMallat_block(params, hvy_block(:,:,:,:,hvyID), wc)

            do neighborhood = 1, 8
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                    ! this neighbor will coarsen
                    if (lgt_block(lgtID_neighbor, IDX_REFINE_STS)==-1) then
                        select case(neighborhood)
                        case(1)
                            wc(1:params%Nwcl, :, :, 1:nc, 2) = 0.0_rk
                            wc(1:params%Nwcl, :, :, 1:nc, 3) = 0.0_rk
                            wc(1:params%Nwcl, :, :, 1:nc, 4) = 0.0_rk
                        case(2)
                            wc(:, ny-params%Nwcr:ny, :, 1:nc, 2) = 0.0_rk
                            wc(:, ny-params%Nwcr:ny, :, 1:nc, 3) = 0.0_rk
                            wc(:, ny-params%Nwcr:ny, :, 1:nc, 4) = 0.0_rk
                        case(3)
                            wc(nx-params%Nwcr:nx, :, :, 1:nc, 2) = 0.0_rk
                            wc(nx-params%Nwcr:nx, :, :, 1:nc, 3) = 0.0_rk
                            wc(nx-params%Nwcr:nx, :, :, 1:nc, 4) = 0.0_rk
                        case(4)
                            wc(:, 1:params%Nwcl, :, 1:nc, 2) = 0.0_rk
                            wc(:, 1:params%Nwcl, :, 1:nc, 3) = 0.0_rk
                            wc(:, 1:params%Nwcl, :, 1:nc, 4) = 0.0_rk
                        case(5)
                            wc(1:params%Nwcl, ny-params%Nwcr:ny, :, 1:nc, 2) = 0.0_rk
                            wc(1:params%Nwcl, ny-params%Nwcr:ny, :, 1:nc, 3) = 0.0_rk
                            wc(1:params%Nwcl, ny-params%Nwcr:ny, :, 1:nc, 4) = 0.0_rk
                        case(6)
                            wc(1:params%Nwcl, 1:params%Nwcl, :, 1:nc, 2) = 0.0_rk
                            wc(1:params%Nwcl, 1:params%Nwcl, :, 1:nc, 3) = 0.0_rk
                            wc(1:params%Nwcl, 1:params%Nwcl, :, 1:nc, 4) = 0.0_rk
                        case(7)
                            wc(nx-params%Nwcr:ny, ny-params%Nwcr:ny, :, 1:nc, 2) = 0.0_rk
                            wc(nx-params%Nwcr:ny, ny-params%Nwcr:ny, :, 1:nc, 3) = 0.0_rk
                            wc(nx-params%Nwcr:ny, ny-params%Nwcr:ny, :, 1:nc, 4) = 0.0_rk
                        case(8)
                            wc(nx-params%Nwcr:ny, 1:params%Nwcl, :, 1:nc, 2) = 0.0_rk
                            wc(nx-params%Nwcr:ny, 1:params%Nwcl, :, 1:nc, 3) = 0.0_rk
                            wc(nx-params%Nwcr:ny, 1:params%Nwcl, :, 1:nc, 4) = 0.0_rk
                        end select
                    endif
                endif
            enddo

            ! repack to Spaghetti-ordering
            call inflatedMallat2spaghetti_block(params, wc(:, :, :, 1:nc, :), hvy_block(:,:, :,1:nc,hvyID))
        endif
    end do


    ! wc are modified - we have to sync
    call sync_ghosts_tree( params, hvy_block, tree_ID )


    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvyID))
    end do


    call sync_ghosts_tree( params, hvy_block, tree_ID )

    ! call executeCoarsening_tree( params, hvy_block, tree_ID, .true. )
    call abort(150623, "we removed the ignore_prefilter, check if that is okay still")

    call updateMetadata_tree(params, tree_ID)

    call saveHDF5_tree("coarsening_332.h5", time, iteration, 1, params, hvy_block, tree_ID )

    ! TODO: update this. CE MUST be done level wise. add wrapper to do this.
    ! call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID), lgt_n(tree_ID),inputDataSynced=.true. )
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
    call sync_ghosts_tree( params, hvy_block, tree_ID )

    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    end do

    ! call sync_ghosts_tree( params, hvy_block, tree_ID )

    ! flag for coarsening
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        ! keep blocks 3001 = 3*1                                  + 1* (2**n_dim)**3
        !             0331 =       3* 2**n_dim + 3* (2**n_dim)**2 + 1* (2**n_dim)**3
        !             2121 = 2*1 + 1* 2**n_dim + 2* (2**n_dim)**2 + 1* (2**n_dim)**3
        coarsen = .true.
        treecode = get_tc(lgt_block(lgtID, IDX_TC_1 : IDX_TC_2))
        treecode1 = 3                                           + 1* 2**(3*params%dim)
        treecode2 =     3* 2**params%dim + 3* 2**(2*params%dim) + 1* 2**(3*params%dim)
        treecode3 = 2 + 1* 2**params%dim + 2* 2**(2*params%dim) + 1* 2**(3*params%dim)
        ! block1 = ( (lgt_block(lgtID,1)==3).and.(lgt_block(lgtID,2)==0).and.(lgt_block(lgtID,3)==0).and.(lgt_block(lgtID,4)==1) )
        ! block2 = ( (lgt_block(lgtID,1)==0).and.(lgt_block(lgtID,2)==3).and.(lgt_block(lgtID,3)==3).and.(lgt_block(lgtID,4)==1) )
        ! block3 = ( (lgt_block(lgtID,1)==2).and.(lgt_block(lgtID,2)==1).and.(lgt_block(lgtID,3)==2).and.(lgt_block(lgtID,4)==1) )
        coarsen = (.not. ((treecode == treecode1) .or. (treecode == treecode2) .or. (treecode == treecode3)) )

        if (coarsen) then
            lgt_block( lgtID, IDX_REFINE_STS ) = -1
        else
            lgt_block( lgtID, IDX_REFINE_STS ) = 0
        endif
    end do

    call ensureGradedness_tree( params, tree_ID )


    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)

    if (.not. allocated(wc  )) allocate(  wc(1:nx, 1:ny, 1:nz, 1:nc, 1:8) )
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )


    ! removal of WC
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)

        ! on blocks that are to be coarsened ...
        if (lgt_block( lgtID, IDX_REFINE_STS ) == -1) then
            ! block will coarsen - delete all WC on this block
            hvy_block( (g+2):(Bs(1)+g):2, (g+1):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+1):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
            hvy_block( (g+2):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :, hvyID) = 0.0_rk
        endif
    end do


    ! wc are modified - we have to sync
    call sync_ghosts_tree( params, hvy_block, tree_ID )


    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvyID))
    end do


    call sync_ghosts_tree( params, hvy_block, tree_ID )

    call executeCoarsening_tree( params, hvy_block, tree_ID )
    call abort(150623, "we removed the ignore_prefilter, check if that is okay still")
    call updateMetadata_tree(params, tree_ID)

    ! do kk = 1, 50
        ! call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID),lgt_n(tree_ID), inputDataSynced=.true. )
    ! enddo
    call saveHDF5_tree("coarsening_334.h5", time, iteration, 1, params, hvy_block, tree_ID )

call refine_tree(params, hvy_block, hvy_tmp, "everywhere", tree_ID)

call saveHDF5_tree("coarsening_335.h5", time, iteration, 1, params, hvy_block, tree_ID )
! call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID),lgt_n(tree_ID), inputDataSynced=.true. )
call saveHDF5_tree("coarsening_336.h5", time, iteration, 1, params, hvy_block, tree_ID )
enddo
! end of loop here



! do k = 1, hvy_n(tree_ID)
!     hvyID = hvy_active(k,tree_ID)
!     call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
!     ! keep blocks 3001,
!     !             0331,
!     !             2121
!     block1 = ( (lgt_block(lgtID,1)==3).and.(lgt_block(lgtID,2)==0).and.(lgt_block(lgtID,3)==0).and.(lgt_block(lgtID,4)==1) )
!     block2 = ( (lgt_block(lgtID,1)==0).and.(lgt_block(lgtID,2)==3).and.(lgt_block(lgtID,3)==3).and.(lgt_block(lgtID,4)==1) )
!     block3 = ( (lgt_block(lgtID,1)==2).and.(lgt_block(lgtID,2)==1).and.(lgt_block(lgtID,3)==2).and.(lgt_block(lgtID,4)==1) )
!
!     if (block1) then
!
!         call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
!
!         call spaghetti2inflatedMallat_block(params, hvy_block(:,:,:,:,hvyID), sc, wcx, wcy, wcxy)
!
!         open(unit=32, file="wcx.csv", status="replace")
!         do ii = g+1, Bs(2)+g
!             write(32,'(48(es12.4,";"))') wcx( g+1:Bs(1)+g, ii, 1, 1)
!         enddo
!         close(32)
!         open(unit=32, file="wcy.csv", status="replace")
!         do ii = g+1, Bs(2)+g
!             write(32,'(48(es12.4,";"))') wcy( g+1:Bs(1)+g, ii, 1, 1)
!         enddo
!         close(32)
!         open(unit=32, file="wcxy.csv", status="replace")
!         do ii = g+1, Bs(2)+g
!             write(32,'(48(es12.4,";"))') wcxy( g+1:Bs(1)+g, ii, 1, 1)
!         enddo
!         close(32)
!     endif
! end do




end subroutine




subroutine waveletVsRHS_timingTest(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_time_step

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: filename
    real(kind=rk)                      :: time, x, y, t0
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g, ix,iy,iz, ic
    integer(kind=ik), dimension(3)     :: Bs
    character(len=2)                   :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_mask(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID, dim

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: nwork, ii, hvy_ID
    real(kind=rk), allocatable :: wc(:,:,:,:,:)
    logical :: coarsen

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    Bs = params%Bs
    g  = params%g

    call get_command_argument( 2, filename )
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, filename )
    call setup_wavelet(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, filename, params%N_mask_components )

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_mask=hvy_mask, hvy_work=hvy_work, nrhs_slots1=1)

    ! N_MAX_COMPONENTS = 4
    params%max_grid_density = 0.10_rk
    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    call createRandomGrid_tree( params, hvy_block, hvy_tmp, 2, .true., 4, tree_ID )


    t0 = MPI_wtime()
    do k = 1, 100
        call RHS_wrapper(0.0_rk, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, tree_ID)
    enddo
    write(*,*) "rhs", MPI_wtime()-t0

    t0 = MPI_wtime()
    do ii = 1, 100
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call WaveDecomposition_dim1( params, hvy_block(:,:,:,:,hvy_id) )
            call WaveReconstruction_dim1( params, hvy_block(:,:,:,:,hvy_id) )
        enddo
    enddo
    write(*,*) "wavelet", MPI_wtime()-t0
end subroutine
