!! Module for MODEL ORDER REDUCTION
!! Current Methods:
!!              * snapshot POD
!-----------------------------------------------------------------

module module_MOR

    use mpi
    use module_globals
    use module_globals
    use module_params
    use module_helpers
    use module_forestMetaData
    use module_mesh

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    !**********************************************************************************************
    ! These are the important routines that are visible to WABBIT:
    !**********************************************************************************************
    PUBLIC :: snapshot_POD, post_POD, post_reconstruct, post_PODerror, compute_tree_L2norm ,post_timecoef_POD
    !**********************************************************************************************
contains


    !##############################################################
    !> Implementation of the multiresolution snapshot POD
    !> Input:
    !>  1.) All Snapshots needed for the POD in the first N_snapshots trees.
    !>  2.) The maximal forest size needs to be at least 2*N_snapshots+1.
    !>  3.) Truncation Rank (optional) is the number of POD MODES computed in the algorithm (maxium is N_snapshots)
    !>      (default is N_snapshots)
    !>  4.) Truncation Error (default: 1e-9) is the smallest squared singularvalue used
    !>      for POD modes
    !> The output is as the following:
    !>  1.) All input snapshots are returned with tree: id 1, 2, ..., N_snapshots
    !>  2.) The POD modes are saved in trees with tree_IDs:
    !>            N_snapshots+1, N_snapshots+2, ..., N_snapshots+truncation_rank
    !>  3.) On output the truncation rank is updated to the acutal number of modes with
    !>      squared singular value bigger then the given or default truncation error
    subroutine snapshot_POD( params, hvy_block, hvy_tmp, truncation_error, truncation_rank, save_all)
        implicit none

        !-----------------------------------------------------------------
        type (type_params), intent(inout)        :: params
        real(kind=rk),  intent(inout)            :: hvy_block(:, :, :, :, :)  !> heavy data array - block data
        !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
        real(kind=rk),  intent(inout)            :: hvy_tmp(:, :, :, :, :)
        integer(kind=ik),optional, intent(inout) :: truncation_rank           !> number of POD modes
        !> Threshold value for truncating POD modes. If the singular value is smaller,tree_ID_dest
        !> then the given threshold we discard the corresponding POD MODE.
        real(kind=rk), optional, intent(in)      :: truncation_error
        !> if true we write out all temporal coefficients and eigenvalues
        !> filenames eigenvalues.txt, acoef.txt
        logical, optional, intent(in)            :: save_all
        !---------------------------------------------------------------
        real(kind=rk) :: C(tree_n,tree_n), V(tree_n,tree_n), work(5*tree_n), &
        eigenvalues(tree_n), alpha(tree_n), max_err, t_elapse, Volume
        integer(kind=ik):: N_snapshots, root, ierr, i, rank, pod_mode_tree_ID, &
        free_tree_ID, tree_ID, max_nr_pod_modes, it
        character(len=cshort):: filename
        character(len=30) :: rowfmt

        ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
        ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
        ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
        ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
        ! to the last subroutine.)  -Thomas


        !---------------------------------------------------------------------------
        ! check inputs and set default values
        !---------------------------------------------------------------------------
        rank= params%rank
        N_snapshots=tree_n
        if (present(truncation_rank)) then
            max_nr_pod_modes=truncation_rank
        else
            max_nr_pod_modes=N_snapshots
        endif

        if (present(truncation_error)) then
            max_err=truncation_error
        else
            max_err=0.0_rk
        endif

        if (rank == 0) then
            write(*, *)
            write(*,'(80("─"))')
            write(*,'(30("#"), " SNAPSHOT POD ", 30("#"))')
            write(*,'(80("─"))')
            write(*,'("Number of SNAPSHOTS used: ",i4)') N_snapshots
            write(*,'("Desired Truncation Rank: ", i4)') max_nr_pod_modes
            write(*,'("Maximal Error in L2 norm: ",es12.4)') max_err
            if (params%adapt_tree) write(*,'("Compression threshold eps: ",es12.4)') params%eps
            write(*,'(80("─"))')
            write(*, *)
        endif

        if ( params%forest_size <= N_snapshots + truncation_rank) call abort(1003191,"Error! Need more Trees. Tip: increase forest_size")

        !---------------------------------------------------------------------------
        ! Covariance Matrix of snapshot matrix X: C = X^T*X
        !---------------------------------------------------------------------------
        if (rank == 0) write(*,*)
        if (rank == 0) write(*,*) "Processing Covariance Matrix C = X^T*X"
        if (rank == 0) write(*,*)
        call compute_covariance_matrix( params, C, hvy_block, hvy_tmp)

        ! Save covariance matrix
        if (save_all .and. rank==0) then
            filename ="covariance_matrix.txt"
            write(*,'( "covariance_matrix saved to file: ", A30 )') filename
            write(*,*)
            write(rowfmt,'(A,I4,A)') '(',N_snapshots,'(es15.8,1x))'
            open(14,file=filename, status='replace')
            do i = 1,N_snapshots
                write(14,FMT=rowfmt) C(i,:)
            enddo
            close(14)
        end if

        !---------------------------------------------------------------------------
        ! eigenvalues of covariance matrix
        !---------------------------------------------------------------------------
        call DSYEV('V', 'U', N_snapshots, C, N_snapshots, eigenvalues, work, 5*N_snapshots, ierr)
        ! as in matlab the eigenvalues are sorted in ascending order...
        ! on output, C now contains the eigenvectors:
        V = C

        if (ierr /= 0) call abort(333,"The eigenvalue solver failed...")


        if (save_all .and. rank==0) then
            ! Save Eigenvalues if requested:
            filename ="eigenvalues.txt"
            write(*,'( "eigenvalues saved to file: ", A30 )') filename
            write(*,*)
            open(14,file=filename,status='replace')
            do i = 1, N_snapshots
                write (14,'(i4,1x,g18.8)') i, eigenvalues(i)
            end do
            close(14)

            ! save eigenvectors
            write(*,*)
            filename = "eigenvectors.txt"
            write(*,'( "Eigenvectors saved to file: ", A30 )') filename
            open(14, file=filename, status='replace')
            do i = 1, N_snapshots
                write(14,FMT=rowfmt) V(i, :)
            enddo
            close(14)
        end if

        !---------------------------------------------------------------------------
        ! construct POD basis functions (modes)
        !---------------------------------------------------------------------------
        call construct_modes( params, hvy_block, hvy_tmp, truncation_error, truncation_rank, V, eigenvalues)

    end subroutine snapshot_POD
    !##############################################################



    !##############################################################
    subroutine construct_modes( params, hvy_block, hvy_tmp, truncation_error, truncation_rank, eigenbasis, eigenvalues)
        implicit none

        !-----------------------------------------------------------------
        type (type_params), intent(inout)        :: params
        real(kind=rk),  intent(inout)            :: hvy_block(:, :, :, :, :)  !> heavy data array - block data
        !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
        real(kind=rk),  intent(inout)            :: hvy_tmp(:, :, :, :, :)
        !> POD eigenvalues and eigenbasis to construct modes
        real(kind=rk),  intent(in)               :: eigenvalues(:), eigenbasis(:, :)
        integer(kind=ik),optional, intent(inout) :: truncation_rank           !> number of POD modes
        !> Threshold value for truncating POD modes. If the singular value is smaller,tree_ID_dest
        !> then the given threshold we discard the corresponding POD MODE.
        real(kind=rk), optional, intent(in)      :: truncation_error
        !---------------------------------------------------------------
        real(kind=rk)                            ::  max_err, t_elapse, Volume
        real(kind=rk) ,allocatable               ::  alpha(:), a_coefs(:,:)
        integer(kind=ik):: N_snapshots, root, ierr, i, rank, pod_mode_tree_ID, &
        tree_ID, N_modes, max_nr_pod_modes, it, j
        character(len=cshort)                        :: filename

        ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
        ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
        ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
        ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
        ! to the last subroutine.)  -Thomas

        !---------------------------------------------------------------------------
        ! check inputs and set default values
        !---------------------------------------------------------------------------
        rank= params%rank
        N_snapshots=size(eigenvalues)
        allocate(alpha(N_snapshots),a_coefs(N_snapshots,N_snapshots))
        alpha=0.0_rk
        a_coefs=0.0_rk

        if (present(truncation_rank)) then
            max_nr_pod_modes = truncation_rank
        else
            max_nr_pod_modes = N_snapshots
        endif

        if (present(truncation_error)) then
            max_err = truncation_error
        else
            max_err = 0.0_rk
        endif

        if (params%forest_size <= N_snapshots + 3) call abort(1003191,"Error! Need more Trees. Tip: increase forest_size")
        if (rank == 0) then
            write(*,*) "----v eigenvalues v-----"
            do i = 1, N_snapshots
                write(*,'(1(es12.4,1x))') eigenvalues(i)
            enddo
            write(*,*) "----^ eigenvalues ^-----"
            write(*,'("sum(eigs) = ", g18.8)') sum(eigenvalues)
            write(*,*)
        endif

        !---------------------------------------------------------------------------
        ! construct POD basis functions (modes)
        !---------------------------------------------------------------------------
        N_modes = 0
        pod_mode_tree_ID= N_snapshots + 1
        if ( rank == 0 ) write(*,*) "Constructing POD modes (X*V)"
        do i = N_snapshots, 1, -1
            ! compute normalized eigenvectors. If eigenvalues are to small discard modes
            ! We loop in an ascending order since the eigenvalues are sorted from smallest to
            ! largest.
            ! If the max nr of modes (=desired POD rank) or the desired precission (smallest eigenvalue=sigma**2)
            ! is reached, we exit the for loop, since we have build enough POD modes.

            if ( eigenvalues(i) < max_err .or. N_modes > max_nr_pod_modes-1 ) exit
            t_elapse = MPI_wtime()

            N_modes = N_modes +1
            alpha = eigenbasis(:,i)/sqrt(dble(N_snapshots)*eigenvalues(i))
            ! calculate pod modes:
            call copy_tree(params, hvy_block, pod_mode_tree_ID, 1)
            call multiply_tree_with_scalar(params, hvy_block, pod_mode_tree_ID, alpha(1))

            ! Linear Combination of the snapshots to build the i-th POD MODE
            do tree_ID = 2, N_snapshots
                call add_two_trees(params, hvy_block, hvy_tmp, pod_mode_tree_ID, tree_ID, b=alpha(tree_ID))
            end do

            if ( params%adapt_tree) then
                call adapt_tree( 0.0_rk, params, hvy_block, pod_mode_tree_ID, params%coarsening_indicator, hvy_tmp )
            endif

            t_elapse = MPI_WTIME() - t_elapse
            if (rank == 0) then
                write(*,'("POD mode ", i4," constructed in t_cpu=",es12.4, "sec [Jmin,Jmax]=[",i2,",",i2,"] Nblocks=",i12)') &
                N_modes, t_elapse,&
                minActiveLevel_tree( pod_mode_tree_ID ), &
                maxActiveLevel_tree( pod_mode_tree_ID ), &
                sum(lgt_n(1:N_snapshots+2))
            endif

            !----------------------------------
            ! Save Modes
            !----------------------------------
            do j = 1, params%n_eqn
                write( filename, '("mode",i1,"_", i12.12, ".h5")') j, N_modes
                call saveHDF5_tree(filename, 0.0_rk, 0_ik, j, params, hvy_block, tree_ID)
            end do

            !---------------------------------------------------------------------------
            ! temporal coefficients
            !---------------------------------------------------------------------------
            t_elapse = MPI_WTIME()
            do it = 1, N_snapshots
                a_coefs(it,N_modes) = scalar_product_two_trees( params, hvy_block, hvy_tmp, it, pod_mode_tree_ID)
            enddo
            t_elapse = MPI_WTIME() - t_elapse
            if ( rank == 0 ) write(*,'("Time Coefficient Mode ",i3," computed in t_cpu=", es12.4, "sec" )') &
            N_modes, t_elapse

        end do

        Volume = product(params%domain_size(1:params%dim))
        ! V is the matrix of eigenvectors
        if (Volume>0.0_rk) then
            a_coefs = a_coefs / Volume
        else
            call abort(012301,"something wrong")
        end if
        ! the truncation_rank can be lower then the default (N_snapshots) or input value,
        ! when the singular values are smaller as the desired presicion! Therefore we update
        ! the truncation rank here.
        truncation_rank = N_modes

        if (rank==0) then
            write(*,*)
            filename = "a_coefs.txt"
            write(*,'( "Temporal coefficients saved to file: ", A30 )') filename
            open(14, file=filename, status='replace')
            do i = 1, N_snapshots
                write(14,'(400(es15.8,1x))') a_coefs(i, 1:truncation_rank)
            enddo
            close(14)
        end if

        if (rank == 0) then
            write(*, *)
            write(*,'(80("─"))')
            write(*,'("Estimated Truncation Rank r= ",i4)') truncation_rank
            if (truncation_rank < N_snapshots ) then
                if ( eigenvalues(N_snapshots-truncation_rank)>0) &
                write(*,'("Error in L2 norm (sigma_{r+1}): ",es12.4)') sqrt(eigenvalues(N_snapshots-truncation_rank))
            else
                write(*,'("Error in L2 norm (sigma_{r+1}): 0")')
            endif
            write(*,'(80("─"))')
            write(*, *)
        endif

    end subroutine
    !##############################################################



    ! subroutine  print_mat(Mat)

    !     real(kind=rk), intent(in)       :: Mat(:,:)
    !     integer                         :: i, j
    !     character(len=16) :: fmt
    !     character(len=3) :: ncols_str

    !     write(*,*) " "
    !     write(ncols_str,'(i3.3)') size(Mat,2)
    !     fmt = '('//ncols_str//'(es12.4,1x))'
    !     do j = 1,size(Mat,1)
    !         write(*,fmt) Mat(j,:)
    !     end do

    ! end subroutine print_mat


    !##############################################################
    !> Main routine of the POD postprocessing procedure.
    !> If post_POD is called from wabbit-post it will read in files
    !> specified on input and use it as snapshot data. After
    !> decomposition the Modes will be safed to a file!
    subroutine post_POD(params)
        use module_globals
        use module_params
        use module_mpi

        implicit none

        !> parameter struct
        !--------------------------------------------
        type (type_params), intent(inout)  :: params
        !--------------------------------------------
        character(len=cshort)      :: file_out, args,order
        character(len=cshort),dimension(:), allocatable :: fname_list,  eigenbasis_files
        character(len=cshort),dimension(:,:), allocatable :: file_in
        real(kind=rk)   , allocatable :: time(:),M(:,:),eigenvectors(:, :), eigenvalues(:,:)
        integer(kind=ik), allocatable :: iteration(:)
        real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
        real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :)
        integer(kind=ik)                        :: level, k, Bs(3), tc_length
        integer(hid_t)                          :: file_id
        real(kind=rk), dimension(3)             :: domain
        integer(hsize_t), dimension(2)          :: dims_treecode
        integer(kind=ik) :: treecode_size, number_dense_blocks, tree_ID, truncation_rank_in = -1
        integer(kind=ik) :: i, N_snapshots, dim, fsize, lgt_n_tmp, truncation_rank = 3
        integer(kind=ik) :: j, n_components=1, io_error
        real(kind=rk) :: truncation_error=0.0_rk, truncation_error_in=-1.0_rk, maxmem=-1.0_rk, &
        eps=-1.0_rk, L2norm, Volume,t_elapse(2)
        logical :: verbose = .false., save_all = .true., start_from_eigenbasis=.false.
        character(len=30) :: rowfmt

        ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
        ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
        ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
        ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
        ! to the last subroutine.)  -Thomas


        call get_command_argument(2, args)
        if ( args== '--help' .or. args == '--h') then
            if ( params%rank==0 ) then
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) "mpi_command -n number_procs ./wabbit-post --POD --components=3 --list=filelist.txt [list_uy.txt] [list_uz.txt]"
                write(*,*) "[--save_all --order=CDF[2|4]0 --nmodes=3 --error=1e-9 --adapt=0.1 --eps-norm=[L2|Linfty] --start_from_eigenbasis='eigenvalues.txt,eigenvectors.txt']"
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) " Wavelet adaptive Snapshot POD "
                write(*,*) " --list               list of files containing all snapshots"
                write(*,*) " --components         number of components in statevector"
                write(*,*) " --save_all           reconstruct all snapshots (default is false)"
                write(*,*) " --order              order of the predictor"
                write(*,*) " --adapt              threshold for wavelet adaptation of modes and snapshot"
                write(*,*) " --eps-norm           normalization of wavelets"
                write(*,*) "--start_from_eigenbasis  when you want to restart the computation of the modes but you"
                write(*,*) "                         have allready computed the diagonlaized covariance matrix."
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            end if
            return
        end if

        !----------------------------------
        ! read parameters
        !----------------------------------
        call get_cmd_arg_dbl( "--adapt", eps, default=-1.0_rk )
        call get_cmd_arg( "--eps-norm", params%eps_norm, default="L2" )
        call get_cmd_arg( "--order", order, default="CDF44" )
        call get_cmd_arg( "--nmodes", truncation_rank_in, default=-1_ik )
        call get_cmd_arg( "--error", truncation_error_in, default=-1.0_rk )
        call get_cmd_arg( "--list", fname_list )
        call get_cmd_arg( "--memory", args, default="-2.0GB")
        read(args(1:len_trim(args)-2),* ) maxmem
        call get_cmd_arg( "--save_all", save_all, default=.true.)
        call get_cmd_arg( "--start_from_eigenbasis", eigenbasis_files)
        call get_cmd_arg( "--components", n_components, default=1_ik)

        !-------------------------------
        ! Set some wabbit specific params
        !-------------------------------
        ! block distirbution:
        params%block_distribution="sfc_hilbert"
        ! no time stepping:
        params%time_step_method="none"
        params%Jmin=1
        ! coarsening indicator
        params%physics_type="POD"
        params%eps_normalized=.True. ! normalize the statevector before thresholding
        params%coarsening_indicator="threshold-state-vector"
        params%threshold_mask=.False.

        call set_wavelet_params(params, order )

        if (eps > 0) then
            ! adapt the mesh if possible
            params%adapt_tree = .True.! .False.!.True.
            params%eps=eps
        else
            params%adapt_tree = .False.
            params%eps = 0.0_rk
        endif

        !-------------------------------
        ! set defaults for POD truncation
        !-------------------------------
        if (truncation_rank_in /=-1) truncation_rank = truncation_rank_in
        if (truncation_error_in >0.0_rk ) truncation_error = truncation_error_in
        if (allocated(eigenbasis_files)) start_from_eigenbasis = .True.

        call count_lines_in_ascii_file_mpi(fname_list(1), N_snapshots, n_header=0)

        !--------------------------
        ! check if files exists:
        !-------------------------------
        ! fname_list is the name of the file which contains a list of files you want to
        ! read in (i.e. rho_000000000.h5 rho_000000001.h5 ... etc)
        if (.not. allocated(fname_list)) call abort(207191,"you must pass at least one file list! Use: --list my_filelist.txt")
        do j = 1, n_components
            call check_file_exists ( fname_list(j) )
            if (params%rank==0) write(*,*) "Reading list of files from "//fname_list(j)
        enddo

        if (start_from_eigenbasis) then
            if (params%rank==0) write(*,*) "Eigenbasis read from: "
            do j = 1, 2
                call check_file_exists ( eigenbasis_files(j) )
                if (params%rank==0) write(*,*) trim(eigenbasis_files(j))
            enddo
        endif
        !-----------------------------------------------------------------------------
        ! read in the file, loop over lines
        !-----------------------------------------------------------------------------
        allocate(params%input_files(n_components))
        allocate(file_in(N_snapshots,n_components))
        allocate(time(N_snapshots))
        allocate(iteration(N_snapshots))
        if (start_from_eigenbasis) allocate(eigenvalues(N_snapshots, 2), &
        eigenvectors(N_snapshots,N_snapshots))
        !-------------------------------------------
        ! check and find common params in all h5-files
        !-------------------------------------------
        ! open all files
        do j = 1, n_components
            open( unit=10+j, file=fname_list(j), action='read', status='old' )
        enddo

        io_error = 0
        i = 1
        do while(i <= N_snapshots)
            do j = 1, n_components
                read (10+j, '(A)', iostat=io_error) file_in(i,j)

                call check_file_exists ( file_in(i,j) )

                if ( i == 1 .and. j == 1 ) then
                    ! read all geometric parameters of grid
                    call read_attributes(file_in(i,j), lgt_n_tmp, time(1), iteration(1), params%domain_size, &
                    params%Bs, params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
                endif
                call read_attributes(file_in(i,j), lgt_n_tmp, time(i), iteration(i), domain, bs, level, dim)

                params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot

                if (any(params%Bs /= Bs)) call abort( 203191, " Block size is not consistent ")

                if (abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")

                if (params%dim /= dim) call abort( 203193, "Dimension is not consistent ")
            end do
            i = i + 1
        end do

        ! now we have all information to allocate the grid and set up the forest:
        if (truncation_rank_in /= -1) then
            fsize = N_snapshots + truncation_rank + 2 !we need some extra fields for storing etc
        else
            fsize = 2*N_snapshots + 2 !we need some extra fields for storing etc
        end if

        params%forest_size = fsize
        number_dense_blocks = 2_ik**(dim*params%Jmax)*fsize
        params%n_eqn = n_components
        allocate(params%threshold_state_vector_component(params%n_eqn))
        params%threshold_state_vector_component(1:params%n_eqn)=1
        if (maxmem < 0.0_rk) then
            params%number_blocks = ceiling( 4.0_rk * N_snapshots * number_dense_blocks / params%number_procs * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk) )
        endif

        call adjust_n_ghosts_for_scalarproduct(params)

        !----------------------------------
        ! allocate data
        !----------------------------------
        call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

        hvy_neighbor = -1_ik
        tree_n= 0_ik ! reset number of trees in forest
        !----------------------------------
        ! READ ALL SNAPSHOTS
        !----------------------------------
        t_elapse(1) = MPI_WTIME()
        do tree_ID = 1, N_snapshots
            call readHDF5vct_tree(file_in(tree_ID,:), params, hvy_block, tree_ID)
            !----------------------------------
            ! Adapt the data to the given eps
            !----------------------------------
            if (params%adapt_tree) then
                ! now, evaluate the refinement criterion on each block, and coarsen the grid where possible.
                ! adapt-mesh also performs neighbor and active lists updates
                call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp )
            endif
        end do
        ! elapsed time for reading and coarsening the data
        t_elapse(1)= MPI_Wtime()-t_elapse(1)

        ! --------------------
        ! Compute the L2 norm
        !---------------------
        L2norm = compute_L2norm(params, hvy_block, hvy_tmp, N_snapshots, verbose)
        Volume = product(domain(1:params%dim))

        if (params%rank==0) then
            write(*,'(80("─"))')
            write(*,*) "WABBIT POD."
            write(*,*) "Snapshot matrix X build from:"
            do i = 1, N_snapshots
                write(*, '("Snapshot=",i3 , " it=",i7,1x," time=",f16.9,1x," Nblocks=", i6," sparsity=(",f5.1,"% / ",f3.1,"%) [Jmin,Jmax]=[",i2,",",i2,"]")')&
                i, iteration(i), time(i), lgt_n(i), &
                100.0*dble(lgt_n(i))/dble( (2**maxActiveLevel_tree(i))**params%dim ), &
                100.0*dble(lgt_n(i))/dble( (2**params%Jmax)**params%dim ), &
                minActiveLevel_tree(i), &
                maxActiveLevel_tree(i)
                write(*, '("Files:")')
                do j = 1, n_components
                    write(*, '(i2, "   ", A)') j, trim(file_in(i,j))
                end do
            end do
            write(*,'(80("─"))')
            write(*,'("L2 norm ||X||_2^2/N_snapshots/Volume: ",g18.8)') L2norm**2/dble(N_snapshots)/Volume
            write(*,'("Data dimension: ",i1,"D")') params%dim
            write(*,'("Domain size is ",3(g12.4,1x))') domain
            write(*,'("NCPU=",i6)') params%number_procs
            write(*,'("Number Trees=",1x,i4)') fsize
            write(*,'("[Jmin,Jmax] =[",i2,",",i2,"]")')params%Jmin, params%Jmax
            write(*,'("Nblocks Available from Memory =",i9)') params%number_blocks
            write(*,'("Nblocks (if all trees dense)=",i9)') number_dense_blocks
            write(*,'("Nblocks used (sparse)=",i6)') sum(lgt_n(1:tree_n))
            write(*,'("Predictor=",A)') params%order_predictor
            write(*,'("wavelet=",A)') params%wavelet
            write(*,'("Number ghosts=",i2)') params%g
            write(*,'("block_distribution=",A)') params%block_distribution
            write(*,'(80("─"))')
        endif

        !----------------------------------
        ! COMPUTE POD Modes
        !----------------------------------
        t_elapse(2) = MPI_WTIME()
        if (start_from_eigenbasis) then
            ! eigenvalues.txt contains the number of the eigenvalue in the first column and its value in the second
            call read_array_from_ascii_file_mpi(eigenbasis_files(1), eigenvalues, n_header=0_ik)
            call read_array_from_ascii_file_mpi(eigenbasis_files(2), eigenvectors, n_header=0_ik)
            call construct_modes(params, hvy_block, hvy_tmp, truncation_error, truncation_rank, eigenvectors, eigenvalues(:,2))
        else
            call snapshot_POD( params, hvy_block, hvy_tmp, truncation_error, truncation_rank, save_all)
        endif
        t_elapse(2) = MPI_Wtime()-t_elapse(2)
        call toc( "post_POD (read + coarse): ", 350, t_elapse(1) )
        call toc( "post_POD (snapshot_POD): ", 351, t_elapse(2) )
        call toc( "post_POD (total wPOD algo): ", 352, t_elapse(1) + t_elapse(2) )

    end subroutine
    !##############################################################



    !##############################################################
    subroutine post_PODerror(params)
        use module_globals
        use module_params
        use module_mpi

        implicit none

        !> parameter struct
        !--------------------------------------------
        type (type_params), intent(inout)  :: params
        !--------------------------------------------
        character(len=cshort) :: fname_acoefs, filename, args, order
        character(len=cshort),dimension(:), allocatable   :: fsnapshot_list, fmode_list
        character(len=100),dimension(:,:), allocatable ::snapshot_in, mode_in
        real(kind=rk)   , allocatable :: L2error(:),time(:)
        integer(kind=ik), allocatable :: iter_list(:)
        real(kind=rk),    allocatable :: hvy_block(:, :, :, :, :)
        real(kind=rk),    allocatable :: hvy_tmp(:, :, :, :, :),a_coefs(:,:)
        integer(kind=ik)                 :: level, k, Bs(3), tc_length
        integer(hid_t)                   :: file_id
        real(kind=rk), dimension(3)      :: domain
        integer(hsize_t), dimension(2)   :: dims_treecode
        integer(kind=ik) :: treecode_size, number_dense_blocks, tree_ID, dF, N_modes, min_lvl
        integer(kind=ik) :: i, n_opt_args, N_snapshots, dim, fsize, lgt_n_tmp, r, iteration=-1
        integer(kind=ik) :: j, n_components=1, io_error, reconst_tree_ID, unused_int
        real(kind=rk)    :: maxmem=-1.0_rk, eps=-1.0_rk, L2norm, L2norm_snapshots, Volume, norm
        real(kind=rk)    :: unused_var
        logical :: verbose = .false., save_all = .false., all_snapshots_dense = .false.

        ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
        ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
        ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
        ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
        ! to the last subroutine.)  -Thomas


        call get_command_argument(2, args)
        if ( args== '--help' .or. args == '--h') then
            if ( params%rank==0 ) then
                write(*,*)
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) " Compute L2-error of wPOD "
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) "Call mpi_command -n number_procs ./wabbit-post --POD-error --time_coefficients=a_coefs.txt --components=3"
                write(*,*) '                                --snapshot-list="list_u1.txt [list_u2.txt] [list_u3.txt]"'
                write(*,*) '                                --mode-list="list_m1.txt [list_m2.txt] [list_m3.txt]"'
                write(*,*) "                                [--order=[CDF40] --adapt=0.1 --memory=2GB --iteration=1]"
                write(*,*) "                                [--save_all ,--eps-norm=[L2|Linfty]]"
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) " --time_coefficients  list of temporal coefficients for all modes"
                write(*,*) " --snapshot-list      list of files containing all snapshots"
                write(*,*) " --mode-list          list of files containing all modes for reconstruction"
                write(*,*) " --components         number of components in statevector"
                write(*,*) " --save_all           reconstruct all snapshots (default is false)"
                write(*,*) " --iteration          reconstruct single snapshot (default is no snapshots)"
                write(*,*) " --order              order of the predictor"
                write(*,*) " --adapt              threshold for wavelet adaptation of modes and snapshot"
                write(*,*) " --eps-norm           normalization of wavelets"
                write(*,*) " --iteration          reconstruct single snapshot"
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            end if
            return
        end if

        !----------------------------------
        ! read predefined params
        !----------------------------------
        call get_cmd_arg_str( "--eps-norm", params%eps_norm, default="L2" )
        call get_cmd_arg_dbl( "--adapt", eps, default=-1.0_rk )
        call get_cmd_arg_str( "--order", order, default="CDF44" )
        call get_cmd_arg_str_vct( "--snapshot-list", fsnapshot_list )
        call get_cmd_arg_str_vct( "--mode-list", fmode_list )
        call get_cmd_arg_str( "--time_coefficients", fname_acoefs, default= "acoefs.txt" )
        call get_cmd_arg( "--memory", args, default="2GB")
        read(args(1:len_trim(args)-2),* ) maxmem
        call get_cmd_arg( "--save_all", save_all, default=.False.)
        call get_cmd_arg( "--components", n_components, default=1_ik)
        call get_cmd_arg( "--iteration", iteration, default=1_ik)

        if ( iteration>0 ) then
            if ( params%rank == 0 ) write(*,*) "Iteration reconstructed: " ,iteration
        endif

        if ( save_all ) then
            if ( params%rank == 0 ) write(*,*) "Reconstructing all snapshots!"
        endif
        !-------------------------------
        ! Set some wabbit specific params
        !-------------------------------
        ! block distirbution:
        params%block_distribution="sfc_hilbert"
        ! no time stepping:
        params%time_step_method="none"
        params%Jmin=1
        ! coarsening indicator
        params%physics_type="POD"
        params%eps_normalized=.True. ! normalize the statevector before thresholding
        params%coarsening_indicator="threshold-state-vector"
        params%threshold_mask=.False.
        params%n_eqn=n_components

        call set_wavelet_params(params, order )

        if ( eps > 0.0_rk) then
            ! adapt the mesh if possible
            params%adapt_tree = .True.! .False.!.True.
            params%eps=eps
        else
            params%adapt_tree = .False.
            params%eps = 0.0_rk
        endif


        !-------------------------------
        ! check if files exists:
        !-------------------------------
        ! fname_list is the name of the file which contains a list of files you want to
        ! read in (i.e. rho_000000000.h5 rho_000000001.h5 ... etc)
        if (.not. allocated(fsnapshot_list)) call abort(207191,"you must pass at least one snapshot list! Use: --snaphot-list my_filelist.txt")
        if (.not. allocated(fmode_list))     call abort(207191,"you must pass at least one mode list! Use: --mode-list my_filelist.txt")
        do j = 1, n_components
            call check_file_exists ( fsnapshot_list(j) )
            call check_file_exists ( fmode_list(j) )
            if (params%rank==0) write(*,*) "Snapshot list: "//trim(fsnapshot_list(j)) //"   Mode list: "// fmode_list(j)
        enddo
        call check_file_exists(fname_acoefs)
        call count_lines_in_ascii_file_mpi(fname_acoefs, N_snapshots, n_header=0_ik)
        call count_lines_in_ascii_file_mpi(fmode_list(1), N_modes, n_header=0_ik)

        !-----------------------------------------------------------------------------
        ! ALLOCATE additional arrays
        !-----------------------------------------------------------------------------
        allocate(a_coefs(1:N_snapshots, 1:N_modes))
        allocate(L2error(1:N_modes))
        allocate(params%input_files(n_components))
        allocate(mode_in(1:N_modes,1:n_components))
        allocate(snapshot_in(N_snapshots,n_components))
        allocate(time(N_snapshots))
        allocate(iter_list(N_snapshots))
        !----------------------------------
        ! READ a_coefficient from file
        !----------------------------------
        call read_array_from_ascii_file_mpi(fname_acoefs, a_coefs, n_header=0_ik)
        call count_lines_in_ascii_file_mpi(fsnapshot_list(1), N_snapshots, n_header=0_ik)

        !-------------------------------------------
        !-------------------------------------------
        ! SNAPSHOTS
        !-------------------------------------------
        !-------------------------------------------
        ! check and find common params in all h5-files
        ! of all snapshots
        !-------------------------------------------
        ! open all files
        do j = 1, n_components
            open( unit=10+j, file=fsnapshot_list(j), action='read', status='old' )
        enddo

        io_error = 0
        i = 1
        do while(i <= N_snapshots)
            do j = 1, n_components
                read (10+j, '(A)', iostat=io_error) snapshot_in(i,j)

                call check_file_exists ( snapshot_in(i,j) )

                if ( i == 1 .and. j == 1 ) then
                    ! read all geometric parameters of grid
                    call read_attributes(snapshot_in(i,j), lgt_n_tmp, time(1), iter_list(1), params%domain_size, &
                    params%Bs, params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
                endif
                call read_attributes(snapshot_in(i,j), lgt_n_tmp, time(i), iter_list(i), domain, bs, level, dim)
                params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot
                if (any(params%Bs .ne. Bs)) call abort( 204191, " Block size is not consistent ")
                if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")
                if (params%dim .ne. dim) call abort( 203193, "Dimension is not consistent ")
            end do
            i = i + 1
        end do

        !--------------------------------
        ! close the list of files
        !--------------------------------
        do j = 1, n_components
            close(10+j)
        enddo
        !-------------------------------------------
        !-------------------------------------------
        ! MODES
        !-------------------------------------------
        !-------------------------------------------
        !--------------------------------
        ! MODES: open all the given txt-files
        !--------------------------------
        do j = 1, n_components
            open( unit=15+j, file=fmode_list(j), action='read', status='old' )
        enddo
        !--------------------------------
        ! in every txt-file there is a list
        ! of h5-files which contain the
        ! statevector components of the
        ! POD modes
        !--------------------------------
        do i = 1, N_modes
            do j = 1, n_components
                ! here we read in each line of file fname_list(j)
                read (15+j, '(A)', iostat=io_error) mode_in(i,j)
                call check_file_exists ( mode_in(i,j) )
                !-------------------------------------------
                ! check and find common params in all h5-files
                !-------------------------------------------
                call read_attributes(mode_in(i,j), lgt_n_tmp, unused_var, unused_int, domain, bs, level, dim)
                params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot
                if (any(params%Bs .ne. Bs)) then
                    write(*,*) "Bs: ",Bs, " vs ",params%Bs;
                    call abort( 203199, " Block size is not consistent ")
                end if
                if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")
                if (params%dim .ne. dim) call abort( 204193, "Dimension is not consistent ")
            end do
        enddo
        !--------------------------------
        ! close the list of files
        !--------------------------------
        do j = 1, n_components
            close(15+j)
        enddo

        ! now we have all information to allocate the grid and set up the forest:
        fsize = N_snapshots + N_modes + 2 !we need some extra fields for storing etc
        params%forest_size = fsize
        number_dense_blocks = 2_ik**(dim*params%Jmax)*fsize
        params%n_eqn = n_components
        allocate(params%threshold_state_vector_component(params%n_eqn))
        params%threshold_state_vector_component(1:params%n_eqn)=1
        if (maxmem < 0.0_rk) then
            params%number_blocks = ceiling( 4.0_rk *fsize * number_dense_blocks / params%number_procs )
        endif

        call adjust_n_ghosts_for_scalarproduct(params)
        !----------------------------------
        ! allocate data
        !----------------------------------
        call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

        call reset_forest(params)


        hvy_neighbor = -1
        lgt_n = 0 ! reset number of acitve light blocks
        tree_n= 0 ! reset number of trees in forest
        !----------------------------------
        ! READ ALL SNAPSHOTS
        !----------------------------------
        min_lvl = 9999
        do tree_ID = 1, N_snapshots
            call readHDF5vct_tree(snapshot_in(tree_ID,:), params, hvy_block, tree_ID)
            min_lvl = min( min_lvl, minActiveLevel_tree(tree_ID) )
        end do

        if (min_lvl==params%Jmax) then
            all_snapshots_dense = .True.
            params%adapt_tree=.False.
            ! it is faster to bring all modes to the same mesh at this point,
            ! since otherwise they will be refined and coarsed when reconstructing
            ! the snapshots. this produces a lot of overhead and makes the calculations slow
            ! The results will stay the same.
            if ( params%rank==0 ) write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            if ( params%rank==0 ) write(*,*) "All modes will be refined to dense grid !!"
            if ( params%rank==0 ) write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        else
            all_snapshots_dense = .False.
        endif
        !----------------------------------
        ! READ ALL MODES
        !----------------------------------
        do j = 1, N_modes
            tree_ID = N_snapshots+j
            call readHDF5vct_tree(mode_in(tree_ID,:), params, hvy_block, tree_ID)
            ! it is faster to bring all modes to the same mesh at this point,
            ! since otherwise they will be refined and coarsed when reconstructing
            ! the snapshots. this produces a lot of overhead and makes the calculations slow
            if ( all_snapshots_dense ) then
                call refineToEquidistant_tree(params, hvy_block, hvy_tmp, tree_ID, target_level=params%Jmax)
            endif
            if ( params%rank == 0 ) then
                write(*,'("Mode", i3," Stored in Tree_id: ",i3)') j, tree_ID
                write(*,'("Number blocks ",i12)') lgt_n(tree_ID)
            endif
        end do

        if ( all_snapshots_dense ) params%Jmin = params%Jmax
        ! --------------------
        ! Compute the L2 norm
        !---------------------
        L2norm_snapshots = compute_L2norm(params, hvy_block, hvy_tmp, N_snapshots, verbose)
        Volume = product(domain(1:params%dim))
        if (params%rank==0) then
            write(*,'(80("─"))')
            write(*,*) "WABBIT POD."
            write(*,*) "Snapshot matrix X build from:"
            do i = 1, N_snapshots
                write(*, '("Snapshot=",i3 , " time=",f16.9,1x," Nblocks=", i6," sparsity=(",f5.1,"% / ",f5.1,"%) [Jmin,Jmax]=[",i2,",",i2,"]")')&
                i, time(i), lgt_n(i), &
                100.0*dble(lgt_n(i))/dble( (2**maxActiveLevel_tree(i))**params%dim ), &
                100.0*dble(lgt_n(i))/dble( (2**params%Jmax)**params%dim ), &
                minActiveLevel_tree(i), &
                maxActiveLevel_tree(i)
            end do
            write(*,*) "Modes used for reconstruction:"
            do i = 1, N_modes
                do j = 1, params%n_eqn
                    write(*, '("file=",A )') mode_in(i,j)
                end do
            end do
            write(*,'(80("─"))')
            write(*,'("L2 norm ||X||_2^2: ",g18.8)') L2norm_snapshots**2
            write(*,'("Data dimension: ",i1,"D")') params%dim
            write(*,'("Domain size is ",3(g12.4,1x))') domain
            write(*,'("NCPU=",i6)') params%number_procs
            write(*,'("Number Trees=",1x,i4)') fsize
            write(*,'("[Jmin,Jmax] =[",i2,",",i2,"]")')params%Jmin, params%Jmax
            write(*,'("Nblocks Available from Memory =",i6)') params%number_blocks
            write(*,'("Nblocks (if all trees dense)=",i6)') number_dense_blocks
            write(*,'("Nblocks used (sparse)=",i6)') sum(lgt_n(1:tree_n))
            write(*,'("Predictor=",A)') params%order_predictor
            write(*,'("wavelet=",A)') params%wavelet
            write(*,'("Number ghosts=",i2)') params%g
            write(*,'("block_distribution=",A)') params%block_distribution
            write(*,'(80("─"))')
        endif

        !---------------------------------------------------------------
        !---------------------------------------------------------------
        !          HERE WE DO THE ACTUAL error COMPUTATION
        !
        !   rel error = (sum_i \norm{u_i(x) - \tilde u_i(x)}_{L2}^2)
        !                        /(sum_i \norm{u_i(x)}_{L2}^2)
        !
        ! 1. compute the reconstructed snapshots \tilde{u}(x,t)
        !   for different numbers of modes r = 1,...,Nmodes
        ! 2. substract original snapshot by the reconstructed one
        !                  \delta u = u_i(x)-\tilde{u}_i(x)
        ! 3. compute the error for each time: err_i = \norm{delta u(x,t_i)}
        ! 4. add all errors sum_i(err_i)
        ! 5. normalize it for the relative error
        !---------------------------------------------------------------
        !---------------------------------------------------------------
        if (tree_n < fsize) then
            reconst_tree_ID = tree_n + 1
        else
            call abort(2808191, "Take care of the trees, they will take care of you.")
        end if
        do r = 1, N_modes
            L2norm = 0.0_rk
            do j = 1, N_snapshots

                call reconstruct_iteration( params, hvy_block, hvy_tmp, &
                a_coefs, j , r, reconst_tree_ID, opt_offset_mode_tree_ID=N_snapshots)

                !---------------------------------
                ! save the reconstructed snapshot
                if (j == iteration .or. save_all) then
                    ! loop over every component
                    do dF = 1, n_components
                        ! filename is constructed from: reconst<dF>-<r>_<time>.h5
                        !                   - componentnumber dF
                        !                   - number of modes used for reconstruction r
                        !                   - time of reconstructed snapshot
                        write( filename, '("reconst",i1,"-",i3.3,"_", i12.12, ".h5")') dF, r, nint(time(j) * 1.0e6_rk)
                        call saveHDF5_tree(filename, time(j), iter_list(j), dF, params, hvy_block, reconst_tree_ID)
                    end do
                end if
                !---------------------------------

                call substract_two_trees(params, hvy_block, hvy_tmp, reconst_tree_ID, j)
                ! compute L2 norm
                norm = compute_tree_L2norm(params, hvy_block, hvy_tmp, reconst_tree_ID, verbose)
                L2norm = L2norm + norm**2
            end do
            L2error(r) = L2norm/L2norm_snapshots**2

            if (params%rank == 0) then
                write(*, *)
                write(*,'(80("─"))')
                write(*,'(30("#"), " Error using Nmodes = ", i4 ,30("#"))') r
                write(*,'(80("─"))')
                write(*,'("relative L2 error ",f12.8)') L2error(r)
                write(*,'(80("─"))')
                write(*, *)
            endif

        end do

        !--------------------------
        ! save L2error
        !--------------------------
        if (params%rank==0 ) then
            write(*,*)
            filename = "L2error.txt"
            write(*,'( "L2error of wPOD is saved to file: ", A30 )') filename
            open(14, file=filename, status='replace')
            do r = 1, N_modes
                write(14,'(es15.8,1x)') L2error(r)
            end do
            close(14)
        end if


    end subroutine
    !###############################################################



    !##############################################################
    !> Compute the L2 norm of the snapshotmatrix
    function compute_L2norm(params, hvy_block, hvy_tmp, N_snapshots, verbosity) &
        result(L2norm)

        implicit none
        !-----------------------------------------------------------------
        ! inputs:
        type (type_params), intent(inout)   :: params
        integer(kind=ik), intent(in)        :: N_snapshots
        real(kind=rk),  intent(inout)       :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
        !> heavy temp data: needed in blockxfer which is called in add_two_trees
        real(kind=rk),   intent(inout)      :: hvy_tmp(:, :, :, :, :)
        logical, intent(in),optional        :: verbosity !< if true: additional information of processing
        !-----------------------------------------------------------------
        ! result
        real(kind=rk)                       :: L2norm
        !-----------------------------------------------------------------
        integer(kind=ik)                    :: lgt_id, hvy_id, ierr, tree_ID = 1
        integer(kind=ik)                    :: k, rank, g, Bs(3), dF=1
        real(kind=rk)                       :: norm, x0(3), dx(3), Volume
        logical                             :: verbose=.false.

        ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
        ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
        ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
        ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
        ! to the last subroutine.)  -Thomas


        if (present(verbosity)) verbose = verbosity
        rank = params%rank
        Bs= params%Bs
        g = params%g
        Volume = product(params%domain_size(1:params%dim))
        L2norm = 0.0_rk
        ! Loop over the active hvy_data
        do tree_ID =1, N_snapshots
            norm = compute_tree_L2norm( params, hvy_block, hvy_tmp, tree_ID, verbose)
            L2norm = L2norm + norm**2
        end do

        L2norm =sqrt(L2norm)
        if (params%rank == 0 .and. verbose ) write(*,*) "sum(eigs)= ", L2norm**2/dble(N_snapshots)/Volume

    end function
    !##############################################################


    !##############################################################
    !> construction of covariance matrix from snapshot data
    subroutine compute_covariance_matrix(params, C, hvy_block, hvy_tmp)
        implicit none

        !-----------------------------------------------------------------
        type (type_params), intent(inout)   :: params
        real(kind=rk),  intent(inout)       :: C(:,:)                     !> Covariance matrix for snapshot POD
        real(kind=rk),  intent(inout)       :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
        !> heavy temp data: needed in blockxfer which is called in add_two_trees
        real(kind=rk),   intent(inout)      :: hvy_tmp(:, :, :, :, :)
        !---------------------------------------------------------------
        integer(kind=ik)                    :: tree_ID1, tree_ID2, N_snapshots, rank
        real(kind=rk)                       :: C_val, Volume, t_elapse, t_inc(2)
        integer(kind=ik)  , allocatable     :: lgt_active_tmp(:,:), lgt_block_tmp(:,:),lgt_n_tmp(:)

        ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
        ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
        ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
        ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
        ! to the last subroutine.)  -Thomas


        allocate(lgt_block_tmp(size(lgt_block,1),size(lgt_block,2)))
        allocate(lgt_n_tmp(size(lgt_n,1)))
        allocate(lgt_active_tmp(size(lgt_active,1), size(lgt_active,2)))

        !-----------------------------------
        ! copy light data:
        ! Is very expensive but we dont do it offen
        lgt_active_tmp= lgt_active
        lgt_block_tmp = lgt_block
        lgt_n_tmp     = lgt_n
        !-----------------------------------

        t_elapse = MPI_wtime()
        N_snapshots = size(C,1)
        rank = params%rank
        Volume = product(params%domain_size(1:params%dim))
        ! We loop over all snapshots X_i, i=1,...,N to calculate the values of the symmetric
        ! covariance matrix C_{j,i} = C_{i,j} = <X_i, X_j>
        do tree_ID1 = 1, N_snapshots
            do tree_ID2 = tree_ID1, N_snapshots
                t_inc(1) = MPI_wtime()
                ! L2 scalarproduct between tree 1 and tree 2
                C_val = scalar_product_two_trees( params, hvy_block, hvy_tmp, tree_ID1, tree_ID2)

                ! Construct symmetric correlation matrix
                C(tree_ID1, tree_ID2) = C_val
                C(tree_ID2, tree_ID1) = C_val
                !
                t_inc(1) = MPI_WTIME() - t_inc(1)
                if (rank == 0) then
                    write(*,'("Matrixelement (i,j)= (", i3,",", i3, ") constructed in t_cpu=",es10.2, "sec Nblocks=",i10)') &
                    tree_ID1, tree_ID2, t_inc(1), sum(lgt_n(1:N_snapshots))
                endif
            end do
        end do
        ! Normalice C to the Volume of the integrated area (L2 scalar product)
        ! and by the number of snapshots
        C = C / Volume
        C = C / real(N_snapshots,kind=rk)

        call toc( "compute_covariance_matrix", 353, MPI_WTIME() - t_elapse )

    end subroutine
    !##############################################################


    !##############################################################
    !> Main routine of the POD postprocessing procedure.
    !> If post_reconstruct is called from wabbit-post it will read in files
    !> specified on input and use it as snapshot data. After
    !> decomposition the Modes will be safed to a file!
    subroutine post_reconstruct(params)
        use module_globals
        use module_params
        use module_mpi

        implicit none

        !> parameter struct
        !--------------------------------------------
        type (type_params), intent(inout)  :: params
        !--------------------------------------------
        character(len=cshort), allocatable :: file_in(:,:), fname_list(:)
        character(len=cshort) :: fname_acoefs, tmp_name, file_out, args, order
        real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
        real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :),a_coefs(:,:)
        integer(kind=ik), allocatable           :: mode_number(:)

        integer(kind=ik)                        :: level, k, Bs(3), tc_length
        integer(hid_t)                          :: file_id
        real(kind=rk), dimension(3)             :: domain
        integer(hsize_t), dimension(2)          :: dims_treecode
        integer(kind=ik) :: N_modes_used=1_ik, max_nr_modes, iteration=-1, n_components, min_lvl
        integer(kind=ik) :: treecode_size,iter, number_dense_blocks, tree_ID, reconst_tree_ID
        integer(kind=ik) :: i,j, n_opt_args, N_snapshots, dim, fsize, lgt_n_tmp, rank, io_error
        real(kind=rk) ::  maxmem=-1.0_rk, eps=-1.0_rk, Volume, tmp_time
        logical :: verbosity = .false., save_all

        ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
        ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
        ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
        ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
        ! to the last subroutine.)  -Thomas


        rank = params%rank
        call get_command_argument(2, args)
        if (args == '--help' .or. args == '--h') then
            if ( params%rank==0 ) then
                write(*,*)
                write(*,*)
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) " RECONSTRUCTION OF SNAPSHOTS WITH wPOD MODES "
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) " CALL: ./wabbit-post --POD-reconstruct --time_coefficients=a_coefs.txt --components=3 "
                write(*,*) "                     --mode-list=POD_mode1.txt [POD_mode2.txt] [POD_mode3.txt]"
                write(*,*) "                     [--save_all --iteration=10 --order=[2|4] --nmodes=3 --adapt=0.1] "
                write(*,*) "                     [--eps-norm=[L2,Linfty]] "
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                write(*,*) " --time_coefficients  list of temporal coefficients for all modes"
                write(*,*) " --components         number of components in statevector"
                write(*,*) " --mode-list          textfiles including the path to a mode computed with --POD for each component"
                write(*,*) " --save_all           reconstruct all snapshots"
                write(*,*) " --iteration          reconstruct single snapshot"
                write(*,*) " --nmodes             number of modes used for reconstruction"
                write(*,*) " --order              order of the predictor"
                write(*,*) " --adapt              threshold for wavelet adaptation of modes and snapshot"
                write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            end if
            return
        end if

        !----------------------------------
        ! read predefined params
        !----------------------------------
        call get_cmd_arg_str( "--eps-norm", params%eps_norm, default="L2" )
        call get_cmd_arg_dbl( "--adapt", eps, default=-1.0_rk )
        call get_cmd_arg_str( "--order", order, default="CDF44" )
        call get_cmd_arg_str_vct( "--mode-list", fname_list )
        call get_cmd_arg_str( "--time_coefficients", fname_acoefs, default= "a_coefs.txt" )
        call get_cmd_arg( "--memory", args, default="2GB")
        read(args(1:len_trim(args)-2),* ) maxmem
        call get_cmd_arg( "--save_all", save_all, default=.true.)
        call get_cmd_arg( "--components", n_components, default=1_ik)
        call get_cmd_arg( "--iteration", iteration, default=-1_ik)
        call get_cmd_arg( "--nmodes", N_modes_used, default=1_ik)

        if ( iteration>0 ) then
            save_all = .False.
            if ( rank == 0 ) write(*,*) "Iteration reconstructed: " ,iteration
        else
            save_all = .True.
            if ( rank == 0 ) write(*,*) "Reconstructing all snapshots!"
        endif
        !----------------------------------
        ! READ a_coefficient from file
        !----------------------------------
        call check_file_exists(fname_acoefs)
        call count_lines_in_ascii_file_mpi(fname_acoefs, N_snapshots, n_header=0_ik)
        call count_cols_in_ascii_file_mpi(fname_acoefs, max_nr_modes, n_header=0_ik)
        if (N_modes_used>max_nr_modes) call abort(270419,"Try to use more modes then available")
            if (N_modes_used< 0) N_modes_used = max_nr_modes
            allocate( a_coefs(1:N_snapshots, 1:N_modes_used))
            allocate(file_in(1:max_nr_modes,1:n_components))
            call read_array_from_ascii_file_mpi(fname_acoefs, a_coefs, n_header=0_ik)
            if (rank==0) then
                write(*,*) "~~~~v      ", fname_acoefs
                do k = 1, N_snapshots
                    write(*,'(400(es15.8,1x))') a_coefs(k, :)
                enddo
                write(*,*) "~~~~^      ", fname_acoefs
            endif

            !----------------------------------
            ! some wabbit params
            !----------------------------------
            call set_wavelet_params(params, order )

            if ( eps > 0.0_rk ) then
                ! adapt the mesh if possible
                params%adapt_tree = .True.! .False.!.True.
                params%eps=eps
            else
                params%adapt_tree = .False.
                params%eps = 0.0_rk
            endif

            !----------------------------------
            ! set addtitional params
            !----------------------------------
            ! number of components is the number of quantities saved in the
            ! statevector. Therefore we set:
            params%n_eqn = n_components
            ! block distirbution:
            params%block_distribution="sfc_hilbert"
            ! no time stepping:
            params%time_step_method="no"
            params%Jmin=1
            ! coarsening indicator
            params%physics_type="POD"
            params%eps_normalized=.True. ! normalize the statevector before thresholding
            params%coarsening_indicator="threshold-state-vector"
            params%threshold_mask=.False.
            ! read ini-file and save parameters in struct
            allocate(params%input_files(params%n_eqn))
            allocate(mode_number(N_modes_used))

            !--------------------------------
            ! open all the given txt-files
            !--------------------------------
            do j = 1, n_components
                open( unit=15+j, file=fname_list(j), action='read', status='old' )
            enddo
            !--------------------------------
            ! in every txt-file there is a list
            ! of h5-files which contain the
            ! statevector components of the
            ! POD modes
            !--------------------------------
            do i = 1, N_modes_used
                do j = 1, n_components
                    ! here we read in each line of file fname_list(j)
                    read (15+j, '(A)', iostat=io_error) file_in(i,j)

                    call check_file_exists ( file_in(i,j) )

                    !-------------------------------------------
                    ! check and find common params in all h5-files
                    !-------------------------------------------
                    if ( i == 1 .and. j == 1 ) then
                        ! read all geometric parameters of grid
                        call read_attributes(file_in(i,j), lgt_n_tmp, tmp_time, mode_number(1), params%domain_size, &
                        params%Bs, params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
                    endif
                    call read_attributes(file_in(i,j), lgt_n_tmp, tmp_time, mode_number(i), domain, bs, level, dim)
                    params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot
                    if (any(params%Bs .ne. Bs)) call abort( 203191, " Block size is not consistent ")
                    if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")
                    if (params%dim .ne. dim) call abort( 203193, "Dimension is not consistent ")
                end do
            enddo

            !--------------------------------
            ! close the list of files
            !--------------------------------
            do j = 1, n_components
                close(10+j)
            enddo
            fsize = N_modes_used +  2 !we need one more additional field for the reconstructed field
            params%forest_size = fsize
            number_dense_blocks = 2_ik**(dim*params%Jmax)*fsize
            allocate(params%threshold_state_vector_component(params%n_eqn))
            params%threshold_state_vector_component=1
            if (maxmem < 0.0_rk) then
                params%number_blocks = ceiling( 4.0_rk * fsize * number_dense_blocks / params%number_procs )
            endif

            call adjust_n_ghosts_for_scalarproduct(params)
            !----------------------------------
            ! allocate data
            !----------------------------------
            call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

            call reset_forest(params)

            hvy_neighbor = -1
            lgt_n = 0 ! reset number of acitve light blocks
            tree_n= 0 ! reset number of trees in forest
            !----------------------------------
            ! READ ALL Modes needed
            !----------------------------------
            min_lvl = 9999
            do tree_ID = 1, N_modes_used
                call readHDF5vct_tree(file_in(tree_ID,:), params, hvy_block, tree_ID)
                min_lvl = min(min_lvl, minActiveLevel_tree(tree_ID) )
            end do

            if (min_lvl == params%Jmax) then
                params%Jmin= params%Jmax
            endif



            if (params%rank==0) then
                write(*,'(80("─"))')
                write(*,*) "WABBIT POD-reconstruction."
                write(*,*) "POD modes for reconstruction:"
                do i = 1, N_Modes_used
                    write(*, '("Mode=",i3, " Nblocks=", i6," sparsity=(",f5.1,"% / ",f3.1,"%) [Jmin,Jmax]=[",i2,",",i2,"]")')&
                    mode_number(i), lgt_n(i), &
                    100.0*dble(lgt_n(i))/dble( (2**maxActiveLevel_tree( i ))**params%dim ), &
                    100.0*dble(lgt_n(i))/dble( (2**params%Jmax)**params%dim ), &
                    minActiveLevel_tree( i ), &
                    maxActiveLevel_tree( i )
                    do j = 1, n_components
                        write(*, '(i2, "   ", A)') j, trim(file_in(i,j))
                    end do
                end do
                write(*,'(80("─"))')
                write(*,'("Temporal coefficients from:",A)') fname_acoefs
                write(*,'("Number modes used/available=",1x,i4,"/",i4)') N_modes_used, max_nr_modes
                write(*,'("Data dimension: ",i1,"D")') params%dim
                write(*,'("Domain size is ",3(g12.4,1x))') domain
                write(*,'("NCPU=",i6)') params%number_procs
                write(*,'("Number Trees=",1x,i4)') fsize
                write(*,'("[Jmin,Jmax] =[",i2,",",i2,"]")')params%Jmin, params%Jmax
                write(*,'("Nblocks Available from Memory =",i6)') params%number_blocks
                write(*,'("Nblocks (if all trees dense)=",i6)') number_dense_blocks
                write(*,'("Nblocks used (sparse)=",i6)') sum(lgt_n(1:tree_n))
                write(*,'("Predictor=",A)') params%order_predictor
                write(*,'("block_distribution=",A)') params%block_distribution
                write(*,'(80("─"))')
            endif

            !----------------------------------
            ! reconstruct iteration
            !----------------------------------
            reconst_tree_ID = tree_n + 1
            if (save_all) then
                do iteration = 1, N_snapshots
                    call reconstruct_iteration( params, hvy_block, hvy_tmp, a_coefs, iteration , N_modes_used, reconst_tree_ID)
                    !----------------------------------
                    ! Save components of reconstructed
                    ! Snapshot
                    !----------------------------------
                    do j = 1, n_components
                        write( file_out, '("reconst",i1,"-",i3.3,"_", i12.12, ".h5")') j, N_modes_used, iteration
                        call saveHDF5_tree(file_out, real(iteration,kind=rk), iteration, j, params, hvy_block, reconst_tree_ID)
                    end do
                end do
            else
                !! reconstruct only single snapshot
                call reconstruct_iteration( params, hvy_block, hvy_tmp, a_coefs, iteration , N_modes_used, reconst_tree_ID)
                !----------------------------------
                ! Save components of reconstructed
                ! Snapshot
                !----------------------------------
                do j = 1, n_components
                    write( file_out, '("reconst",i1,"-",i3.3,"_", i12.12, ".h5")') j, N_modes_used, iteration

                    call saveHDF5_tree(file_out, real(iteration,kind=rk), iteration, j, params, hvy_block, reconst_tree_ID)
                end do
            endif



        end subroutine
        !##############################################################


        !##############################################################
        subroutine reconstruct_iteration( params, hvy_block, hvy_tmp, a_coefs, iteration, N_modes, dest_tree_ID, opt_offset_mode_tree_ID)
            implicit none

            !-----------------------------------------------------------------
            type (type_params), intent(inout)   :: params
            real(kind=rk),  intent(inout)       :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
            !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
            real(kind=rk),  intent(inout)       :: hvy_tmp(:, :, :, :, :)
            integer(kind=ik), intent(inout)     :: dest_tree_ID
            integer(kind=ik), intent(in)        :: iteration, N_modes         !> number of POD modes
            !> Threshold value for truncating POD modes. If the singular value is smaller,
            !> then the given threshold we discard the corresponding POD MODE.
            real(kind=rk),  intent(in)          :: a_coefs(:,:)
            !> optional argument for the tree_ID of the first mode:
            !> this can be useful if you have some trees which are stored at beginning
            !>  tree_1 tree_2 tree_3, tree_mode1, tree_mode2 etc.
            !> default is 0
            integer(kind=ik), optional, intent(in)       :: opt_offset_mode_tree_ID
            !---------------------------------------------------------------
            integer(kind=ik)                    ::  i, rank, free_tree_ID, offset = 0_ik
            real(kind=rk)                       :: a, t_elapse
            character(len=cshort)               :: filename

            ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
            ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
            ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
            ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
            ! to the last subroutine.)  -Thomas

            !---------------------------------------------------------------------------
            ! check inputs and set default values
            !---------------------------------------------------------------------------
            rank= params%rank

            if (present(opt_offset_mode_tree_ID)) then
                offset=opt_offset_mode_tree_ID
            endif

            if (rank == 0) then
                write(*, *)
                write(*,'(80("─"))')
                write(*,'(30("#"), " Reconstruct ", 30("#"))')
                write(*,'(80("─"))')
                write(*,'("Number of Modes used: ",i4)') N_modes
                if (params%adapt_tree) write(*,'("Compression threshold eps: ",es12.4)') params%eps
                write(*,'(80("─"))')
                write(*, *)
            endif

            ! in this algorithm we need at least N_mode trees plus 2 additional fields for
            ! saving the reconstructed field and a temporary field
            if ( params%forest_size < N_modes + 2 ) call abort(1003191,"Error! Need more Trees. Tip: increase forest_size")
            if ( dest_tree_ID<= N_modes) call abort(200419,"Error! Destination tree id overwrites POD Modes!")
            call copy_tree(params, hvy_block, dest_tree_ID, offset+1)
            call multiply_tree_with_scalar(params, hvy_block, dest_tree_ID, a_coefs(iteration,1))
            free_tree_ID = dest_tree_ID + 1
            !---------------------------------------------------------------------------
            ! reconstruct field
            !---------------------------------------------------------------------------
            t_elapse = MPI_wtime()
            do i = 2, N_modes

                a = a_coefs(iteration, i)
                ! calculate pod modes:
                call copy_tree(params, hvy_block, free_tree_ID , offset + i)
                call multiply_tree_with_scalar(params, hvy_block, free_tree_ID, a)
                call add_two_trees(params, hvy_block, hvy_tmp, dest_tree_ID, free_tree_ID)

            end do
            !---------------------------------------------------------------------------
            ! adapted reconstructed field
            !---------------------------------------------------------------------------
            call updateNeighbors_tree( params, dest_tree_ID, search_overlapping=.false.)

            if ( params%adapt_tree) then
                call adapt_tree( 0.0_rk, params, hvy_block, dest_tree_ID, params%coarsening_indicator, hvy_tmp )
            else
                call sync_ghosts_tree( params, hvy_block, dest_tree_ID)
            endif

            t_elapse = MPI_WTIME() - t_elapse
            if (rank == 0) then
                write(*,'("Snapshot ", i4," reconstructed in t_cpu=",es12.4, "sec [Jmin,Jmax]=[",i2,",",i2,"]")') &
                iteration, t_elapse, &
                minActiveLevel_tree( dest_tree_ID ), &
                maxActiveLevel_tree( dest_tree_ID )
            endif

        end subroutine reconstruct_iteration
        !##############################################################


        !##############################################################
        subroutine post_timecoef_POD(params)
            use module_globals
            use module_params
            use module_mpi

            implicit none

            !> parameter struct
            !--------------------------------------------
            type (type_params), intent(inout)  :: params
            !--------------------------------------------
            character(len=cshort) ::  filename, args, order
            character(len=cshort),dimension(:), allocatable   :: fsnapshot_list, fmode_list
            character(len=cshort),dimension(:,:), allocatable ::snapshot_in, mode_in
            real(kind=rk)   , allocatable :: L2error(:),time(:)
            integer(kind=ik), allocatable :: iter_list(:)
            real(kind=rk),    allocatable :: hvy_block(:, :, :, :, :)
            real(kind=rk),    allocatable :: hvy_tmp(:, :, :, :, :),a_coefs(:,:)
            integer(kind=ik)                 :: level, k, Bs(3), tc_length
            integer(hid_t)                   :: file_id
            real(kind=rk), dimension(3)      :: domain
            integer(hsize_t), dimension(2)   :: dims_treecode
            integer(kind=ik) :: treecode_size, number_dense_blocks, tree_ID, dF, N_modes
            integer(kind=ik) :: i,it, N_snapshots, dim, fsize, lgt_n_tmp, r, iteration=-1
            integer(kind=ik) :: j, n_components=1, io_error, free_tree_ID,pod_mode_tree_ID, unused_int
            real(kind=rk)    :: maxmem=-1.0_rk, eps=-1.0_rk, L2norm, L2norm_snapshots, Volume, norm
            real(kind=rk)    :: unused_var
            logical :: verbose = .false., save_all = .false.

            ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
            ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
            ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
            ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
            ! to the last subroutine.)  -Thomas


            call get_command_argument(2, args)
            if ( args== '--help' .or. args == '--h') then
                if ( params%rank==0 ) then
                    write(*,*)
                    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    write(*,*) " Compute Temporal coefficients of POD BASIS"
                    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    write(*,*) "Call mpi_command -n number_procs ./wabbit-post --POD-time --components=3"
                    write(*,*) "                                --snapshot-list='list_u1.txt [list_u2.txt] [list_u3.txt]'"
                    write(*,*) "                                --mode-list='list_m1.txt [list_m2.txt] [list_m3.txt]'"
                    write(*,*) "                                [--order=[2|4] --adapt=0.1 --memory=2GB --iteration=1]"
                    write(*,*) "                                [--save_all, --eps-norm=[L2,Linfty]]"
                    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    write(*,*) " returns:"
                    write(*,*) " a_coef.txt           list of temporal coefficients for all modes"
                    write(*,*) " "
                    write(*,*) " Input:"
                    write(*,*) " --snapshot-list      list of files containing all snapshots"
                    write(*,*) " --mode-list          list of files containing all modes for reconstruction"
                    write(*,*) " --components         number of components in statevector"
                    write(*,*) " --save_all           reconstruct all snapshots (default is false)"
                    write(*,*) " --iteration          reconstruct single snapshot (default is no snapshots)"
                    write(*,*) " --order              order of the predictor"
                    write(*,*) " --adapt              threshold for wavelet adaptation of modes and snapshot"
                    write(*,*) " --iteration          reconstruct single snapshot"
                    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                end if
                return
            end if

            !----------------------------------
            ! read predefined params
            !----------------------------------
            call get_cmd_arg_str( "--eps-norm", params%eps_norm, default="L2" )
            call get_cmd_arg_dbl( "--adapt", eps, default=-1.0_rk )
            call get_cmd_arg_str( "--order", order, default="CDF44" )
            call get_cmd_arg_str_vct( "--snapshot-list", fsnapshot_list )
            call get_cmd_arg_str_vct( "--mode-list", fmode_list )
            call get_cmd_arg( "--memory", args, default="2GB")
            read(args(1:len_trim(args)-2),* ) maxmem
            call get_cmd_arg( "--save_all", save_all, default=.true.)
            call get_cmd_arg( "--components", n_components, default=1_ik)
            call get_cmd_arg( "--iteration", iteration, default=1)

            if ( iteration>0 ) then
                if ( params%rank == 0 ) write(*,*) "Iteration reconstructed: " ,iteration
            endif

            if ( save_all ) then
                if ( params%rank == 0 ) write(*,*) "Reconstructing all snapshots!"
            endif
            !-------------------------------
            ! Set some wabbit specific params
            !-------------------------------
            ! block distirbution:
            params%block_distribution="sfc_hilbert"
            ! no time stepping:
            params%time_step_method="none"
            params%Jmin=1
            ! coarsening indicator
            params%physics_type="POD"
            params%eps_normalized=.True. ! normalize the statevector before thresholding
            params%coarsening_indicator="threshold-state-vector"
            params%threshold_mask=.False.
            params%n_eqn=n_components

            call set_wavelet_params(params, order )

            if ( eps > 0.0_rk) then
                ! adapt the mesh if possible
                params%adapt_tree = .True.! .False.!.True.
                params%eps=eps
            else
                params%adapt_tree = .False.
                params%eps = 0.0_rk
            endif


            !-------------------------------
            ! check if files exists:
            !-------------------------------
            ! fname_list is the name of the file which contains a list of files you want to
            ! read in (i.e. rho_000000000.h5 rho_000000001.h5 ... etc)
            if (.not. allocated(fsnapshot_list)) call abort(207191,"you must pass at least one snapshot list! Use: --snaphot-list my_filelist.txt")
            if (.not. allocated(fmode_list))     call abort(207191,"you must pass at least one mode list! Use: --mode-list my_filelist.txt")
            do j = 1, n_components
                call check_file_exists ( fsnapshot_list(j) )
                call check_file_exists ( fmode_list(j) )
                if (params%rank==0) write(*,*) "Snapshot list: "//trim(fsnapshot_list(j)) //"   Mode list: "// fmode_list(j)
            enddo
            call count_lines_in_ascii_file_mpi(fsnapshot_list(1), N_snapshots, n_header=0_ik)
            call count_lines_in_ascii_file_mpi(fmode_list(1), N_modes, n_header=0_ik)

            !-----------------------------------------------------------------------------
            ! ALLOCATE additional arrays
            !-----------------------------------------------------------------------------
            allocate(a_coefs(1:N_snapshots, 1:N_modes))
            allocate(L2error(1:N_modes))
            allocate(params%input_files(n_components))
            allocate(mode_in(1:N_modes,1:n_components))
            allocate(snapshot_in(N_snapshots,n_components))
            allocate(time(N_snapshots))
            allocate(iter_list(N_snapshots))
            !----------------------------------
            ! READ a_coefficient from file
            !----------------------------------
            call count_lines_in_ascii_file_mpi(fsnapshot_list(1), N_snapshots, n_header=0_ik)

            !-------------------------------------------
            !-------------------------------------------
            ! SNAPSHOTS
            !-------------------------------------------
            !-------------------------------------------
            ! check and find common params in all h5-files
            ! of all snapshots
            !-------------------------------------------
            ! open all files
            do j = 1, n_components
                open( unit=10+j, file=fsnapshot_list(j), action='read', status='old' )
            enddo

            io_error = 0
            i = 1
            do while(i <= N_snapshots)
                do j = 1, n_components
                    read (10+j, '(A)', iostat=io_error) snapshot_in(i,j)

                    call check_file_exists ( snapshot_in(i,j) )

                    if ( i == 1 .and. j == 1 ) then
                        ! read all geometric parameters of grid
                        call read_attributes(snapshot_in(i,j), lgt_n_tmp, time(1), iter_list(1), params%domain_size, &
                        params%Bs, params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
                    endif
                    call read_attributes(snapshot_in(i,j), lgt_n_tmp, time(i), iter_list(i), domain, bs, level, dim)
                    params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot
                    if (any(params%Bs .ne. Bs)) call abort( 204191, " Block size is not consistent ")
                    if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")
                    if (params%dim .ne. dim) call abort( 203193, "Dimension is not consistent ")
                end do
                i = i + 1
            end do

            !--------------------------------
            ! close the list of files
            !--------------------------------
            do j = 1, n_components
                close(10+j)
            enddo
            !-------------------------------------------
            !-------------------------------------------
            ! MODES
            !-------------------------------------------
            !-------------------------------------------
            !--------------------------------
            ! MODES: open all the given txt-files
            !--------------------------------
            do j = 1, n_components
                open( unit=15+j, file=fmode_list(j), action='read', status='old' )
            enddo
            !--------------------------------
            ! in every txt-file there is a list
            ! of h5-files which contain the
            ! statevector components of the
            ! POD modes
            !--------------------------------
            do i = 1, N_modes
                do j = 1, n_components
                    ! here we read in each line of file fname_list(j)
                    read (15+j, '(A)', iostat=io_error) mode_in(i,j)
                    call check_file_exists ( mode_in(i,j) )
                    !-------------------------------------------
                    ! check and find common params in all h5-files
                    !-------------------------------------------
                    call read_attributes(mode_in(i,j), lgt_n_tmp, unused_var, unused_int, domain, bs, level, dim)
                    params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot
                    if (any(params%Bs .ne. Bs)) then
                        write(*,*) "Bs: ",Bs, " vs ",params%Bs;
                        call abort( 203199, " Block size is not consistent ")
                    end if
                    if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")
                    if (params%dim .ne. dim) call abort( 204193, "Dimension is not consistent ")
                end do
            enddo
            !--------------------------------
            ! close the list of files
            !--------------------------------
            do j = 1, n_components
                close(15+j)
            enddo

            ! now we have all information to allocate the grid and set up the forest:
            fsize = N_snapshots + N_modes + 2 !we need some extra fields for storing etc
            params%forest_size = fsize
            number_dense_blocks = 2_ik**(dim*params%Jmax)*fsize
            params%n_eqn = n_components
            allocate(params%threshold_state_vector_component(params%n_eqn))
            params%threshold_state_vector_component(1:params%n_eqn)=1
            if (maxmem < 0.0_rk) then
                params%number_blocks = ceiling( 4.0_rk *fsize * number_dense_blocks / params%number_procs )
            endif

            call adjust_n_ghosts_for_scalarproduct(params)
            !----------------------------------
            ! allocate data
            !----------------------------------
            call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

            call reset_forest(params)


            hvy_neighbor = -1
            lgt_n = 0 ! reset number of acitve light blocks
            tree_n= 0 ! reset number of trees in forest
            !----------------------------------
            ! READ ALL SNAPSHOTS
            !----------------------------------
            do tree_ID = 1, N_snapshots
                call readHDF5vct_tree(snapshot_in(tree_ID,:), params, hvy_block, tree_ID)
            end do
            !----------------------------------
            ! READ ALL MODES
            !----------------------------------
            do j = 1, N_modes
                tree_ID = N_snapshots+j
                call readHDF5vct_tree(mode_in(tree_ID,:), params, hvy_block, tree_ID)
            end do


            ! --------------------
            ! Compute the L2 norm
            !---------------------
            L2norm_snapshots = compute_L2norm(params, hvy_block, hvy_tmp, N_snapshots, verbose)
            Volume = product(domain(1:params%dim))
            if (params%rank==0) then
                write(*,'(80("─"))')
                write(*,*) "WABBIT POD."
                write(*,*) "Snapshot matrix X build from:"
                do i = 1, N_snapshots
                    write(*, '("Snapshot=",i3 , " time=",f16.9,1x," Nblocks=", i6," sparsity=(",f5.1,"% / ",f5.1,"%) [Jmin,Jmax]=[",i2,",",i2,"]")')&
                    i, time(i), lgt_n(i), &
                    100.0*dble(lgt_n(i))/dble( (2**maxActiveLevel_tree( i ))**params%dim ), &
                    100.0*dble(lgt_n(i))/dble( (2**params%Jmax)**params%dim ), &
                    minActiveLevel_tree( i ), &
                    maxActiveLevel_tree( i )
                end do
                write(*,*) "Modes used for reconstruction:"
                do i = 1, N_modes
                    do j = 1, params%n_eqn
                        write(*, '("file=",A )') mode_in(i,j)
                    end do
                end do
                write(*,'(80("─"))')
                write(*,'("L2 norm ||X||_2^2: ",g18.8)') L2norm_snapshots**2
                write(*,'("Data dimension: ",i1,"D")') params%dim
                write(*,'("Domain size is ",3(g12.4,1x))') domain
                write(*,'("NCPU=",i6)') params%number_procs
                write(*,'("Number Trees=",1x,i4)') fsize
                write(*,'("[Jmin,Jmax] =[",i2,",",i2,"]")')params%Jmin, params%Jmax
                write(*,'("Nblocks Available from Memory =",i6)') params%number_blocks
                write(*,'("Nblocks (if all trees dense)=",i6)') number_dense_blocks
                write(*,'("Nblocks used (sparse)=",i6)') sum(lgt_n(1:tree_n))
                write(*,'("Predictor=",A)') params%order_predictor
                write(*,'("wavelet=",A)') params%wavelet
                write(*,'("Number ghosts=",i2)') params%g
                write(*,'("block_distribution=",A)') params%block_distribution
                write(*,'(80("─"))')
            endif

            !---------------------------------------------------------------------------
            ! temporal coefficients
            !---------------------------------------------------------------------------
            if (params%rank ==0) write(*,*) "Computing temporal coefficients a"

            do it = 1, N_snapshots
                do i = 1, N_modes
                    ! scalar product (inner product)
                    pod_mode_tree_ID = N_snapshots + i
                    a_coefs(it,i) = scalar_product_two_trees( params, hvy_block, hvy_tmp, it, pod_mode_tree_ID)
                enddo
            enddo

            Volume = product(params%domain_size(1:params%dim))
            ! V is the matrix of eigenvectors
            a_coefs = a_coefs / Volume

            if (params%rank==0) then
                write(*,*)
                filename = "temporal_coefs.txt"
                write(*,'( "Temporal coefficients saved to file: ", A30 )') filename
                open(14, file=filename, status='replace')
                do i = 1, N_snapshots
                    write(14,'(400(es15.8,1x))') a_coefs(i, 1:N_modes)
                enddo
                close(14)
            end if



        end subroutine
        !###############################################################


        subroutine set_wavelet_params(params, order )
            implicit none
            !> params structure of WABBIT
            type(type_params),intent(inout)  :: params
            !> order of the wavelt
            character(len=cshort), intent(in)    :: order

            ! initialize wavelet transform
            ! also, set number of ghost nodes params%G to minimal value for this wavelet
            params%wavelet = order
            call setup_wavelet(params, params%g)
        end subroutine set_wavelet_params

        subroutine adjust_n_ghosts_for_scalarproduct(params)
            !> This routine increases the number of ghost nodes if it is allowed by the
            !> block size. This is needed in order to make use of the scaling functions
            !> autocorrelation in the scalarproduct. The autocorrelation function has
            !> a support size of 10 points (5 in each direction)
            implicit none
            !> params structure of WABBIT
            type(type_params),intent(inout)  :: params
            integer(kind=ik) :: Bs(3)

            Bs = params%Bs

            if (params%order_predictor == "multiresolution_4th" .and. params%g<5) then
                if (params%dim == 3) then
                    ! note setting n_ghosts to 5 is only allowed if the block is big enough
                    if (5<(Bs(1)+1)/2 .and. 5<(Bs(2)+1)/2 .and. 5<(Bs(3)+1)/2) params%g = 6 ! however we adjust it to 6 because its even
                else
                    if (5<(Bs(1)+1)/2 .and. 5<(Bs(2)+1)/2) params%g = 6
                endif
            endif

        end subroutine

    end module module_MOR
