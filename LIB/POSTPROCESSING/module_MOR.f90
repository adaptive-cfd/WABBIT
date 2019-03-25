!-----------------------------------------------------------------
!> \file
!> \brief
!! Module for MODEL ORDER REDUCTION
!! Current Methods:
!!              * snapshot POD
!> \version 01.3.2019
!> \author P.Krah
!-----------------------------------------------------------------

module module_MOR

  ! modules
  use mpi
  use module_forest
  use module_precision
  use module_params
  use module_IO

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: snapshot_POD, post_POD
  !**********************************************************************************************
contains


  !##############################################################
  !> Implementation of the multiresolution snapshot POD
  !> Input: 
  !>  1.) All Snapshots needed for the POD in the first N_snapshots trees.
  !>  2.) The maximal forest size needs to be at least 2*N_snapshots.
  !>  3.) Truncation Rank (optional) is the number of POD MODES computed in the algorithm (maxium is N_snapshots)
  !>      (default is N_snapshots)
  !>  4.) Truncation Error (default: 1e-9) is the smallest squared singularvalue used
  !>      for POD modes
  !> The output is as the following:
  !>  1.) All input snapshots are returned with tree: id 1, 2, ..., N_snapshots
  !>  2.) The POD modes are saved in trees with tree_ids:
  !>            N_snapshots+1, N_snapshots+2, ..., N_snapshots+truncation_rank
  !>  3.) On output the truncation rank is updated to the acutal number of modes with
  !>      squared singular value bigger then the given or default truncation error
  subroutine snapshot_POD( params, lgt_block,  lgt_active, lgt_n, lgt_sortednumlist, &
                       hvy_block, hvy_neighbor, hvy_active, hvy_tmp, hvy_n, tree_n, &
                       truncation_error, truncation_rank)
    implicit none

    !-----------------------------------------------------------------
    !> user defined parameter structure
    type (type_params), intent(inout)     :: params
    !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    !> size of active lists
    integer(kind=ik),  intent(inout)        :: lgt_n(:),tree_n, hvy_n
    !> heavy data array - block data
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk),  intent(inout)           :: hvy_tmp(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), intent(inout)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)          :: lgt_active(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)       :: lgt_sortednumlist(:,:)
    !> number of POD modes
    integer(kind=ik),optional, intent(inout)     :: truncation_rank
    !> Threshold value for truncating POD modes. If the singular value is smaller,
    !> then the given treshold we discard the corresponding POD MODE.
    real(kind=rk), optional, intent(in)     :: truncation_error 
    !---------------------------------------------------------------
    real(kind=rk) :: C(tree_n,tree_n), V(tree_n,tree_n), work(5*tree_n), &
                    eigenvalues(tree_n), alpha(tree_n), max_err  
    integer(kind=ik):: N_snapshots, root, ierr, i, rank, pod_mode_tree_id, &
                      free_tree_id,tree_id, N_modes, max_nr_pod_modes

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
      max_err=1e-9_rk
    endif

    if (rank == 0) then
      write(*, *)
      write(*,'(80("-"))')
      write(*,'(30("#"), " SNAPSHOT POD ", 30("#"))')
      write(*,'(80("-"))')
      write(*,'("Number of SNAPSHOTS used: ",i4)') N_snapshots
      write(*,'("Desired Truncation Rank: ", i4)') max_nr_pod_modes
      write(*,'("Maximal Error in L2 norm: ",es12.4)') max_err
      write(*,'(80("-"))')
      write(*, *)
    endif

    if ( params%forest_size <= 2*N_snapshots) call abort(1003191,"Error! Need more Trees. Tip: increase forest_size")

    !---------------------------------------------------------------------------
    ! Covariance Matrix of snapshot matrix X: C = X^T*X 
    !---------------------------------------------------------------------------
    if (rank == 0) write(*,*)
    if (rank == 0) write(*,*) "Processing Covariance Matrix C = X^T*X"
    if (rank == 0) write(*,*) 
    call compute_covariance_matrix( params, C, tree_n, &
                       lgt_block,  lgt_active, lgt_n, lgt_sortednumlist, &
                       hvy_block, hvy_neighbor, hvy_active, hvy_n, hvy_tmp)

    !---------------------------------------------------------------------------
    ! eigenvalues of covariance matrix
    !---------------------------------------------------------------------------
    call DSYEV('V', 'U', N_snapshots, C, N_snapshots, eigenvalues, work, 5*N_snapshots, ierr)
    ! as in matlab the eigenvalues are sorted in ascending order...
    ! on output, C now contains the eigenvectors:
    V = C

    if (ierr /= 0) call abort(333,"The eigenvalue solver failed...")
  
    if (rank == 0) then
      write(*,*) "----v eigenvalues v-----"
      do i = 1, N_snapshots
        write(*,'(1(es12.4,1x))') eigenvalues(i)
      enddo
      write(*,*) "----^ eigenvalues ^-----"
    endif

    
  !---------------------------------------------------------------------------
  ! construct POD basis functions (modes)
  !---------------------------------------------------------------------------
  N_modes = 0
  pod_mode_tree_id= N_snapshots
  free_tree_id = N_snapshots + 1
  if ( rank == 0 ) write(*,*) "Constructing POD modes (X*V)"
  do i = N_snapshots, 1, -1
  ! compute normalized eigenvectors. If eigenvalues are to small discard modes
  ! We loop in an ascending order since the eigenvalues are sorted from largest to 
  ! smallest.
  ! If the max nr of modes (=desired POD rank) or the desired precission (smallest eigenvalue=sigma**2)
  ! is reached, we exit the for loop, since we have build enough POD modes.
  
    if ( eigenvalues(i) < max_err .or. N_modes > max_nr_pod_modes-1 ) exit
    N_modes = N_modes +1
    if (rank == 0) write(*,*) "constructing POD mode", N_modes
    alpha = V(:,i)/sqrt(N_snapshots*eigenvalues(i))
    ! calculate pod modes:
    pod_mode_tree_id = pod_mode_tree_id + 1
    call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_neighbor, pod_mode_tree_id, 1) 
    call scalar_multiplication_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_neighbor, pod_mode_tree_id, alpha(1))

    free_tree_id = free_tree_id + 1

    do tree_id = 2, N_snapshots
      call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_neighbor, free_tree_id, tree_id) 

      call scalar_multiplication_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_neighbor, free_tree_id, alpha(tree_id))
           
      call add_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, pod_mode_tree_id, free_tree_id)    
    end do
  end do
  ! the truncation_rank can be lower then the default (N_snapshots) or input value,
  ! when the singular values are smaller as the desired presicion! Therefore we update
  ! the truncation rank here.
  truncation_rank = N_modes

  if (rank == 0) then
      write(*, *)
      write(*,'(80("-"))')
      write(*,'("Estimated Truncation Rank r= ",i4)') truncation_rank 
      if (truncation_rank < N_snapshots) then
        write(*,'("Error in L2 norm (sigma_{r+1}): ",es12.4)') sqrt(eigenvalues(truncation_rank+1))
      else
        write(*,'("Error in L2 norm (sigma_{r+1}): 0")') 
      endif
      write(*,'(80("-"))')
      write(*, *)
  endif

  end subroutine snapshot_POD
  !##############################################################
  

    subroutine  print_mat(Mat)

        real(kind=rk), intent(in)       :: Mat(:,:)
        integer                         :: i, j
        character(len=16) :: fmt
        character(len=3) :: ncols_str

        write(*,*) " "
        write(ncols_str,'(i3.3)') size(Mat,2)
        fmt = '('//ncols_str//'(es12.4,1x))'
        do j = 1,size(Mat,1)
          write(*,fmt) Mat(j,:)
        end do

    end subroutine print_mat


  !##############################################################
  !> Main routine of the POD postprocessing procedure.
  !> If post_POD is called from wabbit-post it will read in files
  !> specified on input and use it as snapshot data. After 
  !> decomposition the Modes will be safed to a file!
  subroutine post_POD(params)
    use module_precision
    use module_params
    use module_forest
    use module_mpi

    implicit none

    !> parameter struct
    !--------------------------------------------
    type (type_params), intent(inout)  :: params
    !--------------------------------------------
    character(len=80)      :: file_in, file_out, args
    real(kind=rk)   , allocatable :: time(:)
    integer(kind=ik), allocatable :: iteration(:)
    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:), hvy_active(:)
    integer(kind=ik), allocatable           :: lgt_active(:,:), lgt_n(:), tree_n
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(kind=ik)                        :: hvy_n,  max_neighbors, level, k, Bs(3), tc_length
    integer(hid_t)                          :: file_id
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik) :: treecode_size, number_dense_blocks, tree_id, truncation_rank_in = -1
    integer(kind=ik) :: i, n_opt_args, N_snapshots, dim, fsize, lgt_n_tmp, truncation_rank = 3
    real(kind=rk) :: truncation_error=1e-13_rk, truncation_error_in=-1.0_rk
    character(len=80) :: tmp_name
    character(len=2)  :: order

    call get_command_argument(2, file_in)
    if (file_in == '--help' .or. file_in == '--h') then
        if ( params%rank==0 ) then
            write(*,*) "postprocessing subroutine to refine/coarse mesh to a uniform grid (up and downsampling ensured). command line:"
            write(*,*) "mpi_command -n number_procs ./wabbit-post --POD [--order=[2|4] --nmodes=3 --error=1e-9] sources_*.h5 " 
        end if
        return
    end if

    !----------------------------------
    ! read predefined params 
    !----------------------------------
    n_opt_args = 1 ! counting all extra arguments, which are not *h5 files

    do i = 1, command_argument_count()
      !-------------------------------
      ! order of predictor
      if ( index(args,"--order=")==1 ) then
        read(args(9:len_trim(args)),* ) order 
        n_opt_args = n_opt_args + 1
      end if
      !-------------------------------
      ! TRUNCATION RANK
      call get_command_argument(i,args)
      if ( index(args,"--nmodes=")==1 ) then
        read(args(10:len_trim(args)),* ) truncation_rank_in
        n_opt_args = n_opt_args + 1
      end if
      !-------------------------------
      ! TRUNCATION ERROR
      call get_command_argument(i,args)
      if ( index(args,"--error=")==1 ) then
        read(args(9:len_trim(args)),* ) truncation_error_in
        n_opt_args = n_opt_args + 1
      end if
    end do

    if (order == "4") then
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik
    else 
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik
    end if

    if (truncation_rank_in /=-1) truncation_rank = truncation_rank_in
    if (truncation_error_in >1e-16_rk ) truncation_error = truncation_error_in

    !----------------------------------
    ! set addtitional params
    !----------------------------------
    ! block distirbution:
    params%block_distribution="sfc_hilbert"
    ! no time stepping:
    params%time_step_method="no" 
    
    ! read ini-file and save parameters in struct
    N_snapshots = command_argument_count()
    N_snapshots = N_snapshots - n_opt_args ! because of the first nargs arguments which are not files
    allocate(params%input_files(N_snapshots))
    allocate(time(N_snapshots)) 
    allocate(iteration(N_snapshots)) 
    !----------------------------------
    ! check and find common params
    !----------------------------------
    call get_command_argument(n_opt_args+1, params%input_files(1))
    call read_attributes(params%input_files(1), lgt_n_tmp, time(1), iteration(1), params%domain_size, &
                         params%Bs,params%max_treelevel, params%dim)
    do i = 2, N_snapshots
      call get_command_argument(n_opt_args+i, file_in)
      params%input_files(i) = file_in
      call read_attributes(file_in, lgt_n_tmp, time(i), iteration(i), domain, bs, level, dim)
      params%max_treelevel = max(params%max_treelevel, level) ! find the maximal level of all snapshot
      if (any(params%Bs .ne. Bs)) call abort( 203191, " Block size is not consistent ")
      if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")
      if (params%dim .ne. dim) call abort( 203193, "Dimension is not consistent ")
    end do
    fsize = 2*N_snapshots + 1 !we need some extra fields for storing etc
    params%forest_size = fsize
    number_dense_blocks = 2_ik**(dim*level)
    params%number_blocks = ceiling( 4.0_rk * N_snapshots * number_dense_blocks / params%number_procs )
    params%n_eqn = 1
    !----------------------------------
    ! allocate data
    !----------------------------------
    call allocate_hvy_lgt_data(params, lgt_block, hvy_block, hvy_neighbor, &
              lgt_active, lgt_n, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)
    call reset_lgt_data(lgt_block, lgt_active(:, fsize+1), &
              params%max_treelevel, lgt_n(fsize+1), lgt_sortednumlist)
    hvy_neighbor = -1
    lgt_n = 0 ! reset number of acitve light blocks
    tree_n= 0 ! reset number of trees in forest
    !----------------------------------
    ! READ ALL SNAPSHOTS
    !----------------------------------
    do i = 1, N_snapshots
      call read_field2tree(params, params%input_files(i), 1, i, &
                  tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, hvy_block, &
                  hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    end do

    if (params%rank==0) then
        write(*,'(80("-"))')
        write(*,*) "WABBIT POD." 
        write(*,*) "Snapshot matrix build from:"
        do i = 1, N_snapshots
          write(*, '("Snapshot=",i3 ," File=", A, " it=",i7,1x," time=",f16.9,1x, " Nblocks=", i6)') &
          i, trim(params%input_files(i)), iteration(i), time(i), lgt_n(i)
        end do
        write(*,'(A20,1x,A80)') "Predictor used:", params%order_predictor
        write(*,'(A20,1x,i3," using ",i9," Blocks")') "Jmax: ", level, number_dense_blocks
        write(*,'(80("-"))')
        write(*,'("Data dimension: ",i1,"D")') params%dim
        write(*,'("Domain size is ",3(g12.4,1x))') domain
        write(*,'("NCPU=",i6)') params%number_procs
        write(*,'("Memory Nb=",i6)') params%number_blocks
        write(*,'("Dense  Nb=",i6)') number_dense_blocks
        write(*,'(80("-"))')
        write(*,'("Predictor=",A)') params%order_predictor
        write(*,'("block_distribution=",A)') params%block_distribution
        
    endif

    !----------------------------------
    ! COMPUTE POD Modes
    !----------------------------------
    call snapshot_POD( params, lgt_block,  lgt_active, lgt_n, lgt_sortednumlist, &
                       hvy_block, hvy_neighbor, hvy_active, hvy_tmp, hvy_n, tree_n, &
                       truncation_error, truncation_rank)
    !----------------------------------
    ! Save Modes
    !----------------------------------
    do tree_id = N_snapshots+1, N_snapshots + truncation_rank
      i = tree_id - N_snapshots
      tmp_name = "PODmode"
      write( file_out, '(a, "_", i12.12, ".h5")') trim(adjustl(tmp_name)), i 

      call write_tree_field(file_out, params, lgt_block, lgt_active, hvy_block, &
          lgt_n, hvy_n, hvy_active, params%n_eqn, tree_id , 0.0_rk , tree_id - N_snapshots ) 

    end do

  end subroutine 
  !##############################################################




  !##############################################################
  !> construction of covariance matrix from snapshot data
  subroutine compute_covariance_matrix( params, C, tree_n, &
                       lgt_block,  lgt_active, lgt_n, lgt_sortednumlist, &
                       hvy_block, hvy_neighbor, hvy_active, hvy_n, hvy_tmp)
    implicit none

    !-----------------------------------------------------------------
    !> user defined parameter structure
    type (type_params), intent(in)  :: params
    !> Covariance matrix for snapshot POD
    real(kind=rk),  intent(inout)   :: C(:,:) 
    !> light data array
    integer(kind=ik),  intent(inout):: lgt_block(:, :)
    !> size of active lists
    integer(kind=ik),  intent(inout):: lgt_n(:), tree_n, hvy_n
    !> heavy data array - block data
    real(kind=rk),  intent(inout)   :: hvy_block(:, :, :, :, :)
    !> heavy temp data: needed in blockxfer which is called in add_tree 
    real(kind=rk),   intent(inout)  :: hvy_tmp(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), intent(inout) :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout) :: lgt_active(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout) :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)       :: lgt_sortednumlist(:,:)
    !---------------------------------------------------------------
    integer(kind=ik) :: tree_id1, tree_id2, free_tree_id, Jmax, Bs(3), g, &
                        N_snapshots, N, k, lgt_id, hvy_id, rank, i, mpierr
    real(kind=rk) :: C_val, Volume
    real(kind=rk) :: x0(3), dx(3)
    
    N_snapshots = size(C,1)
    N = params%number_blocks
    rank = params%rank
    Jmax = params%max_treelevel
    g = params%n_ghosts
    Bs= params%Bs
    Volume = product(params%domain_size(1:params%dim))
    if (params%forest_size < N_snapshots + 1 ) call abort(030319,"Forest is to small")
    free_tree_id = N_snapshots + 1
    ! We loop over all snapshots X_i, i=1,...,N to calculate the values of the symmetric
    ! covariance matrix C_{j,i} = C_{i,j} = <X_i, X_j>
    do tree_id1 = 1, N_snapshots
      do tree_id2 = tree_id1, N_snapshots
        !---------------------------
        ! copy tree_id1 -> free_tree_id
        !----------------------------
        ! copy tree with tree_id1 to tree with free_tree_id
        ! note: this routine deletes lgt_data of the "free_tree_id" before copying 
        !       the tree
        call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_neighbor, free_tree_id, tree_id1) 

        !---------------------------------------------------
        ! multiply tree_id2 * free_tree_id -> free_tree_id
        !---------------------------------------------------
        call multiply_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, free_tree_id, tree_id2)

        !----------------------------------------------------
        ! sum over all elements of the tree with free_tree_id
        !-----------------------------------------------------
        C_val = 0.0_rk
        do k = 1, hvy_n
          hvy_id = hvy_active(k)
          call hvy_id_to_lgt_id(lgt_id, hvy_id, rank, N )
          ! first we have to find out if the hvy data 
          ! belongs to the tree we want to sum over
          if ( lgt_block(lgt_id, Jmax + idx_tree_id) .ne. free_tree_id) cycle
          ! calculate the lattice spacing.
          ! It is needed to perform the L2 inner product
          call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
          if ( params%dim == 3 ) then
              C_val = C_val + dx(1)*dx(2)*dx(3)* sum( hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, :, hvy_id))
          else
              C_val = C_val + dx(1)*dx(2)*sum( hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, :, hvy_id))
          endif
        end do 

        ! Construct correlation matrix
        C(tree_id1, tree_id2) = C_val 
        C(tree_id2, tree_id1) = C_val
      end do
    end do
    ! sum over all Procs
    call MPI_ALLREDUCE(MPI_IN_PLACE, C, N_snapshots**2, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,WABBIT_COMM, mpierr)
    ! Normalice C to the Volume of the integrated area (L2 scalar product)
    ! and by the number of snapshots
    C = C / Volume
    C = C / dble(N_snapshots)
  end subroutine 
  !##############################################################
  





end module module_MOR

