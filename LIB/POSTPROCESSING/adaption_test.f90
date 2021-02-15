subroutine adaption_test(params)
  use module_precision
  use module_params
  use module_forest
  use module_IO
  use module_mpi
  use module_globals
  use module_mesh

  implicit none

  !> parameter struct
  !--------------------------------------------
  type (type_params), intent(inout)  :: params
  !--------------------------------------------
  character(len=80):: file_out, order, args
  real(kind=rk), allocatable      :: error(:)
  integer(kind=ik), allocatable   :: lgt_block(:, :),Nb_adapt(:)
  character(len=80),allocatable   :: eps_str_list(:)
  real(kind=rk), allocatable      :: hvy_block(:, :, :, :, :)
  real(kind=rk), allocatable      :: hvy_tmp(:, :, :, :, :)
  integer(kind=ik), allocatable   :: hvy_neighbor(:,:), hvy_active(:, :)
  integer(kind=ik), allocatable   :: lgt_active(:,:), lgt_n(:), hvy_n(:)
  integer(kind=tsize), allocatable:: lgt_sortednumlist(:,:,:)
  integer(kind=ik)                :: max_neighbors, level, k, Bs(3), tc_length
  real(kind=rk), dimension(3)     :: domain
  integer(kind=ik) :: treecode_size, number_dense_blocks, tree_id_input, tree_id_adapt
  integer(kind=ik) :: i, dim, fsize, n_eps, rank, iteration
  integer(kind=ik) :: j, n_components=1, tree_n, lgt_n_tmp,Jmin, Jmax
  real(kind=rk) :: maxmem=-1.0_rk, eps=-1.0_rk, L2norm, Volume, t_elapse(2), time
  logical :: verbose = .false., save_all = .true.
  character(len=30) :: rowfmt

  call get_command_argument(2, args)
  if ( args== '--help' .or. args == '--h') then
      if ( params%rank==0 ) then
          write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          write(*,*) "mpi_command -n number_procs ./wabbit-post --adaption_test --eps-list='1e-5 1e-4' --list=filelist.txt [list_uy.txt] [list_uz.txt]"
          write(*,*) "[--order=CDF[22|44|40] --eps-norm=[L2|Linfty] --save_all]"
          write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          write(*,*) " Wavelet adaption test "
          write(*,*) " --list               list of files containing all snapshots"
          write(*,*) " --components         number of components in statevector"
          write(*,*) " --eps-list           eps used"
          write(*,*) " --order              order of the predictor"
          write(*,*) " --adapt              threshold for wavelet adaptation of modes and snapshot"
          write(*,*) " --eps-norm           normalization of wavelets"
          write(*,*) " --save_all           saves adapted snapshots"
          write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      end if
      return
  end if

  !----------------------------------
  ! read parameters
  !----------------------------------
  call get_cmd_arg( "--save_all", save_all, default=.false.)
  call get_cmd_arg_str( "--eps-norm", params%eps_norm, default="L2" )
  call get_cmd_arg_str( "--order", order, default="CDF44" )
  call get_cmd_arg( "--list", params%input_files )
  call get_cmd_arg( "--eps-list", eps_str_list )
  call get_cmd_arg( "--memory", args, default="2GB")
  read(args(1:len_trim(args)-2),* ) maxmem


  !-------------------------------
  ! Set some wabbit specific params
  !-------------------------------
  rank = params%rank
  n_components = size(params%input_files)
  n_eps = size(eps_str_list)
  params%block_distribution="sfc_hilbert"
  params%time_step_method="none"
  params%min_treelevel=1
  params%physics_type="POD"
  params%eps_normalized=.True. ! normalize the statevector before thresholding
  params%adapt_mesh = .True.! .False.!.True.
  params%coarsening_indicator="threshold-state-vector"
  params%threshold_mask=.False.
  params%n_eqn = n_components
  params%forest_size = 3

  ! Check parameters for correct inputs:
  if (order == "CDF20" .or. order == "CDF2,0") then
      params%harten_multiresolution = .true.
      params%order_predictor = "multiresolution_2nd"
      params%wavelet='CDF2,0'
      params%n_ghosts = 2_ik
  elseif (order == "CDF40" .or. order == "CDF4,0") then
      params%harten_multiresolution = .true.
      params%order_predictor = "multiresolution_4th"
      params%wavelet='CDF4,0'
      params%n_ghosts = 4_ik
  elseif (order == "CDF44" .or. order == "CDF4,4") then
      params%harten_multiresolution = .false.
      params%wavelet_transform_type = 'biorthogonal'
      params%order_predictor = "multiresolution_4th"
      params%wavelet='CDF4,4'
      params%n_ghosts = 6_ik
  else
      call abort(20030202, "The --order parameter is not correctly set [CDF40, CDF20, CDF44]")
  end if

  allocate(error(n_eps))
  allocate(Nb_adapt(n_eps))

  !-------------------------------
  ! check if files exists:
  !-------------------------------
  if (.not. allocated(params%input_files)) call abort(207191,"you must pass at least one file list! Use: --list my_filelist.txt")
  do j = 1, n_components
     call check_file_exists ( params%input_files(j) )
     if (params%rank==0) write(*,*) "Reading list of files from "//params%input_files(j)
  enddo
  !-----------------------------------------------------------------------------
  ! read in the file, loop over lines
  !-----------------------------------------------------------------------------
  do j = 1, n_components
      call read_attributes(params%input_files(j), lgt_n_tmp, time, iteration, params%domain_size, &
                       params%Bs, params%max_treelevel, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
  end do

  number_dense_blocks = 2_ik**(dim*params%max_treelevel)*fsize
  allocate(params%threshold_state_vector_component(params%n_eqn))
  params%threshold_state_vector_component(1:params%n_eqn)=.True.
  if (maxmem < 0.0_rk) then
    params%number_blocks = ceiling( 4.0_rk * params%forest_size * number_dense_blocks / params%number_procs )
  endif
  !----------------------------------
  ! allocate data
  !----------------------------------
  call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
  hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)

  call reset_forest(params, lgt_block, lgt_active, lgt_n,hvy_active, hvy_n, &
  lgt_sortednumlist,tree_n)


  hvy_neighbor = -1_ik
  tree_n= 0_ik ! reset number of trees in forest
  tree_id_input = 1
  tree_id_adapt = 2
  call read_field2tree(params, params%input_files , params%n_eqn, tree_id_input, &
                tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, hvy_block, &
                hvy_active, hvy_n, hvy_tmp, hvy_neighbor,  verbosity=.false.)

  !###########################################################
  ! actual error calculation:
  ! error_epsilon=|| u_input - u_epsilon||_2 / ||u_input||_2
  !###########################################################
  ! need the L2 norm of the input for the relative error
  L2norm = compute_tree_L2norm( params, tree_n, lgt_block,  lgt_active, &
               lgt_n, lgt_sortednumlist, hvy_block, hvy_neighbor,&
               hvy_active, hvy_n, hvy_tmp, tree_id_input, verbose )

  do i = 1, n_eps
    ! copy form original data
    call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_adapt, tree_id_input)
    ! adapt to given eps:
    read (unit=eps_str_list(i),fmt=*) params%eps
    call adapt_mesh( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_id_adapt), &
        lgt_n(tree_id_adapt), lgt_sortednumlist(:,:,tree_id_adapt), hvy_active(:,tree_id_adapt), &
        hvy_n(tree_id_adapt), tree_id_adapt, params%coarsening_indicator, hvy_tmp )

    lgt_n_tmp = lgt_n(tree_id_adapt)
    Jmin = min_active_level( lgt_block, lgt_active(:,tree_id_adapt), lgt_n(tree_id_adapt) )
    Jmax = max_active_level( lgt_block, lgt_active(:,tree_id_adapt), lgt_n(tree_id_adapt) )

    if (save_all) then
      do j = 1, n_components
          write( file_out, '("u",i1,"-eps", A,"_",i12.12 ,".h5")') j, trim(adjustl(eps_str_list(i))), nint(time * 1.0e6_rk)

          call write_tree_field(file_out, params, lgt_block, lgt_active, hvy_block, &
          lgt_n, hvy_n, hvy_active, j, tree_id_adapt , time , iteration )
      end do
    end if

    ! compare to original data:
    call substract_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id_adapt, tree_id_input)
    ! compute L2 norm
    error(i) = compute_tree_L2norm( params, tree_n, lgt_block,  lgt_active, &
                                lgt_n, lgt_sortednumlist, hvy_block, hvy_neighbor,&
                                 hvy_active, hvy_n, hvy_tmp, tree_id_adapt, verbose )
    error(i) = error(i)/L2norm
    Nb_adapt(i) = lgt_n_tmp

    ! give some output
    if (rank == 0) then
       write(*,'("Field adapted to eps=",es10.2," rel err=", es10.2," Nblocks=", i6," Nb_adapt/Nb_dense=",f6.1,"% Nb_adapt/Nb_input=",f5.1,"% [Jmin,Jmax]=[",i2,",",i2,"]")') &
              params%eps, error(i), lgt_n_tmp, &
              100.0*dble(lgt_n_tmp)/dble( (2**params%max_treelevel)**params%dim ), &
              100.0*dble(lgt_n_tmp)/dble(lgt_n(tree_id_input)), Jmin, Jmax
    endif
  end do
  ! elapsed time for reading and coarsening the data
  t_elapse(1)= MPI_Wtime()-t_elapse(1)

  !  Save values
  if (save_all .and. params%rank==0) then
    file_out ="compression_error.txt"
    write(*,'( "compresion error saved to: ", A30 )') file_out
    write(*,*)
    open(14,file=file_out, status='replace')
    write(14,'("eps", 1x, " L2error", 1x ,"Nb blocks", 1x , "Nb_adapt/Nb_input")')
    do i = 1, n_eps
      write(14,FMT='(A," ",es15.8," ",i7," ",es10.3)') &
      trim(adjustl(eps_str_list(i))),error(i),Nb_adapt(i), 1.0_rk*dble(Nb_adapt(i))/dble(lgt_n(tree_id_input))
    enddo
    close(14)
  end if


end subroutine
!##############################################################
