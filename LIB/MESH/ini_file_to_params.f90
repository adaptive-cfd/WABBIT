! initialize params struct from INI file
subroutine ini_file_to_params( params, filename )
   use module_t_files, only : flush_frequency
   implicit none

   !> user defined parameter structure
   type (type_params), intent(inout)               :: params
   character(len=*), intent(in)                    :: filename
   ! process rank
   integer(kind=ik)                                :: rank
   ! number of processes
   integer(kind=ik)                                :: number_procs
   ! inifile structure
   type(inifile)                                   :: FILE
   ! maximum memory available on all cpus
   real(kind=rk)                                   :: maxmem, mem_per_block, nstages
   ! string read from command line call
   character(len=cshort)                           :: memstring
   integer(kind=ik)                                :: d,i, Nblocks_Jmax, g, Neqn, Nrk, g_RHS_min, diff_L, diff_R, Bs(1:3)

   rank         = params%rank
   number_procs = params%number_procs


   ! read the file, only process 0 should create output on screen
   call set_lattice_spacing_mpi(1.0_rk)
   call read_ini_file_mpi(FILE, filename, .true.)

   ! which physics module is used? (note that the initialization of different parameters takes
   ! place in those modules, i.e., they are not read here.)
   call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )

   ! what wavelet to use?
   ! (check that here as default for number ghost nodes g depends on it)
   call read_param_mpi(FILE, 'Wavelet', 'wavelet', params%wavelet, 'CDF40')

   call ini_domain(params, FILE)
   call ini_blocks(params, FILE)
   call ini_time(params, FILE)

   allocate(params%symmetry_vector_component(1:params%n_eqn))
   params%symmetry_vector_component = "0"
   call read_param_mpi(FILE, 'Domain', 'symmetry_vector_component', params%symmetry_vector_component, params%symmetry_vector_component )

   !**************************************************************************
   ! read INITIAL CONDITION parameters

   ! which physics module is used? (note that the initialization of different parameters takes
   ! place in those modules, i.e., they are not read here.)
   call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )

   ! if the initial condition is read from file, it is handled by wabbit itself, i.e. not
   ! by the physics modules. the pyhsics modules cannot do this, because they just see 'blocks'
   ! and never the entire grid as such.
   call read_param_mpi(FILE, 'Physics', 'read_from_files', params%read_from_files, .false. )

   ! sometimes you want to save each iteration, therefore it is better to use the
   ! iteration number as an identifier of the file
   call read_param_mpi(FILE, 'Physics', 'use_iteration_as_fileid', params%use_iteration_as_fileid, .false. )


   if (params%read_from_files ) then
      ! read variable names
      allocate( params%input_files( params%n_eqn ) )

      params%input_files = "---"
      call read_param_mpi(FILE, 'Physics', 'input_files', params%input_files, params%input_files, check_file_exists=.true. )
   end if

   ! wabbit does need to know how many fiels are written to disk when saving is triggered.
   ! e.g. saving ux, uy and p would mean 3. The names of these files as well as their contents
   ! are defined by the physics modules.
   call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params%N_fields_saved, 3 )

   !***************************************************************************
   ! read DISCRETIZATION parameters
   !
   ! discretization order
   call read_param_mpi(FILE, 'Discretization', 'order_discretization', params%order_discretization, "---" )
   ! poisson order
   call read_param_mpi(FILE, 'Discretization', 'poisson_order', params%poisson_order, "FD_4th_comp_1_3")
   call read_param_mpi(FILE, 'Discretization', 'poisson_cycle_end_criteria', params%poisson_cycle_end_criteria, "fixed_iterations")
   call read_param_mpi(FILE, 'Discretization', 'poisson_cycle_it', params%poisson_cycle_it, 6)
   call read_param_mpi(FILE, 'Discretization', 'poisson_cycle_tol', params%poisson_cycle_tol, 1.0e-6_rk)
   call read_param_mpi(FILE, 'Discretization', 'poisson_cycle_max_it', params%poisson_cycle_max_it, 100)
   call read_param_mpi(FILE, 'Discretization', 'poisson_GS_it', params%poisson_GS_it, 8)
   call read_param_mpi(FILE, 'Discretization', 'poisson_Sync_it', params%poisson_Sync_it, 2)
   call read_param_mpi(FILE, 'Discretization', 'poisson_coarsest', params%poisson_coarsest, "FFT")
   call read_param_mpi(FILE, 'Discretization', 'nprojection_NSI', params%nprojection_NSI, 1)
   call read_param_mpi(FILE, 'Discretization', 'FFT_accuracy', params%FFT_accuracy, "spectral")

   ! filter frequency
   call read_param_mpi(FILE, 'Discretization', 'filter_type', params%filter_type, "no_filter" )
   call read_param_mpi(FILE, 'Discretization', 'filter_only_maxlevel', params%filter_only_maxlevel, .false. )
   call read_param_mpi(FILE, 'Discretization', 'filter_all_except_maxlevel', params%filter_all_except_maxlevel, .false. )

   if (params%filter_type /= "no_filter") then
      call read_param_mpi(FILE, 'Discretization', 'filter_freq', params%filter_freq, -1 )
   endif

   !***************************************************************************
   ! read statistics parameters
   call read_param_mpi(FILE, 'Statistics', 'nsave_stats', params%nsave_stats, 99999999_ik )
   call read_param_mpi(FILE, 'Statistics', 'tsave_stats', params%tsave_stats, 9999999.9_rk )

   !***************************************************************************
   ! read time statistics parameters
   call read_param_mpi(FILE, 'Time-Statistics', 'time_statistics', params%time_statistics, .false.)
   if (params%time_statistics) then
      call read_param_mpi(FILE, 'Time-Statistics', 'N_time_statistics', params%N_time_statistics, 1)
      allocate( params%time_statistics_names(1:params%N_time_statistics) )
      call read_param_mpi(FILE, 'Time-Statistics', 'time_statistics_names', params%time_statistics_names, (/ "none" /))
   endif

   !***************************************************************************
   ! WABBIT needs to know about the mask function (if penalization is used): does it contain
   ! a time-dependent-part (e.g. moving obstacles, time-dependent forcing)? does it contain
   ! a time-independent part (fixed walls, homogeneous forcing)? or both? WABBIT needs to know
   ! that since we try to create the time-independent mask function only once, but the time-dependent
   ! part of course in every time step.
   call read_param_mpi(FILE, 'VPM', 'penalization', params%penalization, .false.)
   call read_param_mpi(FILE, 'VPM', 'mask_time_dependent_part', params%mask_time_dependent_part, .true.)
   call read_param_mpi(FILE, 'VPM', 'mask_time_independent_part', params%mask_time_independent_part, .true.)
   call read_param_mpi(FILE, 'VPM', 'dont_use_pruned_tree_mask', params%dont_use_pruned_tree_mask, .false.)



   !***************************************************************************
   ! read DEBUG parameters
   !
   ! unit test treecode flag
   call read_param_mpi(FILE, 'Debug', 'test_treecode', params%test_treecode, .false.)
   call read_param_mpi(FILE, 'Debug', 'test_ghost_nodes_synch', params%test_ghost_nodes_synch, .true.)
   call read_param_mpi(FILE, 'Debug', 'test_wavelet_decomposition', params%test_wavelet_decomposition, .true.)

   ! Hack.
   ! Small ascii files are written with the module_t_files, which is just a buffered wrapper.
   ! Instead of directly dumping the files to disk, it collects data and flushes after "flush_frequency"
   ! samples. In 2D, the code generally runs fast and does many time steps, hence
   ! we flush more rarely. In 3D, we can flush more often, because time steps take longer
   call read_param_mpi(FILE, 'Debug', 'flush_frequency', flush_frequency, -1)
   if (flush_frequency <= 0) then
      if (params%dim == 2) then
         flush_frequency = 50
      else
         flush_frequency = 10
      endif
   endif
   if (params%rank==0) write(*, '(A, i0, A)') "INIT: Flushing t-files every ", flush_frequency, " appends."

   !***************************************************************************
   ! read MPI parameters
   !
   ! data exchange method
   call ini_MPI(params, FILE )

   ! clean up
   if (params%rank==0) write(*,'("INIT: cleaning ini file")')
   call clean_ini_file_mpi(FILE)

   ! initialize wavelet (needs to be after reading parameters because it also checks the discretization)
    call setup_wavelet(params)


   ! check ghost nodes number
   if (params%rank==0) write(*,'("INIT: checking if g and predictor work together")')
   if ( (params%g < 6 .and. params%order_predictor == 'multiresolution_12th') .or. &
        (params%g < 5 .and. params%order_predictor == 'multiresolution_10th') .or. &
        (params%g < 4 .and. params%order_predictor == 'multiresolution_8th') .or. &
        (params%g < 3 .and. params%order_predictor == 'multiresolution_6th') .or. &
        (params%g < 2 .and. params%order_predictor == 'multiresolution_4th') .or. &
        (params%g < 1 .and. params%order_predictor == 'multiresolution_2nd') ) then
      call abort("ERROR: need more ghost nodes for order of supplied refinement interpolatior")
   end if
   if ( (params%g < 3 .and. (params%order_discretization == 'FD_4th_central_optimized' .or. params%order_discretization == 'FD_6th_central')) .or. &
        (params%g < 2 .and. params%order_discretization == 'FD_4th_central') .or. &
        (params%g < 1 .and. params%order_discretization == 'FD_2th_central') ) then
      call abort("ERROR: need more ghost nodes for order of supplied finite distance scheme")
   end if

   ! alter g_RHS if necessary, CDF4Y wavelets need only g_RHS=2 for example
   ! g_RHS is also dependent on the wavelet due to how we synch each stage:
   ! the third stage (prediction) needs points from the boundary to correctly interpolate values
   if ( params%order_discretization == 'FD_4th_central_optimized' .or. params%order_discretization == 'FD_6th_central' ) g_RHS_min = 3
   if ( params%order_discretization == 'FD_4th_central' ) g_RHS_min = 2
   if ( params%order_discretization == 'FD_2th_central' ) g_RHS_min = 1
   if (params%g_RHS < g_RHS_min) then
      if (params%rank==0) then
         write(*,  '(A, i0, A, i0, A)') "Warning!! 'number_ghost_nodes_rhs' was set smaller as required for FD scheme, adapting it from ", params%g_RHS, " to ", g_RHS_min, " (ignore this if it was not explicitly set)"
      endif
      params%g_RHS = g_RHS_min
   endif

   ! JB ToDo - correct once this is more settled
   ! ! we want to know the stencil size for laplacian schemes usually, so lets save this here
   ! if (params%poisson_order == 'CFD_2nd') then
   !    params%poisson_stencil_size = 1
   ! elseif (params%poisson_order == 'CFD_4th') then
   !    params%poisson_stencil_size = 2
   ! elseif (params%poisson_order == 'CFD_6th') then
   !    params%poisson_stencil_size = 3
   ! elseif (params%poisson_order == 'CFD_8th') then
   !    params%poisson_stencil_size = 4
   ! elseif (params%poisson_order == 'MST_6th') then
   !    params%poisson_stencil_size = 1
   ! else
   !    call abort(1234567, "Error: no laplacian order specified or not supported!")
   ! endif
   ! if ( params%g < params%poisson_stencil_size ) then
   !    call abort("ERROR: need more ghost nodes for order of supplied finite difference laplacian scheme")
   ! end if

 
end subroutine ini_file_to_params




!> @brief     reads parameters for initializing a bridge from file
subroutine ini_MPI(params, FILE )

   implicit none
   !> pointer to inifile
   type(inifile) ,intent(inout)     :: FILE
   !> params structure of WABBIT
   type(type_params),intent(inout)  :: params

   if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: MPI Communication!"
      write(*,'(" ----------------------------")')
   endif

   ! READ Bridge Parameter
   ! ----------------------
   ! decide if we need a bridge
   call read_param_mpi(FILE, 'BRIDGE', 'connect_with_bridge', params%bridge_exists, .false.)
   if (params%bridge_exists) then               ! if a bridge structure is required
      call read_param_mpi(FILE, 'BRIDGE', 'bridgeCommonMPI', params%bridgeCommonMPI, .false. )
      call read_param_mpi(FILE, 'BRIDGE', 'bridgeFluidMaster', params%bridgeFluidMaster, .false. )
      call read_param_mpi(FILE, 'BRIDGE', 'particleCommand', params%particleCommand, "---" )
   endif

end subroutine ini_MPI

!-------------------------------------------------------------------------!!!!


!> @brief     reads parameters for initializing grid parameters
subroutine ini_domain(params, FILE )
   implicit none
   !> pointer to inifile
   type(inifile) ,intent(inout)     :: FILE
   !> params structure of WABBIT
   type(type_params),intent(inout)  :: params

   integer :: i

   if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Domain"
      write(*,'(" -----------------")')
   endif


   call read_param_mpi(FILE, 'Domain', 'dim', params%dim, 2 )
   if ( .not.(params%dim==2 .or. params%dim==3) ) then
      call abort(234534,"Hawking: ERROR! The idea of 10 dimensions might sound exciting, &
      & but they would cause real problems if you forget where you parked your car. Tip: &
      & Try dim=2 or dim=3 ")
   endif

   params%domain_size=(/ 1.0_rk, 1.0_rk, 0.0_rk /) !default
   call read_param_mpi(FILE, 'Domain', 'domain_size', params%domain_size(1:params%dim), &
      params%domain_size(1:params%dim) )

   params%periodic_BC = .true.
   call read_param_mpi(FILE, 'Domain', 'periodic_BC', params%periodic_BC, params%periodic_BC )

   if (params%rank==0 .and. .not. all(params%periodic_BC)) then
      write(*, '(A)') "Symmetric BC are currently an experimental feature, you should know what you are doing!"
   endif

   params%symmetry_BC = .not. params%periodic_BC
   call read_param_mpi(FILE, 'Domain', 'symmetry_BC', params%symmetry_BC, params%symmetry_BC )

   do i = 1, 3
      if (.not. (params%periodic_BC(i) .or. params%symmetry_BC(i)) ) then
         call abort(92841123, "Get your crowbar: the arrays for periodic_BC and symmetry_BC are incompatible.")
      endif
   enddo
end subroutine ini_domain


!> @brief     reads parameters for initializing grid parameters
subroutine ini_blocks(params, FILE )
   implicit none
   !> pointer to inifile
   type(inifile) ,intent(inout)     :: FILE
   !> params structure of WABBIT
   type(type_params),intent(inout)  :: params
   integer(kind=ik) :: i, g_default, g_RHS_default
   real(kind=rk), dimension(:), allocatable  :: tmp
   logical :: lifted_wavelet

   if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Blocks"
      write(*,'(" -----------------")')
   endif

   ! read number_block_nodes
   params%Bs =read_Bs(FILE, 'Blocks', 'number_block_nodes', params%Bs,params%dim, params%rank)

   ! annoyingly, this is redundant with the definition in setup_wavelet, but as setup_wavelet is part of module_wavelets
   ! which use module_params it cannot be called here (circular dependency....)
   ! compute the number of ghost points, which for CDFXY is X for unlifted and X+Y-1 for lifted wavelets
   ! check if it is lifted or unlifted, which depends on Y (0 for unlifted and >0 for lifted wavelets)

   ! check for X in CDFXY
   ! g_RHS is also dependent on the wavelet due to how we synch each stage:
   ! the third stage (prediction) needs points from the boundary to correctly interpolate values
   if (params%wavelet(4:4) == "2") then
      g_default = 1
      g_RHS_default = 1
   elseif (params%wavelet(4:4) == "4") then
      g_default = 3
      g_RHS_default = 2
   elseif (params%wavelet(4:4) == "6") then
      g_default = 5
      g_RHS_default = 3
   else 
      call abort(2320242, "no default specified for this wavelet...")
   endif

   ! check for Y in CDFXY
   lifted_wavelet = params%wavelet(5:5) /= "0"
   if (params%wavelet(5:5) == "0") then
      ! g stays the same
   elseif (params%wavelet(5:5) == "2") then
      g_default = g_default + 1
   elseif (params%wavelet(5:5) == "4") then
      g_default = g_default + 3
   elseif (params%wavelet(5:5) == "6") then
      g_default = g_default + 5
   elseif (params%wavelet(5:5) == "8") then
      g_default = g_default + 7
   else
      call abort(2320243, "no default specified for this wavelet...")
   endif

   call read_param_mpi(FILE, 'Blocks', 'max_forest_size', params%forest_size, 3 )
   call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%g, g_default )
   call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes_rhs', params%g_RHS, g_RHS_default )  ! might be overwritten later if larger is needed
   call read_param_mpi(FILE, 'Blocks', 'number_blocks', params%number_blocks, -1 )
   call read_param_mpi(FILE, 'Blocks', 'number_equations', params%n_eqn, 1 )
   call read_param_mpi(FILE, 'Blocks', 'number_equations_rhs', params%n_eqn_rhs, params%n_eqn )
   call read_param_mpi(FILE, 'Blocks', 'eps', params%eps, 1e-3_rk )
   call read_param_mpi(FILE, 'Blocks', 'eps_normalized', params%eps_normalized, .false. )
   call read_param_mpi(FILE, 'Blocks', 'eps_norm', params%eps_norm, "Linfty" )
   call read_param_mpi(FILE, 'Blocks', 'azzalini_iterations', params%azzalini_iterations, 1 )
   call read_param_mpi(FILE, 'Blocks', 'threshold_wc', params%threshold_wc, .false. )
   call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params%Jmax, 5 )
   call read_param_mpi(FILE, 'Blocks', 'min_treelevel', params%Jmin, 1 )
   call read_param_mpi(FILE, 'Blocks', 'ini_treelevel', params%Jini, params%Jmin )

   if (params%g_RHS < g_RHS_default) then
      if (params%rank==0) then
         write(*,  '(A, i0, A, i0, A)') "Warning!! 'number_ghost_nodes_rhs' was set smaller as required for Wavelet, adapting it from ", params%g_RHS, " to ", g_RHS_default, " (ignore this if it was not explicitly set)"
      endif
      params%g_RHS = g_RHS_default
   endif

   if (params%g_RHS > params%g) then
      if (params%rank==0) then
         write(*,  '(A, i0, A, i0, A)') "Warning!! 'number_ghost_nodes_rhs' was explicitly set larger than number_ghost_nodes, adapting number_ghost_nodes from ", params%g, " to ", params%g_RHS, " (ignore this if you know what you are doing)"
      endif
      params%g = params%g_RHS
   endif

   if ( params%Jmax < params%Jmin ) then
      call abort(2609181,"Error: Minimal Treelevel cant be larger then Max Treelevel! ")
   end if

   if ( (params%dim==3 .and. params%Jmax > 21) .or. (params%dim==2 .and. params%Jmax > 32) ) then
      ! as we internally convert the treecode to a single integer number, the number of digits is
      ! limited by what we can encode with 64 bits. dim=3 needs 3 bits while dim=2 needs 2 bits resulting in 21 or 32 levels possible
      if (params%dim == 3) then
         call abort(170619,"Error: Max treelevel cannot be larger 21 (64bit long integer problem) ")
      else
         call abort(170619,"Error: Max treelevel cannot be larger 32 (64bit long integer problem) ")
      endif
   end if

   ! read switch to turn on|off mesh refinement
   call read_param_mpi(FILE, 'Blocks', 'adapt_tree', params%adapt_tree, .true. )
   call read_param_mpi(FILE, 'Blocks', 'adapt_inicond', params%adapt_inicond, params%adapt_tree )
   call read_param_mpi(FILE, 'Blocks', 'inicond_refinements', params%inicond_refinements, 0 )
   call read_param_mpi(FILE, 'Blocks', 'inicond_grid_from_file', params%inicond_grid_from_file, "no" )
   call read_param_mpi(FILE, 'Blocks', 'block_dist', params%block_distribution, "sfc_hilbert" )
   call read_param_mpi(FILE, 'Blocks', 'refinement_indicator', params%refinement_indicator, "everywhere" )
   call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params%coarsening_indicator, "threshold-state-vector" )
   call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator_inicond', params%coarsening_indicator_inicond, params%coarsening_indicator )
   call read_param_mpi(FILE, 'Blocks', 'threshold_mask', params%threshold_mask, .false. )
   call read_param_mpi(FILE, 'Blocks', 'force_maxlevel_dealiasing', params%force_maxlevel_dealiasing, .false. )
   call read_param_mpi(FILE, 'Blocks', 'N_dt_per_grid', params%N_dt_per_grid, 1_ik )

   ! for unlifted wavelets we need coarse extension if any WC are changed and reconstruction is requested
   call read_param_mpi(FILE, 'Blocks', 'useCoarseExtension', params%useCoarseExtension, lifted_wavelet &
      .or. params%coarsening_indicator == "threshold-cvs" .or. params%coarsening_indicator=="threshold-image-denoise" )
   call read_param_mpi(FILE, 'Blocks', 'useSecurityZone', params%useSecurityZone, lifted_wavelet &
      .or. params%coarsening_indicator == "threshold-cvs" .or. params%coarsening_indicator=="threshold-image-denoise" )

   ! Which components of the state vector (if indicator is "threshold-state-vector") shall we use?
   ! in ACM, it can be good NOT to apply it to the pressure (=0) or treat the velocity together (>1)
   allocate(tmp(1:params%n_eqn))
   allocate(params%threshold_state_vector_component(1:params%n_eqn))
   ! as default, use ones (all components used for indicator)
   tmp = 1.0_rk
   call read_param_mpi(FILE, 'Blocks', 'threshold_state_vector_component',  tmp, tmp )
   do i = 1, params%n_eqn
      params%threshold_state_vector_component(i) = nint(tmp(i))
   enddo
   deallocate(tmp)

   if (params%Jini > params%Jmax) then
      call abort(1021024, "You cannot set ini_treelevel > max_treelevel!")
   endif

end subroutine ini_blocks

 !-------------------------------------------------------------------------!!!!


 !> @brief     reads parameters for time stepping
subroutine ini_time(params, FILE )
   implicit none
   !> pointer to inifile
   type(inifile) ,intent(inout)     :: FILE
   !> params structure of WABBIT
   type(type_params),intent(inout)  :: params
   real(kind=rk) :: butcher_RK4(1:5,1:5)

   if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Time"
      write(*,'(" --------------")')
   endif

   ! time to reach in simulation
   call read_param_mpi(FILE, 'Time', 'time_max', params%time_max, 1.0_rk )
   ! maximum walltime before ending job
   call read_param_mpi(FILE, 'Time', 'walltime_max', params%walltime_max, 24.0_rk*7-0_rk )
   ! number of time steps to be performed. default value is very large, so if not set
   ! the limit will not be reached
   call read_param_mpi(FILE, 'Time', 'nt', params%nt, 99999999_ik )
   call read_param_mpi(FILE, 'Time', 'time_step_method', params%time_step_method, "RungeKuttaGeneric" )
   call read_param_mpi(FILE, 'Time', 'M_krylov', params%M_krylov, 12 )
   call read_param_mpi(FILE, 'Time', 'krylov_err_threshold', params%krylov_err_threshold, 1.0e-3_rk )
   call read_param_mpi(FILE, 'Time', 'krylov_subspace_dimension', params%krylov_subspace_dimension, "fixed" )
   ! read output write method
   call read_param_mpi(FILE, 'Time', 'write_method', params%write_method, "fixed_freq" )
   ! read output write frequency
   call read_param_mpi(FILE, 'Time', 'write_freq', params%write_freq, 25 )
   ! read output write interval
   call read_param_mpi(FILE, 'Time', 'write_time', params%write_time, 1.0_rk )
   call read_param_mpi(FILE, 'Time', 'write_time_first', params%write_time_first, 0.0_rk )
   ! read output write frequency
   call read_param_mpi(FILE, 'Time', 'walltime_write', params%walltime_write, 99999.9_rk )
   ! read value of fixed time step
   call read_param_mpi(FILE, 'Time', 'dt_fixed', params%dt_fixed, 0.0_rk )
   ! read value of fixed time step
   call read_param_mpi(FILE, 'Time', 'dt_max', params%dt_max, 0.0_rk )
   ! read CFL number
   call read_param_mpi(FILE, 'Time', 'CFL', params%CFL, 0.5_rk )
   ! number of stages "s" for the RungeKuttaChebychev method. Memory is always 6 registers
   ! independent of stages.
   call read_param_mpi(FILE, 'Time', 's', params%s, 4 )

   call read_param_mpi(FILE, 'Time', 'RKC_custom_scheme', params%RKC_custom_scheme, .false. )
   if (params%RKC_custom_scheme) then
      call read_param_mpi(FILE, 'Time', 'RKC_mu', params%RKC_mu(1:params%s) )
      call read_param_mpi(FILE, 'Time', 'RKC_mu_tilde', params%RKC_mu_tilde(1:params%s) )
      call read_param_mpi(FILE, 'Time', 'RKC_nu', params%RKC_nu(1:params%s) )
      call read_param_mpi(FILE, 'Time', 'RKC_gamma_tilde', params%RKC_gamma_tilde(1:params%s) )
      call read_param_mpi(FILE, 'Time', 'RKC_c', params%RKC_c(1:params%s) )
   endif


   ! read butcher tableau (set default value to RK4)
   butcher_RK4(1,1:5) = (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/)
   butcher_RK4(2,1:5) = (/0.5_rk, 0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk/)
   butcher_RK4(3,1:5) = (/0.5_rk, 0.0_rk, 0.5_rk, 0.0_rk, 0.0_rk/)
   butcher_RK4(4,1:5) = (/1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk/)
   butcher_RK4(5,1:5) = (/0.0_rk, 1.0_rk/6.0_rk, 1.0_rk/3.0_rk, 1.0_rk/3.0_rk, 1.0_rk/6.0_rk/)

   call read_param_mpi(FILE, 'Time', 'butcher_tableau', params%butcher_tableau, butcher_RK4)
end subroutine ini_time