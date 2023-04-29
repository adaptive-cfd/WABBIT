! initialize params struct from INI file
subroutine ini_file_to_params( params, filename )
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
    real(kind=rk)                                   :: maxmem, mem_per_block, max_neighbors, nstages
    ! string read from command line call
    character(len=cshort)                               :: memstring
    !
    integer(kind=ik)                                :: d,i, Nblocks_Jmax, g, Neqn, Nrk
    integer(kind=ik), dimension(3)                  :: Bs

    rank         = params%rank
    number_procs = params%number_procs


    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    ! which physics module is used? (note that the initialization of different parameters takes
    ! place in those modules, i.e., they are not read here.)
    call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )


    call ini_domain(params, FILE )
    call ini_blocks(params,FILE)
    call ini_time(params,FILE)

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
        call read_param_mpi(FILE, 'Physics', 'input_files', params%input_files, params%input_files)
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
    ! WABBIT needs to know about the mask function (if penalization is used): does it contain
    ! a time-dependent-part (e.g. moving obstacles, time-dependent forcing)? does it contain
    ! a time-independent part (fixed walls, homogeneous forcing)? or both? WABBIT needs to know
    ! that since we try to create the time-independent mask function only once, but the time-dependent
    ! part of course in every time step.
    call read_param_mpi(FILE, 'VPM', 'penalization', params%penalization, .false.)
    call read_param_mpi(FILE, 'VPM', 'mask_time_dependent_part', params%mask_time_dependent_part, .true.)
    call read_param_mpi(FILE, 'VPM', 'mask_time_independent_part', params%mask_time_independent_part, .true.)
    call read_param_mpi(FILE, 'VPM', 'dont_use_pruned_tree_mask', params%dont_use_pruned_tree_mask, .false.)

    ! decide if we use hartens point value multiresolution transform, which uses a coarsening operator
    ! that just takes every 2nd grid point or biorthogonal wavlets, which apply a smoothing filter (lowpass)
    ! prior to downsampling.
    call read_param_mpi(FILE, 'Wavelet', 'wavelet', params%wavelet, 'CDF40')

    !***************************************************************************
    ! read DEBUG parameters
    !
    ! unit test treecode flag
    call read_param_mpi(FILE, 'Debug', 'test_treecode', params%test_treecode, .false.)
    call read_param_mpi(FILE, 'Debug', 'test_ghost_nodes_synch', params%test_ghost_nodes_synch, .true.)
    call read_param_mpi(FILE, 'Debug', 'test_wavelet_decomposition', params%test_wavelet_decomposition, .true.)

    !***************************************************************************
    ! read MPI parameters
    !
    ! data exchange method
    call ini_MPI(params, FILE )

    ! clean up
    if (params%rank==0) write(*,'("INIT: cleaning ini file")')
    call clean_ini_file_mpi(FILE)


    ! check ghost nodes number
    if (params%rank==0) write(*,'("INIT: checking if g and predictor work together")')
    if ( (params%g < 2) .and. (params%order_predictor == 'multiresolution_4th') ) then
        call abort("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%g < 6) .and. (params%wavelet=='CDF44') )  then
        call abort(050920194, "ERROR: for CDF44 wavelet, 6 ghost nodes are required")
    end if
    if ( (params%g < 1) .and. (params%order_predictor == 'multiresolution_2nd') ) then
        call abort("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%g < 3) .and. (params%order_discretization == 'FD_4th_central_optimized') ) then
        call abort("ERROR: need more ghost nodes for given derivative order")
    end if
    if ( (params%g < 2) .and. (params%order_discretization == 'FD_4th_central') ) then
        call abort("ERROR: need more ghost nodes for given derivative order")
    end if

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
    !> power used for dimensionality (d=2 or d=3)
    integer(kind=ik) :: i
    real(kind=rk), dimension(:), allocatable  :: tmp
    if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Blocks"
      write(*,'(" -----------------")')
    endif

    ! read number_block_nodes
    params%Bs =read_Bs(FILE, 'Blocks', 'number_block_nodes', params%Bs,params%dim)

    call read_param_mpi(FILE, 'Blocks', 'max_forest_size', params%forest_size, 3 )
    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%g, 1 )
    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes_rhs', params%g_rhs, params%g )
    call read_param_mpi(FILE, 'Blocks', 'number_blocks', params%number_blocks, -1 )
    call read_param_mpi(FILE, 'Blocks', 'number_equations', params%n_eqn, 1 )
    call read_param_mpi(FILE, 'Blocks', 'eps', params%eps, 1e-3_rk )
    call read_param_mpi(FILE, 'Blocks', 'eps_normalized', params%eps_normalized, .false. )
    call read_param_mpi(FILE, 'Blocks', 'eps_norm', params%eps_norm, "Linfty" )
    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params%Jmax, 5 )
    call read_param_mpi(FILE, 'Blocks', 'min_treelevel', params%Jmin, 1 )

    if ( params%Jmax < params%Jmin ) then
        call abort(2609181,"Error: Minimal Treelevel cant be larger then Max Treelevel! ")
    end if

    if ( params%Jmax > 18 ) then
        ! as we internally convert the treecode to a single integer number, the number of digits is
        ! limited by that type. The largest 64-bit integer is 9 223 372 036 854 775 807
        ! which is 19 digits, but the 18th digit cannot be arbitrarily set. Therefore, 18 refinement levels
        ! are the maximum this code can currently perform.
        call abort(170619,"Error: Max treelevel cannot be larger 18 (64bit long integer problem) ")
    end if

    ! read switch to turn on|off mesh refinement
    call read_param_mpi(FILE, 'Blocks', 'adapt_tree', params%adapt_tree, .true. )
    call read_param_mpi(FILE, 'Blocks', 'adapt_inicond', params%adapt_inicond, params%adapt_tree )
    call read_param_mpi(FILE, 'Blocks', 'inicond_refinements', params%inicond_refinements, 0 )
    call read_param_mpi(FILE, 'Blocks', 'inicond_grid_from_file', params%inicond_grid_from_file, "no" )
    call read_param_mpi(FILE, 'Blocks', 'block_dist', params%block_distribution, "sfc_hilbert" )
    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params%coarsening_indicator, "threshold-state-vector" )
    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator_inicond', params%coarsening_indicator_inicond, params%coarsening_indicator )
    call read_param_mpi(FILE, 'Blocks', 'threshold_mask', params%threshold_mask, .false. )
    call read_param_mpi(FILE, 'Blocks', 'force_maxlevel_dealiasing', params%force_maxlevel_dealiasing, .false. )
    call read_param_mpi(FILE, 'Blocks', 'N_dt_per_grid', params%N_dt_per_grid, 1_ik )

    ! Which components of the state vector (if indicator is "threshold-state-vector") shall we
    ! use? in ACM, it can be good NOT to apply it to the pressure.
    allocate(tmp(1:params%n_eqn))
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    ! as default, use ones (all components used for indicator)
    tmp = 1.0_rk
    call read_param_mpi(FILE, 'Blocks', 'threshold_state_vector_component',  tmp, tmp )
    do i = 1, params%n_eqn
        if (tmp(i)>0.0_rk) then
            params%threshold_state_vector_component(i) = .true.
        else
            params%threshold_state_vector_component(i) = .false.
        endif
    enddo
    deallocate(tmp)

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



    !-------------------------------------------------------------------------!
    !> @brief Read Bs from inifile for unknown number of Bs in inifile
    function read_Bs(FILE, section, keyword, default_Bs, dims) result(Bs)
        type(inifile) ,intent(inout)     :: FILE
        character(len=*), intent(in)    :: section ! What section do you look for? for example [Resolution]
        character(len=*), intent(in)    :: keyword ! what keyword do you
        integer(kind=ik), intent(in)    :: default_Bs(:)
        integer(kind=ik), intent(in)    :: dims !number of dimensions
        integer(kind=ik):: Bs(3)
        integer(kind=ik):: i, n_entries
        character(len=cshort):: output

        Bs = 1
        ! read number_block_nodes
        call read_param_mpi(FILE, section, keyword, output, "empty")

        if (trim(output) .eq. "empty") then
            write(*,'("Warning!! ", A, "[",A,"] is empty! Using default! ")') keyword, section
            Bs = default_Bs

        else
            call count_entries(output, n_entries)

            ! check if the number of entries is valid
            if (n_entries > dims) call abort(10519,"Dimensions and number of Bs entries dissagree!")

            ! Cast the output string into the integer
            read(output,*) (Bs(i), i=1,n_entries)

            ! If only one Bs is given in the ini file, we duplicate it
            ! for the rest of the Bs array:
            if (n_entries==1) then
                Bs(1:dims) = Bs(1)
            endif
        endif

        do i = 1, dims
            if (mod(Bs(i), 2) /= 0) then
                write(*,*) "Bs=", Bs
                call abort(202392929, "Block-size must be EVEN number")
            end if
        end do
    end function
