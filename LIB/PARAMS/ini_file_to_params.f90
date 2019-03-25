!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name ini_file_to_params.f90
!> \version 0.5
!> \author msr
!
!> \brief distribute blocks at start => create light data array
!
!>
!! input:    - filename \n
!! output:   - filled parameter struct \n
!!
!!
!! = log ======================================================================================
!! \n
!! 25/01/17    - create \n
!! 29/01/17    - add filter parameter \n
!! 30/01/17    - add automatic memory management
!
! ********************************************************************************************

subroutine ini_file_to_params( params, filename )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)               :: params

    character(len=*), intent(in)                   :: filename

    ! process rank
    integer(kind=ik)                                :: rank
    ! number of processes
    integer(kind=ik)                                :: number_procs
    ! inifile structure
    type(inifile)                                   :: FILE
    ! maximum memory available on all cpus
    real(kind=rk)                                   :: maxmem, mem_per_block, max_neighbors, nstages
    ! string read from command line call
    character(len=80)                               :: memstring
    !
    integer(kind=ik)                                :: d,i, Nblocks_Jmax, g, Neqn, Nrk
    integer(kind=ik), dimension(3)                  :: Bs

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

!---------------------------------------------------------------------------------------------
! main body

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    call ini_domain(params, FILE )
    call ini_blocks(params,FILE)
    call ini_time(params,FILE)


    !**************************************************************************
    ! read INITIAL CONDITION parameters

    ! which physics module is used? (note that the initialization of different parameters takes
    ! place in those modules, i.e., they are not read here.)
    call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )

    ! if the initial condition is read from file, it is handled by wabbit itself, i.e. not
    ! by the physics modules. the pyhsics modules cannot do this, because they just see 'blocks'
    ! and never the entire grid as such.
    call read_param_mpi(FILE, 'Physics', 'read_from_files', params%read_from_files, .false. )

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
    ! order of predictor for refinement
    call read_param_mpi(FILE, 'Discretization', 'order_predictor', params%order_predictor, "---" )
    ! filter frequency
    call read_param_mpi(FILE, 'Discretization', 'filter_type', params%filter_type, "no_filter" )
    if (params%filter_type /= "no_filter") then
        call read_param_mpi(FILE, 'Discretization', 'filter_freq', params%filter_freq, -1 )
    endif

    !***************************************************************************
    ! read statistics parameters
    !
    ! data exchange method
    call read_param_mpi(FILE, 'Statistics', 'nsave_stats', params%nsave_stats, 99999999_ik )
    call read_param_mpi(FILE, 'Statistics', 'tsave_stats', params%tsave_stats, 9999999.9_rk )
    !> assume start at time 0.0 /todo change if start with reloaded data
    params%next_stats_time = 0.0_rk + params%tsave_stats


    !***************************************************************************
    ! read DEBUG parameters
    !
    ! timing flag
    !call read_param_mpi(FILE, 'Debug', 'timing', params%timing, .true. )
    call read_param_mpi(FILE, 'Timing', 'write_individual_timings', params%write_individual_timings, .false. )
    ! unit test treecode flag
    call read_param_mpi(FILE, 'Debug', 'test_treecode', params%test_treecode, .false.)
    call read_param_mpi(FILE, 'Debug', 'test_ghost_nodes_synch', params%test_ghost_nodes_synch, .false.)
    call read_param_mpi(FILE, 'Debug', 'check_redundant_nodes', params%check_redundant_nodes, .false.)

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
    if ( (params%n_ghosts < 4) .and. (params%order_predictor == 'multiresolution_4th') ) then
        call abort("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%n_ghosts < 2) .and. (params%order_predictor == 'multiresolution_2nd') ) then
        call abort("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%n_ghosts < 2) .and. (params%order_discretization == 'FD_4th_central_optimized') ) then
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
    call read_param_mpi(FILE, 'Domain', 'periodic_BC', params%periodic_BC(1:params%dim), &
                                                       params%periodic_BC(1:params%dim) )
  end subroutine ini_domain


!> @brief     reads parameters for initializing grid parameters
  subroutine ini_blocks(params, FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    !> params structure of WABBIT
    type(type_params),intent(inout)  :: params
    !> power used for dimensionality (d=2 or d=3)
    integer(kind=ik) :: i, n_entries=0
    real(kind=rk), dimension(:), allocatable  :: tmp
    character(len=80) :: Bs_str, Bs_conc
    character(:), allocatable :: Bs_short
    real(kind=rk), dimension(3) :: Bs_real
    if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Blocks"
      write(*,'(" -----------------")')
    endif

    ! read number_block_nodes
    call read_param_mpi(FILE, 'Blocks', 'number_block_nodes', Bs_str, "empty")
    call merge_blancs(Bs_str)
    Bs_short=trim(Bs_str)
    call count_entries(Bs_short, " ", n_entries)
    if (Bs_str .eq. "empty") then
      Bs_conc="17 17 17"
    elseif (n_entries==1) then
      if (params%dim==3) then
        Bs_conc=Bs_short // " " // Bs_short // " " // Bs_short
      elseif (params%dim==2) then
        Bs_conc=Bs_short//" "//Bs_short//" 1"
      endif
    elseif (n_entries==2) then
      if (params%dim==3) then
          call abort(231191737,"ERROR: You only gave two values for Bs, but want three to be read...")
      elseif (params%dim==2) then
        Bs_conc=Bs_short//" 1"
      endif
    elseif (n_entries==3) then
      if (params%dim==2) then
        call abort(231191738,"ERROR: You gave three values for Bs, but only want two to be read...")
      elseif (params%dim==3) then
        Bs_conc=trim(adjustl(Bs_str))
      endif
    elseif (n_entries .gt. 3) then
      call abort(231191752,"ERROR: You gave too many arguments for Bs...")
    endif
    read(Bs_conc, *) Bs_real
    params%Bs = int(Bs_real)

    call read_param_mpi(FILE, 'Blocks', 'max_forest_size', params%forest_size, 1 )
    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%n_ghosts, 1 )
    call read_param_mpi(FILE, 'Blocks', 'number_blocks', params%number_blocks, -1 )
    call read_param_mpi(FILE, 'Blocks', 'number_equations', params%n_eqn, 1 )
    call read_param_mpi(FILE, 'Blocks', 'eps', params%eps, 1e-3_rk )
    call read_param_mpi(FILE, 'Blocks', 'eps_normalized', params%eps_normalized, .false. )
    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params%max_treelevel, 5 )
    call read_param_mpi(FILE, 'Blocks', 'min_treelevel', params%min_treelevel, 1 )
    if ( params%max_treelevel < params%min_treelevel ) then
      call abort(2609181,"Error: Minimal Treelevel cant be larger then Max Treelevel! ")
    end if
    ! read switch to turn on|off mesh refinement
    call read_param_mpi(FILE, 'Blocks', 'adapt_mesh', params%adapt_mesh, .true. )
    call read_param_mpi(FILE, 'Blocks', 'adapt_inicond', params%adapt_inicond, params%adapt_mesh )
    call read_param_mpi(FILE, 'Blocks', 'inicond_refinements', params%inicond_refinements, 0 )
    call read_param_mpi(FILE, 'Blocks', 'block_dist', params%block_distribution, "---" )
    call read_param_mpi(FILE, 'Blocks', 'loadbalancing_freq', params%loadbalancing_freq, 1 )
    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params%coarsening_indicator, "threshold-state-vector" )
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
      ! read output write frequency
      call read_param_mpi(FILE, 'Time', 'write_time', params%write_time, 1.0_rk )
      ! assume start at time 0.0
      params%next_write_time = 0.0_rk + params%write_time
      ! read value of fixed time step
      call read_param_mpi(FILE, 'Time', 'dt_fixed', params%dt_fixed, 0.0_rk )
      ! read value of fixed time step
      call read_param_mpi(FILE, 'Time', 'dt_max', params%dt_max, 0.0_rk )
      ! read CFL number
      call read_param_mpi(FILE, 'Time', 'CFL', params%CFL, 0.5_rk )

      ! read butcher tableau (set default value to RK4)
      butcher_RK4(1,1:5) = (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/)
      butcher_RK4(2,1:5) = (/0.5_rk, 0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk/)
      butcher_RK4(3,1:5) = (/0.5_rk, 0.0_rk, 0.5_rk, 0.0_rk, 0.0_rk/)
      butcher_RK4(4,1:5) = (/1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk/)
      butcher_RK4(5,1:5) = (/0.0_rk, 1.0_rk/6.0_rk, 1.0_rk/3.0_rk, 1.0_rk/3.0_rk, 1.0_rk/6.0_rk/)

      call read_param_mpi(FILE, 'Time', 'butcher_tableau', params%butcher_tableau, butcher_RK4)
    end subroutine ini_time

    !-------------------------------------------------------------------------!
    !> @brief remove (multiple) blancs as seperators in a string
    subroutine merge_blancs(string_merge)
      ! this routine removes blanks at the beginning and end of an string
      ! and multiple blanks which are right next to each other

      implicit none
      character(len=*), intent(inout) :: string_merge
      integer(kind=ik) :: i, j, len_str, count

      len_str = len(string_merge)
      count = 0

      string_merge = string_merge
      do i=1,len_str-1
        if (string_merge(i:i)==" " .and. string_merge(i+1:i+1)==" ") then
          count = count + 1
          string_merge(i+1:len_str-1) = string_merge(i+2:len_str)
        end if
      end do

      string_merge = adjustl(string_merge)

    end subroutine merge_blancs

    !-------------------------------------------------------------------------!
    !> @brief count number of vector elements in a string
    subroutine count_entries(string_cnt, seperator, n_entries)
      ! only to be used after merged blaks
      ! this routine counts the seperators and gives back this value +1

      implicit none
      character(len=1), intent(in) :: seperator
      character(len=*), intent(in) :: string_cnt
      integer(kind=ik), intent(out) :: n_entries
      integer(kind=ik) :: count_seperator, i, l_string
      character(len=Len_trim(adjustl(string_cnt))):: string_trim

      count_seperator = 0
      string_trim = trim(adjustl(string_cnt))
      l_string = LEN(string_trim)

      do i=1,l_string
        if (string_trim(i:i)==seperator) then
          count_seperator = count_seperator + 1
        end if
      end do

      n_entries = count_seperator + 1

    end subroutine count_entries
