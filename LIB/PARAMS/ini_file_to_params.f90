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
    !> inifile name
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
    integer(kind=ik)                                :: d,i, Nblocks_Jmax, Bs, g, Neqn, Nrk

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


    call read_param_mpi(FILE, 'Dimensionality', 'dim', params%dim, 2 )
    ! domain size
    call read_param_mpi(FILE, 'DomainSize', 'Lx', params%Lx, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', params%Ly, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Lz', params%Lz, 1.0_rk )
    ! saving options.
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params%N_fields_saved, 3 )

    call ini_blocks(params,FILE)
    call ini_time(params,FILE)


    !**************************************************************************
    ! read INITIAL CONDITION parameters

    ! initial condition. NOTE: nowadays, there is only two distinct ones: read_from_files
    ! which, well, reads some files. This is the same for all physics modules. The only
    ! other initial condition is "physics-module", which means the modules decide what inicond
    ! is set.

    call read_param_mpi(FILE, 'Physics', 'initial_cond', params%initial_cond, "---" )

    if (params%initial_cond == 'read_from_files') then
        ! read variable names
        allocate( params%input_files( params%number_data_fields ) )

        params%input_files = "---"
        call read_param_mpi(FILE, 'Physics', 'input_files', params%input_files, params%input_files)
    end if

    !***************************************************************************
    ! read DISCRETIZATION parameters
    !
    ! discretization order
    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params%order_discretization, "---" )
    ! order of predictor for refinement
    call read_param_mpi(FILE, 'Discretization', 'order_predictor', params%order_predictor, "---" )
    ! boundary condition
    call read_param_mpi(FILE, 'Discretization', 'boundary_cond', params%boundary_cond, "---" )

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
    ! debug flag
    call read_param_mpi(FILE, 'Debug', 'debug', params%debug, .true. )
    ! unit test time_stepper flag
    call read_param_mpi(FILE, 'Debug', 'test_time_stepper', params%test_time_stepper, .false.)
    ! unit test spatial flag
    call read_param_mpi(FILE, 'Debug', 'test_spatial', params%test_spatial, .false.)
    ! unit test wavelet compression flag
    call read_param_mpi(FILE, 'Debug', 'test_wavelet_comp', params%test_wavelet_comp, .false.)
    ! unit test treecode flag
    call read_param_mpi(FILE, 'Debug', 'test_treecode', params%test_treecode, .false.)
    call read_param_mpi(FILE, 'Debug', 'test_ghost_nodes_synch', params%test_ghost_nodes_synch, .false.)
    call read_param_mpi(FILE, 'Debug', 'check_redundant_nodes', params%check_redundant_nodes, .true.)

    !***************************************************************************
    ! read PHYSICS parameters
    !
    ! first: read physics type
    call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )
    call read_param_mpi(FILE, 'Physics', 'initial_cond', params%initial_cond, "physics-module" )
    !***************************************************************************
    ! read MPI parameters
    !
    ! data exchange method
    call ini_MPI(params, FILE )


    ! NOTE: this routine initializes WABBIT AND NOT THE PHYSICS MODULES THEMSELVES!

    !---------------------------------------------------------------------------
    ! Automatic memory management. If specified --memory=0.3GB in the call line,
    ! wabbit will automatically select the number of blocks per rank to be allocated
    ! to consume this amount of memory. helpful on local machines.
    !---------------------------------------------------------------------------
    ! loop over all arguments until you find the string "--memory" in them (or not)
    ! this ensures compatibility with future versions, as the argument may be anywhere in the call.
    do i = 1, command_argument_count()
        call get_command_argument(i, memstring)
        ! is memory limitation used?
        if ( index(memstring,"--memory=")==1 ) then
            if (params%rank==0) write(*,'(80("-"))')
            if (params%rank==0) write(*,'("INIT: automatic selection of blocks per rank is active!")')
            if (params%rank==0) write(*,'(80("-"))')

            ! read memory from command line (in GB)
            read(memstring(10:len_trim(memstring)-2),* ) maxmem
            ! memory per MPIRANK (in GB)
            maxmem = maxmem / dble(params%number_procs)

            if (params%rank==0) write(*,'("INIT: memory-per-rank: ",f9.4,"GB")') maxmem

            if ( params%threeD_case ) then
                d = 3
                max_neighbors = 56.0
            else
                d = 2
                max_neighbors = 12.0
            endif

            Bs      = params%number_block_nodes
            g       = params%number_ghost_nodes
            Neqn    = params%number_data_fields
            Nrk     = max( size(params%butcher_tableau,1)-1, params%N_fields_saved ) + 2
            nstages = 2.0

            mem_per_block = real(Neqn) * real(Nrk) * (real(Bs+2*g))**d & ! hvy_work+hvy_block
            + 2.0 * nstages * real(Neqn) * real((Bs+2*g)**d - Bs**d) &  ! real buffer ghosts
            + 2.0 * nstages * max_neighbors * 3 * real(params%number_procs) / 2.0 ! int bufer (4byte hence /2)

            ! in GB:
            mem_per_block = mem_per_block * 8.0e-9

            params%number_blocks = nint( maxmem / mem_per_block)

            if (params%rank==0) then
                write(*,'("INIT: for the desired memory we can allocate ",i8," blocks per rank")') params%number_blocks
                write(*,'("INIT: we allocated ",i8," blocks per rank (total: ",i8," blocks) ")') params%number_blocks, &
                params%number_blocks*params%number_procs
            endif
        endif
    end do

    ! clean up
    if (params%rank==0) write(*,'("INIT: cleaning ini file")')
    call clean_ini_file_mpi(FILE)


    ! check ghost nodes number
    if (params%rank==0) write(*,'("INIT: checking if g and predictor work together")')
    if ( (params%number_ghost_nodes < 4) .and. (params%order_predictor == 'multiresolution_4th') ) then
        call abort("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%number_ghost_nodes < 2) .and. (params%order_predictor == 'multiresolution_2nd') ) then
        call abort("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%number_ghost_nodes < 2) .and. (params%order_discretization == 'FD_4th_central_optimized') ) then
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
  subroutine ini_blocks(params, FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    !> params structure of WABBIT
    type(type_params),intent(inout)  :: params
    !> power used for dimensionality (d=2 or d=3)
    integer(kind=ik)                                :: i
    real(kind=rk), dimension(:), allocatable        :: tmp

    if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Blocks"
      write(*,'(" -----------------")')
    endif

    ! read number_block_nodes
    call read_param_mpi(FILE, 'Blocks', 'number_block_nodes', params%number_block_nodes, 1 )
    ! read number_ghost_nodes
    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%number_ghost_nodes, 1 )
    ! read number_blocks
    call read_param_mpi(FILE, 'Blocks', 'number_blocks', params%number_blocks, 1 )

    ! read number_data_fields
    call read_param_mpi(FILE, 'Blocks', 'number_data_fields', params%number_data_fields, 1 )
    ! set number of fields in heavy work data
    ! every datafield has 5 additional fields: old, k1, k2, k3, k4
    params%number_fields = params%number_data_fields*5
    ! read threshold value
    call read_param_mpi(FILE, 'Blocks', 'eps', params%eps, 1e-3_rk )
    call read_param_mpi(FILE, 'Blocks', 'eps_normalized', params%eps_normalized, .false. )
    ! read treelevel bounds
    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params%max_treelevel, 5 )
    call read_param_mpi(FILE, 'Blocks', 'min_treelevel', params%min_treelevel, 1 )
    ! read switch to turn on|off mesh refinement
    call read_param_mpi(FILE, 'Blocks', 'adapt_mesh', params%adapt_mesh, .true. )
    call read_param_mpi(FILE, 'Blocks', 'adapt_inicond', params%adapt_inicond, params%adapt_mesh )
    call read_param_mpi(FILE, 'Blocks', 'inicond_refinements', params%inicond_refinements, 0 )
    ! block distribution
    call read_param_mpi(FILE, 'Blocks', 'block_dist', params%block_distribution, "---" )
    call read_param_mpi(FILE, 'Blocks', 'loadbalancing_freq', params%loadbalancing_freq, 1 )
    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params%coarsening_indicator, "threshold-state-vector" )
    call read_param_mpi(FILE, 'Blocks', 'force_maxlevel_dealiasing', params%force_maxlevel_dealiasing, .false. )

    ! Which components of the state vector (if indicator is "threshold-state-vector") shall we
    ! use? in ACM, it can be good NOT to apply it to the pressure.
    allocate(tmp(1:params%number_data_fields))
    allocate(params%threshold_state_vector_component(1:params%number_data_fields))
    ! as default, use ones (all components used for indicator)
    tmp = 1.0_rk
    call read_param_mpi(FILE, 'Blocks', 'threshold_state_vector_component',  tmp, tmp )
    do i = 1, params%number_data_fields
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
      call read_param_mpi(FILE, 'Time', 'butcher_tableau', params%butcher_tableau, &
      reshape((/ 0.0_rk, 0.5_rk, 0.5_rk, 1.0_rk, 0.0_rk, &
      0.0_rk, 0.5_rk, 0.0_rk, 0.0_rk, 1.0_rk/6.0_rk, &
      0.0_rk, 0.0_rk, 0.5_rk, 0.0_rk, 1.0_rk/3.0_rk,&
      0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk/3.0_rk,&
      0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/6.0_rk /), (/ 5,5 /)))


    end subroutine ini_time
