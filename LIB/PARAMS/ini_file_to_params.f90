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
    real(kind=rk)                                   :: maxmem, byte_int=4.0_rk, &
    byte_real=8.0_rk, Ncpu, Nfriend, Bs, g, max_neighbors, Nrk4, Neqn,Jmax
    ! power used for dimensionality (d=2 or d=3)
    integer(kind=ik)                                :: d,i, Nblocks_Jmax
    real(kind=rk), dimension(:), allocatable        :: tmp
    ! string read from command line call
    character(len=80)                               :: memstring

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

    !***************************************************************************
    ! read BLOCK parameters
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
    ! use non-uniform mesh correction
    call read_param_mpi(FILE, 'Blocks', 'non_uniform_mesh_correction', params%non_uniform_mesh_correction, .true. )
    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params%coarsening_indicator, "threshold-state-vector" )

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

    ! domain size
    call read_param_mpi(FILE, 'DomainSize', 'Lx', params%Lx, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', params%Ly, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Lz', params%Lz, 1.0_rk )

    ! saving options.
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params%N_fields_saved, 3 )


    !***************************************************************************
    ! read TIME parameters
    !
    ! time to reach in simulation
    call read_param_mpi(FILE, 'Time', 'time_max', params%time_max, 1.0_rk )
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
    call init_MPI(params, FILE )


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

            ! read memory from command line and convert to bytes
            read(memstring(10:len_trim(memstring)-2),* ) maxmem
            ! how much memory is reserved, in bytes?
            maxmem = maxmem * 1000.0 * 1000.0 * 1000.0 ! in bytes
            if ( params%threeD_case ) then
                d = 3
                max_neighbors = 56.0
            else
                max_neighbors = 12.0
                d = 2
            endif

            Ncpu = real(params%number_procs)
            Neqn = real(params%number_data_fields)
            Bs = real(params%number_block_nodes)
            g = real(params%number_ghost_nodes)
            ! we assume sort of a worst-case scenario and assume 20 friends.
            ! in 2D that is too much (11 in practice), in 3d we dont know yet
            Nfriend = min(Ncpu, 20.0_rk)
            Jmax = real(params%max_treelevel)
            Nrk4 = real(size(params%butcher_tableau, 1)+1)

            maxmem = 0.90_rk * maxmem / Ncpu

            if ( d == 2 ) then
                ! in 2d, we try to be precise and take all fat arrays into account. this worked well
                params%number_blocks = ceiling( maxmem/( 2.0*max_neighbors*3.0*Nfriend*byte_int &
                + 2.0*max_neighbors*(Bs+g+1.0)*((g+1.)**(d-1))*Neqn*Nfriend*byte_real &
                + (Nrk4+1.)*((Bs+2.0*g)**d)*Neqn*byte_real &
                + max_neighbors*byte_int &
                + Ncpu*(Jmax+2.0+2.0+1.0)*byte_int ) )
            else
                ! in 3d, use a simpler formula which takes only hvy and lgt data into account.
                ! the buffers are negligible in front of it. I still do not understand it fully. 21/06/2018
                params%number_blocks = ceiling( maxmem/( (Nrk4+1.)*((Bs+2.0*g)**d)*Neqn*byte_real + Ncpu*(Jmax+2.0+2.0+1.0)*byte_int ) )
            endif

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

    if (params%rank==0) then
        write(*,'("INIT: resetting some *.t files.")')

        open (44, file='dt.t', status='replace')
        close(44)
        open (44, file='timesteps_info.t', status='replace')
        close(44)
        open (44, file='blocks_per_mpirank.t', status='replace')
        close(44)
        open (44, file='blocks_per_mpirank_rhs.t', status='replace')
        close(44)
        open (44, file='eps_norm.t', status='replace')
        close(44)
    endif
end subroutine ini_file_to_params




!> @brief     reads parameters for initializing a bridge from file
  subroutine init_MPI(params, FILE )
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

    call read_param_mpi(FILE, 'MPI', 'mpi_data_exchange', params%mpi_data_exchange, "---" )

    ! READ Bridge Parameter
    ! ----------------------
    ! decide if we need a bridge
    call read_param_mpi(FILE, 'BRIDGE', 'connect_with_bridge', params%bridge_exists, .false.)
    if (params%bridge_exists) then               ! if a bridge structure is required
      call read_param_mpi(FILE, 'BRIDGE', 'bridgeCommonMPI', params%bridgeCommonMPI, .false. )
      call read_param_mpi(FILE, 'BRIDGE', 'bridgeFluidMaster', params%bridgeFluidMaster, .false. )
      call read_param_mpi(FILE, 'BRIDGE', 'particleCommand', params%particleCommand, "---" )
    endif

  end subroutine init_MPI

!!!!-------------------------------------------------------------------------!!!!
