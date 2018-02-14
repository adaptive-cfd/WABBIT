!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name    main.f90
!> \version 0.5
!> \author  msr
!
!> \brief main program, init all data, start time loop, output on screen during program run
!
!>
!!
!! = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4 \n
!! 23/11/16 - use computing time array for simple performance tests \n
!! 07/12/16 - now uses heavy work data array \n
!! 25/01/17 - switch to 3D, v0.5
!! 30/01/18 - start postprocessing/unittesting
! ********************************************************************************************
!> \image html rhs.svg width=600
!> \image html rhs.eps

program main_post

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1
    ! if help mode is on we print some output on screen
    logical                             :: help

    type (type_params)                  :: params
    character(len=80)                   :: mode
    character(len=80)                   :: filename

!---------------------------------------------------------------------------------------------
! main body

    ! init mpi
    call MPI_Init(ierr)
    ! determine process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    params%rank         = rank
    ! determine process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)
    params%number_procs = number_procs
    ! output MPI status
    if (rank==0) then
        write(*,'(80("_"))')
        write(*, '("MPI: using ", i5, " processes")') params%number_procs
    end if


    ! start time
    sub_t0 = MPI_Wtime()
    call cpu_time(t0)

    ! are we running in 2D or 3D mode? Check that from the command line call.
    call decide_if_running_2D_or_3D(params)

    !---------------------------------------------------------------------------
    ! Initialize parameters and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    ! get the second command line argument: this should be the ini-file name
    call get_command_argument( 2, mode )

    if (mode=="--h") then
        help=.true.
    else
        if (rank==0) write(*,'("Starting postprocessing in ", a20, "mode")'), mode
        help=.false.
    end if

    select case(mode)
    case("--sparse-to-dense")
        call sparse_to_dense(help, params)
    case("--vorticity")
        call compute_vorticity_post(help, params) 
    case("--keyvalues")
        call get_command_argument(3,filename)
        call keyvalues(filename, params)
    case default
    if (params%rank==0) then
        write(*,*) "Available Postprocessing tools are:"
        write(*,*) "--sparse-to-dense"
        write(*,*) "--vorticity"
        write(*,*) "--keyvalues"
        write(*,*) "Your postprocessing option is "// trim(adjustl(mode)) //", which I don't know"
    end if
    end select

    ! end mpi
    call MPI_Finalize(ierr)

end program main_post
