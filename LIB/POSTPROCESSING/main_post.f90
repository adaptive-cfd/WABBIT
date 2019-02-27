!> \file
! WABBIT
!> \name    main.f90
!> \version 0.5
!> \author  sm
!
!> \brief main postprocessing program. get command argument to decide which postprocessing mode to run
!
! = log ======================================================================================
!
!> \version 30/1/2018 - create hashcode: commit 13cb3d25ab12e20cb38e5b87b9a1e27a8fe387e8
! ********************************************************************************************

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

    type (type_params)                  :: params
    character(len=80)                   :: mode
    character(len=80)                   :: filename, key1, key2

!---------------------------------------------------------------------------------------------
! main body

    ! init mpi
    call MPI_Init(ierr)
    ! determine process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    params%rank = rank

    ! determine process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)
    params%number_procs = number_procs
    WABBIT_COMM         = MPI_COMM_WORLD
    ! output MPI status
    if (rank==0) then
        write(*,'(40("*"),A,40("*"))') "STARTING wabbit-post"
        write(*, '("MPI: using ", i5, " processes")') params%number_procs
#ifdef BLOCKINGSENDRECV
        write(*,'("MPI: code build with blocking send/recv in transfer (block_xfer_blocking.f90)")')
#else
        write(*,'("MPI: code build with NON-blocking send/recv in transfer (block_xfer_nonblocking.f90)")')
#endif
    end if

    !---------------------------------------------------------------------------
    ! Initialize parameters and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    call get_command_argument( 1, mode )
    if (rank==0) write(*,'("Starting postprocessing in ", a20, "mode")') mode

    select case(mode)
    case("--compute-rhs")
        call post_rhs(params)

    case("--mult-mask")
        call mult_mask(params)

    case("--mean")
        call post_mean(params)

    case("--sparse-to-dense")
        call sparse_to_dense(params)

    case("--vorticity", "--divergence", "--vor-abs")
        call compute_vorticity_post(params)

    case("--keyvalues")
        call get_command_argument(2,filename)
        call keyvalues(filename, params)

    case ("--compare-keys")
        if (rank == 0) then
            call get_command_argument(2,key1)
            call get_command_argument(3,key2)
            call compare_keys(key1,key2)
        end if

    case ("--flusi-to-wabbit")
        call flusi_to_wabbit(params)

    case default

        if (params%rank==0) then
            write(*,*) "Available Postprocessing tools are:"
            write(*,*) "--sparse-to-dense"
            write(*,*) "--mean"
            write(*,*) "--vorticity"
            write(*,*) "--keyvalues"
            write(*,*) "--compare-keys"
            write(*,*) "--flusi-to-wabbit"
            if (mode=="--h" .or. mode=="--help") then
                write(*,*) "To get more information about each postprocessing tool type: wabbit-post --[one of the listed tools] --help"
            else
                write(*,*) "Your postprocessing option is "// trim(adjustl(mode)) //", which I don't know"
            end if
        end if
    end select

    if (rank==0) then
        write(*,'(40("*"),A,40("*"))') "(regular) EXIT wabbit-post"
    endif

    ! end mpi
    call MPI_Finalize(ierr)

end program main_post
