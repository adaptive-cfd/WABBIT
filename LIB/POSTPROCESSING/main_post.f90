!> \brief main postprocessing program. get command argument to decide which postprocessing mode to run
! ********************************************************************************************

program main_post
    use mpi
    use module_params
    use module_MOR, only : post_POD, post_reconstruct, post_PODerror, post_timecoef_POD
    use module_timing

    implicit none

    integer(kind=ik)                    :: ierr                   ! MPI error variable
    integer(kind=ik)                    :: rank                   ! process rank
    integer(kind=ik)                    :: number_procs           ! number of processes
    type (type_params)                  :: params
    character(len=cshort)               :: mode
    character(len=clong)                :: filename, key1, key2
    real(kind=rk)                       :: elapsed_time

    call MPI_Init(ierr)                                           ! init mpi
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)                ! determine process rank
    params%rank = rank

    ! determine process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)
    params%number_procs = number_procs
    WABBIT_COMM         = MPI_COMM_WORLD
    elapsed_time        = MPI_wtime()

    ! output MPI status
    if (rank==0) then
        write(*,'(20("─"), A26, 34("─"))') "   STARTING wabbit-post   "
        write(*,'("MPI: using ", i5, " processes")') params%number_procs
        write(*,'("MPI: code build with NON-blocking send/recv in transfer (block_xfer_nonblocking.f90)")')
    end if

    ! this is used to document a bit (one often forgets to write down the params in the command line call)
    call print_command_line_arguments()

    !---------------------------------------------------------------------------
    ! Initialize parameters and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    call get_command_argument( 1, mode )
    if (rank==0) write(*,'("Starting postprocessing in ", a20, " mode")') mode

    select case(mode)
    case ("--extract-slice")
        call post_extract_slice(params)

    case ("--evaluate-wavelet-thresholding")
        call post_evaluate_thresholding(params)

    case ("--wavelet-decompose", "--wavelet-reconstruct")
        call post_wavelet_transform(params)

    case ("--dump-neighbors")
        call post_dump_neighbors(params)

    case ("--OP")
        call operator_reconstruction(params)

    case ("--OP-rhs")
        call rhs_operator_reconstruction(params)

    case ("--prune-tree")
        call post_prune_tree(params)

    case ("--superstl")
        call post_superstl(params)

    case ("--add-two-masks", "--add", "--subtract", "--multiply", "--test_operations", "--grid1-to-grid2", "--noise-like-grid1")
        call post_add_two_masks(params)

    case ("--stl2dist")
        call post_stl2dist(params)

    case("--compute-rhs")
        call post_rhs(params)

    case("--mult-mask", "--mult-mask-direct", "--mult-mask-inverse")
        call mult_mask(params)

        !mean of a given field sum(q_ij)/size(q_ij) the result is a scalar
    case("--mean")
        call post_mean(params)

        ! average of multiple snapshots the result is the averaged snapshot
    case("--average")
        call post_average_snapshots(params)

    case("--sparse-to-dense", "--refine-everywhere", "--refine-everywhere-forced", "--coarsen-everywhere")
        call sparse_to_dense(params)

    case("--refine-coarsen-test", "--ghost-nodes-test","--wavelet-decomposition-unit-test", "--wavelet-decomposition-invertibility-test", "--sync-test", "--treecode-test")
        call post_unit_test(params)

    case ("--performance-test")
        call performance_test(params)

    case ("--adaption-test")
        call adaption_test(params)

    case ("--compression-unit-test")
        call post_compression_unit_test(params)

    case("--dense-to-sparse")
        call dense_to_sparse(params)

    case("--dry-run")
        call post_dry_run()

    case("--vorticity", "--divergence", "--vor-abs", "--Q", "--copy")
        call compute_vorticity_post(params)

    case("--derivative")
        call post_derivative(params)

    case("--gradient")
        call compute_scalar_field_post(params)

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

    case ("--POD")
        call post_POD(params)

    case ("--filter")
        call post_filtertest(params)

    case ("--POD-reconstruct")
        call post_reconstruct(params)

    case ("--POD-error")
        call post_PODerror(params)

    case ("--POD-time")
        call post_timecoef_POD(params)

    case ("--generate_forest")
        call post_generate_forest(params)

    case ("--denoise")
        call post_denoising(params)

    case ("--denoising-test")
        call post_denoising_test(params)

    case ("--cvs-invertibility-test")
        call post_cvs_invertibility_test(params)

    case ("--proto-GS-multigrid")
        call proto_GS_multigrid(params)
        
    case default

        if (params%rank==0) then
            write(*, '(A)') "Available Postprocessing tools are:"
            write(*, '(A)') "--sparse-to-dense"
            write(*, '(A)') "--dense-to-sparse"
            write(*, '(A)') "--mean"
            write(*, '(A)') "--vorticity"
            write(*, '(A)') "--vor-abs"
            write(*, '(A)') "--divergence"
            write(*, '(A)') "--Q"
            write(*, '(A)') "--OP-rhs"
            write(*, '(A)') "--OP"
            write(*, '(A)') "--gradient"
            write(*, '(A)') "--keyvalues"
            write(*, '(A)') "--dry-run"
            write(*, '(A)') "--compare-keys"
            write(*, '(A)') "--flusi-to-wabbit"
            write(*, '(A)') "--POD"
            write(*, '(A)') "--POD-reconstruct"
            write(*, '(A)') "--POD-error"
            write(*, '(A)') "--POD-time"
            write(*, '(A)') "--stl2dist"
            write(*, '(A)') "--add-two-masks"
            write(*, '(A)') "--mult-mask"
            write(*, '(A)') "--mult-mask-direct"
            write(*, '(A)') "--mult-mask-inverse"
            write(*, '(A)') "--post_rhs"
            write(*, '(A)') "--average"
            write(*, '(A)') "--generate_forest"
            write(*, '(A)') "--evaluate-wavelet-thresholding"
            write(*, '(A)') "--refine-everywhere"
            write(*, '(A)') "--denoise"
            ! tests
            write(*, '(A)') "--compression-unit-test"
            write(*, '(A)') "--performance_test"
            write(*, '(A)') "--adaption-test"
            write(*, '(A)') "--refine-coarsen-test"
            write(*, '(A)') "--ghost-nodes-test"
            write(*, '(A)') "--wavelet-decomposition-unit-test"
            write(*, '(A)') "--wavelet-decomposition-invertibility-test"
            write(*, '(A)') "--sync-test"
            write(*, '(A)') "--treecode-test"

            if (mode=="--h" .or. mode=="--help") then
                write(*, '(A)') "To get more information about each postprocessing tool type: wabbit-post --[one of the listed tools] --help"
            else
                write(*, '(A)') "Your postprocessing option is "// trim(adjustl(mode)) //", which I don't know"
            end if
        end if
    end select

    ! close and flush all existings *.t files
    call close_all_t_files()

    elapsed_time = MPI_wtime() - elapsed_time
    if (rank==0) then
        write(*, '(A)')
        write(*,'("Elapsed time:", f16.4, " s")') elapsed_time
        write(*,'(20("─"),A,28("─"))') "   (regular) EXIT wabbit-post   "
    endif
    ! MPI Barrier before program ends
    call MPI_Barrier(WABBIT_COMM, ierr)

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    if ( mode(:5) == "--POD") call summarize_profiling( WABBIT_COMM )

    ! end mpi
    call MPI_Finalize(ierr)

end program main_post