!> \brief postprocessing routine for applying a filter to the data in a file
!-----------------------------------------------------------------------------------------------------

subroutine post_filter(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_time_step

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_in, file_out, filter_type
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs
    character(len=cshort)                   :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: nwork
    logical          :: help1, help2, help3

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    !-----------------------------------------------------------------------------------------------------
    call get_cmd_arg( "--help", help1, default=.false. )
    call get_cmd_arg( "--h", help2, default=.false. )
    call get_cmd_arg( "-h", help3, default=.false. )
    ! does the user need help?
    if (help1 .or. help2 .or. help3 .and. params%rank==0) then
        if (params%rank==0) then
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " Wabbit postprocessing: filter application"
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " Applies a filter operation to the data in the input files, and saves the result in a new file."
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " ./wabbit-post --filter --filter_type=explicit_3pt -i=source.h5 -o=target.h5 --g=2"
            write(*,'(A)') ""
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " --filter_type         - type of filter to apply, choose from:"
            ! write(*,'(A)') "    CDFXY_H - apply wavelet filters. X is interpolation order, Y is reconstruction order, H=HD for decomposition scaling filter, H=HR for reconstruction scaling filter, H=GD for decomposition wavelet filter, H=GR for reconstruction wavelet filter"
            write(*,'(A)') "    explicit_[POINT]pt          - superviscosity filter with [POINT] points, being 3,5,7,...,21. Needs g=[POINT-1/2] ghost points to work properly."
            write(*,'(A)') " -i                    - input file"
            write(*,'(A)') " -o                    - output file, defaults to filter-[input file]"
            write(*,'(A)') " --g                   - set amount of ghost points to use for the filter operation, as this is not automized"
            write(*,'(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    ! get values from command line (filename and level for interpolation)
    call get_cmd_arg("--filter_type", filter_type, default="explicit_3pt")
    call get_cmd_arg("-i", file_in, default="source.h5")
    call get_cmd_arg("-o", file_out, default="---")
    call get_cmd_arg("--g", g, default=-1)

    call check_file_exists(trim(file_in))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_in, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    if (filter_type(1:3) == "CDF") then
        ! we need to set everything up to "_" as wavelet
        params%wavelet = filter_type(1:index(filter_type,"_")-1)
    else
        ! CDF20 is sufficient
        params%wavelet = "CDF20"
    endif

    params%Jmax = tc_length
    params%n_eqn = 1
    params%domain_size(1:3) = domain(1:3)
    params%Bs = Bs

    allocate(params%filter_component(params%n_eqn))
    params%filter_component(:) = .true.

    call setup_wavelet(params, params%g)

    Bs = params%Bs
    g  = max(params%g, g)
    params%g = g

    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )

    nwork = 1

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=nwork)

    ! read input data
    call readHDF5vct_tree( (/file_in/), params, hvy_block, tree_ID)
    call sync_ghosts_tree( params, hvy_block, tree_ID)

    ! now apply filter
    if (filter_type(1:3) == "CDF") then
        ! ToDo? I actually don't know if this makes sense, we can only apply filter in one direction after all
    else
        ! just call filter wrapper
        params%filter_type = filter_type
        call filter_wrapper(time, params, hvy_block, tree_ID_flow)
    end if

    ! save filtered data
    if (file_out == "---") then
        file_out = "filter-"//trim(file_in)
    end if
    call saveHDF5_tree(file_out, time, iteration, 1, params, hvy_block, tree_ID )

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine post_filter
