subroutine post_mean(params)
    use module_globals
    use module_params
    use module_operators
    use module_mesh
    use mpi
    use module_forestMetaData

    implicit none
    character(len=cshort)                   :: fname, fname_out                 !> name of the file
    type (type_params), intent(inout)       :: params                           !> parameter struct

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1
    integer(kind=ik)                   :: rank, iteration
    real(kind=rk)                      :: time
    real(kind=rk)                      :: norms(1:5)
    logical                            :: exist

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    !-----------------------------------------------------------------------------------------------------
    rank = params%rank
    !-----------------------------------------------------------------------------------------------------
    call get_command_argument(2,fname)
    if (fname =='--help' .or. fname=='--h') then
        if (rank==0) then
            write(*,*) "WABBIT postprocessing: compute mean value of field. We output all values to a file."
            write(*,*) "mpi_command -n number_procs ./wabbit-post --mean filename.h5 result.txt"
        end if
        return
    endif

    if (rank==0) write (*,*) "Computing spatial mean of file: "//trim(adjustl(fname))
    call check_file_exists( fname )

    ! add some parameters from the file
    call read_attributes(fname, lgt_n(tree_ID), time, iteration, params%domain_size, params%Bs, params%Jmax, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%n_eqn = 1
    params%g = 2_ik
    params%order_predictor = "multiresolution_2nd"
    params%number_blocks = ceiling( real(lgt_n(tree_ID)) / real(params%number_procs) )
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1

    call allocate_forest(params, hvy_block)

    ! read input data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID)

    ! use functions from module_operators
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "Mean", norms(1:1), threshold_state_vector=.false.)
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L1", norms(2:2), threshold_state_vector=.false.)
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norms(3:3), threshold_state_vector=.false.)
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "Linfty", norms(4:4), threshold_state_vector=.false.)
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "H1", norms(5:5), threshold_state_vector=.false.)

    if (rank == 0) then
        write(*,'(A, es16.8)') "Mean:        ", norms(1)
        write(*,'(A, es16.8)') "L1 norm:     ", norms(2)
        write(*,'(A, es16.8)') "L2 norm:     ", norms(3)
        write(*,'(A, es16.8)') "LInfty norm: ", norms(4)
        write(*,'(A, es16.8)') "H1 norm:     ", norms(5)

        ! write volume integral to disk
        call get_command_argument(3,fname_out)
        inquire ( file=fname_out, exist=exist )
        if (exist) then
            open(14,file=fname_out, status='replace')
            write(14,'(A)') "# Mean, L1 norm, L2 norm, LInfty norm, H1 norm" 
            write(14,'(5(es16.8, ", "))') norms(1), norms(2), norms(3), norms(4), norms(5)
            close(14)
        else
            write(*,'(A)') "Output file "//trim(adjustl(fname_out))//" does not exist. Not writing output. Use e.g. touch"
        endif
    endif
end subroutine post_mean
