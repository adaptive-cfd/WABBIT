subroutine post_mean(params)
    use module_globals
    use module_params
    use module_operators
    use module_mesh
    use mpi
    use module_forestMetaData

    implicit none
    character(len=cshort)                   :: fname, fname_out, option         !> name of the file
    type (type_params), intent(inout)       :: params                           !> parameter struct

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1
    integer(kind=ik)                   :: rank, iteration, lgtID, hvyID, k
    real(kind=rk)                      :: time
    real(kind=rk)                      :: norms(1:5)
    logical                            :: exist, subtract_mean
    integer(kind=ik), dimension(3)     :: Bs

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
            write(*,*) ""
            write(*,*) "To subtract the mean from the file and overwrite it:"
            write(*,*) "mpi_command -n number_procs ./wabbit-post --mean filename.h5 --subtract"
        end if
        return
    endif
    
    ! Check for --subtract option
    subtract_mean = .false.
    call get_command_argument(3, option)
    if (trim(adjustl(option)) == '--subtract') then
        subtract_mean = .true.
    endif

    if (rank==0) write (*,*) "Computing spatial mean of file: "//trim(adjustl(fname))
    call check_file_exists( fname )

    ! add some parameters from the file
    call read_attributes(fname, lgt_n(tree_ID), time, iteration, params%domain_size, params%Bs, params%Jmax, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%n_eqn = 1
    params%g = 3_ik
    params%order_predictor = "multiresolution_4th"
    params%number_blocks = ceiling( real(lgt_n(tree_ID)) / real(params%number_procs) )
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1

    call allocate_forest(params, hvy_block)

    ! read input data, syncing data not needed
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID, synchronize_ghosts=.false.)

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
    endif
    
    ! If --subtract option is set, subtract the mean and overwrite the file
    if (subtract_mean) then
        if (rank == 0) write(*,'(A, es16.8)') "Subtracting mean value from file: ", norms(1)
        
        ! Loop over all active blocks and subtract the mean
        do k = 1, hvy_n(tree_ID)
            hvyID = hvy_active(k, tree_ID)
            call hvy2lgt(lgtID, hvyID, rank, params%number_blocks)
            Bs = params%Bs
            
            ! Subtract the mean from the block data (excluding ghost nodes, these are later synched anyways)
            if (params%dim == 2) then
                hvy_block(params%g+1:Bs(1)+params%g, params%g+1:Bs(2)+params%g, 1, 1, hvyID) = &
                hvy_block(params%g+1:Bs(1)+params%g, params%g+1:Bs(2)+params%g, 1, 1, hvyID) - norms(1)
            else
                ! 3D case
                hvy_block(params%g+1:Bs(1)+params%g, params%g+1:Bs(2)+params%g, params%g+1:Bs(3)+params%g, 1, hvyID) = &
                hvy_block(params%g+1:Bs(1)+params%g, params%g+1:Bs(2)+params%g, params%g+1:Bs(3)+params%g, 1, hvyID) - norms(1)
            endif
        enddo

        ! for adaptive grids we need wavelet setup for synching for writing
        params%wavelet = "CDF40"
        call setup_wavelet(params)
        
        ! Overwrite the original file
        if (rank == 0) write(*,'(A)') "Overwriting file with mean-subtracted data: "//trim(adjustl(fname))
        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_block, tree_ID)
        
        if (rank == 0) write(*,'(A)') "File successfully overwritten."
    else
        ! Original behavior: write norms to disk
        if (rank == 0) then
            ! write volume integral to disk
            call get_command_argument(3,fname_out)
            inquire ( file=fname_out, exist=exist )
            if (exist) then
                open(14,file=fname_out, status='replace')
                write(14,'(A)') "# Mean, L1 norm, L2 norm, LInfty norm, H1 norm" 
                write(14,'(4(es16.8, " "), es16.8)') norms(1), norms(2), norms(3), norms(4), norms(5)
                close(14)
            else
                write(*,'(A)') "Output file "//trim(adjustl(fname_out))//" does not exist. Not writing output. Use e.g. touch"
            endif
        endif
    endif
end subroutine post_mean
