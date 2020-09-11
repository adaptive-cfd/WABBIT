subroutine post_remesh(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_operators

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=80)      :: infile2, outfile, wavelet

    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration
    integer(kind=ik)       :: k
    integer(kind=ik)       :: tc_length, lgt_n_tmp

    integer(kind=ik) :: lgt_n
    integer(kind=ik) :: hvy_n
    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:, :)
    logical :: help1, help2


    ! call get_cmd_arg( "--help", help1, default=.false. )
    ! call get_cmd_arg( "-h", help2, default=.false. )
    !
    ! if ((help1 .or. help2) .and. (params%rank==0)) then
    !     write(*,*) "-----------------------------------------------------------"
    !     write(*,*) " Wabbit postprocessing: remeshing. Read in a file, apply thresholding with given eps."
    !     write(*,*) " Note: if the data is already sparse, setting an eps<eps_file makes not much sense."
    !     write(*,*) "-----------------------------------------------------------"
    !     write(*,*) "./wabbit-post --remesh --input=mask_0000.h5 --output=mask2_0000.h5 --eps=1.0e-4 --wavelet=CDF40"
    !     write(*,*) " --wavelet: CDF40, CDF44, CDF20, CDF22"
    !     write(*,*) " --memory=2.0gb"
    !     write(*,*) "-----------------------------------------------------------"
    ! end if
    !
    ! call get_cmd_arg( "--eps", params%eps, default=1.0_rk )
    ! call get_cmd_arg( "--input", infile2, default="none" )
    ! call get_cmd_arg( "--output", outfile, default="none" )
    ! call get_cmd_arg( "--wavelet", wavelet, default="CDF40" )
    !
    ! call check_file_exists(infile2)
    !
    !
    ! select case(wavelet)
    ! case ('CDF40')
    !     params%n_ghosts = 3_ik
    !     params%order_predictor = "multiresolution_4th"
    !     params%wavelet = "CDF40"
    !     params%wavelet_transform_type = "harten-multiresolution"
    ! case ('CDF44')
    !     params%n_ghosts = 6_ik
    !     params%order_predictor = "multiresolution_4th"
    !     params%wavelet = "CDF4,4"
    !     params%wavelet_transform_type = "biorthogonal"
    !
    ! case default
    !     call abort(20200906, "unknown wavelet")
    !
    ! end select
    !
    !
    ! params%max_treelevel = tc_length
    ! params%n_eqn = 1
    ! params%coarsening_indicator = "threshold-state-vector"
    ! params%threshold_state_vector_component = (/.true./)
    ! allocate(params%butcher_tableau(1,1))
    !
    ! call read_attributes(infile2, lgt_n_tmp, time, iteration, params%domain_size, &
    ! params%Bs, tc_length, params%dim)
    !
    ! call allocate_tree(params, lgt_block, hvy_block, hvy_neighbor, &
    ! lgt_active, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)
    !
    ! ! read mesh and field
    ! call read_mesh(infile2, params, lgt_n, hvy_n, lgt_block )
    ! call read_field(infile2, 1, params, hvy_block, hvy_n)
    !
    ! ! create lists of active blocks (light and heavy data)
    ! ! update list of sorted nunmerical treecodes, used for finding blocks
    ! call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    ! lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)
    !
    ! ! update neighbor relations
    ! call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    ! lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
    !
    ! call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    !
    ! call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,1), &
    ! lgt_n(1), lgt_sortednumlist(:,:,1), hvy_active(:,1), hvy_n(1), 1, params%coarsening_indicator, hvy_tmp )
    !
    ! call write_field(outfile, time, iteration, 1, params, lgt_block,&
    ! hvy_tmp, lgt_active, lgt_n, hvy_n, hvy_active )
end subroutine
