!> \file
! WABBIT
!> \name sparse_to_dense.f90
!> \version 0.5
!> \author sm
!
!> \brief postprocessing routine for interpolation of a given field to the desired level
!
! = log ======================================================================================
!> \date  31/01/18 - create hashcode: commit 13cb3d25ab12e20cb38e5b87b9a1e27a8fe387e8
!-----------------------------------------------------------------------------------------------------

subroutine sparse_to_dense(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_globals

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: file_in
    character(len=80)      :: file_out
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration

    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(kind=ik)                        :: hvy_n, lgt_n, max_neighbors, level, k, tc_length
    integer(kind=ik), dimension(3)          :: Bs
    integer(hid_t)                          :: file_id
    character(len=2)                        :: level_in, order
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: treecode_size, number_dense_blocks
!-----------------------------------------------------------------------------------------------------

    call get_command_argument(2, file_in)
    call get_command_argument(3, file_out)

    if (file_in == '--help' .or. file_in == '--h') then
        if ( params%rank==0 ) then
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "                SPARSE to DENSE "
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "postprocessing subroutine to refine/coarse mesh to a uniform"
            write(*,*) "grid (up and downsampling ensured)."
            write(*,*) "Command:"
            write(*,*) "mpi_command -n number_procs ./wabbit-post --sparse-to-dense "
            write(*,*) "source.h5 target.h5 [target_treelevel order-predictor(2 or 4)]"
            write(*,*) "-------------------------------------------------------------"
            write(*,*) "Optional Inputs: "
            write(*,*) "  1. target_treelevel = number specifying the desired treelevel"
            write(*,*) "  (default is the max treelevel of the source file) " 
            write(*,*) "  2. order-predictor = consistency order or the predictor stencil"
            write(*,*) "  (default is preditor order 4) " 
            write(*,*)
            write(*,*)
        end if
        return
    end if

    ! get values from command line (filename and level for interpolation)
    call check_file_exists(trim(file_in))
    call read_attributes(file_in, lgt_n, time, iteration, &
                         domain, Bs,tc_length, params%dim)

    if (len_trim(file_out)==0) then
      call abort(0909191,"You must specify a name for the target! See --sparse-to-dense --help")
    endif
    
    ! check if optional arguments are specified
    if (command_argument_count()<4) then
      ! set defaults
      order = "2"
      level = tc_length
    else
      call get_command_argument(4, level_in)
      read(level_in,*) level
      call get_command_argument(5, order)
    endif

    if (order == "4") then
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik
    elseif (order == "2") then
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik
    else
        call abort(392,"ERROR: chosen predictor order invalid or not (yet) implemented. choose between 4 (multiresolution_4th) and 2 (multiresolution_2nd)")
    end if

    ! in postprocessing, it is important to be sure that the parameter struct is correctly filled:
    ! most variables are unfortunately not automatically set to reasonable values. In simulations,
    ! the ini files parser takes care of that (by the passed default arguments). But in postprocessing
    ! we do not read an ini file, so defaults may not be set.
    allocate(params%butcher_tableau(1,1))
    ! we read only one datafield in this routine
    params%n_eqn  = 1
    params%block_distribution="sfc_hilbert"

     if (params%dim==3) then
        ! how many blocks do we need for the desired level?
        number_dense_blocks = 8_ik**level
        max_neighbors = 74
    else
        number_dense_blocks = 4_ik**level
        max_neighbors = 12
    end if

    if (params%rank==0) then
        write(*,'(80("-"))')
        write(*,*) "Wabbit sparse-to-dense. Will read a wabbit field and return a"
        write(*,*) "full grid with all blocks at the chosen level."
        write(*,'(A20,1x,A80)') "Reading file:", file_in
        write(*,'(A20,1x,A80)') "Writing to file:", file_out
        write(*,'(A20,1x,A80)') "Predictor used:", params%order_predictor
        write(*,'(A20,1x,i3," => ",i9," Blocks")') "Target level:", level, number_dense_blocks
        write(*,'(80("-"))')
    endif

    ! set max_treelevel for allocation of hvy_block
    params%max_treelevel = max(level, tc_length)
    params%min_treelevel = level
    params%Bs = Bs
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)

    ! is lgt_n > number_dense_blocks (downsampling)? if true, allocate lgt_n blocks
    !> \todo change that for 3d case
    params%number_blocks = ceiling( 4.0*dble(max(lgt_n, number_dense_blocks)) / dble(params%number_procs) )

    if (params%rank==0) then
        write(*,'("Data dimension: ",i1,"D")') params%dim
        write(*,'("File contains Nb=",i6," blocks of size Bs=",i4," x ",i4," x ",i4)') lgt_n, Bs(1),Bs(2),Bs(3)
        write(*,'("Domain size is ",3(g12.4,1x))') domain
        write(*,'("Time=",g12.4," it=",i9)') time, iteration
        write(*,'("Length of treecodes in file=",i3," in memory=",i3)') tc_length, params%max_treelevel
        write(*,'("   NCPU=",i6)') params%number_procs
        write(*,'("File   Nb=",i6," blocks")') lgt_n
        write(*,'("Memory Nb=",i6)') params%number_blocks
        write(*,'("Dense  Nb=",i6)') number_dense_blocks
    endif

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
        hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)

    ! read field
    call read_mesh(file_in, params, lgt_n, hvy_n, lgt_block)
    call read_field(file_in, 1, params, hvy_block, hvy_n)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n, tree_ID=1)

    ! balance the load
    call balance_load(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID=1)

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! refine/coarse to attain desired level, respectively
    !coarsen
    do while (max_active_level( lgt_block, lgt_active, lgt_n )>level)
        ! check where coarsening is actually needed and set refinement status to -1 (coarsen)
        do k = 1, lgt_n
            if (treecode_size(lgt_block(lgt_active(k),:), params%max_treelevel) > level)&
                lgt_block(lgt_active(k), params%max_treelevel + IDX_REFINE_STS) = -1
        end do
        ! this might not be necessary since we start from an admissible grid
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n )

        call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist, &
        hvy_active, hvy_n, tree_ID=1)

        call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, tree_ID=1)
    end do
    ! refine
    do while (min_active_level( lgt_block, lgt_active, lgt_n )<level)
        ! check where refinement is actually needed
        do k = 1, lgt_n
            if (treecode_size(lgt_block(lgt_active(k),:), params%max_treelevel) < level)&
                lgt_block(lgt_active(k), params%max_treelevel + IDX_REFINE_STS) = 1
        end do
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n )
        if ( params%dim == 3 ) then
            ! 3D:
            call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
        else
            ! 2D:
            call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:),&
                hvy_active, hvy_n )
        end if

        call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, tree_ID=1)

        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    end do

    call balance_load( params, lgt_block, hvy_block, &
        hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID=1 )

    call write_field(file_out, time, iteration, 1, params, lgt_block, &
        hvy_block, lgt_active, lgt_n, hvy_n, hvy_active)

    if (params%rank==0 ) then
        write(*,'("Wrote data of input-file: ",A," now on uniform grid (level",i3, ") to file: ",A)') &
            trim(adjustl(file_in)), level, trim(adjustl(file_out))
         write(*,'("Minlevel:", i3," Maxlevel:", i3, " (should be identical now)")') &
             min_active_level( lgt_block, lgt_active, lgt_n ),&
             max_active_level( lgt_block, lgt_active, lgt_n )
    end if

    call deallocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
        hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp)
end subroutine sparse_to_dense
