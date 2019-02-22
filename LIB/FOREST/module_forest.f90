
!-----------------------------------------------------------------
!> \file
!> \brief
!! Module for implementing multiple trees in a forest
!> \version 10.1.2019
!> \author P.Krah
!-----------------------------------------------------------------

module module_forest

  ! modules
  use mpi
  use module_precision
  use module_ini_files_parser_mpi, only : read_param_mpi
  use module_params
  use module_mesh
  use module_IO

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: read_field2tree, read_forest_params, write_tree_field, &
            add_tree, count_tree_hvy_n
  !**********************************************************************************************
contains






  !##############################################################
  !> read in all parameters needed for initializing trees in
  !> the forest
  subroutine read_forest_params( params, filename )

      implicit none
      !-----------------------------------------------
      type (type_params), intent(inout)  :: params
      character(len=*), intent(in)       :: filename
      !-----------------------------------------------
      type(inifile) :: FILE
      ! read the file, only process 0 should create output on screen
      call set_lattice_spacing_mpi(1.0d0)
      call read_ini_file_mpi(FILE, filename, .true.)
      call ini_MPI(params, FILE ) ! read MPI information for bridge interface
      call ini_discretization( params, FILE ) ! reads order predictor etc.

      call read_param_mpi(FILE, 'Blocks', 'max_forest_size', params%forest_size, 1 )
      call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%n_ghosts, 1 )
      call read_param_mpi(FILE, 'Blocks', 'number_blocks', params%number_blocks, 1 )
      call read_param_mpi(FILE, 'Blocks', 'block_dist', params%block_distribution , "equal" )
      ! check parameter if needed
      if (params%number_blocks < 1) call abort(11219, "number of blocks is to small")
  end subroutine read_forest_params
  !##############################################################













  !##############################################################
  !> This function reads a set of data into a specified tree.
  !> The set of data can involve multiple files. Each file should
  !> contain a quantity at the same snapshot (i.e. same iteration/time).
  !> If no data has been read in before. This function will allocate
  !> the heavy data for you.
  subroutine read_field2tree(params, fnames, N_files, tree_id, tree_n, &
                            lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                            hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(inout) :: params           !< params structure
      integer(kind=ik), intent(in)      :: N_files     !< number of fields/quantities to read
      character(len=*), intent(in)      :: fnames(N_files)  !< file names
      integer(kind=ik), intent(inout)   :: hvy_n, tree_n !< number of heavy and light active blocks
      integer(kind=ik), intent(in)      :: tree_id     !< number of the tree
      integer(kind=ik), allocatable, intent(inout)   :: lgt_n(:)  !< number of trees and active blocks
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_active(:,:), hvy_active(:) !< active lists
      integer(kind=tsize), ALLOCATABLE, intent(inout):: lgt_sortednumlist(:,:)
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
      !-------------------------------- ---------------------------------
      integer(kind=ik)  :: iteration, dF, tc_length, dim, i, lgt_n_tmp,  &
                           rank, level, Bs, fsize
      real(kind=rk)     :: time, domain(3)
      
      !N =jparams%number_blocks
      ! set MPI parameter
      rank  = params%rank
      fsize = params%forest_size
      ! grid parameter
      call read_attributes(fnames(1), lgt_n_tmp, time, iteration, domain, Bs, tc_length, dim)
      ! Check if one tree already exists in the forest:
      ! If it doesnt initialice some improtant parameters like the Block size Bs
      ! and the spatial dimension of the data.
      ! If it does check if the parameters are suitable for the forest
      if (allocated(hvy_block)) then
          if ( params%Bs.ne.Bs) call abort(100119,"Block size should not change);!")
          if ( params%max_treelevel < tc_length) call abort(10119,"maximal treelevel incompatible")
      else
          params%max_treelevel = tc_length
          params%dim = dim
          params%threeD_case = (dim == 3)
          params%n_eqn = N_files
          params%N_fields_saved = N_files
          params%time_step_method="no" ! is necessary for
          params%Bs = Bs
          params%number_blocks = calculate_max_nr_blocks ( params )
          params%domain_size = domain
          ! we have to allocate grid if this routine is called for the first time
          call allocate_hvy_lgt_data(params, lgt_block, hvy_block, hvy_neighbor, &
                        lgt_active, lgt_n, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)
          call reset_lgt_data(lgt_block, lgt_active(:, fsize+1), &
                            params%max_treelevel, lgt_n(fsize+1), lgt_sortednumlist)
          hvy_neighbor = -1
          lgt_n = 0 ! reset number of acitve light blocks
          tree_n= 0 ! reset number of trees in forest
      endif
      ! read treecode from first input file
      call read_tree(fnames, N_files, params, lgt_n_tmp, hvy_n, lgt_block, &
                        hvy_block, hvy_tmp, tree_id)
     
      ! * We save the number of all active lgt_ids in lgt_n_tmp, which is actually
      ! not needed here, since
      ! create_active_and_sorted_lists is calculating sum(lgt_n) again.
      ! However lgt_n_tmp is used to double check the calculation of the lgt_n's!
      lgt_n_tmp = lgt_n( fsize + 1) + lgt_n_tmp

      call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. , tree_n)
      
      ! check the calculation of create_active and sorted list
      if (lgt_n_tmp .ne. lgt_n( fsize + 1 )) call abort(132191,"There is something wrong with the light data")
     
      call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,fsize + 1),&
                            lgt_n(fsize + 1), lgt_sortednumlist, hvy_active, hvy_n )

      call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

      !> USED FOR DEBUGGING \TODO : REMOVE when you read this. thank you(:
      !write(*,*)"sum over all tree_n=", tree_n
      !do i =1,2*lgt_n( fsize + 1)
         !write(*,'("lgt_data=",7(i9))') lgt_block(i,1:params%max_treelevel+idx_tree_id) 
      !end do          


      !write(*,*)"all lgt data for tree_id", tree_n
      !do i =1,lgt_n( tree_id )
         !write(*,'("lgt_data=",7(i9))') lgt_block(lgt_active(i,tree_id),1:params%max_treelevel+idx_tree_id) 
      !end do          

  end subroutine read_field2tree
  !##############################################################


  !##############################################################
  !> This function takes to trees and refines their treestructures
  !> until they are refined to the same max level. 
  !> The resulting trees are identical in treestructure. Hence
  !> they have the same treecodes and number of blocks.
  !> Remark: You may want to balance the load after using this routine.
  subroutine refine_trees2same_lvl(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, & 
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2)

      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(inout) :: params   !< params structure
      integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
      integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
      integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
      integer(kind=ik), allocatable, intent(inout)   :: lgt_n(:) !< number of light active blocks
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
      integer(kind=tsize), ALLOCATABLE, intent(inout):: lgt_sortednumlist(:,:)
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
      !-----------------------------------------------------------------
      integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
      integer(kind=ik)    :: k1, k2, Nblocks_refined, level_min 
      integer(kind=tsize) :: treecode1, treecode2 
      
      Jmax = params%max_treelevel ! max treelevel
      fsize= params%forest_size   ! maximal number of trees in forest
      
      write(*,'("Refining trees to same level: ",i9,",",x,i9)') tree_id1, tree_id2
      ! The Trees can only be added when their grids are identical. At present we
      ! try to keep the finest levels of all trees. This means we refine
      ! all blocks which are not on the same level.
      do while( .true. ) ! loop until both trees have the same grid structure
                         ! meaning that no block has to be refined anymore
          Nblocks_refined = 0
          ! Loop over the active light data of both trees and compare the meshlevels
          do k1 = 1, lgt_n(tree_id1)
             !--------
             ! TREE 1
             !--------
             lgt_id1  = lgt_active(k1, tree_id1)
             level1   = lgt_block(lgt_id1, Jmax + idx_mesh_lvl)
             
             do k2 = 1, lgt_n(tree_id2)
             !--------
             ! TREE 2
             !--------
                lgt_id2  = lgt_active(k2, tree_id2)
                ! Take this block only if it is not already marked for refinement
                if ( lgt_block(lgt_id2, Jmax + idx_refine_sts) .ne. 0 ) cycle
                level2   = lgt_block(lgt_id2, Jmax + idx_mesh_lvl)
                
                ! The treecodes can only be compared if they have the same size
                ! (i.e. treecode length). Therefore we have to find the minimum of both 
                ! max levels.
                level_min= min(level1, level2)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level_min))
                treecode2= treecode2int(lgt_block(lgt_id2, 1 : level_min))
                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                ! Comparison of the trees:
                !     * check if treecodes on the trees are identical up to level1 (level of tree1)
                !     * yes:
                !             - if they are not on the same level
                !               refine the coarser block on the corresponding tree 
                !     * no: move to next treecode and compare again 
                if (treecode1 == treecode2 .and. level1 .ne. level2 )then
                    Nblocks_refined = Nblocks_refined + 1
                    if (level1 > level2) then
                         ! mark the block on tree 2 for refinement 
                         lgt_block(lgt_id2, Jmax + idx_refine_sts) = +1
                    else
                         ! mark the block on tree 1 for refinement 
                         lgt_block(lgt_id1, Jmax + idx_refine_sts) = +1
                    end if
                    !###########################
                    exit  ! exit the inner loop
                    !###########################
                end if
                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             end do ! loop over tree2
          end do ! loop over tree1

          ! Decide if trees have the same treestructure or not (i.e. exit loop or refine)
          if (Nblocks_refined == 0) then
              exit   ! EXIT the loop when nothing has to be refined anymore
          else
              !----------------------------
              ! refine the marked blocks
              !----------------------------
              write(*,'("Number of blocks marked for refinement: ",i9)') Nblocks_refined 
              ! 1) check gradedness of the grid (meshlevel of ajacent blocks should not differe by more than 1
              call ensure_gradedness( params, lgt_block, hvy_neighbor, &
              lgt_active(:, fsize+1), lgt_n(fsize+1), lgt_sortednumlist, hvy_active, hvy_n )
              ! 2) refine blocks
              if ( params%threeD_case ) then
                ! 3D:
                call refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
              else
                ! 2D:
                call refinement_execute_2D( params, lgt_block, hvy_block(:,:,1,:,:),&
                hvy_active, hvy_n )
              end if
              ! since lgt_block was synced we have to create the active lists again 
              call create_active_and_sorted_lists( params, lgt_block, lgt_active,&
              lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true., tree_n )
              ! update neighbor relations
              call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:, fsize+1), &
              lgt_n(fsize+1), lgt_sortednumlist, hvy_active, hvy_n )
              ! synchronice 
              call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
          endif    
      end do
  end subroutine 
  !##############################################################




  !##############################################################
  subroutine add_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, & 
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2)

      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(inout) :: params   !< params structure
      integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
      integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
      integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
      integer(kind=ik), allocatable, intent(inout)   :: lgt_n(:) !< number of light active blocks
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
      integer(kind=tsize), ALLOCATABLE, intent(inout):: lgt_sortednumlist(:,:)
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
      !-----------------------------------------------------------------
      integer(kind=ik)    :: my_rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
      integer(kind=ik)    :: k1, k2, N, level_min, rank1, rank2 
      integer(kind=tsize) :: treecode1, treecode2 
      
      Jmax = params%max_treelevel ! max treelevel
      fsize= params%forest_size   ! maximal number of trees in forest
      N    = params%number_blocks ! number of blocks per rank
      my_rank = params%rank  ! proc rank
      write(*,'("Adding trees: ",i9,",",i9)') tree_id1, tree_id2
      ! The Trees can only be added when their grids are identical. At present we
      ! try to keep the finest levels of all trees. This means we refine
      ! all blocks which are not on the same level.
      
      call refine_trees2same_lvl(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, & 
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2)

         return
         
      !call balance_load( params, lgt_block, hvy_block,  hvy_neighbor, &
          !lgt_active(:, fsize + 1), lgt_n(fsize + 1), lgt_sortednumlist, hvy_active, hvy_n, hvy_tmp )
         
          do k1 = 1, lgt_n(tree_id1)
             !--------
             ! TREE 1
             !--------
             lgt_id1  = lgt_active(k1, tree_id1)
             level1   = lgt_block(lgt_id1, Jmax + idx_mesh_lvl)
             treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))

             do k2 = 1, lgt_n(tree_id2)
             !--------
             ! TREE 2
             !--------
                lgt_id2  = lgt_active(k2, tree_id2)
                level2   = lgt_block(lgt_id2, Jmax + idx_mesh_lvl)
                treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                if (treecode1 == treecode2 )then
                   call lgt_id_to_proc_rank( rank1, lgt_id1, N)
                   call lgt_id_to_proc_rank( rank2, lgt_id2, N)
                end if
                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             end do ! loop over tree2
          end do ! loop over tree1

    do k1 = 1, hvy_n
       call hvy_id_to_lgt_id( lgt_id1, hvy_active(k1), my_rank, N) 

    end do
  end subroutine add_tree
  !##############################################################








  !##############################################################
  subroutine read_tree(fnames, N_files, params, lgt_n, hvy_n, lgt_block, &
                        hvy_block, hvy_tmp, tree_id_optional)

    implicit none

      !-------------------------------- ---------------------------------
      integer(kind=ik), intent(in)   :: N_files!< number of files
      character(len=*), intent(in)      :: fnames(N_files) !< name of file
      type (type_params), intent(in)    :: params !< user defined parameter structure
      integer(kind=ik), intent(inout)   :: hvy_n, lgt_n !< number of active blocks (heavy and light data)
      integer(kind=ik), intent(inout)   :: lgt_block(:,:) !< light data array
      integer(kind=ik), optional, intent(in)   :: tree_id_optional !< index of the tree you want to save the field in
      real(kind=rk), intent(inout):: hvy_block(:, :, :, :, :) !< heavy data array - block data
      real(kind=rk), intent(inout):: hvy_tmp(:, :, :, :, :) !< heavy data array - block data
      !-------------------------------- ---------------------------------
      integer(kind=ik) :: k, N, rank, number_procs, ierr, treecode_size, tree_id, Bs, g, dF
      integer(kind=ik) :: ubounds(2), lbounds(2), blocks_per_rank_list(0:params%number_procs-1)
      integer(kind=ik) :: free_hvy_id, free_lgt_id, free_hvy_n 
      integer(hsize_t)      :: dims_treecode(2)
      integer(kind=ik), dimension(:,:), allocatable :: block_treecode
      integer(hid_t)        :: file_id
      character(len = 80) :: fname
      if (present(tree_id_optional)) then
        tree_id=tree_id_optional
      else
        tree_id=1
      endif
      ! set MPI parameters
      rank         = params%rank
      number_procs = params%number_procs
      N            = params%number_blocks
      fname        = fnames(1) ! take the first of the N_files to read the lgt_data from
      Bs = params%Bs
      g  = params%n_ghosts
      
      ! open the file
      call check_file_exists(fname)
      call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

      if ( rank == 0 ) then
        write(*,'(80("_"))')
        write(*,'(A)') "READING: initializing grid from file..."
        write(*,'("Filename: ",A)') trim(adjustl(fname))
        write(*,'("Tree_id: ",i9)') tree_id
        ! NOTE: we pass the routine lgt_n, which must be previsouly read from the
        ! file. This is done using read_attributes. This can possibly be merged here
        ! and mustnt be done in the caller
        write(*,'("Expected Nblocks=",i6," (on all cpus)")') lgt_n
      end if
      ! Nblocks per CPU
      ! this list contains (on each mpirank) the number of blocks for each mpirank. note
      ! zero indexing as required by MPI
      ! set list to the average value
      blocks_per_rank_list(:) = lgt_n / number_procs

      ! as this does not necessarily work out, distribute remaining blocks on the first CPUs
      if (mod(lgt_n, number_procs) > 0) then
        blocks_per_rank_list(0:mod(lgt_n, number_procs)-1) = &
        blocks_per_rank_list(0:mod(lgt_n, number_procs)-1) + 1
      end if
      ! some error control -> did we loose blocks? should never happen.
      if ( sum(blocks_per_rank_list) /= lgt_n) then
        call abort(1028,"ERROR: while reading from file, we seem to have gained/lost some blocks during distribution...")
      end if

      ! number of active blocks on my process
      free_hvy_n = blocks_per_rank_list(rank)

      ! check what Jmax was saved in file (check dimensions of treecode in file)
      call get_size_datafield(2, file_id, "block_treecode", dims_treecode)

      ! compare dimensions
      if (dims_treecode(1)>params%max_treelevel) then
        ! treecode in input file is greater than the new one, abort and output on screen
        ! NOTE this can be made working if not all levels in the file are actually used (e.g. level_max=17
        ! but active level=4). On the other han11219d, that appears to be rare.
        call abort(73947887, "ERROR: Treecode in file is longer than what is set in INI file.")
      end if

      allocate(block_treecode(1:dims_treecode(1), 1:free_hvy_n))
      block_treecode = -1

      ! tell the hdf5 wrapper what part of the global [ n_active x max_treelevel + idx_refine_sts]
      ! array we want to hold, so that all CPU can read from the same file simultaneously
      ! (note zero-based offset):
      lbounds = (/0, sum(blocks_per_rank_list(0:rank-1))/)
      ubounds = (/int(dims_treecode(1),4)-1, lbounds(2) + free_hvy_n - 1/)
      call read_dset_mpi_hdf5_2D(file_id, "block_treecode", lbounds, ubounds, block_treecode)

      ! close file and HDF5 library
      call close_file_hdf5(file_id)

      ! read datafields from files into hvy_tmp
      do dF = 1, N_files
        ! read data from file
        call check_file_exists(trim(fnames(dF)))
        call read_field(fnames(dF), dF, params, hvy_tmp, free_hvy_n )
      end do

      do k = 1, free_hvy_n
        call get_free_local_light_id( params, rank, lgt_block, free_lgt_id)
        
        call lgt_id_to_hvy_id( free_hvy_id, free_lgt_id, rank, N )
        ! copy treecode
        lgt_block(free_lgt_id, 1:dims_treecode(1)) = block_treecode(1:dims_treecode(1), k)
        ! set mesh level
        lgt_block(free_lgt_id, params%max_treelevel+idx_mesh_lvl) = treecode_size(block_treecode(:,k), size(block_treecode,1))
        ! set refinement status
        lgt_block(free_lgt_id, params%max_treelevel+idx_refine_sts) = 0
        ! set number of the tree
        lgt_block(free_lgt_id, params%max_treelevel+idx_tree_id) = tree_id
        do dF = 1, params%n_eqn
             hvy_block( :, :, :, df, free_hvy_id ) = hvy_tmp(:, :, :, df, k)
        end do
     end do

    ! synchronize light data. This is necessary as all CPUs above created their blocks locally.
    ! As they all pass the same do loops, the counter array blocks_per_rank_list does not have to
    ! be synced. However, the light data has to.
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

    deallocate(block_treecode)

    ! it is useful to print out the information on active levels in the file
    ! to get an idea how it looks like and if the desired dense level is larger
    ! or smaller
    if (params%rank==0) then
        write(*,'("In the file we just read, Jmin=",i3," Jmax=",i3)') min_active_level( lgt_block ), &
        max_active_level( lgt_block )
    endif
   
    ! add the amount of hvy_n to the new data
    hvy_n = hvy_n + free_hvy_n

end subroutine read_tree
!##############################################




  !##############################################################
  subroutine allocate_hvy_lgt_data(params, lgt_block, hvy_block, hvy_neighbor, &
                  lgt_active, lgt_n, hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp)
    implicit none

    !-----------------------------------------------------------------
    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :) 
    !> light data size
    integer(kind=ik), allocatable, intent(out)          :: lgt_n(:) 
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: 00, k1, k2 etc)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:)
    ! local shortcuts:
    integer(kind=ik)                                    :: Bs, g, Neqn, number_blocks,&
                                                      rank, number_procs,number_trees
    integer(kind=ik)    :: rk_steps
    real(kind=rk)       :: effective_memory
    integer             :: status, nrhs_slots, nwork
    !-----------------------------------------------------------------
    ! variables initialization
    ! set parameters for readability
    rank            = params%rank
    number_blocks   = params%number_blocks
    Bs              = params%Bs
    g               = params%n_ghosts
    Neqn            = params%n_eqn
    number_procs    = params%number_procs

    ! 19 oct 2018: The work array hvy_work is modified to be used in "register-form"
    ! that means one rhs is stored in a 5D subset of a 6D array.
    ! Hence, nrhs_slots is number of slots for RHS saving:
    if (params%time_step_method == "RungeKuttaGeneric") then
        nrhs_slots = size(params%butcher_tableau,1)

    elseif (params%time_step_method == "Krylov") then
        nrhs_slots = params%M_krylov +3
    elseif (params%time_step_method == "no") then
        nrhs_slots = 0 ! no time stepping
    else
        call abort(191018161, "time_step_method is unkown: "//trim(adjustl(params%time_step_method)))
    endif

    nwork = max( 2*Neqn, params%N_fields_saved)

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'(A)') "INIT: Beginning memory allocation and initialization."
        write(*,'("INIT: mpisize=",i6)') params%number_procs
        write(*,'("INIT: nwork=",i6)') nwork
        write(*,'("INIT: Bs=",i7," blocks-per-rank=",i7," total blocks=", i7)') Bs, number_blocks, number_blocks*number_procs
    endif
    
    ! allocate memory
    if ( params%dim==3 ) then
        ! 3D:
        if (rank==0) write(*,'("INIT: Allocating a 3D case.")')

        !---------------------------------------------------------------------------
        allocate( hvy_block( Bs+2*g, Bs+2*g, Bs+2*g, Neqn, number_blocks ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_block", product(real(shape(hvy_block)))*8.0e-9, shape(hvy_block)
        endif

        !---------------------------------------------------------------------------
        ! work data (Runge-Kutta substeps and old time level)
        if (present(hvy_work)) then
            allocate( hvy_work( Bs+2*g, Bs+2*g, Bs+2*g, Neqn, number_blocks, nrhs_slots ) )
            if (rank==0) then
                write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
                "hvy_work", product(real(shape(hvy_work)))*8.0e-9, shape(hvy_work)
            endif
        end if

        if ( present(hvy_tmp) ) then
        allocate( hvy_tmp( Bs+2*g, Bs+2*g, Bs+2*g, nwork, number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_tmp", product(real(shape(hvy_tmp)))*8.0e-9, shape(hvy_tmp)
        endif
        endif

        !---------------------------------------------------------------------------
        ! 3D: maximal 74 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 74 ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_neighbor", product(real(shape(hvy_neighbor)))*8.0e-9, shape(hvy_neighbor)
        endif

    else

        ! 2D:
        if (rank==0) write(*,'("INIT: Allocating a 2D case.")')

        !---------------------------------------------------------------------------
        allocate( hvy_block( Bs+2*g, Bs+2*g, 1, Neqn, number_blocks ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_block", product(real(shape(hvy_block)))*8.0e-9, shape(hvy_block)
        endif

        !---------------------------------------------------------------------------
        ! work data (Runge-Kutta substeps and old time level)
        if ( present(hvy_work) ) then
            allocate( hvy_work( Bs+2*g, Bs+2*g, 1, Neqn, number_blocks, nrhs_slots ) )
            if (rank==0) then
                write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",6(i9,1x))') &
                "hvy_work", product(real(shape(hvy_work)))*8.0e-9, shape(hvy_work)
            endif
        end if

        if ( present(hvy_tmp) ) then
        allocate( hvy_tmp( Bs+2*g, Bs+2*g, 1, nwork, number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_tmp", product(real(shape(hvy_tmp)))*8.0e-9, shape(hvy_tmp)
        endif
        endif

        !---------------------------------------------------------------------------
        ! 2D: maximal 16 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 16 ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_neighbor", product(real(shape(hvy_neighbor)))*8.0e-9, shape(hvy_neighbor)
        endif

    end if

    !---------------------------------------------------------------------------)
    allocate( lgt_block( number_procs*number_blocks, params%max_treelevel+extra_lgt_fields) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_block", product(real(shape(lgt_block)))*4.0e-9, shape(lgt_block)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_sortednumlist( size(lgt_block,1), 2) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_sortednumlist", product(real(shape(lgt_sortednumlist)))*4.0e-9, shape(lgt_sortednumlist)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_active( size(lgt_block, 1), params%forest_size + 1 ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_active", product(real(shape(lgt_active)))*4.0e-9, shape(lgt_active)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_n( params%forest_size + 1 ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_n", product(real(shape(lgt_n)))*4.0e-9, shape(lgt_n)
    endif

    !---------------------------------------------------------------------------
    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_active", product(real(shape(hvy_active)))*4.0e-9, shape(hvy_active)
    endif


    if (rank == 0) then
        write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') number_blocks, Neqn
        write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*number_blocks
        write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') number_blocks

        effective_memory = (dble(size(hvy_block)) + & ! real data
        dble(size(lgt_block)+size(lgt_sortednumlist)+size(hvy_neighbor)+size(lgt_active)+size(hvy_active))/2.0 & ! integer (hence /2)
        )*8.0e-9 ! in GB

        if (present(hvy_tmp)) effective_memory = effective_memory + dble(size(hvy_tmp))*8.0e-9 ! in GB
        if (present(hvy_work)) effective_memory = effective_memory + dble(size(hvy_work))*8.0e-9 ! in GB

        ! note we currently use 8byte per real and integer by default, so all the same bytes per point
        write(*,'("INIT: Measured (true) local (on 1 cpu) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB per mpirank")') &
        effective_memory

        write(*,'("INIT-GLOBAL: Measured (true) TOTAL (on all CPU) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB")') &
        effective_memory*dble(params%number_procs)
    end if

  end subroutine allocate_hvy_lgt_data
  !##############################################################



  !##############################################################
  !> Counts the number of hvy_blocks owned by the tree with the 
  !> specified tree_id
  function count_tree_hvy_n( params, lgt_block, hvy_active, hvy_n, tree_id) &
           result(tree_hvy_n)

    implicit none

    !-----------------------------------------------------------------
    integer(kind=ik),  intent(in)  :: hvy_active(:) !< list of active blocks on the current rank
    integer(kind=ik), intent(in)  :: lgt_block(:,:) !< light data array
    integer(kind=ik),  intent(in)  :: tree_id, hvy_n !< index of the tree and number of active hvy_blocks 
    type (type_params), intent(in) :: params !< params structure of wabbit
    !-----------------------------------------------------------------
    integer (kind=ik) :: k, rank, lgt_id, tree_hvy_n, idxtree, N

    N = params%number_blocks
    rank = params%rank
    idxtree = params%max_treelevel + idx_tree_id
    tree_hvy_n = 0

    do k = 1, hvy_n
        ! calculate lgtblock id corresponding to hvy id
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )
        if (lgt_block(lgt_id, idxtree) == tree_id) then
            tree_hvy_n = tree_hvy_n + 1
        end if
    end do
  end function count_tree_hvy_n
  !##############################################################


  !##############################################################
  !> save a specific field with specified tree_id
  !> \TODO This function will replace write_field in the future
  subroutine write_tree_field(fname, params, lgt_block, lgt_active, hvy_block, &
                    lgt_n, hvy_n, hvy_active, dF, tree_id, time, iteration ) 

    implicit none

    !-----------------------------------------------------------------
    character(len=*), intent(in) :: fname  !< filename
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in)   :: lgt_block(:, :) !< ligh block data
    integer(kind=ik), intent(in)   :: lgt_active(:,:) !< list of active blocks for each tree
    integer(kind=ik), intent(in)   :: lgt_n(:), hvy_n !< length of active lists
    integer(kind=ik), intent(in)   :: hvy_active(:) !< list of active hvy blocks
    real(kind=rk), intent(in)      :: hvy_block(:, :, :, :, :)
    !--------------------
    ! optional parameter
    !--------------------
    integer(kind=ik), optional, intent(in) :: tree_id !< id of the tree (default: 1)
    real(kind=rk), optional, intent(in)    :: time !< time loop parameters (default: 0.0)
    integer(kind=ik), optional, intent(in) :: iteration !< iteration of the solver (default: 1)
    integer(kind=ik), optional,intent(in)  :: dF !< datafield number (default: 0)
    !-----------------------------------------------------------------
    integer (kind=ik) ::tree_hvy_n, treeid, it, dataField
    real (kind=rk) :: t

    ! check the optional params
    if (present(dF)) then
        dataField = dF 
    else
        dataField = 1
    endif
    if (present(tree_id)) then
        treeid = tree_id 
    else
        treeid = 1
    endif
    if (present(time)) then
        t = time 
    else
        t = 0.0_rk
    endif
    if (present(iteration)) then
       it = iteration
    else
        it = 0
    endif

    ! calculate how many hvy_n belong to the tree
    tree_hvy_n = count_tree_hvy_n(params, lgt_block, hvy_active, hvy_n, treeid) 
    ! write the data
    call write_field(fname, t, it, dataField, params, lgt_block, hvy_block, &
                    lgt_active(:,treeid), lgt_n(treeid), tree_hvy_n, hvy_active)
  end subroutine
  !##############################################################





  !##############################################################





end module module_forest
