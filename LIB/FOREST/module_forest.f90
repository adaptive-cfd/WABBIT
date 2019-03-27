
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
  PUBLIC :: read_field2tree,  write_tree_field, &
            add_tree, count_tree_hvy_n, allocate_hvy_lgt_data, delete_tree, &
            copy_tree, multiply_tree, scalar_multiplication_tree
  !**********************************************************************************************
contains




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
                           rank, level, Bs(3), fsize, number_dense_blocks
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
          if ( any(params%Bs.ne.Bs )) call abort(100119,"Block size should not change);!")
          if ( params%max_treelevel < tc_length) call abort(10119,"maximal treelevel incompatible")
      else
          params%max_treelevel = tc_length
          params%dim = dim
          params%n_eqn = N_files
          params%N_fields_saved = N_files
          params%Bs = Bs
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

  end subroutine read_field2tree
  !##############################################################

  
  !##############################################################
  !> deletes the lgt_block data of tree with given tree_id and creates new
  !> active and sorted lists
  subroutine delete_tree(params, tree_n, lgt_block, lgt_active, lgt_n, &
                        lgt_sortednumlist, hvy_active, hvy_n, tree_id)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in)      :: params    !< user defined parameter structure
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)!< light data array
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)!< number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)!< list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: hvy_n!< number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: tree_n!< highest tree id
    integer(kind=ik), intent(in)        :: tree_id!< highest tree id
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)!< sorted light data with numerical treecodes
    !-----------------------------------------------------------------
    integer(kind=ik)                    :: k, lgt_id

    ! loop over active list of tree 
    do k = 1, lgt_n(tree_id) 
        ! block is active
        lgt_id = lgt_active(k, tree_id)
        lgt_block(lgt_id,:) = -1_ik
    end do
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. , tree_n)


  end subroutine delete_tree
  !##############################################################



  !##############################################################
  !> copy hvy and lgt data of tree_id2 to tree_id1
  !> after copy the active and sorted lists are updated, as well 
  !> as neighbors and ghost nodes
  subroutine copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id1, tree_id2)

      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(in) :: params   !< params structure
      integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
      integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
      integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< all data from tree_id2 gets copied to tree_id1
      integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
      integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
      integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
      !-----------------------------------------------------------------
      integer(kind=ik)    :: Jmax, lgt_id1, lgt_id2, hvy_id1, hvy_id2, fsize
      integer(kind=ik)    :: k2, N, rank

      Jmax = params%max_treelevel ! max treelevel
      fsize= params%forest_size   ! maximal number of trees in forest
      rank = params%rank
      N = params%number_blocks

      !if (params%rank == 0 ) write(*,'("Copy tree: ",i3,"<-",i3)') tree_id1, tree_id2
      if (tree_id1 > fsize) call abort(0403191,"tree_id you are asking for does not exist!")

      ! first we delete tree_id1 if it is allocated
      if (tree_id1 <= tree_n)then ! tree_id1 only exists if it is in the list of active trees
        call delete_tree(params, tree_n, lgt_block, lgt_active, lgt_n, &
                        lgt_sortednumlist, hvy_active, hvy_n, tree_id1)
      end if

      ! Loop over the active hvy_data
      do k2 = 1, hvy_n
        hvy_id2 = hvy_active(k2)
        call hvy_id_to_lgt_id(lgt_id2, hvy_id2, rank, N )
        ! first we have to find out if the hvy data belongs to the tree we want to copy
        if ( lgt_block(lgt_id2, Jmax + idx_tree_id) .ne. tree_id2) cycle
        call get_free_local_light_id( params, rank, lgt_block, lgt_id1)
        call lgt_id_to_hvy_id( hvy_id1, lgt_id1, rank, N )
        !--------------------
        ! Light DATA
        !--------------------
        lgt_block(lgt_id1, :) = lgt_block(lgt_id2, :)  ! copy light data
        lgt_block(lgt_id1, Jmax + idx_tree_id) = tree_id1  ! asign new tree_id
        !--------------------
        ! Heavy DATA
        !--------------------
        hvy_block( :, :, :, :, hvy_id1) = hvy_block( :, :, :, :, hvy_id2)
      end do ! loop over tree2

      ! always synchronice lgt_data when you have changed lgt_block locally 
      ! (i.e. when looping over lgt_ids which originate from a hvy_id)
      call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )
      call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. , tree_n)
      call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,fsize + 1),&
                            lgt_n(fsize + 1), lgt_sortednumlist, hvy_active, hvy_n )
      call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

   
  end subroutine
  !##############################################################

  !##############################################################
  !> multiply every element of the tree with a given value alpha
  subroutine scalar_multiplication_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id, alpha, verbosity)

      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(in) :: params   !< params structure
      integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
      integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
      integer(kind=ik), intent(in)      :: tree_id !< all data from tree_id2 gets copied to tree_id1
      integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
      integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      real(kind=rk), intent(in)      :: alpha !< heavy data array - block data
      integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
      integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
      logical, intent(in),optional      :: verbosity !< if true aditional stdout is printed
      !-----------------------------------------------------------------
      integer(kind=ik)    :: Jmax, lgt_id, hvy_id, fsize
      integer(kind=ik)    :: k, N, rank
      logical :: verbose=.false.
      if (present(verbosity)) verbose=verbosity
      Jmax = params%max_treelevel ! max treelevel
      fsize= params%forest_size   ! maximal number of trees in forest
      rank = params%rank
      N = params%number_blocks

      if (params%rank == 0 .and. verbose ) write(*,'("scalar multiplication tree: ",i3)') tree_id

      ! Loop over the active hvy_data
      do k = 1, hvy_n
        hvy_id = hvy_active(k)
        call hvy_id_to_lgt_id(lgt_id, hvy_id, rank, N )
        ! first we have to find out if the hvy data belongs to the tree we want to copy
        if ( lgt_block(lgt_id, Jmax + idx_tree_id) .ne. tree_id) cycle
        !--------------------
        ! Heavy DATA
        !--------------------
        hvy_block( :, :, :, :, hvy_id) = alpha * hvy_block( :, :, :, :, hvy_id)
      end do ! loop over tree2

   
  end subroutine
  !##############################################################



  !##############################################################
  !> This function takes to trees and refines their treestructures
  !> until they are refined to the same max level.
  !> The resulting trees are identical in treestructure. Hence
  !> they have the same treecodes and number of blocks.
  !> Remark: You may want to balance the load after using this routine.
  subroutine refine_trees2same_lvl(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, verbosity)

      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(in) :: params   !< params structure
      integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
      integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
      integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
      integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
      integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
      integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
      real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
      logical, intent(in),optional      :: verbosity
      !-----------------------------------------------------------------
      integer(kind=ik)    :: rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
      integer(kind=ik)    :: k1, k2, Nblocks_refined, level_min
      integer(kind=tsize) :: treecode1, treecode2
      logical :: verbose=.false.

      Jmax = params%max_treelevel ! max treelevel
      fsize= params%forest_size   ! maximal number of trees in forest

      if (present(verbosity)) verbose=verbosity
      if (params%rank == 0 .and. verbose ) write(*,'("Refining trees to same level: ",i9,",",i9)') tree_id1, tree_id2
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
              if (params%rank == 0) write(*,'("Number of blocks marked for refinement: ",i9)') Nblocks_refined
              ! 1) check gradedness of the grid (meshlevel of ajacent blocks should not differe by more than 1
              call ensure_gradedness( params, lgt_block, hvy_neighbor, &
              lgt_active(:, fsize+1), lgt_n(fsize+1), lgt_sortednumlist, hvy_active, hvy_n )
              ! 2) refine blocks
              if ( params%dim==3 ) then
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



!#############################################################
!> Perform pointwise operations (+, -, /, *) for two given trees.
!> The trees will be refined to same treestructure and then the
!> operation will be executed on each block pointwise.
!> Its like the Matlabs .* or ./ for trees  
!                                                      @@@
!                       @                              @@@                               @
!         @             @    @ @                                          @             @    @ @
!           @    @@      @ @                     @ @ @ @ @ @ @ @             @    @@      @ @
!           @  @          @        @@                                       @  @          @        @@
!      @     @@        @  @      @@@                   @@@              @     @@        @  @      @@@
!      @ @     @@@ @    @ @     @    @                 @@@             @ @     @@@ @    @ @     @    @
!       @@        @@     @     @@@@@@@                                  @@        @@     @     @@@@@@@
!         @@@      @@  @@    @@                         @                 @@@      @@  @@    @@
!           @@@    @@@@@ @@@@@@@         @@             @                   @@@    @@@@@ @@@@@@@         @@
!    @@@@@@          @@@@       @@@@@@@@                @            @@@@@@          @@@@       @@@@@@@@
!        @           @@@         @               @ @ @ @ @ @ @ @         @           @@@         @
!      @             @@@                                @              @             @@@
!                   @@@                                 @                           @@@
!                   @@@                                 @                           @@@
!                   @@@                                                             @@@
!                   @@@                             @       @                       @@@
!                   @@@                              @     @                        @@@
!                   @@@                               @   @                         @@@
!                   @@@@                                @                           @@@@
!                    @@@@                        @ @ @ @ @ @ @ @                    @@@@
!                  @@@@@@@@    @                        @                         @@@@@@@@    @
!             @@@     @    @@@                        @   @                   @@@     @    @@@
!          @@@  @      @       @                     @     @               @@@  @      @       @
!                                                   @       @
! (done with png to ascii)
  subroutine tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, &
             operation)

      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(in) :: params   !< params structure
      integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
      integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
      integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
      integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
      integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
      integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
      real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
      character (len=*)                 :: operation !< which arithmetical operation (+,-,*,/) which is applied
      !-----------------------------------------------------------------
      integer(kind=ik)    :: my_rank, level1, level2, Jmax, lgt_id1, lgt_id2, fsize
      integer(kind=ik)    :: k1, k2, N, level_min, rank1, rank2
      integer(kind=tsize) :: treecode1, treecode2, hvy_id1, hvy_id2

      Jmax = params%max_treelevel ! max treelevel
      fsize= params%forest_size   ! maximal number of trees in forest
      N    = params%number_blocks ! number of blocks per rank
      my_rank = params%rank  ! proc rank

      !=============================================
      ! Prepare the trees for pointwise arithmentic
      !=============================================
      ! The Trees can only be added when their grids are identical. At present we
      ! try to keep the finest levels of all trees. This means we refine
      ! all blocks which are not on the same level.
      call refine_trees2same_lvl(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2)
      ! because all trees have the same treestructrue thier hilbertcurve is identical
      ! and therefore balance the load will try to distribute blocks with the same
      ! treecode (but different trees) at the same rank.
      call balance_load( params, lgt_block, hvy_block,  hvy_neighbor, &
          lgt_active(:, fsize + 1), lgt_n(fsize + 1), lgt_sortednumlist, hvy_active, hvy_n, hvy_tmp )
      ! after balance_load the active and sorted lists have to be allways updated.
      call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. , tree_n)
      ! because it may happen that after balance_load blocks with the same treecode
      ! are not on the same rank, we have to make sure that they are
      ! before we can add them!
      call same_block_distribution(params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                                hvy_block, hvy_active, hvy_n, hvy_tmp, tree_n, tree_id1, tree_id2)

      !=================================================
      ! Decide which pointwice arithmetic shell be used
      !=================================================
!      _
!     / \    from here we assume that trees have the same grid structure
!    /   \   and the blocks of the same treecode are on the same rank !
!   /     \
!  /caution\
! +--------+
    select case(operation)
        case("+")  
             !         # 
             !         #
             !    ###########           ADDITION
             !         #
             !         #
             do k1 = 1, hvy_n
                hvy_id1 = hvy_active(k1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, my_rank, N )
                ! we want to add everything to tree1
                if (lgt_block(lgt_id1, Jmax + idx_tree_id) .ne. tree_id1) cycle
                level1   = lgt_block(lgt_id1, Jmax + idx_mesh_lvl)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n
                    hvy_id2 = hvy_active(k2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, my_rank, N )
                    if (lgt_block(lgt_id2, Jmax + idx_tree_id) .ne. tree_id2) cycle
                    level2   = lgt_block(lgt_id2, Jmax + idx_mesh_lvl)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        !> \todo make this one faster by looping only over inner grid points
                        hvy_block(:,:,:,:,hvy_id1) = hvy_block(:,:,:,:,hvy_id1) + &
                                                     hvy_block(:,:,:,:,hvy_id2)
                        exit
                    end if
                end do
             end do
         case("-")
             
             !          
             !         
             !    ###########           SUBSTRACTION
             !         
             !         
             do k1 = 1, hvy_n
                hvy_id1 = hvy_active(k1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, my_rank, N )
                ! we want to add everything to tree1
                if (lgt_block(lgt_id1, Jmax + idx_tree_id) .ne. tree_id1) cycle
                level1   = lgt_block(lgt_id1, Jmax + idx_mesh_lvl)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n
                    hvy_id2 = hvy_active(k2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, my_rank, N )
                    if (lgt_block(lgt_id2, Jmax + idx_tree_id) .ne. tree_id2) cycle
                    level2   = lgt_block(lgt_id2, Jmax + idx_mesh_lvl)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        !> \todo make this one faster by looping only over inner grid points
                        hvy_block(:,:,:,:,hvy_id1) = hvy_block(:,:,:,:,hvy_id1) - &
                                                     hvy_block(:,:,:,:,hvy_id2)
                        exit
                    end if
                end do
             end do
         case("/")
             !         # 
             !         
             !    ###########          Division 
             !         
             !         # 
             do k1 = 1, hvy_n
                hvy_id1 = hvy_active(k1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, my_rank, N )
                ! we want to add everything to tree1
                if (lgt_block(lgt_id1, Jmax + idx_tree_id) .ne. tree_id1) cycle
                level1   = lgt_block(lgt_id1, Jmax + idx_mesh_lvl)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n
                    hvy_id2 = hvy_active(k2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, my_rank, N )
                    if (lgt_block(lgt_id2, Jmax + idx_tree_id) .ne. tree_id2) cycle
                    level2   = lgt_block(lgt_id2, Jmax + idx_mesh_lvl)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        !> \todo make this one faster by looping only over inner grid points
                        hvy_block(:,:,:,:,hvy_id1) = hvy_block(:,:,:,:,hvy_id1) / &
                                                     hvy_block(:,:,:,:,hvy_id2)
                        exit
                    end if
                end do
             end do
         case("*")
             !        # # 
             !         # 
             !     #########        Multiply 
             !         # 
             !        # # 
             do k1 = 1, hvy_n
                hvy_id1 = hvy_active(k1)
                call hvy_id_to_lgt_id(lgt_id1, hvy_id1, my_rank, N )
                ! we want to add everything to tree1
                if (lgt_block(lgt_id1, Jmax + idx_tree_id) .ne. tree_id1) cycle
                level1   = lgt_block(lgt_id1, Jmax + idx_mesh_lvl)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level1))
                do k2 = 1, hvy_n
                    hvy_id2 = hvy_active(k2)
                    call hvy_id_to_lgt_id(lgt_id2, hvy_id2, my_rank, N )
                    if (lgt_block(lgt_id2, Jmax + idx_tree_id) .ne. tree_id2) cycle
                    level2   = lgt_block(lgt_id2, Jmax + idx_mesh_lvl)
                    treecode2= treecode2int(lgt_block(lgt_id2, 1 : level2))
                    if (treecode1 .ne. treecode2 ) then
                        cycle
                    else
                        !> \todo make this one faster by looping only over inner grid points
                        hvy_block(:,:,:,:,hvy_id1) = hvy_block(:,:,:,:,hvy_id1) * &
                                                     hvy_block(:,:,:,:,hvy_id2)
                        exit
                    end if
                end do
             end do

         case default
             call abort(135,"Operation unknown")
    end select
  end subroutine 
  !##############################################################



!##############################################################
subroutine add_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params   !< params structure
    integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
    integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
    integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
    integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
    integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
    integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
    integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    logical :: verbose=.false.

    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose) write(*,'("Adding trees: ",i4,",",i4)') tree_id1, tree_id2
    call tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
         hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2,"+")

end subroutine
!##############################################################



!##############################################################
subroutine multiply_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
             hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2, verbosity)

    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(in) :: params   !< params structure
    integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
    integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
    integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
    integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
    integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
    real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
    integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
    integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
    real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) !< used for saving, filtering, and helper qtys
    logical, intent(in),optional      :: verbosity !< if true: additional information of processing
    !-----------------------------------------------------------------
    logical :: verbose=.false.

    if (present(verbosity)) verbose=verbosity
    if (params%rank == 0 .and. verbose ) write(*,'("Multiply trees: ",i4,",",i4)') tree_id1, tree_id2
    call tree_pointwise_arithmetic(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
         hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1, tree_id2,"*")

end subroutine
!##############################################################



  !##############################################################
  !> This function compares the tree structure and block distribution of two trees.
  !> If treecodes are identical but the processor ranks not, then we redistribute
  !> one of the corresponding blocks, such that blocks with same treecode are on the
  !> same rank.
  subroutine same_block_distribution(params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                                hvy_block, hvy_active, hvy_n, hvy_tmp, &
                                tree_n, tree_id1, tree_id2)

      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(in) :: params   !< params structure
      integer(kind=ik), intent(inout)   :: hvy_n    !< number of active heavy blocks
      integer(kind=ik), intent(inout)   :: tree_n   !< number of trees in forest
      integer(kind=ik), intent(in)      :: tree_id1, tree_id2 !< number of the tree
      integer(kind=ik), intent(inout)   :: lgt_n(:) !< number of light active blocks
      integer(kind=ik), intent(inout)   :: lgt_block(:, : )  !< light data array
      real(kind=rk), intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), intent(inout)   :: lgt_active(:, :), hvy_active(:) !< active lists
      integer(kind=tsize), intent(inout):: lgt_sortednumlist(:,:)
      real(kind=rk), intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
      !-----------------------------------------------------------------
      integer(kind=ik)    :: level1, level2, lgt_id1, lgt_id2, N, fsize
      integer(kind=ik)    :: k1, k2, rank1, rank2, level_min, n_comm, Jmax
      integer(kind=tsize) :: treecode1, treecode2
      integer(kind=ik), allocatable, save :: comm_list(:,:)

      N = params%number_blocks
      fsize = params%forest_size
      Jmax = params%max_treelevel
      if (.not.allocated(comm_list)) allocate( comm_list( params%number_procs*N, 3 ) )

      !! Loop over all treecodes of both trees and check if they are identical.
      !! If the treecodes are the same but the blocks are not on the same rank,
      !! we have to send a block.
      !! The comm_list is created to save all of the block transfers.
      n_comm = 0 ! number of communications
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
                level2   = lgt_block(lgt_id2, Jmax + idx_mesh_lvl)
                level_min= min(level1,level2)
                treecode1= treecode2int( lgt_block(lgt_id1, 1 : level_min))
                treecode2= treecode2int(lgt_block(lgt_id2, 1 : level_min))
                if (treecode1 == treecode2 )then
                   ! check treestructure
                   if (level1 .ne. level2) then
                       call abort(270219, "Trees need same treestructure to fetch processor block distribution")
                   endif
                   ! check the block distribution (leaves are on same rank)
                   call lgt_id_to_proc_rank( rank1, lgt_id1, N)
                   call lgt_id_to_proc_rank( rank2, lgt_id2, N)
                   if (rank1 .ne. rank2) then
                     n_comm = n_comm + 1
                     comm_list(n_comm, 1) = rank1   ! sender mpirank
                     comm_list(n_comm, 2) = rank2   ! receiver mpirank
                     comm_list(n_comm, 3) = lgt_id1 ! block lgt_id to send
                   end if
                end if
             end do ! loop over tree2
        end do ! loop over tree1
        !if (params%rank == 0) write(*,'("nr blocks redistributed: ",i9)') n_comm
        ! Transfer the blocks
        call block_xfer( params, comm_list, n_comm, lgt_block, hvy_block, &
                 lgt_active(:, fsize+1), lgt_n(fsize+1), lgt_sortednumlist, hvy_tmp )
        ! after block transfer we have to create new lists
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
           lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. , tree_n)
  end subroutine
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
      integer(kind=ik) :: k, N, rank, number_procs, ierr, treecode_size, tree_id, Bs(3), g, dF
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
                  lgt_active, lgt_n, hvy_active, lgt_sortednumlist, hvy_work, &
                  hvy_gridQ, hvy_tmp)
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
    !> grid qtys (qtys that do depend only on the grid, i.e. mask funciton and geometry factors)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_gridQ(:, :, :, :, :)
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
    integer(kind=ik)                                    :: Bs(3), g, Neqn, number_blocks,&
                                                      rank, number_procs,number_trees
    integer(kind=ik)    :: rk_steps, i
    real(kind=rk)       :: effective_memory, mem_per_block, maxmem
    real(kind=rk), parameter ::  nstages = 2.0_rk ! stages for ghost node synching
    character(len=80)   :: memstring
    integer             :: status, nrhs_slots, nwork, nx, ny, nz, max_neighbors
    !-----------------------------------------------------------------
    ! variables initialization
    ! set parameters for readability
    rank            = params%rank
    number_blocks   = params%number_blocks
    Bs              = params%Bs
    g               = params%n_ghosts
    Neqn            = params%n_eqn
    number_procs    = params%number_procs
    nx = Bs(1)+2*g
    ny = Bs(2)+2*g
    nz = Bs(3)+2*g

    if (params%dim==2) then
        nz = 1
        max_neighbors = 16
    else
      max_neighbors   = 74
    endif


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
        write(*,'("INIT: Bs=",3(i7)," blocks-per-rank=",i7," total blocks=", i7)') Bs, number_blocks, number_blocks*number_procs
    endif

   !Automatic memory management. If specified --memory=0.3GB in the call line,
    if (params%number_blocks < 1) then
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
              if (params%rank==0) write(*,'("INIT: automatic selection of blocks per rank is active!")')

              ! read memory from command line (in GB)
              read(memstring(10:len_trim(memstring)-2),* ) maxmem

              if (params%rank==0) write(*,'("INIT: total memory: ",f9.4,"GB")') maxmem

              ! memory per MPIRANK (in GB)
              maxmem = maxmem / dble(params%number_procs)

              if (params%rank==0) write(*,'("INIT: memory-per-rank: ",f9.4,"GB")') maxmem

              if (params%dim==3) then
                mem_per_block = real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g)) & ! hvy_block
                + real(params%n_gridQ) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g)) & ! hvy_gridQ
                + real(max( 2*Neqn, params%N_fields_saved)) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g)) & ! hvy_tmp
                + real(Neqn) * real(nrhs_slots) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g)) & ! hvy_work
                + 2.0 * nstages * real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g) - ((Bs(1))*(Bs(2))*(Bs(3)))) &  ! real buffer ghosts
                + 2.0 * nstages * real(max_neighbors) * 5 / 2.0 ! int bufer (4byte hence /2)
              else
                mem_per_block = real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)) & ! hvy_block
                + real(params%n_gridQ) * real((Bs(1)+2*g)*(Bs(2)+2*g)) & ! hvy_gridQ
                + real(max(2*Neqn, params%N_fields_saved)) * real((Bs(1)+2*g)*(Bs(2)+2*g)) & ! hvy_tmp
                + real(Neqn) * real(nrhs_slots) * real((Bs(1)+2*g)*(Bs(2)+2*g)) & ! hvy_work
                + 2.0 * nstages * real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g) - (Bs(1)*Bs(2))) &  ! real buffer ghosts
                + 2.0 * nstages * real(max_neighbors) * 5 / 2.0 ! int bufer (4byte hence /2)
              endif

              ! in GB:
              mem_per_block = mem_per_block * 8.0e-9
              params%number_blocks = nint( maxmem / mem_per_block)
              number_blocks = params%number_blocks
              if (params%rank==0) then
                  write(*,'("INIT: for the desired memory we can allocate ",i8," blocks per rank")') params%number_blocks
                  write(*,'("INIT: we allocated ",i8," blocks per rank (total: ",i8," blocks) ")') params%number_blocks, &
                  params%number_blocks*params%number_procs
              endif

              if ((params%adapt_mesh .or. params%adapt_inicond) .and. (params%number_blocks<2**params%dim)) then
                  call abort(1909181740,"[ini_file_to_params.f90]: The number of blocks for the --memory specification is lower than 2^d&
                  & and comp is adaptive: we cannot fetch all 2^d blocks on one CPU for merging. Use more memory or less cores.")
              endif
          endif
      end do
    endif


    !==========================
    ! hvy_data
    !==========================
    allocate( hvy_block( nx, ny, nz, Neqn, params%number_blocks ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_block", product(real(shape(hvy_block)))*8.0e-9, shape(hvy_block)
    endif

    !---------------------------------------------------------------------------
    ! work data (Runge-Kutta substeps and old time level)
    if (present(hvy_work)) then
        allocate( hvy_work( nx, ny, nz, Neqn, params%number_blocks, nrhs_slots ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_work", product(real(shape(hvy_work)))*8.0e-9, shape(hvy_work)
        endif
    end if

    !---------------------------------------------------------------------------
    if ( present(hvy_tmp) ) then
        allocate( hvy_tmp( nx, ny, nz, nwork, params%number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_tmp", product(real(shape(hvy_tmp)))*8.0e-9, shape(hvy_tmp)
        endif
    endif

    !---------------------------------------------------------------------------
    if ( present(hvy_gridQ) .and. params%n_gridQ > 0 ) then
        allocate( hvy_gridQ( nx, ny, nz, params%n_gridQ, params%number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_gridQ", product(real(shape(hvy_gridQ)))*8.0e-9, shape(hvy_gridQ)
        endif
    endif

    !---------------------------------------------------------------------------
    allocate( hvy_neighbor( params%number_blocks, max_neighbors ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_neighbor", product(real(shape(hvy_neighbor)))*8.0e-9, shape(hvy_neighbor)
    endif


    !---------------------------------------------------------------------------
    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_active", product(real(shape(hvy_active)))*4.0e-9, shape(hvy_active)
    endif

    !==========================
    ! lgt_data
    !==========================
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






end module module_forest
