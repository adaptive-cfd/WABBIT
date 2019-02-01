
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
  PUBLIC :: read_field2tree, read_forest_params
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
      call read_param_mpi(FILE, 'Blocks', 'max_forest_size', params%max_forest_size, 1 )
      call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%n_ghosts, 1 )

  end subroutine read_forest_params
  !##############################################################

  !##############################################################
  !> This function reads a set of data into a specified tree.
  !> The set of data can involve multiple files. Each file should
  !> contain a quantity at the same snapshot (i.e. same iteration/time).
  !> If no data has been read in before. This function will allocate
  !> the heavy data for you.
  subroutine read_field2tree(params, fnames, N_files, tree_id, &
                            lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
                            hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
      implicit none
      !-----------------------------------------------------------------
      type (type_params), intent(inout) :: params           !< params structure
      integer(kind=ik), intent(in)      :: N_files     !< number of fields/quantities to read
      character(len=*), intent(in)      :: fnames(N_files)  !< file names
      integer(kind=ik), intent(inout)   :: hvy_n, lgt_n     !< number of heavy and light active blocks
      integer(kind=ik), intent(in)      :: tree_id     !< number of the tree
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_block(:,:,:)  !< light data array
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
      integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_active(:,:), hvy_active(:) !< active lists
      integer(kind=tsize), ALLOCATABLE, intent(inout):: lgt_sortednumlist(:,:)
      real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
      !-------------------------------- ---------------------------------
      integer(kind=ik)       :: iteration, dF, Bs, tc_length, dim
      real(kind=rk):: time, domain(3)

      call read_attributes(fnames(1), lgt_n, time, iteration, domain, Bs, tc_length, dim)
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
          params%n_eqn = N_files
          params%N_fields_saved = N_files
          params%time_step_method="no" ! is necessary for
          params%Bs = Bs
          params%number_blocks = calculate_max_nr_blocks ( params )
          ! we have to allocate grid if this routine is called for the first time
          call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
                        lgt_active, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)
      endif
      ! read treecode from first input file
      call read_mesh(fnames(1), params, lgt_n, hvy_n, lgt_block, tree_id)
      ! read datafields from files into hvy_block array
      do dF = 1, N_files
          call check_file_exists(trim(fnames(dF)))
          call read_field(fnames(dF), dF, params, hvy_block, hvy_n )
      end do

  end subroutine read_field2tree
  !##############################################################





end module module_forest
