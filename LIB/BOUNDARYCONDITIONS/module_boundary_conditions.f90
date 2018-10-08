!----------------------------------------------------------------
!> This Module provides generic functions for the implementation
!> of the boundary conditions
!> \details
!> \date 10.08.2018 - creation
!> \author Philipp Krah
!----------------------------------------------------------------
module module_boundary_conditions

    use module_globals
    use module_params
    use module_treelib,only: adjacent_block_2D
    implicit none

    ! make all functions privat if not explicit marked as public
    PRIVATE

    ! functions and routines acessible outside of this module
    !**********************************************************************************************
    PUBLIC :: get_adjacent_boundary_surface_normal,block_is_adjacent_to_boundary
    !**********************************************************************************************

    character(len=3), parameter, dimension(16)  :: DIR_2D  = (/ '__N', '__E', &
                                                          '__S', '__W', '_NE', '_NW', '_SE', &
                                                          '_SW', 'NNE', 'NNW', 'SSE', 'SSW', &
                                                          'ENE', 'ESE', 'WNW', 'WSW' /)

    !


contains

 !-----------------------------------------------------------------------------
 !> \brief Computes the surface normal of the boundary, if the current block is adjacent to the boundary.\n
 !> \details
 !>   - The surface normal is 0 if the block is not adjacent to the boundary or the boundaries are periodicaly
 !> continued.
 !>   - The surface normal is computed from the treecode
 subroutine get_adjacent_boundary_surface_normal(params,lgt_id,lgt_block,max_treelevel,surface)
   implicit none
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> user defined parameter structure
   type (type_params), intent(in) :: params
   !> light data
   integer(kind=ik)  , intent(in) :: lgt_id,lgt_block(:,:),max_treelevel
   !> surface normal of the adjacent boundary wall
   integer(kind=2)   ,intent(out) :: surface(3)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=ik) :: neighbor(max_treelevel),my_treecode(max_treelevel),my_level,i

   surface      =0
   my_treecode  =lgt_block(lgt_id,1:max_treelevel)
   my_level     =lgt_block(lgt_id,max_treelevel+idx_mesh_lvl)

   ! reset surface normal
   surface=0

   do i=1,4  ! loop over the four directions N E S W
     ! 1. calculate the neighbor treecode in this direction
     call adjacent_block_2D( my_treecode, neighbor, DIR_2D(i), my_level,max_treelevel)
     ! 2. find out if the neigbor is on the other side of the domain. If this is the case,
     !    then mark the side as boundary of the domain
     !    The boundary is crossed when the treecode changes on the coarsest level
     if (block_is_adjacent_to_boundary(params,DIR_2D(i),my_treecode,neighbor,my_level,max_treelevel) ) then
       select case ( DIR_2D(i))
       case( '__N' )
         surface(1) = -1 ! surface normal in the negative x direction n=e_1*(-1)
       case( '__W' )
         surface(2) = -1 ! surface normal in the negative y direction n=e_2*(-1)
       case( '__S' )
         surface(1) = 1 ! surface normal in the positive x direction n=e_1*(1)
       case( '__E' )
         surface(2) =  1 ! surface normal in the positive y direction n=e_2*1
       case default
         call abort(2808181, "Could not determine a surface normal for the boundary block!")
       end select
     end if
   end do

   ! write(*,('("my=",2(i1)," surf=",3(I2))')) my_treecode,surface

 end subroutine get_adjacent_boundary_surface_normal


!> is the boarder between my block and its neighbor corssing a domain boundary?
recursive logical function block_is_adjacent_to_boundary(params,dir,my_treecode,neighbor,my_level,max_treelevel) result(is_adjacent)
   implicit none
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   type (type_params), intent(in) :: params !< user defined parameter structure
   integer(kind=ik), intent(in)   :: max_treelevel,my_level
   !> treecode of the block and the adjacent block
   integer(kind=ik), intent(in) :: neighbor(max_treelevel), my_treecode(max_treelevel)
   !> direction of the block boundary
   character(len=3), intent(in) :: dir
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=ik),save, allocatable :: tmp_neighbor(:)

   if (.not. allocated(tmp_neighbor)) allocate(tmp_neighbor(max_treelevel))

   is_adjacent =.false.
   ! To find out if the boarder between two blocks is a domain boundary we
   ! compare the treecodes on the first level.
   ! Treecodes in the +-x direction (West/East) differe by +-1
   ! Treecodes in the +-y direction (North/South) differe by +-2
   !TREE CODES
   !      |              |              |
   !  3   |      2       |      3       |   2
   ! -------------------------------------------
   !      | 00   :    01 |      :   11  |
   !      |      :       |      :       |
   !  1   |......0.......|......1.......|   0
   !      | 02   :    03 |      :       |
   !      |      :       |      :       |
   ! ------------------------------------------ E
   !      |      :       |      :       |
   !      |      :       |      :       |
   !   3  |......2.......|......3.......|  2
   !      |      :       |      :       |
   !      | 22   :       |      :  33   |
   ! ------------------------------------------
   !      |              |              |
   !  1   |       0      |      1       |    0

   ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! check boundary in north, south, east, west
   ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if (.not. params%periodic_BC(2)) then! non periodic BC in y direction          ! NOTE:
    if((neighbor(1)-my_treecode(1) ==-1 .and. dir == DIR_2D(2)) .or. & !east      ! apparently the block synchronization of the
       (neighbor(1)-my_treecode(1) == 1 .and. dir == DIR_2D(4))) then  !west      ! block ghost node layer in east/west and north/south
              is_adjacent =.true.                                                 ! seams to be switched,
              return                                                              ! such that east/west is in x direction and north/south is in y direction
    endif
   endif

   if ( .not. params%periodic_BC(1) ) then ! non periodic BC in x direction
    if((neighbor(1)-my_treecode(1) == 2 .and. dir == DIR_2D(1)) .or. & !north
       (neighbor(1)-my_treecode(1) ==-2 .and. dir == DIR_2D(3))) then !south
           is_adjacent =.true.
           return
    endif
   endif

   ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! check boundary in NE,NW,SE,SW
   ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! For the tilted directions (NW/SE/NE/SW) we check if there is a adjacent
   ! boundary in any of the two directions:
   ! Exmaple:
   ! North East: 1. Check North: Is the my block adjacent to its northern boundary?
   !              if yes: return true
   !              if not:
   !             2. Check East!
   if (dir==DIR_2D(5)) then     ! north east
      ! check if a boundary is in the north
      call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(1), my_level,max_treelevel)
      if(block_is_adjacent_to_boundary(params,DIR_2D(1),my_treecode,tmp_neighbor,my_level,max_treelevel)) then
        is_adjacent = .true.
        return
      endif
      ! check if a boundary is in the east
      call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(2), my_level,max_treelevel)
      if( block_is_adjacent_to_boundary(params,DIR_2D(2),my_treecode,tmp_neighbor,my_level,max_treelevel)) then
            is_adjacent = .true.
            return
      endif

    endif

    if (dir==DIR_2D(6)) then    ! north west
      ! is my block adjacent to a domain boundary in the north?
      call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(1), my_level,max_treelevel)
      if ( block_is_adjacent_to_boundary(params,DIR_2D(1),my_treecode,tmp_neighbor,my_level,max_treelevel)) then
        is_adjacent = .true.
        return
      endif
      ! is my block adjacent to a domain boundary in the west?
      call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(4), my_level,max_treelevel)
      if( block_is_adjacent_to_boundary(params,DIR_2D(4),my_treecode,tmp_neighbor,my_level,max_treelevel)) then
            is_adjacent = .true.
            return
      endif

    endif

    if (dir==DIR_2D(7)) then     ! south east

      call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(2), my_level,max_treelevel)
      if (block_is_adjacent_to_boundary(params,DIR_2D(2),my_treecode,tmp_neighbor,my_level,max_treelevel)) then
        is_adjacent = .true.
        return
      endif
      call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(3), my_level,max_treelevel)
      if(block_is_adjacent_to_boundary(params,DIR_2D(3),my_treecode,tmp_neighbor,my_level,max_treelevel)) then
           is_adjacent = .true.
           return
      endif

    endif

    if (dir==DIR_2D(8)) then    ! south west

      call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(3), my_level,max_treelevel)
      if (block_is_adjacent_to_boundary(params,DIR_2D(3),my_treecode,tmp_neighbor,my_level,max_treelevel))then
        is_adjacent = .true.
        return
    endif
    call adjacent_block_2D( my_treecode, tmp_neighbor, dir_2d(4), my_level,max_treelevel)
    if(block_is_adjacent_to_boundary(params,DIR_2D(4),my_treecode,tmp_neighbor,my_level,max_treelevel)) then
          is_adjacent = .true.
          return
      endif

    endif


 end function
 !

end module module_boundary_conditions
