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

    character(len=3), parameter, dimension(16)  :: dir_2D  = (/ '__N', '__E', &
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
 subroutine get_adjacent_boundary_surface_normal(lgt_id,lgt_block,max_treelevel,surface)
   implicit none
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> light data
   integer(kind=ik)       , intent(in) :: lgt_id,lgt_block(:,:),max_treelevel
   !> surface normal of the adjacent boundary wall
   integer(kind=2)          , intent(out):: surface(3)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=ik)                   :: neighbor(max_treelevel),my_treecode(max_treelevel),level,i

   !integer(kind=2)                    :: normal_vec
   surface      =0
   my_treecode  =lgt_block(lgt_id,1:max_treelevel)
   level        =lgt_block(lgt_id,max_treelevel+idx_mesh_lvl)

   do i=1,4  ! loop over the four directions N E S W

     ! 1. calculate the neighbor treecode in this direction
     call adjacent_block_2D( my_treecode, neighbor, dir_2d(i), level, max_treelevel)

     ! 2. find out if the neigbor is on the other side of the domain. If this is the case,
     !    then mark the side as boundary of the domain
     !    The boundary is crossed when the treecode changes on the coarsest level
     if (    block_is_adjacent_to_boundary(dir_2D(i),my_treecode,neighbor,max_treelevel) ) then
       surface(2) =  1 ! surface normal in the positive y direction n=e_2*1
     elseif ( block_is_adjacent_to_boundary(dir_2D(i),my_treecode,neighbor,max_treelevel)) then
       surface(1) =  1 ! surface normal in the positive x direction n=e_1*(1)
     elseif ( block_is_adjacent_to_boundary(dir_2D(i),my_treecode,neighbor,max_treelevel)) then
       surface(2) = -1 ! surface normal in the negative y direction n=e_2*(-1)
     elseif ( block_is_adjacent_to_boundary(dir_2D(i),my_treecode,neighbor,max_treelevel)) then
       surface(1) = -1 ! surface normal in the negative x direction n=e_1*(-1)
     else
      ! do nothing
     endif

   end do


 end subroutine get_adjacent_boundary_surface_normal



recursive logical function block_is_adjacent_to_boundary(dir,my_treecode,neighbor,max_treelevel) result(is_adjacent)
   implicit none
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=ik), intent(in) :: max_treelevel
   integer(kind=ik), intent(in) :: neighbor(max_treelevel), my_treecode(max_treelevel)
   character(len=3), intent(in) :: dir
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if (  (neighbor(1)-my_treecode(1) == 2 .and. dir == dir_2D(1)) .or. & !north
         (neighbor(1)-my_treecode(1) ==-1 .and. dir == dir_2D(2)) .or. & !east
         (neighbor(1)-my_treecode(1) ==-2 .and. dir == dir_2D(3)) .or. & !south
         (neighbor(1)-my_treecode(1) == 1 .and. dir == dir_2D(4))) then  !west

            is_adjacent =.true.

   elseif (dir==dir_2D(5) .and. &     ! north east
      (block_is_adjacent_to_boundary(dir_2D(1),my_treecode,neighbor,max_treelevel) .or. &
      block_is_adjacent_to_boundary(dir_2D(2),my_treecode,neighbor,max_treelevel))) then

            is_adjacent = .true.

  elseif (dir==dir_2D(6) .and. &     ! north west
      (block_is_adjacent_to_boundary(dir_2D(1),my_treecode,neighbor,max_treelevel) .or. &
      block_is_adjacent_to_boundary(dir_2D(4),my_treecode,neighbor,max_treelevel))) then

            is_adjacent = .true.

  elseif (dir==dir_2D(7) .and. &     ! south east
      (block_is_adjacent_to_boundary(dir_2D(2),my_treecode,neighbor,max_treelevel) .or. &
      block_is_adjacent_to_boundary(dir_2D(3),my_treecode,neighbor,max_treelevel))) then

           is_adjacent = .true.

   elseif (dir==dir_2D(8) .and. &     ! south west
      (block_is_adjacent_to_boundary(dir_2D(3),my_treecode,neighbor,max_treelevel) .or. &
      block_is_adjacent_to_boundary(dir_2D(4),my_treecode,neighbor,max_treelevel))) then

          is_adjacent = .true.

   else

          is_adjacent =.false.
   endif


 end function
 !

end module module_boundary_conditions
