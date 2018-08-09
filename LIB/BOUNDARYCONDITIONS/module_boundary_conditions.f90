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
    PUBLIC :: get_adjacent_boundary_surface_normal
    !**********************************************************************************************

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
   integer(kind=1)          , intent(out):: surface(3)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=ik)                   :: neighbor(max_treelevel),i_bc,my_treecode(max_treelevel),level,i
   character(len=3),  parameter,dimension(16)     :: dir  &
   = [ character(len=3):: '__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW']

   !integer(kind=1)                    :: normal_vec
   surface      =0
   my_treecode  =lgt_block(lgt_id,1:max_treelevel)
   level        =lgt_block(lgt_id,max_treelevel+idx_mesh_lvl)

   do i=1,4  ! loop over the four directions

     ! 1. calculate the neighbor treecode in this directory
     call adjacent_block_2D( my_treecode, neighbor, dir(i), level, max_treelevel)

     ! 2. find out if the neigbor is on the other side of the domain. If this is the case,
     !    then mark the side as boundary of the domain
     !    The boundary is crossed when the treecode changes on the coarsest level.
     if (     neighbor(1)-my_treecode(1) == 2 .and. dir(i) == dir(1) ) then
       surface(2)=1 ! surface normal in the positive y direction n=e_2*1
     elseif ( neighbor(1)-my_treecode(1) ==-1 .and. dir(i) == dir(2)) then
       surface(1)=1 ! surface normal in the positive x direction n=e_1*(1)
     elseif ( neighbor(1)-my_treecode(1) ==-2 .and. dir(i) == dir(3)) then
       surface(2)=-1 ! surface normal in the negative y direction n=e_2*(-1)
     elseif ( neighbor(1)-my_treecode(1) ==1 .and. dir(i) == dir(4)) then
       surface(1)=1 ! surface normal in the negative x direction n=e_1*(-1)
     else
      ! do nothing
     endif

   end do


 end subroutine get_adjacent_boundary_surface_normal


 !
 !  function block_is_adjacent_to_boundary(surface,lgt_id,lgt_block,hvy_neighbor,max_treelevel)
 !   implicit none
 !   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !   integer(kind=ik)       , intent(in) :: lgt_id,lgt_block(:,:),max_treelevel
 !   !++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !   logical                             ::block_is
 !   integer(kind=ik)                    ::lvl, neighbor(1:maxtreelevel),level
 !
 !   my_treecode=lgt_block(lgt_id,1:max_treelevel)
 !   level   =lgt_block(lgt_id,max_treelevel+idx_mesh_lvl)
 !   call adjacent_block_2D( my_treecode, neighbor, '__N', level, max_treelevel)
 !   !! positiv y direction
 !   !! if my treecode differs from neighbor treecode by two on every digit, then we have
 !   !! crossed the upper domain boundary
 !   if (abs(my_treecode(1)-neighbor(1))==2) then
 !     surface(2)=1
 !   endif
 !
 ! end function
 !

end module module_boundary_conditions
