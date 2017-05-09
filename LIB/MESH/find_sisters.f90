!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> name: find_sisters.f90
!> version: 0.5
!> author: engels
!
!> \brief To a given light id "my_id", find the 3 (2D) or 7 (3D) sister block that have a common mother
!! block. They are returned in the sisters array.
!
!> \details
!! input:    - light data array \n
!! output:   - light data array
!!
!!
!! = log ======================================================================================
!! \n
!! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine find_sisters( params, lgt_my_id, lgt_sisters_id, lgt_block, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! modules


!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> this is the block whose sisters we look for
    integer(kind=ik), intent(in)        :: lgt_my_id
    !> here we will return the sisters
    integer(kind=ik), intent(inout)     :: lgt_sisters_id(:)
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! loop variables
    integer(kind=ik)                    :: i
    ! treecode variable
    integer(kind=ik), allocatable       :: all_treecodes(:,:)
    ! block level
    integer(kind=ik)                    :: N_sisters
    integer(kind=ik)                    :: mother_level, my_level
    logical                             :: exists


!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization
  ! check out how many sisters we look for. The number can be 4 or 8 in 2D or 3D. Note the
  ! block whose sisters we look for is returned as well, if 4 or 8 is called. if 3 or 7 is called,
  ! omit returning this block.
  N_sisters = size(lgt_sisters_id)
  if ( N_sisters /= 4 .and. N_sisters /= 8 .and. N_sisters /= 3 .and. N_sisters /= 7 ) then
    call error_msg("find_sisters: you don't ask for a valid number of sisters")
  endif

  ! allocate an array for all treecodes (including all 4/8 sisters)
  allocate( all_treecodes(1:N_sisters,1:params%max_treelevel) )
  ! initialize array as -1, since we do not use all of it, possibly (if we do not happen to
  ! be on the highest level)
  all_treecodes = -1

  my_level = lgt_block( lgt_my_id, params%max_treelevel+1 )

  lgt_sisters_id = -1
!---------------------------------------------------------------------------------------------
! main body

  ! Find sisters. The sister blocks have the same mother, that means their treecode
  ! is idential up to the last entry
  mother_level = my_level - 1

  if (N_sisters==3 .or. N_sisters==7) then
    ! fill the Array of sisters treecodes, skip block in question
    do i = 1, N_sisters
      if ( i-1 /= lgt_block( lgt_my_id, my_level )) then
        ! copy the idential mother level
        all_treecodes(i,1:mother_level) = lgt_block( lgt_my_id, 1:mother_level )
        ! the last index is (0..3) or (0..7)
        all_treecodes(i,mother_level+1) = i-1
        ! look for the sisters in the list of blocks (light data), store their ID if found
        ! (-1 otherwise)
        call does_block_exist( all_treecodes(i,:), lgt_block, params%max_treelevel, exists, lgt_sisters_id(i), lgt_active, lgt_n)
      end if
    end do

  else
    ! fill the Array of sisters treecodes, including block in question
    do i = 1, N_sisters
      ! copy the idential mother level
      all_treecodes(i,1:mother_level) = lgt_block( lgt_my_id, 1:mother_level )
      ! the last index is (0..3) or (0..7)
      all_treecodes(i,mother_level+1) = i-1
      ! look for the sisters in the list of blocks (light data), store their ID if found
      ! (-1 otherwise)
      call does_block_exist( all_treecodes(i,:), lgt_block, params%max_treelevel, exists, lgt_sisters_id(i), lgt_active, lgt_n)
    end do

  end if

  deallocate( all_treecodes )

end subroutine find_sisters
