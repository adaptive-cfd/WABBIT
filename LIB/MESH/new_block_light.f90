! ********************************
! WABBIT
! --------------------------------
!
! create new block, light data
!
! name: new_block_light.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine new_block_light(k, treecode, treeN)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                :: k, treeN
    integer(kind=ik), dimension(treeN), intent(in)              :: treecode

    integer(kind=ik)                                            :: treecode_size

    ! error handling
    if ( (k <= 0) .or. (k > blocks_params%number_max_blocks) ) then

      write(*,*) "ERROR! You try to create a block outside of the list"
      write(*,'("your id: ",i8," N_max_blocks=",i8)') k, blocks_params%number_max_blocks
      stop

    endif

    ! set block status
    blocks(k)%active                = .true.

    ! save treecode
    blocks(k)%treecode(1:treeN)     = treecode
    blocks(k)%level                 = treecode_size(treecode)

    ! refinement status
    blocks(k)%refinement            = 0

    ! reset neighbor status
    blocks(k)%neighbor_id           = -1
    blocks(k)%neighbor_treecode     = -1

end subroutine new_block_light
