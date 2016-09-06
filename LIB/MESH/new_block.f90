! ********************************
! 2D AMR prototype
! --------------------------------
!
! create new block
!
! name: new_block.f90
! date: 12.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine new_block(k, treecode, treeN, data_, ix, iy, Bs, g, dF)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                :: k, Bs, g, treeN, dF
    integer(kind=ik), dimension(treeN), intent(in)              :: treecode

    real(kind=rk), dimension(Bs+2*g, Bs+2*g), intent(in)        :: data_
    real(kind=rk), dimension(Bs), intent(in)                    :: ix, iy

    integer(kind=ik)                                            :: treecode_size

    ! error handling
    if ( (k <= 0) .or. (k > blocks_params%number_max_blocks) ) then

      write(*,*) "ERROR! You try to create a block outside of the list"
      write(*,'("your id: ",i8," N_max_blocks=",i8)') k, blocks_params%number_max_blocks
      stop

    endif

    ! save data
    blocks(k)%data_fields(dF)%data_ = data_

    ! set block status
    blocks(k)%active                = .true.

    ! save coordinates
    blocks(k)%coord_x               = ix
    blocks(k)%coord_y               = iy

    ! update spacing
    blocks(k)%dx                    = abs(ix(2) - ix(1))
    blocks(k)%dy                    = abs(iy(2) - iy(1))

    ! save treecode
    blocks(k)%treecode(1:treeN)     = treecode
    blocks(k)%level                 = treecode_size(treecode)

    ! refinement status
    blocks(k)%refinement            = 0

    ! reset neighbor status
    blocks(k)%neighbor_id           = -1
    blocks(k)%neighbor2_id          = -1
    blocks(k)%neighbor_treecode     = -1
    blocks(k)%neighbor2_treecode    = -1

end subroutine new_block
