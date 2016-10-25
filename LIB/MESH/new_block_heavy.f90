! ********************************
! WABBIT
! --------------------------------
!
! create new block, heavy data
!
! name: new_block_heavy.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine new_block_heavy(proc_k, k, data_, ix, iy, Bs, g, dF)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                :: proc_k, k, Bs, g, dF

    real(kind=rk), dimension(Bs+2*g, Bs+2*g), intent(in)        :: data_
    real(kind=rk), dimension(Bs), intent(in)                    :: ix, iy

    ! error handling
    if ( (k <= 0) .or. (k > blocks_params%number_max_blocks) ) then

      write(*,*) "ERROR! You try to create a block outside of the list"
      write(*,'("your id: ",i8," N_max_blocks=",i8)') k, blocks_params%number_max_blocks
      stop

    endif

    ! save data
    blocks_data(proc_k)%data_fields(dF)%data_       = data_

    ! save coordinates
    blocks_data(proc_k)%coord_x                     = ix
    blocks_data(proc_k)%coord_y                     = iy

    ! update spacing
    blocks_data(proc_k)%dx                          = abs(ix(2) - ix(1))
    blocks_data(proc_k)%dy                          = abs(iy(2) - iy(1))

    ! save light-data id
    blocks_data(proc_k)%block_id                    = k

end subroutine new_block_heavy
