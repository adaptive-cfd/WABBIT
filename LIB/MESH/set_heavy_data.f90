! ********************************
! WABBIT
! --------------------------------
!
! set heavy block data
!
! name: set_heavy_data.f90
! date: 02.11.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine set_heavy_data(k, data_, Bs, g, dF)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                :: k, Bs, g, dF

    real(kind=rk), dimension(Bs+2*g, Bs+2*g), intent(in)        :: data_

    ! error handling
    if ( (k <= 0) .or. (k > blocks_params%number_max_blocks) ) then

      write(*,*) "ERROR! You try to create a block outside of the list"
      write(*,'("your id: ",i8," N_max_blocks=",i8)') k, blocks_params%number_max_blocks
      stop

    endif

    ! save data
    blocks_data(k)%data_fields(dF)%data_       = data_

end subroutine set_heavy_data
