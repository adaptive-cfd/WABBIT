!> \file
!> \callgraph
subroutine allocate_init_debugging(params)
  implicit none
    ! user defined parameter structure
    type (type_params), intent(in)          :: params
  ! ------------------------------------------------------------------------------------------------------
  ! init debug data
  ! note: fix size of time measurements array
  if ( params%debug ) then

      ! allocate array for time measurements - data
      allocate( debug%comp_time( 50, 4 )  )

      ! reset times
      debug%comp_time = 0.0_rk

      ! allocate array for time measurements - names
      allocate( debug%name_comp_time( 320 )  )

      ! reset names
      debug%name_comp_time = "---"

  end if

end subroutine
