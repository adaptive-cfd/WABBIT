! ********************************
! WABBIT
! --------------------------------
!
! delete light block data
!
! name: delete_block_light.f90
! date: 27.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine delete_block_light(light_id)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)    :: light_id

    ! delete light data
    blocks(light_id)%treecode                      = -1
    blocks(light_id)%level                         = -1
    blocks(light_id)%refinement                    = 0
    blocks(light_id)%neighbor_id                   = -1
    blocks(light_id)%neighbor_dir                  = ""
    blocks(light_id)%neighbor_treecode             = -1
    blocks(light_id)%proc_rank                     = -1
    blocks(light_id)%proc_data_id                  = -1
    blocks(light_id)%active                        = .false.

end subroutine delete_block_light
