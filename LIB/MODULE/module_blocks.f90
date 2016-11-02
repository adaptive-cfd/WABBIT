! ********************************
! WABBIT
! --------------------------------
!
! blocks data structure
!
! name: module_blocks.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

module module_blocks

    implicit none

    ! data precision parameters
    integer, parameter, public   :: sngl_prec=selected_real_kind(4)
    integer, parameter, public	 :: dble_prec=selected_real_kind(8)

    integer, parameter, public   :: int_prec=selected_int_kind(8)

    integer, parameter, public   :: rk=dble_prec
    integer, parameter, public   :: ik=int_prec

    ! datatype for block parameter
    type type_blocks_params

        ! start data fields
        real(kind=rk), dimension(:,:), allocatable      :: phi

        ! grid parameter
        integer(kind=ik)                                :: size_domain
        integer(kind=ik)                                :: size_block
        integer(kind=ik)                                :: number_ghost_nodes

        ! number of allocated blocks, data fields and active blocks
        integer(kind=ik)                                :: number_max_blocks
        integer(kind=ik)                                :: number_data_fields

        ! switch for mesh adaption
        logical                                         :: adapt_mesh

        ! number of heavy data fields per process
        integer(kind=ik)                                :: number_max_blocks_data

    end type type_blocks_params

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (type_blocks_params), save                     :: blocks_params
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type type_data_fields

        ! data fields
        real(kind=rk), dimension(:,:), allocatable      :: data_, data_old

        ! data fields for runge kutta stages
        real(kind=rk), dimension(:,:), allocatable      :: k1, k2, k3, k4

        ! wavelet coefficent for datafield dF
        real(kind=rk)                                   :: detail

    end type type_data_fields

    ! datatype for light block data
    type type_blocks_light

        ! integer data
        ! ------------
        ! treecode, level and neighbor information
        integer(kind=ik), dimension(:), allocatable     :: treecode
        integer(kind=ik)                                :: level
        integer(kind=ik), dimension(:,:), allocatable   :: neighbor_treecode
        integer(kind=ik), dimension(16)                 :: neighbor_id

        ! refinement flag
        integer(kind=ik)                                :: refinement

        ! block is assigned to one process at one specific data-id
        integer(kind=ik)                                :: proc_rank
        integer(kind=ik)                                :: proc_data_id

        ! character data
        ! --------------
        ! neighbor dirs
        character(len=3), dimension(16)                  :: neighbor_dir

        ! logical data
        ! ------------
        ! active flag
        logical                                          :: active

    end type type_blocks_light

    ! datatype for heavy block data
    type type_blocks_heavy

        ! data fields
        type(type_data_fields), dimension(:), allocatable    :: data_fields

        ! coordinates, spacing
        real(kind=rk), dimension(:), allocatable        :: coord_x, coord_y
        real(kind=rk)                                   :: dx, dy

        ! block id at light block data
        integer(kind=rk)                                :: block_id

    end type type_blocks_heavy

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (type_blocks_light), save, dimension(:), allocatable :: blocks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (type_blocks_heavy), save, dimension(:), allocatable :: blocks_data
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module module_blocks
