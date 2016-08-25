! ********************************
! 2D AMR prototype
! --------------------------------
! 
! blocks data structure
!
! name: module_blocks.f90
! date: 02.08.2016
! author: msr
! version: 0.1
! 
! ********************************

module module_blocks

    implicit none

    integer, parameter, public   :: sngl_prec=selected_real_kind(4)
    integer, parameter, public	 :: dble_prec=selected_real_kind(8)

    integer, parameter, public   :: int_prec=selected_int_kind(8)

    integer, parameter, public   :: rk=dble_prec
    integer, parameter, public   :: ik=int_prec

    ! datatype for block parameter
    type type_blocks_params

        real(kind=rk), dimension(:,:), allocatable      :: phi, D1, D2

        integer(kind=ik), dimension(:), allocatable     :: active_list
        integer(kind=ik)                                :: size_domain
        integer(kind=ik)                                :: size_block
        integer(kind=ik)                                :: number_max_blocks
        integer(kind=ik)                                :: number_ghost_nodes

    end type type_blocks_params

    type (type_blocks_params), save                     :: blocks_params

    ! datatype for block data
    type type_blocks

        real(kind=rk), dimension(:,:), allocatable      :: data1, data2, data_old, k1, k2, k3, k4
        real(kind=rk), dimension(:), allocatable        :: coord_x, coord_y
        real(kind=rk)                                   :: dx, dy

        integer(kind=ik), dimension(:), allocatable     :: treecode
        integer(kind=ik)                                :: level
        integer(kind=ik), dimension(:,:), allocatable   :: neighbor_treecode
        character(len=2), dimension(8)                  :: neighbor_dir
        integer(kind=ik), dimension(8)                  :: neighbor_id
        integer(kind=ik), dimension(:,:), allocatable   :: neighbor2_treecode
        character(len=2), dimension(4)                  :: neighbor2_dir
        integer(kind=ik), dimension(4)                  :: neighbor2_id
        integer(kind=ik), dimension(8)                  :: neighbor_number

        integer(kind=ik)                                :: refinement

        logical                                         :: active

    end type type_blocks

    type (type_blocks), save, dimension(:), allocatable :: blocks

end module module_blocks
