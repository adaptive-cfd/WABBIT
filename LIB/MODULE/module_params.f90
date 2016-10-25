! ********************************
! WABBIT
! --------------------------------
!
! params data structure
!
! name: module_params.f90
! date: 29.09.2016
! author: msr
! version: 0.2
!
! ********************************

module module_params

    use module_blocks

    implicit none

    type type_params

        real(kind=rk) 			        :: time_max
        real(kind=rk) 			        :: CFL
        real(kind=rk), dimension(2)	    :: u0
        real(kind=rk)			        :: nu
        real(kind=rk)			        :: Lx, Ly
        integer(kind=ik)			    :: write_freq
        real(kind=rk)                   :: eps
        integer(kind=ik)                :: max_treelevel
        integer(kind=ik)                :: min_treelevel
        character(len=80)               :: order_discretization, order_predictor, inicond, boundary

    end type type_params

    type (type_params), save :: params

end module module_params
