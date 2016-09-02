! ********************************
! 2D AMR prototype
! --------------------------------
!
! params data structure
!
! name: module_params.f90
! date: 02.08.2016
! author: msr
! version: 0.1
!
! ********************************

module module_params

    use module_blocks

    implicit none

    type type_params

        real(kind=rk) 			    :: time_max
        real(kind=rk) 			    :: CFL
        real(kind=rk), dimension(2)	:: u0
        real(kind=rk)			    :: nu
        real(kind=rk)			    :: Lx, Ly
        integer(kind=ik)			:: write_freq
        real(kind=rk)               :: eps_coarsen
        real(kind=rk)               :: eps_refine
        integer(kind=ik)            :: max_treelevel
        integer(kind=ik)            :: min_treelevel
        character(len=80)           :: order_discretization, order_predictor, inicond

    end type type_params

    type (type_params), save :: params

end module module_params
