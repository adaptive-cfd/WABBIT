! ********************************
! WABBIT
! --------------------------------
!
! params data structure
!
! name: module_params.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

module module_params

    use module_blocks

    implicit none

    type type_params

        real(kind=rk) 			                    :: time_max
        real(kind=rk) 			                    :: CFL
        ! each data field can use separat velocity and diffusion coefficient
        real(kind=rk), dimension(:,:), allocatable	:: u0
        real(kind=rk), dimension(:), allocatable    :: nu
        real(kind=rk)			                    :: Lx, Ly
        integer(kind=ik)			                :: write_freq
        real(kind=rk)                               :: eps
        integer(kind=ik)                            :: max_treelevel
        integer(kind=ik)                            :: min_treelevel
        character(len=80)                           :: order_discretization, order_predictor, inicond, boundary

    end type type_params

    type (type_params), save :: params

end module module_params
