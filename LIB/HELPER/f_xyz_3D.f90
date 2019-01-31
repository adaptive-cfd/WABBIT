!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name f_xyz_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief calculate periodic f(x,y) for 3D field
!
!>
!! input:    - params, coordinate arrays \n
!! output:   - f(x,y) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 22/03/17 - create
!
! ********************************************************************************************

subroutine f_xyz_3D( x, y, z, f, Bs, g, Lx, Ly, Lz, frequ )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> grid parameter
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)               :: Lx, Ly, Lz, frequ

    !> coordinate arrays
    real(kind=rk), intent(in)               :: x(Bs(1)+2*g), y(Bs(2)+2*g), z(Bs(3)+2*g)

    !> function array
    real(kind=rk), intent(inout)            :: f(Bs(1)+2*g,Bs(2)+2*g,Bs(3)+2*g)

    ! loop variables
    integer(kind=ik)                        :: i, j, k

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    do i = 1, Bs(1)+2*g
        do j = 1, Bs(2)+2*g
            do k = 1, Bs(3)+2*g

                ! calculate function
                f(i,j,k) = sin( frequ* x(i)/Lx * 2*pi ) * sin( frequ*y(j)/Ly * 2*pi ) * sin( frequ*z(k)/Lz * 2*pi )
                !f(i,j,k) = z(k)

            end do
        end do
    end do

end subroutine f_xyz_3D
