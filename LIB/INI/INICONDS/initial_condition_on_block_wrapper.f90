!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name initial_condition_on_block_wrapper.f90
!> \version 0.5
!> \author engels
!
!> \brief module for all init subroutines
!
!>
!! = log ======================================================================================
!! \n
!! 03 Apr 2017 - create
! ********************************************************************************************

! This routine sets the initial condition on a single, arbitrary block
subroutine initial_condition_on_block_wrapper( params, u, x0, dx, inicond )
  implicit none
  !> user defined parameter structure
  type (type_params), intent(inout)    :: params
  !> actual block data
  real(kind=rk), intent(inout) :: u(:,:,:,:)
  !> spacing and origin of block
  real(kind=rk), intent(in) :: x0(1:3),dx(1:3)
  !> what function to use
  character(len=*), intent(in) :: inicond
  integer(kind=ik)              :: Bs, g

  Bs    = params%number_block_nodes
  g     = params%number_ghost_nodes

  select case( inicond )
  case ("sinus_2d","sinus2d","sin2d")
    call inicond_sinus_2D( params, u, x0, dx )

  case ("zeros")
    u = 0.0_rk

  case ("gauss-blob","gauss_blob","ns_pressure_blob")
    call inicond_gauss_blob( params, u, x0, dx )

  case ("3D_sphere")
    call inicond_sphere( params, u, x0, dx )

  case("constant_acm")
    call inicond_constant_acm( u )

  case ("shear_layer")
    call inicond_shear_layer( params, u, x0, dx )

  case default
    call error_msg("the initial condition is unkown: "//trim(adjustl(inicond)))

  end select


end subroutine initial_condition_on_block_wrapper
