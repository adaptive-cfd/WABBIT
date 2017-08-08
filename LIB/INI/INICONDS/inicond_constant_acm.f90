!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_constant_acm.f90
!> \version 0.5
!> \author sm
!
!> \brief initialize field for acm
!
!>
!! input:    - phi \n
!! output:   - phi \n
!!
!!
!! = log ======================================================================================
!! \n
!! 
!
! ********************************************************************************************

subroutine inicond_constant_acm(phi)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! initial data field
    real(kind=rk), intent(inout)            :: phi(:, :, :, :)


!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! create field phi
    phi = 0.0_rk

    ! it sometimes causes bizarre effects not to delete extremely small numbers:
    ! so we do that now.
    where ( phi<1.0e-13_rk )
       phi = 0.0_rk
    end where

    phi(:,:,:,1) = 1.0_rk

end subroutine inicond_constant_acm
