! ********************************
! 2D AMR prototype
! --------------------------------
!
! this is a wrapper for different
! initial conditions
!
! name: initial_condition_dense_field.f90
! date: 31.08.2016
! author: engels, msr
! version: 0.1
!
! ********************************

subroutine initial_condition_dense_field()

      use module_params
      use module_blocks

      implicit none

      select case( params%inicond )

          case ("gauss-blob","gauss_blob")
              call inicond_gauss_blob()

          case ("sinus")
              call inicond_sinus()

          case default
              write(*,*) "params%inicond is unkown"
              write(*,*) params%inicond
              stop

      end select

end subroutine initial_condition_dense_field
