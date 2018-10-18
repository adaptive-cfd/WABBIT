!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_operators.f90
!> \version 0.5
!> \author sm
!
!> \brief module for all operator routines
!
!>
!! = log ======================================================================================
!! \n
!! 28/7/17 - create
! *********************************************************************************************

module module_operators

!---------------------------------------------------------------------------------------------
! modules

use mpi
! global parameters
use module_params
! debug module
use module_debug
! use mesh module, since we want to compute origin/spacing of blocks
! use module_mesh, only : get_block_spacing_origin
implicit none


PRIVATE

!**********************************************************************************************
! These are the important routines that are visible to WABBIT:
!**********************************************************************************************
PUBLIC :: compute_vorticity, divergence

!####################
! Derivative Matrices
!####################
! in x-direction
! In the beginning of the program we derive 2 different kinds of
! derivative Matices:
!           * Derivatives which only include inner grid cells
real(kind=rk), private, ALLOCATABLE, save :: Dx(:,:)
!           * Derivatives which include boundary and inner grid cells
real(kind=rk), private, ALLOCATABLE, save :: Dxplus(:,:)
real(kind=rk), private, ALLOCATABLE, save :: Dxminus(:,:)
! in y-direction
! In the beginning of the program we derive 2 different kinds of
! derivative Matices:
!           * Derivatives which only include inner grid cells
real(kind=rk), private, ALLOCATABLE, save :: Dy(:,:)
!           * Derivatives which include boundary and inner grid cells
real(kind=rk), private, ALLOCATABLE, save :: Dyplus(:,:)
real(kind=rk), private, ALLOCATABLE, save :: Dyminus(:,:)

contains

    ! include "volume_integral.f90"
    include "compute_vorticity.f90"
    include "divergence.f90"

    !> this routine inits the matrices used for derivatives
    subroutine initialice_derivatives(boundary_type,Bs,g)
      implicit none
      !-----------------------------------------------
      CHARACTER(LEN=80), INTENT(IN) :: boundary_type
      INTEGER(KIND=ik), INTENT(IN) :: Bs,g
      !-----------------------------------------------
      INTEGER(KIND=ik) :: Q

      Q=Bs+2*g

      if ( .not. ALLOCATED(Dx) ) then
        allocate(Dx(Q**2,Q**2))
      else
        call abort(1110181, "EHHHH: Have allready initialiced derivatives")
      end if




    end subroutine initialice_derivatives


    !> derivative in x direction
    subroutine  D_x(u_x,dx, u)
        external dgemv

        real(kind=rk), intent(in)       :: dx
        real(kind=rk), intent(in)       :: u(:,:)
        real(kind=rk), intent(out)      :: u_x(:,:)

        integer :: i, j ,Q, M, N
        real(kind=rk) :: dx_inv
        real(kind=rk),allocatable,save :: u_vec(:)

        N=size(u,1)
        M=size(u,2)
        Q=M*N

        if (.not. allocated(u_vec)) then
          allocate(u_vec(Q))
        endif



        ! vectorice
        u_vec=reshape(u,(/Q/))
        !invert lattice spacing
        dx_inv = 1.0_rk/dx
        !-------------------------------------------------------------
        ! Compute the derivative:
        !           u_vec= dx_inv Dc*u_vec+0*u_vec
        ! for this step we use the fast blas Implementation
        ! Explanation of the arguments from left to right:
        !         1.'N'= normal Matrix (not "T" transposed)
        !         2. Number of rows of Matrix Dx
        !         3. Number of columns of matrix Dx
        !         4. scaling factor of the matrix vector product
        !         5. Matrix of dimension (LDA,Ncols)
        !         6. LDA is the allocated memory of the 1 dimension of the Matrix
        !         7. vector
        !         8. increment between vector elements
        !         9. scaling factor of the second term
        !         10. Vector of the second term
        !         11. increment of the vector
        call dgemv('N', Q,   Q,   dx_inv,  Dx, M*N,  u_vec, 1,   0.0, u_vec, 1)
        !-------------------------------------------------------------
        u_x=reshape(u_vec,(/N,M/))
    end subroutine D_x

    !> derivative in y direction
    subroutine  D_y(u_y,dy, u)
        external dgemv

        real(kind=rk), intent(in)       :: dy
        real(kind=rk), intent(in)       :: u(:,:)
        real(kind=rk), intent(out)      :: u_y(:,:)

        integer :: i, j ,Q, M, N
        real(kind=rk) :: dy_inv
        real(kind=rk),allocatable,save :: u_vec(:)

        N=size(u,1)
        M=size(u,2)
        Q=M*N

        if (.not. allocated(u_vec)) then
          allocate(u_vec(Q))
        endif
        ! vectorice
        u_vec=reshape(u,(/N*M/))
        !invert lattice spacing
        dy_inv = 1.0_rk/dy
        !-------------------------------------------------------------
        ! Compute the derivative:
        !           u_vec= dy_inv Dc*u_vec+0*u_vec
        ! for this step we use the fast blas Implementation
        ! Explanation of the arguments from left to right:
        !         1.'N'= normal Matrix (not "T" transposed)
        !         2. Number of rows of Matrix Dx
        !         3. Number of columns of matrix Dx
        !         4. scaling factor of the matrix vector product
        !         5. Matrix of dimension (LDA,Ncols)
        !         6. LDA is the allocated memory of the 1 dimension of the Matrix
        !         7. vector
        !         8. increment between vector elements
        !         9. scaling factor of the second term
        !         10. Vector of the second term
        !         11. increment of the vector
        call dgemv('N', Q,   Q,   dy_inv,  Dy, M*N,  u_vec, 1,   0.0, u_vec, 1)
        !-------------------------------------------------------------
        u_y=reshape(u_vec,(/N,M/))
    end subroutine D_y


    subroutine  print_mat(Mat)

        real(kind=rk), intent(in)       :: Mat(:,:)
        integer                         :: i, j
        character(len=16) :: fmt
        character(len=3) :: ncols_str

        write(ncols_str,'(i3.3)') size(Mat,2)
        fmt = '('//ncols_str//'(f9.3,1x))'
        do j = 1,size(Mat,1)
          write(*,fmt) Mat(j,:)
        end do

    end subroutine print_mat

    function  diag( diag_vals, offset) result(D)
        !-----------------------------------------------------
        integer(kind=ik), intent(in),optional  :: offset
        real(kind=rk), intent(in)              :: diag_vals(:)
        real(kind=rk),ALLOCATABLE :: D(:,:)
        !-----------------------------------------------------
        integer(kind=ik) :: k,i,Nd

        k=0
        if ( PRESENT(offset) ) then
          k=offset
        end if

        ! size of the matrix
        Nd=size(diag_vals)+abs(k)
        allocate(D(Nd,Nd))
        D=0.0_rk

        if ( k<0 ) then
          forall(i=1:Nd+k)
            D(i-k,i)=diag_vals(i)
          end forall
        else
          forall(i=1:Nd-k)
            D(i,i+k)=diag_vals(i)
          end forall
        end if

    end function diag



    subroutine  example( Bs, g, dx, u, dudx)
        external dgemm
        external sgeprt

        integer(kind=ik), intent(in)    :: g, Bs
        real(kind=rk), intent(in)       :: dx
        real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g)
        real(kind=rk), intent(out)      :: dudx(Bs+2*g, Bs+2*g)

        integer                         :: i, j
        real(kind=rk)                   :: dx_inv
        real(kind=rk),ALLOCATABLE         :: D(:,:)
        real(kind=rk),ALLOCATABLE         :: A(:,:)
        real(kind=rk)        :: C(5,5)
        real::C_(5,5),A_(5,5),D_(5,5)

        character(len=16) :: fmt
        character(len=3) :: ncols_str
        A=diag((/2.0_rk,2.0_rk, 2.0_rk, 2.0_rk,2.0_rk/))
        D=diag((/2.0_rk, 1.0_rk, 3.0_rk,3.0_rk/),-1)+diag((/2.0_rk, 1.0_rk, 3.0_rk,1.0_rk/),1)
      !  tfm  tfm  rowA colB K  alpha a lda  b  ldb beta c  ldc */
        ! A_=real(A,8)
        ! D_=real(D,8)
        ! C_=real(C,8)
        call dgemm('N', 'N', 5,   5,   5, 1.0_rk,D  ,  5,A, 5,  0.0_rk, C, 5)
        ! A=A_
        ! D=D_
        ! C=C_
        !call sgeprt(5, 5, D, 'a= ')
        call print_mat(A)
        call print_mat(D)
        call print_mat(C)



    end subroutine example



    subroutine  mat_vec_mul( A, v, b)

        real(kind=rk), intent(in)       :: A(:,:)
        real(kind=rk), intent(out)      :: v(:)
        real(kind=rk), intent(out)      :: b(size(A,1))
        integer                         :: i, j,N,M
        real(kind=rk)                   :: tmp
        N =size(A,1)
        M =size(A,2)

        if ( M.ne.size(v) ) then
          call abort(2344,"Error: Dimension of vector and matrix do not fit")
        end if

        do j = 1,N
            tmp=0
            do i = 1,M
              tmp=tmp+A(j,i)*v(i)
            end do
            b(j)=tmp
        end do

    end subroutine mat_vec_mul



end module module_operators
