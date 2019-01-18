!-----------------------------------------------------------------
!>
!> \details
!> \version 23.2.2018
!> \author P.Krah
!-----------------------------------------------------------------

module module_sparse_operators

use mpi
! global parameters
use module_params

#ifdef SBLAS
use blas_sparse

implicit none


PRIVATE
!**********************************************************************************************
! These are the important routines that are visible to WABBIT:
!**********************************************************************************************
PUBLIC :: initialice_derivatives, DUS_Dx, DUS_Dy
!####################
! Derivative Matrices
!####################
! in x-direction
! In the beginning of the program we derive 2 different kinds of
! derivative Matices:
!           * Derivatives which only include inner grid cells
integer, private, save :: diffx_handl
!           * Derivatives which include boundary and inner grid cells
integer, private, save :: diffxplus_handl
integer, private, save :: diffxminus_handl
! in y-direction
! In the beginning of the program we derive 2 different kinds of
! derivative Matices:
!           * Derivatives which only include inner grid cells
integer, private, save :: diffy_handl
!           * Derivatives which include boundary and inner grid cells
integer, private, save :: diffyplus_handl
integer, private, save :: diffyminus_handl

!###############################################
! Interfaces for sparse and full representation
!###############################################
!> namingconvention like it is done in blas libary:
!>        + D  = Double
!>        + US = Unstructured sparse
!>        + GE = General
interface diag
    module procedure DGE_DIAG, DUS_DIAG
end interface diag




contains


    !> this routine inits the matrices used for derivatives
    subroutine initialice_derivatives(boundary_type,Bs,g)
      implicit none
      !-----------------------------------------------
      CHARACTER(LEN=80),INTENT(IN),dimension(3) :: boundary_type
      INTEGER(KIND=ik), INTENT(IN) :: Bs,g
      !-----------------------------------------------
      INTEGER :: Q,ierr,i
      real(kind=rk), allocatable:: D_VAL(:),  EYE_VAL(:)
      integer, allocatable      :: D_INDX(:), EYE_INDX(:)
      integer, allocatable      :: D_JNDX(:), EYE_JNDX(:)
      integer                   :: D_SHAPE(2), EYE_SHAPE(2), Dx_SHAPE(2), Dy_SHAPE(2)
      real(kind=rk), allocatable:: Dx_VAL(:),  Dxminus_VAL(:),  Dxplus_VAL(:)
      integer, allocatable      :: Dx_INDX(:), Dxminus_INDX(:), Dxplus_INDX(:)
      integer, allocatable      :: Dx_JNDX(:), Dxminus_JNDX(:), Dxplus_JNDX(:)
      real(kind=rk), allocatable:: Dy_VAL(:),  Dyminus_VAL(:),  Dyplus_VAL(:)
      integer, allocatable      :: Dy_INDX(:), Dyminus_INDX(:), Dyplus_INDX(:)
      integer, allocatable      :: Dy_JNDX(:), Dyminus_JNDX(:), Dyplus_JNDX(:)

      !--------------------------------
      ! DERIVATIVE in x direction
      !--------------------------------
      Q=Bs+2*g
      allocate(EYE_VAL(Q),EYE_INDX(Q),EYE_JNDX(Q))
      EYE_VAL=1



      select case (boundary_type(1))
      case ('periodic')

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                        D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                        Dxminus_VAL,Dxminus_INDX,Dxminus_JNDX,Dx_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      Dx_VAL,Dx_INDX,Dx_JNDX,Dx_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                        D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                        Dxplus_VAL,Dxplus_INDX,Dxplus_JNDX,Dx_SHAPE)
      case ('symmetric-open')
        ! The matrix of the full Block with 10 points in x direction would
        ! look like:
        !           /-8     9    -1     0     0     0     0     0     0     0 \
        !          | -7     0     8    -1     0     0     0     0     0     0 |
        !          |  1    -8     0     8    -1     0     0     0     0     0 |
        !          |  0     1    -8     0     8    -1     0     0     0     0 |
        !          |  0     0     1    -8     0     8    -1     0     0     0 |
        ! Dx =1/12 |  0     0     0     1    -8     0     8    -1     0     0 |
        !          |  0     0     0     0     1    -8     0     8    -1     0 |
        !          |  0     0     0     0     0     1    -8     0     8    -1 |
        !          |  0     0     0     0     0    -1     6   -18    10     3 |
        !          \  0     0     0     0     0     3   -16    36   -48    25 /
        ! We split it up in 3 different matrix:
        !     * Dxminus: containing boundary and inner grid points
        !     * Dx     : containing  ONLY        inner grid points
        !     * Dxplus : containing boundary and inner grid points

        !                /- 8     9   - 1     0     0     0  \
        !                | -7     0     8    -1     0     0  |
        ! Dxminus= 1/12  |  1    -8     0     8    -1     0  |
        !                |  0     1    -8     0     8    -1  |
        !                |  0     0     1    -8     0     8  |
        !                \  0     0     0     1    -8     0  /
        !          / 0     8    -1     0     0     0     0\
        !          |-8     0     8    -1     0     0     0|
        ! Dx= 1/12 | 1    -8     0     8    -1     0     0|
        !          | 0     1    -8     0     8    -1     0|
        !          | 0     0     1    -8     0     8    -1|
        !          | 0     0     0     1    -8     0     8|
        !          \ 0     0     0     0     1    -8     0/

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        Dxminus_VAL,Dxminus_INDX,Dxminus_JNDX,Dx_SHAPE)

        !          / 0     8    -1     0     0     0     0\
        !          |-8     0     8    -1     0     0     0|
        ! Dx= 1/12 | 1    -8     0     8    -1     0     0|
        !          | 0     1    -8     0     8    -1     0|
        !          | 0     0     1    -8     0     8    -1|
        !          | 0     0     0     1    -8     0     8|
        !          \ 0     0     0     0     1    -8     0/

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        Dx_VAL,Dx_INDX,Dx_JNDX,Dx_SHAPE)


        !                /  0    8    -1     0     0     0     0     0    0   0\
        !                | -8    0     8    -1     0     0     0     0    0   0|
        !                |  1   -8     0     8    -1     0     0     0    0   0|
        !                |  0    1    -8     0     8    -1     0     0    0   0|
        ! Dxplus= 1/12   |  0    0     1    -8     0     8    -1     0    0   0|
        !                |  0    0     0     1    -8     0     8    -1    0   0|
        !                |  0    0     0    -1     6   -18    10     3    0   0|
        !                |  0    0     0     3   -16    36   -48    25    0   0|
        !                |  0    0     8     0     0     0     0     0    0   0|
        !                \  0    0     8     0     0     0     0     0    0   0/
        call D_openFourthOrdplus(g,Bs,D_VAL,D_INDX,D_JNDX,D_SHAPE)
      !  call print_mat(DUS_full(d_val,d_indx,d_jndx,d_shape))
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      Dxplus_VAL,Dxplus_INDX,Dxplus_JNDX,Dx_SHAPE)

      case ('open')
        ! The matrix of the full Block with for example 10 points in x direction would
        ! look like:
        !           /-25    48   -36    16    -3     0     0     0     0     0 \
        !          |  -3   -10    18    -6     1     0     0     0     0     0 |
        !          |   1    -8     0     8    -1     0     0     0     0     0 |
        !          |   0     1    -8     0     8    -1     0     0     0     0 |
        !          |   0     0     1    -8     0     8    -1     0     0     0 |
        ! Dx =1/12 |   0     0     0     1    -8     0     8    -1     0     0 |
        !          |   0     0     0     0     1    -8     0     8    -1     0 |
        !          |   0     0     0     0     0     1    -8     0     8    -1 |
        !          |   0     0     0     0     0    -1     6   -18    10     3 |
        !          \   0     0     0     0     0     3   -16    36   -48    25 /
        ! We split it up in 3 different matrix:
        !     * Dxminus: containing boundary and inner grid points
        !     * Dx     : containing  ONLY        inner grid points
        !     * Dxplus : containing boundary and inner grid points

        !                /  0    0     0     0     0     0     0     0  \
        !                |  0    0     0     0     0     0     0     0  |
        !                |  0    0   -25    48   -36    16    -3     0  \
        !                |  0    0    -3   -10    18    -6     1     0  |   example for g=2 and Bs=4
        ! Dxminus= 1/12  |  0    0     1    -8     0     8    -1     0  |
        !                |  0    0     0     1    -8     0     8    -1  |
        !                |  0    0     0     0     1    -8     0     8  |
        !                \  0    0     0     0     0     1    -8     0  /
        call D_openFourthOrdminus(g,Bs,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        !call print_mat(DUS_full(d_val,d_indx,d_jndx,d_shape))
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        Dxminus_VAL,Dxminus_INDX,Dxminus_JNDX,Dx_SHAPE)

        !          / 0     8    -1     0     0     0     0\
        !          |-8     0     8    -1     0     0     0|
        ! Dx= 1/12 | 1    -8     0     8    -1     0     0|
        !          | 0     1    -8     0     8    -1     0|
        !          | 0     0     1    -8     0     8    -1|
        !          | 0     0     0     1    -8     0     8|
        !          \ 0     0     0     0     1    -8     0/

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        Dx_VAL,Dx_INDX,Dx_JNDX,Dx_SHAPE)
        !                /  0    8    -1     0     0     0     0     0    0   0\
        !                | -8    0     8    -1     0     0     0     0    0   0|
        !                |  1   -8     0     8    -1     0     0     0    0   0|
        !                |  0    1    -8     0     8    -1     0     0    0   0|
        ! Dxplus= 1/12   |  0    0     1    -8     0     8    -1     0    0   0|
        !                |  0    0     0     1    -8     0     8    -1    0   0|
        !                |  0    0     0    -1     6   -18    10     3    0   0|
        !                |  0    0     0     3   -16    36   -48    25    0   0|
        !                |  0    0     8     0     0     0     0     0    0   0|
        !                \  0    0     8     0     0     0     0     0    0   0/
        call D_openFourthOrdplus(g,Bs,D_VAL,D_INDX,D_JNDX,D_SHAPE)
      !  call print_mat(DUS_full(d_val,d_indx,d_jndx,d_shape))
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      Dxplus_VAL,Dxplus_INDX,Dxplus_JNDX,Dx_SHAPE)
      case default
        call abort(2510182,"Commander, there is a leak in the fuel tanks!! Emergency Landing! ")
      end select
      ! this function creates the handle to the different matrices
      diffxminus_handl  = create_blas_handl(Dxminus_VAL,Dxminus_INDX,Dxminus_JNDX,Dx_SHAPE)
      diffx_handl       = create_blas_handl(Dx_VAL,Dx_INDX,Dx_JNDX,Dx_SHAPE)
      diffxplus_handl   = create_blas_handl(Dxplus_VAL,Dxplus_INDX,Dxplus_JNDX,Dx_SHAPE)


      !--------------------------------
      ! DERIVATIVE in y direction
      !--------------------------------
      select case (boundary_type(2))
      case ('periodic')

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      Dyminus_VAL,Dyminus_INDX,Dyminus_JNDX,Dx_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      Dy_VAL,Dy_INDX,Dy_JNDX,Dy_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      Dyplus_VAL,Dyplus_INDX,Dyplus_JNDX,Dy_SHAPE)
      case ('symmetric-open')
        ! The matrix of the full Block with 10 points in x direction would
        ! look like:
        !           /-8     9    -1     0     0     0     0     0     0     0 \
        !          | -7     0     8    -1     0     0     0     0     0     0 |
        !          |  1    -8     0     8    -1     0     0     0     0     0 |
        !          |  0     1    -8     0     8    -1     0     0     0     0 |
        !          |  0     0     1    -8     0     8    -1     0     0     0 |
        ! Dy =1/12 |  0     0     0     1    -8     0     8    -1     0     0 |
        !          |  0     0     0     0     1    -8     0     8    -1     0 |
        !          |  0     0     0     0     0     1    -8     0     8    -1 |
        !          |  0     0     0     0     0    -1     6   -18    10     3 |
        !          \  0     0     0     0     0     3   -16    36   -48    25 /
        ! We split it up in 3 different matrix:
        !     * Dyminus: containing boundary and inner grid points
        !     * Dy     : containing  ONLY        inner grid points
        !     * Dyplus : containing boundary and inner grid points

        !                /- 8     9   - 1     0     0     0  \
        !                | -7     0     8    -1     0     0  |
        ! Dyminus= 1/12  |  1    -8     0     8    -1     0  |
        !                |  0     1    -8     0     8    -1  |
        !                |  0     0     1    -8     0     8  |
        !                \  0     0     0     1    -8     0  /
        !          / 0     8    -1     0     0     0     0\
        !          |-8     0     8    -1     0     0     0|
        ! Dy= 1/12 | 1    -8     0     8    -1     0     0|
        !          | 0     1    -8     0     8    -1     0|
        !          | 0     0     1    -8     0     8    -1|
        !          | 0     0     0     1    -8     0     8|
        !          \ 0     0     0     0     1    -8     0/

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        Dyminus_VAL,Dyminus_INDX,Dyminus_JNDX,Dy_SHAPE)

        !          / 0     8    -1     0     0     0     0\
        !          |-8     0     8    -1     0     0     0|
        ! Dy= 1/12 | 1    -8     0     8    -1     0     0|
        !          | 0     1    -8     0     8    -1     0|
        !          | 0     0     1    -8     0     8    -1|
        !          | 0     0     0     1    -8     0     8|
        !          \ 0     0     0     0     1    -8     0/

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        Dy_VAL,Dy_INDX,Dy_JNDX,Dy_SHAPE)


        !                /  0    8    -1     0     0     0     0     0    0   0\
        !                | -8    0     8    -1     0     0     0     0    0   0|
        !                |  1   -8     0     8    -1     0     0     0    0   0|
        !                |  0    1    -8     0     8    -1     0     0    0   0|
        ! Dyplus= 1/12   |  0    0     1    -8     0     8    -1     0    0   0|
        !                |  0    0     0     1    -8     0     8    -1    0   0|
        !                |  0    0     0    -1     6   -18    10     3    0   0|
        !                |  0    0     0     3   -16    36   -48    25    0   0|
        !                |  0    0     8     0     0     0     0     0    0   0|
        !                \  0    0     8     0     0     0     0     0    0   0/
        call D_openFourthOrdplus(g,Bs,D_VAL,D_INDX,D_JNDX,D_SHAPE)
      !  call print_mat(DUS_full(d_val,d_indx,d_jndx,d_shape))
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      Dyplus_VAL,Dyplus_INDX,Dyplus_JNDX,Dy_SHAPE)

      case ('open', 'symmetryAxis-wall')
        ! The matrix of the full Block with for example 10 points in x direction would
        ! look like:
        !           /-25    48   -36    16    -3     0     0     0     0     0 \
        !          |  -3   -10    18    -6     1     0     0     0     0     0 |
        !          |   1    -8     0     8    -1     0     0     0     0     0 |
        !          |   0     1    -8     0     8    -1     0     0     0     0 |
        !          |   0     0     1    -8     0     8    -1     0     0     0 |
        ! Dy =1/12 |   0     0     0     1    -8     0     8    -1     0     0 |
        !          |   0     0     0     0     1    -8     0     8    -1     0 |
        !          |   0     0     0     0     0     1    -8     0     8    -1 |
        !          |   0     0     0     0     0    -1     6   -18    10     3 |
        !          \   0     0     0     0     0     3   -16    36   -48    25 /
        ! We split it up in 3 different matrix:
        !     * Dyminus: containing boundary and inner grid points
        !     * Dy     : containing  ONLY        inner grid points
        !     * Dyplus : containing boundary and inner grid points

        !                /  0    0     0     0     0     0     0     0  \
        !                |  0    0     0     0     0     0     0     0  |
        !                |  0    0   -25    48   -36    16    -3     0  \
        !                |  0    0    -3   -10    18    -6     1     0  |   example for g=2 and Bs=4
        ! Dyminus= 1/12  |  0    0     1    -8     0     8    -1     0  |
        !                |  0    0     0     1    -8     0     8    -1  |
        !                |  0    0     0     0     1    -8     0     8  |
        !                \  0    0     0     0     0     1    -8     0  /
        call D_openFourthOrdminus(g,Bs,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        !call print_mat(DUS_full(d_val,d_indx,d_jndx,d_shape))
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        Dyminus_VAL,Dyminus_INDX,Dyminus_JNDX,Dy_SHAPE)

        !          / 0     8    -1     0     0     0     0\
        !          |-8     0     8    -1     0     0     0|
        ! Dy= 1/12 | 1    -8     0     8    -1     0     0|
        !          | 0     1    -8     0     8    -1     0|
        !          | 0     0     1    -8     0     8    -1|
        !          | 0     0     0     1    -8     0     8|
        !          \ 0     0     0     0     1    -8     0/

        call D_FourthOrd(Q,D_VAL,D_INDX,D_JNDX,D_SHAPE)
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
        EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
        Dy_VAL,Dy_INDX,Dy_JNDX,Dy_SHAPE)
        !                /  0    8    -1     0     0     0     0     0    0   0\
        !                | -8    0     8    -1     0     0     0     0    0   0|
        !                |  1   -8     0     8    -1     0     0     0    0   0|
        !                |  0    1    -8     0     8    -1     0     0    0   0|
        ! Dyplus= 1/12   |  0    0     1    -8     0     8    -1     0    0   0|
        !                |  0    0     0     1    -8     0     8    -1    0   0|
        !                |  0    0     0    -1     6   -18    10     3    0   0|
        !                |  0    0     0     3   -16    36   -48    25    0   0|
        !                |  0    0     8     0     0     0     0     0    0   0|
        !                \  0    0     8     0     0     0     0     0    0   0/
        call D_openFourthOrdplus(g,Bs,D_VAL,D_INDX,D_JNDX,D_SHAPE)
      !  call print_mat(DUS_full(d_val,d_indx,d_jndx,d_shape))
        call diag(EYE_VAL,EYE_INDX, EYE_JNDX, EYE_SHAPE)
        call DUS_kron(D_VAL,D_INDX,D_JNDX,D_SHAPE,&
                      EYE_VAL,EYE_INDX,EYE_JNDX,EYE_SHAPE,&
                      Dyplus_VAL,Dyplus_INDX,Dyplus_JNDX,Dy_SHAPE)
      case default
        call abort(2510182,"Commander, there is a leak in the fuel tanks!! Emergency Landing! ")
      end select
      !  call print_mat(DUS_full(dx_val,dx_indx,dx_jndx,dx_shape))

      diffyminus_handl = create_blas_handl(Dyminus_VAL,Dyminus_INDX,Dyminus_JNDX,Dy_SHAPE)
      diffy_handl = create_blas_handl(Dy_VAL,Dy_INDX,Dy_JNDX,Dy_SHAPE)
      diffyplus_handl = create_blas_handl(Dyplus_VAL,Dyplus_INDX,Dyplus_JNDX,Dy_SHAPE)





    end subroutine initialice_derivatives



    !> this routine calls the (CR)eate functions of BLAS
    function create_blas_handl(D_VAL,D_INDX,D_JNDX,D_SHAPE) result(D_handl)
      implicit none
      !-----------------------------------------------
      real(kind=rk),intent(in):: D_VAL(:)         !< sparse values
      integer, intent(in)  :: D_INDX(:),D_JNDX(:) !< row and col index of sparse value
      integer, intent(in)  :: D_SHAPE(2)          !< array of matrix shape: nr of rows and cols
      integer :: D_handl                          !< blas handl
      !-----------------------------------------------
      integer :: ierr
!     ----------------------------------
!     Step 1:  Create Sparse BLAS handle
!     ----------------------------------
      CALL DUSCR_BEGIN( D_SHAPE(1), D_SHAPE(2), D_handl, ierr)
!     -----------------------------------
!     Step 2:  Insert entries all at once
!     -----------------------------------
      CALL USCR_INSERT_ENTRIES(D_handl,D_VAL,D_INDX,D_JNDX, ierr)
!     -----------------------------------------------
!     Step 3:  Complete construction of sparse matriX
!     -----------------------------------------------
      CALL USCR_END(D_handl, ierr)

    end function create_blas_handl




    !> This function creates the sparse matrix stencil for periodic BC
    !> It should be only called once in the beginning of the programm to construct
    !> the matrix! Note: It is not optimized for calling it in every iteration!
    !> q_x *delta_x= 1/12 * ( u(ix-2) - 8u(ix-1) + 8u(ix+1) - u(ix+2) )
    subroutine  D_FourthOrd(N,D_VAL,D_INDX,D_JNDX,D_SHAPE)

      integer(kind=ik), intent(in)    :: N
      real(kind=rk), allocatable, intent(out):: D_VAL(:)
      integer, allocatable, intent(out)      :: D_INDX(:)
      integer, allocatable, intent(out)      :: D_JNDX(:)
      integer, intent(out)      :: D_SHAPE(2)

      integer :: nmax,k,i,j

      ! Calculate the matrix size:
      D_SHAPE=(/N, N/)

      ! First we have to clarify how many non zero elements do we have:
      ! It is a tridiagonal matrix and we can estimate the number
      ! of non zero elements by hand.
      ! We have to off diagonals on the upper and lower side:
      nmax= 2 * ( N - 1 ) + 2 * ( N - 2 )

      ! allocate the sparse matrix
      allocate(D_val(nmax),D_INDX(nmax),D_JNDX(nmax))



      ! we now loop over all the rows and columns of this N times N Matrix
      ! to set the desired values for the stencil
      k=1
      do i= 1, N
        do j= 1, N
          if ( j == i + 1 ) then
            D_INDX(k)=i
            D_JNDX(k)=i+1
            D_VAL(k)=8.0_rk/12.0_rk
            k=k+1
          elseif ( j == i - 1 ) then
            D_INDX(k)=i
            D_JNDX(k)=i-1
            D_VAL(k)=-8.0_rk/12.0_rk
            k=k+1
          elseif ( j == i + 2) then
            D_INDX(k)=i
            D_JNDX(k)=i + 2
            D_VAL(k)=-1/12.0_rk
            k=k+1
          elseif ( j == i - 2 ) then
            D_INDX(k)=i
            D_JNDX(k)=i - 2
            D_VAL(k)=1/12.0_rk
            k=k+1
          else
            ! D_val(i)=0
            ! matrix doesn`t save zeros in sparse representation
          endif
        end do
      end do

      if ( k-1 .ne. nmax ) then
        call abort(241020181,"Dear passengers, we need to do an emergency landing! ")
      end if

    end subroutine D_FourthOrd

    !> This function creates the sparse matrix stencil for open boundary CONDITIONS
    !> It should be only called once in the beginning of the programm to construct
    !> the matrix! Note: It is not optimized for calling it in every iteration!
    !> Example Matrix structure for g=3 and Bs=4
    !                /0   0    0     0     0     0     0     0     0  \
    !                |0   0    0     0     0     0     0     0     0  |
    !                |0   0    0     0     0     0     0     0     0  |
    !                |0   0    0   -25    48   -36    16    -3     0  \
    !                |0   0    0    -3   -10    18    -6     1     0  |
    ! Dxminus= 1/12  |0   0    0     1    -8     0     8    -1     0  |
    !                |0   0    0     0     1    -8     0     8    -1  |
    !                |0   0    0     0     0     1    -8     0     8  |
    !                \0   0    0     0     0     0     1    -8     0  /
    subroutine  D_openFourthOrdminus(g,Bs,D_VAL,D_INDX,D_JNDX,D_SHAPE)

      integer(kind=ik), intent(in)    :: g,Bs
      real(kind=rk), allocatable, intent(out):: D_VAL(:)
      integer, allocatable, intent(out)      :: D_INDX(:)
      integer, allocatable, intent(out)      :: D_JNDX(:)
      integer, intent(out)      :: D_SHAPE(2)

      real(kind=rk), ALLOCATABLE:: val_tmp(:)
      integer, ALLOCATABLE      :: indx_tmp(:)
      integer, ALLOCATABLE      :: jndx_tmp(:)
      integer :: nmax,k,i,j,N

      N=Bs+2*g
      ! Calculate the matrix size:
      D_SHAPE=(/N, N/)

      nmax= N*N
      ! allocate the sparse matrix
      allocate(val_tmp(nmax),indx_tmp(nmax),jndx_tmp(nmax))

      indx_tmp=1
      jndx_tmp=1
      val_tmp=0

      ! we now loop over all the rows and columns of this N times N Matrix
      ! to set the desired values for the stencil
      k=1
      do i= g+3, N! skip the first two rows in the inner block grid
                  ! this is necessary to adopt the stencil to a right sided stencil (see below)
        do j= g+1, N

            if ( j == i + 1 ) then
              indx_tmp(k)=i
              jndx_tmp(k)=i+1
              val_tmp(k)=8.0_rk/12.0_rk
              k=k+1
            elseif ( j == i - 1 ) then
              indx_tmp(k)=i
              jndx_tmp(k)=i-1
              val_tmp(k)=-8.0_rk/12.0_rk
              k=k+1
            elseif ( j == i + 2) then
              indx_tmp(k)=i
              jndx_tmp(k)=i + 2
              val_tmp(k)=-1/12.0_rk
              k=k+1
            elseif ( j == i - 2 ) then
              indx_tmp(k)=i
              jndx_tmp(k)=i - 2
              val_tmp(k)=1/12.0_rk
              k=k+1
            else
              ! D_val(i)=0
              ! matrix doesn`t save zeros in sparse representation
            endif
        end do
      end do

      ! initialice the first two rows of the matrix for the right sided stencil
      val_tmp(k:k+9) = 1/12.0_rk * (/-25,  48, -36,  16, -3, &
                      -3,  -10,  18,  -6, 1/)
      indx_tmp(k:k+9)=(/g+1, g+1, g+1, g+1, g+1, &
                       g+2, g+2, g+2, g+2, g+2/)
      jndx_tmp(k:k+9)=(/g+1, g+2, g+3, g+4, g+5, &
                       g+1, g+2, g+3, g+4, g+5/)
      k=k+9

      allocate(D_VAL(k),D_INDX(k),D_JNDX(k))

      D_VAL=val_tmp(1:k)
      D_INDX=indx_tmp(1:k)
      D_JNDX=jndx_tmp(1:k)

    end subroutine D_openFourthOrdminus





    !> This function creates the sparse matrix stencil for open boundary CONDITIONS
    !> It should be only called once in the beginning of the programm to construct
    !> the matrix! Note: It is not optimized for calling it in every iteration!
    !> Example Matrix structure for g=2 and Bs=4
    !                /  0    8    -1     0     0     0     0     0    0   0\
    !                | -8    0     8    -1     0     0     0     0    0   0|
    !                |  1   -8     0     8    -1     0     0     0    0   0|
    !                |  0    1    -8     0     8    -1     0     0    0   0|
    ! Dxplus= 1/12   |  0    0     1    -8     0     8    -1     0    0   0|
    !                |  0    0     0     1    -8     0     8    -1    0   0|
    !                |  0    0     0    -1     6   -18    10     3    0   0|
    !                |  0    0     0     3   -16    36   -48    25    0   0|
    !                |  0    0     8     0     0     0     0     0    0   0|
    !                \  0    0     8     0     0     0     0     0    0   0/

    subroutine  D_openFourthOrdplus(g,Bs,D_val,D_INDX,D_JNDX,D_SHAPE)

      integer(kind=ik), intent(in)    :: g,Bs
      real(kind=rk), allocatable, intent(out):: D_VAL(:)
      integer, allocatable, intent(out)      :: D_INDX(:)
      integer, allocatable, intent(out)      :: D_JNDX(:)
      integer, intent(out)      :: D_SHAPE(2)

      real(kind=rk), ALLOCATABLE:: val_tmp(:)
      integer, ALLOCATABLE      :: indx_tmp(:)
      integer, ALLOCATABLE      :: jndx_tmp(:)
      integer :: nmax,k,i,j,N

      N=Bs+2*g

      if ( g<=2 ) then
        call abort(211181,"Error: number of ghost nodes needs to be increased")
      end if
      ! Calculate the matrix size:
      D_SHAPE=(/N, N/)

      nmax= N*N
      ! allocate the sparse matrix
      allocate(val_tmp(nmax),indx_tmp(nmax),jndx_tmp(nmax))


      indx_tmp=1
      jndx_tmp=1
      val_tmp=0
      ! we now loop over all the rows and columns of this N times N Matrix
      ! to set the desired values for the stencil
      k=1
      do i= g-2, Bs+g-2 !skip the last two rows before the ghostlayer. They are initialiced underneath
        do j= g-2, Bs + g
            if ( j == i + 1 ) then
              indx_tmp(k)=i
              jndx_tmp(k)=i+1
              val_tmp(k)=8.0_rk/12.0_rk
              k=k+1
            elseif ( j == i - 1 ) then
              indx_tmp(k)=i
              jndx_tmp(k)=i-1
              val_tmp(k)=-8.0_rk/12.0_rk
              k=k+1
            elseif ( j == i + 2) then
              indx_tmp(k)=i
              jndx_tmp(k)=i + 2
              val_tmp(k)=-1/12.0_rk
              k=k+1
            elseif ( j == i - 2 ) then
              indx_tmp(k)=i
              jndx_tmp(k)=i - 2
              val_tmp(k)=1/12.0_rk
              k=k+1
            else
              ! D_val(i)=0
              ! matrix doesn`t save zeros in sparse representation
            endif
        end do
      end do

      val_tmp(k:k+9) = -1/12.0_rk * (/-25,  48, -36,  16, -3, &
                      -3,  -10,  18,  -6, 1/)
      indx_tmp(k:k+9)=(/Bs+g, Bs+g, Bs+g, Bs+g, Bs+g, &
                       Bs+g-1, Bs+g-1, Bs+g-1, Bs+g-1, Bs+g-1/)
      jndx_tmp(k:k+9)=(/Bs+g, Bs+g-1, Bs+g-2, Bs+g-3, Bs+g-4, &
                        Bs+g, Bs+g-1, Bs+g-2, Bs+g-3, Bs+g-4/)
      k=k+9

      allocate(D_VAL(k),D_INDX(k),D_JNDX(k))

      D_VAL=val_tmp(1:k)
      D_INDX=indx_tmp(1:k)
      D_JNDX=jndx_tmp(1:k)
    end subroutine D_openFourthOrdplus








    !> derivative in x direction in sparse blas representation
    subroutine  DUS_Dx(u_x, dx, u, boundary)

        real(kind=rk), intent(in)       :: dx
        real(kind=rk), intent(in)       :: u(:,:)
        real(kind=rk), intent(out)      :: u_x(:,:)
        integer(kind=2), intent(in)     :: boundary

        integer :: i, j ,Q, M, N,ierr
        real(kind=rk) :: dx_inv
        real(kind=rk),allocatable:: u_vec(:),u_x_vec(:)

        N=size(u,1)
        M=size(u,2)
        Q=M*N


        if (.not. allocated(u_vec)) then
          allocate(u_vec(Q))
        endif
        if (.not. allocated(u_x_vec)) then
          allocate(u_x_vec(Q))
        endif
        ! vectorice
        u_vec=reshape(u,(/Q/))
        u_x_vec=0.0_rk

        !invert lattice spacing
        dx_inv = 1.0_rk/dx
        !-------------------------------------------------------------
        ! Compute the derivative:
        !           u_vec= dx_inv Dc*u_vec+0*u_vec
        if ( boundary==-1 ) then
          call usmv( diffxminus_handl, u_vec, u_x_vec, ierr, alpha=dx_inv)
        elseif (boundary==1) then
          call usmv( diffxplus_handl, u_vec, u_x_vec, ierr, alpha=dx_inv)
        else
          call usmv( diffx_handl, u_vec, u_x_vec, ierr, alpha=dx_inv)
        end if
        !-------------------------------------------------------------
        u_x=reshape(u_x_vec,(/N,M/))

        !write(*,*) "norm=",SQRT(DOT_PRODUCT(u_x_vec,u_x_vec)),SQRT(DOT_PRODUCT(u_vec,u_vec))


    end subroutine DUS_Dx


    !> derivative in y direction in sparse blas representation
    subroutine  DUS_Dy(u_y, dy, u, boundary)

        real(kind=rk), intent(in)       :: dy
        real(kind=rk), intent(in)       :: u(:,:)
        real(kind=rk), intent(out)      :: u_y(:,:)
        integer(kind=2), intent(in)     :: boundary


        integer :: i, j ,Q, M, N,ierr
        real(kind=rk) :: dy_inv
        real(kind=rk),allocatable:: u_vec(:),u_y_vec(:)

        N=size(u,1)
        M=size(u,2)
        Q=M*N


        if (.not. allocated(u_vec)) then
          allocate(u_vec(Q))
        endif
        if (.not. allocated(u_y_vec)) then
          allocate(u_y_vec(Q))
        endif
        ! vectorice
        u_vec=reshape(u,(/Q/))
        u_y_vec=0.0_rk

        !invert lattice spacing
        dy_inv = 1.0_rk/dy
        !-------------------------------------------------------------
        ! Compute the derivative:
        !           u_vec= dy_inv Dc*u_vec+0*u_vec
        if ( boundary==-1 ) then
          call usmv( diffyminus_handl, u_vec, u_y_vec, ierr, alpha=dy_inv)
        elseif (boundary==1) then
          call usmv( diffyplus_handl, u_vec, u_y_vec, ierr, alpha=dy_inv)
        else
          call usmv( diffy_handl, u_vec, u_y_vec, ierr, alpha=dy_inv)
        end if
        !-------------------------------------------------------------
        u_y=reshape(u_y_vec,(/N,M/))

    end subroutine DUS_Dy

    !> derivative in x direction
    ! subroutine  DGE_Dx(u_x,dx, u)
    !     external dgbmv
    !
    !     real(kind=rk), intent(in)       :: dx
    !     real(kind=rk), intent(in)       :: u(:,:)
    !     real(kind=rk), intent(out)      :: u_x(:,:)
    !
    !     integer :: i, j ,Q, M, N
    !     real(kind=rk) :: dx_inv
    !     real(kind=rk),allocatable,save :: u_vec(:)
    !
    !     N=size(u,1)
    !     M=size(u,2)
    !     Q=M*N
    !
    !     if (.not. allocated(u_vec)) then
    !       allocate(u_vec(Q))
    !     endif
    !
    !     ! vectorice
    !     u_vec=reshape(u,(/Q/))
    !     !invert lattice spacing
    !     dx_inv = 1.0_rk/dx
    !     !-------------------------------------------------------------
    !     ! Compute the derivative:
    !     !           u_vec= dx_inv Dc*u_vec+0*u_vec
    !     ! for this step we use the fast blas Implementation
    !     ! Explanation of the arguments from left to right:
    !     !         1.'N'= normal Matrix (not "T" transposed)
    !     !         2. Number of rows of Matrix Dx
    !     !         3. Number of columns of matrix Dx
    !     !         4. scaling factor of the matrix vector product
    !     !         5. Matrix of dimension (LDA,Ncols)
    !     !         6. LDA is the allocated memory of the 1 dimension of the Matrix
    !     !         7. vector
    !     !         8. increment between vector elements
    !     !         9. scaling factor of the second term
    !     !         10. Vector of the second term
    !     !         11. increment of the vector
    !     call dgbmv('N', Q,   Q,   dx_inv,  Dx, M*N,  u_vec, 1,   0.0, u_vec, 1)
    !     !-------------------------------------------------------------
    !     u_x=reshape(u_vec,(/N,M/))
    ! end subroutine DGE_Dx

    !> derivative in y direction
    ! subroutine  D_y(u_y,dy, u)
    !     external dgemv
    !
    !     real(kind=rk), intent(in)       :: dy
    !     real(kind=rk), intent(in)       :: u(:,:)
    !     real(kind=rk), intent(out)      :: u_y(:,:)
    !
    !     integer :: i, j ,Q, M, N
    !     real(kind=rk) :: dy_inv
    !     real(kind=rk),allocatable,save :: u_vec(:)
    !
    !     N=size(u,1)
    !     M=size(u,2)
    !     Q=M*N
    !
    !     if (.not. allocated(u_vec)) then
    !       allocate(u_vec(Q))
    !     endif
    !     ! vectorice
    !     u_vec=reshape(u,(/N*M/))
    !     !invert lattice spacing
    !     dy_inv = 1.0_rk/dy
    !     !-------------------------------------------------------------
    !     ! Compute the derivative:
    !     !           u_vec= dy_inv Dc*u_vec+0*u_vec
    !     ! for this step we use the fast blas Implementation
    !     ! Explanation of the arguments from left to right:
    !     !         1.'N'= normal Matrix (not "T" transposed)
    !     !         2. Number of rows of Matrix Dx
    !     !         3. Number of columns of matrix Dx
    !     !         4. scaling factor of the matrix vector product
    !     !         5. Matrix of dimension (LDA,Ncols)
    !     !         6. LDA is the allocated memory of the 1 dimension of the Matrix
    !     !         7. vector
    !     !         8. increment between vector elements
    !     !         9. scaling factor of the second term
    !     !         10. Vector of the second term
    !     !         11. increment of the vector
    !     call dgemv('N', Q,   Q,   dy_inv,  Dy, M*N,  u_vec, 1,   0.0, u_vec, 1)
    !     !-------------------------------------------------------------
    !     u_y=reshape(u_vec,(/N,M/))
    ! end subroutine D_y


    subroutine  print_mat(Mat)

        real(kind=rk), intent(in)       :: Mat(:,:)
        integer                         :: i, j
        character(len=16) :: fmt
        character(len=3) :: ncols_str

        write(*,*) " "
        write(ncols_str,'(i3.3)') size(Mat,2)
        fmt = '('//ncols_str//'(f9.3,1x))'
        do j = 1,size(Mat,1)
          write(*,fmt) Mat(j,:)
        end do

    end subroutine print_mat

!> create diagonal matrix with diagonal values specified by input
    subroutine  DGE_diag( diag_vals, offset, D)
      !> \details
      !> namingconvention like it is done in blas:
      !>        + D  = Double
      !>        + GE = General
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

    end subroutine DGE_diag

    !===============================================================================
    !> create a diagonal matrix in the sparse form
    subroutine DUS_diag(D_VAL,D_INDX,D_JNDX,D_SHAPE,offset)
      !> \details
      !> namingconvention like it is done in blas:
      !>        + D  = Double
      !>        + US = Unstructured sparse
      !---------------------------------------------------------------
      real(kind=rk),     intent(in)  :: D_VAL(:)
      integer, optional, intent(in)  :: offset
      integer,           intent(out) ::D_INDX(size(D_VAL)),D_JNDX(size(D_VAL))
      integer,           intent(out) ::D_SHAPE(2)
      !---------------------------------------------------------------
      integer ::k, i,N

      k=0
      if ( PRESENT(offset) ) then
        k=offset
      end if

      N= UBOUND(D_val,DIM=1)

      D_SHAPE(1:2)=N+abs(k)

      if ( k<0 ) then
        forall(i=1:N)
          D_INDX(i)=i-k
          D_JNDX(i)=i
        end forall
      else
        forall(i=1:N)
          D_INDX(i)=i
          D_JNDX(i)=i+k
        end forall
      end if


    end subroutine  DUS_diag
    !===============================================================================




    ! subroutine  example( Bs, g, dx, u, dudx)
    !     external dgemm
    !     external sgeprt
    !
    !     integer(kind=ik), intent(in)    :: g, Bs
    !     real(kind=rk), intent(in)       :: dx
    !     real(kind=rk), intent(in)       :: u(Bs+2*g, Bs+2*g)
    !     real(kind=rk), intent(out)      :: dudx(Bs+2*g, Bs+2*g)
    !
    !     integer                         :: i, j
    !     real(kind=rk)                   :: dx_inv
    !     real(kind=rk),ALLOCATABLE         :: D(:,:)
    !     real(kind=rk),ALLOCATABLE         :: A(:,:)
    !     real(kind=rk)        :: C(5,5)
    !     real::C_(5,5),A_(5,5),D_(5,5)
    !
    !     character(len=16) :: fmt
    !     character(len=3) :: ncols_str
    !     A=diag((/2.0_rk,2.0_rk, 2.0_rk, 2.0_rk,2.0_rk/))
    !     D=diag((/2.0_rk, 1.0_rk, 3.0_rk,3.0_rk/),-1)+diag((/2.0_rk, 1.0_rk, 3.0_rk,1.0_rk/),1)
    !   !  tfm  tfm  rowA colB K  alpha a lda  b  ldb beta c  ldc */
    !     ! A_=real(A,8)
    !     ! D_=real(D,8)
    !     ! C_=real(C,8)
    !     call dgemm('N', 'N', 5,   5,   5, 1.0_rk,D  ,  5,A, 5,  0.0_rk, C, 5)
    !     ! A=A_allocatable
    !     ! D=D_
    !     ! C=C_
    !     !call sgeprt(5, 5, D, 'a= ')
    !     call print_mat(A)
    !     call print_mat(D)
    !     call print_mat(C)
    !
    !
    !
    ! end subroutine example

!
!
!     subroutine  mat_vec_mul( A, v, b)
!
!         real(kind=rk), intent(in)       :: A(:,:)
!         real(kind=rk), intent(out)      :: v(:)
!         real(kind=rk), intent(out)      :: b(size(A,1))
!         integer                         :: i, j,N,M
!         real(kind=rk)                   :: tmp
!         N =size(A,1)
!         M =size(A,2)
!
!         if ( M.ne.size(v) ) then
!           call abort(2344,"Error: Dimension of vector and matrix do not fit")
!         end if
!
!         do j = 1,N
!             tmp=0
!             do i = 1,M
!               tmp=tmp+A(j,i)*v(i)
!             end do
!             b(j)=tmp
!         end do
!
!     end subroutine mat_vec_mul
!
!
!
!
!         subroutine  mat_vec_mul( A, v, b)
!
!             real(kind=rk), intent(in)       :: A(:,:)
!             real(kind=rk), intent(out)      :: v(:)
!             real(kind=rk), intent(out)      :: b(size(A,1))
!             integer                         :: i, j,N,M
!             real(kind=rk)                   :: tmp
!             N =size(A,1)
!             M =size(A,2)
!
!             if ( M.ne.size(v) ) then
!               call abort(2344,"Error: Dimension of vector and matrix do not fit")
!             end if
!
!             do j = 1,N
!                 tmp=0
!                 do i = 1,M
!                   tmp=tmp+A(j,i)*v(i)
!                 end do
!                 b(j)=tmp
!             end do
!
!         end subroutine mat_vec_mul
!
!
subroutine kron(A,B,AB)
!AB = Kronecker product of A and B, both two-dimensional arrays.
!Considers the arrays to be addressed as A(row,column), despite any storage order arrangements.        .
!Creating array AB to fit here, adjusting the caller's array AB, may not work on some compilers.
INTEGER A(:,:),B(:,:)!Two-dimensional arrays, lower bound one.
INTEGER, ALLOCATABLE:: AB(:,:)!To be created to fit.
INTEGER R,RA,RB,C,CA,CB,I,J!Assistants.
 RA = UBOUND(A,DIM = 1)!Ascertain the upper bounds of the incoming arrays.
 CA = UBOUND(A,DIM = 2)!Their lower bounds will be deemed one,
 RB = UBOUND(B,DIM = 1)!And the upper bound as reported will correspond.
 CB = UBOUND(B,DIM = 2)!UBOUND(A) would give an array of two values, RA and CA, more for higher dimensionality.
 IF (ALLOCATED(AB)) DEALLOCATE(AB)  !Discard any lingering storage.
 ALLOCATE (AB(RA*RB,CA*CB))!Obtain the exact desired size.
 R = 0!Syncopation: start the row offset.
 do I = 1,RA!Step down the rows of A.
   C = 0!For each row, start the column offset.
   do J = 1,CA!Step along the columns of A.
     AB(R + 1:R + RB,C + 1:C + CB) = A(I,J)*B!Place a block of B values.
     C = C + CB!Advance a block of columns.
   end do!On to the next column of A.
   R = R + RB!Advance a block of rows.
 end do!On to the next row of A.
end subroutine  kron!No tests for bad parameters, or lack of storage.



!===============================================================================
!> AB = Kronecker product of A and B, both two-dimensional arrays in the sparse form
subroutine DUS_kron(A_VAL,A_INDX,A_JNDX,A_SHAPE,B_VAL,B_INDX,B_JNDX,B_SHAPE,AB_VAL,AB_INDX,AB_JNDX,AB_SHAPE)
!>  \details sparse impelementation
!>  Naming convention is as for the blas library
  !---------------------------------------------------------------
  real(kind=rk),intent(in):: A_VAL(:), B_VAL(:)
  integer, intent(in) ::A_INDX(:),A_JNDX(:)
  integer, intent(in) ::B_INDX(:),B_JNDX(:)
  integer, intent(in) ::A_SHAPE(2),B_SHAPE(2)
  real(kind=rk),intent(inout), ALLOCATABLE:: AB_VAL(:)!To be created to fit.
  integer,intent(inout), ALLOCATABLE:: AB_INDX(:)!To be created to fit.
  integer,intent(inout), ALLOCATABLE:: AB_JNDX(:)!To be created to fit.
  integer, intent(out) ::AB_SHAPE(2)
  !---------------------------------------------------------------
  integer :: n_A,n_B,ia,ib,iab

  n_A= UBOUND(A_val,DIM = 1)!Ascertain the upper bounds of the incoming arrays.
  n_B= UBOUND(B_val,DIM = 1)!And the upper bound as reported will correspond.

  IF (allocated(AB_VAL)) deallocate(AB_VAL)
  IF (allocated(AB_INDX)) deallocate(AB_INDX)
  IF (allocated(AB_JNDX)) deallocate(AB_JNDX)

  allocate (AB_VAL(n_A*n_B))
  allocate (AB_INDX(n_A*n_B))
  allocate (AB_JNDX(n_A*n_B))

  AB_SHAPE=A_SHAPE*B_SHAPE

  iab=1
  do ia = 1, n_A
    do ib = 1, n_B
      AB_VAL(iab)=A_VAL(ia)*B_VAL(ib)
      AB_INDX(iab)=((A_INDX(ia)-1)*B_SHAPE(1))+B_INDX(ib)
      AB_JNDX(iab)=((A_JNDX(ia)-1)*B_SHAPE(2))+B_JNDX(ib)
      iab=iab+1
    end do
  end do

end subroutine  DUS_kron
!===============================================================================


function  DUS_full(a_val, a_indx, a_jndx, a_shape) result(a)
    !-----------------------------------------------------
    integer,intent(in):: a_indx(:)
    integer,intent(in):: a_jndx(:)
    integer,intent(in):: a_shape(:)
    real(kind=rk), intent(in) :: a_val(:)
    real(kind=rk) :: a(a_shape(1),a_shape(2))
    !-----------------------------------------------------
    integer(kind=ik) :: ia

    a=0.0_rk
    do ia = 1, UBOUND(a_val, DIM = 1)
        a(a_indx(ia),a_jndx(ia))=a_val(ia)
    end do

end function DUS_full


subroutine  test_sparse_kron()
  external dgbmv

  real(kind=rk):: A_VAL(2), B_VAL(4)
  integer ::A_INDX(2),A_JNDX(2)
  integer ::B_INDX(4),B_JNDX(4)
  integer ::A_SHAPE(2),B_SHAPE(2)
  real(kind=rk),allocatable:: AB_VAL(:)!To be created to fit.
  integer,allocatable:: AB_INDX(:)!To be created to fit.
  integer,allocatable:: AB_JNDX(:)!To be created to fit.
  integer::AB_SHAPE(2)

  a_val=(/2 ,5 /)
  a_INDX=(/1,2/)
  a_JNDX=(/2,1/)
  a_shape=(/2,2/)
  call print_mat(DUS_full(a_val,a_indx,a_jndx,a_shape))
  ! A_val=
  ! 0.000     2.000
  ! 5.000     0.000


  b_val=(/1,2,3,4/)
  b_INDX=(/1,1,2,2/)
  b_JNDX=(/1,2,1,2/)
  b_shape=(/2,2/)
  call print_mat(DUS_full(b_val,b_indx,b_jndx,b_shape))
  !B_val=    1.000     2.000
  !          3.000     4.000


  call DUS_kron(A_VAL,A_INDX,A_JNDX,A_SHAPE,B_VAL,B_INDX,B_JNDX,B_SHAPE,AB_VAL,AB_INDX,AB_JNDX,AB_SHAPE)

  call print_mat(DUS_full(ab_val,ab_indx,ab_jndx,ab_shape))

  !AB_val=
  !  0.000     0.000     2.000     4.000
  !  0.000     0.000     6.000     8.000
  !  5.000    10.000     0.000     0.000
  ! 15.000    20.000     0.000     0.000
  !


end subroutine test_sparse_kron

#endif


end module module_sparse_operators
