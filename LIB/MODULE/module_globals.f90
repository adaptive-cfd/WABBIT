!! Module of global parameters and functions
!!    * index positions of lgt_data (IDX_TREE_ID, IDX_MESH_LVL, IDX_REFINE_STS etc)
!!    * functions abort
!!    * global prints (todo)
!!    * global MPI communicator
!-----------------------------------------------------------------
module module_globals

    use mpi

    implicit none
    ! index of light data quantities lgt_data(block_index,max_treeleve+idx_<quantity>)
    ! we define the global indices idx_<quantity> here, because it enables us to choose the order of indexing
    ! if necessary
    integer, parameter, public  :: IDX_MESH_LVL        = 1 ! current mesh level of the block
    integer, parameter, public  :: IDX_REFINE_STS      = 2 ! refinement status
    integer, parameter, public  :: IDX_TREE_ID         = 3 ! number of the tree in the forest
    integer, parameter, public  :: EXTRA_LGT_FIELDS    = 3 ! number of data fields additionaly to treecode
    integer, parameter, public  :: tree_ID_flow = 1, tree_ID_mask = 2, tree_ID_mask_coarser = 3
    ! this parameter is a hack. in most parts of the code, a block has n_eqn component entries.
    ! universality dictates that we can also use a different number of components, for example
    ! when syn'ing the mask function (which in many cases has six entries.)
    ! New in 06/2021: the hack continues. We now set this parameter at different places
    ! to save on memory. That can be params%n_eqn (default in simulations), 6 (if mask is synced). The new default is 3,
    ! for postprocessing.
    integer, public  :: N_MAX_COMPONENTS = 3
    ! define data precision parameters
    integer, parameter, public :: sngl_prec=selected_real_kind(4)
    integer, parameter, public :: dble_prec=selected_real_kind(8)
    ! default length of strings (short and long character)
    integer, parameter, public :: cshort=80
    integer, parameter, public :: clong =120
    integer, parameter, public :: int_prec=selected_int_kind(8)
    integer, parameter, public :: maxdigits = 16
    integer, parameter, public :: tsize = selected_int_kind(maxdigits)
    integer, parameter, public :: rk=dble_prec
    integer, parameter, public :: ik=int_prec
    real(kind=rk),parameter, public :: pi  = 4.0_rk * atan(1.0_rk)
    !> global communicator for WABBIT! Dont use MPI_COMM_WORLD!!!!!
    integer(kind=ik) ::  WABBIT_COMM

    !subroutines of this module
    interface abort
        module procedure abort1, abort2, abort3
    end interface

    interface areArraysSameSize
        module procedure areArraysSameSize_1d, areArraysSameSize_2d, areArraysSameSize_3d, areArraysSameSize_4d, areArraysSameSize_5d
    end interface


    contains

        function areArraysSameSize_1d(a,b)
            implicit none
            real(kind=rk),dimension(1:),intent(in) :: a,b
            logical :: areArraysSameSize_1d

            areArraysSameSize_1d = (size(a,1)==size(b,1))
        end function

        function areArraysSameSize_2d(a,b)
            implicit none
            real(kind=rk),dimension(1:,1:),intent(in) :: a,b
            logical :: areArraysSameSize_2d

            areArraysSameSize_2d = (size(a,1)==size(b,1)).and.(size(a,2)==size(b,2))
        end function

        function areArraysSameSize_3d(a,b)
            implicit none
            real(kind=rk),dimension(1:,1:,1:),intent(in) :: a,b
            logical :: areArraysSameSize_3d

            areArraysSameSize_3d = (size(a,1)==size(b,1)).and.(size(a,2)==size(b,2)).and.(size(a,3)==size(b,3))
        end function

        function areArraysSameSize_4d(a,b)
            implicit none
            real(kind=rk),dimension(1:,1:,1:,1:),intent(in) :: a,b
            logical :: areArraysSameSize_4d

            areArraysSameSize_4d = (size(a,1)==size(b,1)).and.(size(a,2)==size(b,2)).and.(size(a,3)==size(b,3))&
            .and.(size(a,4)==size(b,4))
        end function

        function areArraysSameSize_5d(a,b)
            implicit none
            real(kind=rk),dimension(1:,1:,1:,1:,1:),intent(in) :: a,b
            logical :: areArraysSameSize_5d

            areArraysSameSize_5d = (size(a,1)==size(b,1)).and.(size(a,2)==size(b,2)).and.(size(a,3)==size(b,3))&
            .and.(size(a,4)==size(b,4)).and.(size(a,5)==size(b,5))
        end function


        ! helper routine to print when entering a routine. DEBUG only!
        ! all MPIRANKS do print
        subroutine debug_header_barrier(msg)
            use mpi
            implicit none
            character(len=*), intent(in) :: msg
            integer :: mpirank, mpicode

            ! fetch my process id
            call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

            call MPI_BARRIER(WABBIT_COMM, mpicode)
            write(*,'(A,i5,1x,A)') "DEBUG_BARRIER Entry, rank=", mpirank, msg
            call MPI_BARRIER(WABBIT_COMM, mpicode)

        end subroutine

        !> use the abort function instead of stop!
        !> this is necessary since stop, does not kill
        !> the other processes
        subroutine abort1(code,msg)
            implicit none
            character(len=*), intent(in) :: msg
            integer(kind=ik), intent(in) :: code
            integer(kind=ik) :: mpierr

            ! you may be tempted to use (if mprank=0) to have nicer output: do not do that.
            ! if root does not come to the error, then you will not see it. better all ranks
            ! display it.
            write(*,*) msg
            call MPI_ABORT( WABBIT_COMM, code, mpierr)
        end subroutine

        subroutine abort2(code)
            implicit none
            integer(kind=ik), intent(in) :: code
            integer(kind=ik) :: mpierr

            call MPI_ABORT( WABBIT_COMM, code, mpierr)
        end subroutine

        subroutine abort3(msg)
            implicit none
            character(len=*), intent(in) :: msg
            integer(kind=ik) :: mpierr

            ! you may be tempted to use (if mprank=0) to have nicer output: do not do that.
            ! if root does not come to the error, then you will not see it. better all ranks
            ! display it.
            write(*,*) msg
            call MPI_ABORT( WABBIT_COMM, 666, mpierr)
        end subroutine


        !-----------------------------------------------------------------------------
        ! convert degree to radiant
        !-----------------------------------------------------------------------------
        real(kind=rk) function deg2rad(deg)
            implicit none
            real(kind=rk), intent(in) :: deg
            real(kind=rk),parameter :: pi2 = 3.1415926535897932384626433832795028841971693993751058209749445923078164d0
            deg2rad = deg*pi2/180.d0
            return
        end function

        !-----------------------------------------------------------------------------
        ! radiant to degree
        !-----------------------------------------------------------------------------
        real(kind=rk) function rad2deg(deg)
            implicit none
            real(kind=rk), intent(in) :: deg
            real(kind=rk),parameter :: pi2 = 3.1415926535897932384626433832795028841971693993751058209749445923078164d0
            rad2deg = deg*180.d0/pi2
            return
        end function

        !-----------------------------------------------------------------------------
        ! cross product of two vectors
        !-----------------------------------------------------------------------------
        function cross(a,b)
            implicit none
            real(kind=rk),dimension(1:3),intent(in) :: a,b
            real(kind=rk),dimension(1:3) :: cross
            cross(1) = a(2)*b(3)-a(3)*b(2)
            cross(2) = a(3)*b(1)-a(1)*b(3)
            cross(3) = a(1)*b(2)-a(2)*b(1)
        end function

        !-----------------------------------------------------------------------------
        ! angle between two vectors
        !-----------------------------------------------------------------------------
        function angle_between_vectors(a,b)
            implicit none
            real(kind=rk),dimension(1:3),intent(in) :: a,b
            real(kind=rk) :: angle_between_vectors

            angle_between_vectors = acos(sum(a*b) / (norm2(a)*norm2(b)))

        end function

        !-----------------------------------------------------------------------------
        ! turn vector into unit vector
        !-----------------------------------------------------------------------------
        function unit_vector(a)
            implicit none
            real(kind=rk),dimension(1:3),intent(in) :: a
            real(kind=rk),dimension(1:3) :: unit_vector
            unit_vector = a / norm2(a)
        end function

        !-----------------------------------------------------------------------------
        ! check if two vectors are identical
        !-----------------------------------------------------------------------------
        function is_same_vector(a,b)
            implicit none
            real(kind=rk),dimension(1:),intent(in) :: a,b
            logical :: is_same_vector
            integer :: i

            is_same_vector = .true.

            do i = 1, size(a)
                if ( abs(a(i)-b(i)) >= 1.0e-11 ) then
                    is_same_vector = .false.
                    return
                endif
            enddo
        end function

        !-----------------------------------------------------------------------------
        ! 2-norm length of vectors
        !-----------------------------------------------------------------------------
        function norm2(a)
            implicit none
            real(kind=rk),dimension(1:3),intent(in) :: a
            real(kind=rk) :: norm2
            norm2 = sqrt( a(1)*a(1) + a(2)*a(2) + a(3)*a(3) )
        end function

        !-----------------------------------------------------------------------------
        ! given a point x, check if it lies in the computational domain centered at zero
        ! (note: we assume [-xl/2...+xl/2] size this is useful for insects )
        !-----------------------------------------------------------------------------
        function periodize_coordinate(x_glob, box)
            real(kind=rk),intent(in) :: x_glob(1:3), box(1:3)
            real(kind=rk),dimension(1:3) :: periodize_coordinate

            periodize_coordinate = x_glob

            if (x_glob(1)<-box(1)*0.5_rk) periodize_coordinate(1)=x_glob(1)+box(1)
            if (x_glob(2)<-box(2)*0.5_rk) periodize_coordinate(2)=x_glob(2)+box(2)
            if (x_glob(3)<-box(3)*0.5_rk) periodize_coordinate(3)=x_glob(3)+box(3)

            if (x_glob(1)>box(1)*0.5_rk) periodize_coordinate(1)=x_glob(1)-box(1)
            if (x_glob(2)>box(2)*0.5_rk) periodize_coordinate(2)=x_glob(2)-box(2)
            if (x_glob(3)>box(3)*0.5_rk) periodize_coordinate(3)=x_glob(3)-box(3)

        end function

        !---------------------------------------------------------------------------
        ! wrapper for random number generator (this may be compiler dependent)
        !---------------------------------------------------------------------------
        real(kind=rk) function rand_nbr()
            implicit none
            call random_number( rand_nbr )
        end function

    end module module_globals
