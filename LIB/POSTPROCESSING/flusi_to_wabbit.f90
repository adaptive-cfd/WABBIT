!> \file
! WABBIT
!> \name flusi_to_wabbit.f90
!> \version 0.5
!> \author sm
!
!> \brief postprocessing routine that generates a WABBIT-readable .h5 file (a field composed in blocks) 
!! from a .h5 file where all data is stored in one block
!
! = log ======================================================================================
!> \date  07/03/18 - create hashcode: commit 8f4858f429c6c3f537190f48a8e8a931154a01d5
!-----------------------------------------------------------------------------------------------------
!
subroutine flusi_to_wabbit(help, params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use mpi

    implicit none

    !> help flag
    logical, intent(in)                :: help
    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: file_in
    character(len=80)      :: file_out
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration

    integer(kind=ik), allocatable     :: lgt_block(:, :)
    real(kind=rk), allocatable        :: hvy_block(:, :, :, :, :)
    integer(kind=ik), allocatable     :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable     :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable  :: lgt_sortednumlist(:,:)
    integer(kind=ik)                  :: hvy_n, lgt_n,level, Bs
    real(kind=rk), dimension(3)       :: domain
    character(len=3)                  :: Bs_read
    integer(kind=ik), dimension(3)    :: nxyz
    real(kind=rk)                     :: level_tmp

!-----------------------------------------------------------------------------------------------------
    params%number_data_fields  = 1
!----------------------------------------------
    if (help .eqv. .true.) then
        if (params%rank==0) then
            write(*,*) "postprocessing subroutine to read a file from flusi and convert it to wabbit format. command line:"
            write(*,*) "mpi_command -n 1 ./wabbit-post 2D --flusi-to-wabbit source.h5 target.h5 target_blocksize"
        end if
    else
        ! get values from command line (filename and desired blocksize)
        call get_command_argument(3, file_in)
        call check_file_exists(trim(file_in))
        call get_command_argument(4, file_out)
        call get_command_argument(5, Bs_read)
        read(Bs_read,*) Bs
        if (mod(Bs,2)==0) call abort(7844, "ERROR: For WABBIT we need an odd blocksize!")

        ! read attributes such as number of discretisation points, time, domain size
        call get_attributes_flusi(file_in, nxyz, time, domain)
        if (.not. params%threeD_case .and. nxyz(1)/=1) &
            call abort(8714, "ERROR: saved datafield is 3D, WABBIT expects 2D")

        level_tmp = log(dble(nxyz(2))/dble(Bs-1)) / log(2.0_rk)
        level = int(level_tmp)

        ! check the input
        if (nxyz(2)/=nxyz(3)) call abort(8724, "ERROR: nx and ny differ. This is not possible for WABBIT")
        if (mod(nxyz(2),2)/=0) call abort(8324, "ERROR: nx and ny need to be even!")
        if (mod(nxyz(2),(Bs-1))/=0 .or. abs(level-level_tmp)>1.e-14)&
            call abort(2948, "ERROR: I'm afraid your saved blocksize does not match for WABBIT")

        params%max_treelevel=level
        params%Lx = domain(2)
        params%Ly = domain(3)
        params%number_block_nodes = Bs
        params%number_ghost_nodes = 0_ik
        if (params%threeD_case) then
            lgt_n = 8_ik**params%max_treelevel
            params%Lz = domain(1)
        else
            lgt_n = 4_ik**params%max_treelevel
        end if
        hvy_n = lgt_n
        params%number_blocks = lgt_n

        call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
            hvy_active, lgt_sortednumlist, .false.)

        call create_equidistant_grid( params, lgt_block, hvy_block, hvy_neighbor,&
            lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, &
            params%max_treelevel, .true.)

        call read_field_flusi(file_in, hvy_block, lgt_block, hvy_n,&
            hvy_active, params, nxyz(2))

        iteration = 0
        call write_field(file_out, time, iteration, 1, params, lgt_block,&
            hvy_block, lgt_active, lgt_n, hvy_n)
    end if

end subroutine flusi_to_wabbit
