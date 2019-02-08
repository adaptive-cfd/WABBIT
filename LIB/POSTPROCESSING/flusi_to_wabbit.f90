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
!> \date 20/7/2018 - add read_field_flusi_MPI for parallel reading, commit 03b933e706b988828a0b0321baedb2dc9f76773d
!-----------------------------------------------------------------------------------------------------
!
subroutine flusi_to_wabbit(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use mpi
    use module_MPI

    implicit none

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
    integer(kind=ik)                  :: hvy_n, lgt_n, k
    integer(kind=ik)                  :: refinement_status
    integer(kind=ik) , dimension(3)   :: Bs
    real(kind=rk), dimension(3)       :: x0, dx
    real(kind=rk), dimension(3)       :: domain
    character(len=3)                  :: dummy
    integer(kind=ik), dimension(3)    :: nxyz
    real(kind=rk), dimension(3)       :: level_tmp
    integer(kind=ik), dimension(3)    :: level
    integer(kind=ik)                  :: status, start_x, start_y, start_z

!-----------------------------------------------------------------------------------------------------
    params%n_eqn  = 1
!----------------------------------------------
    call get_command_argument(2, file_in)
    ! does the user need help?
    if (file_in=='--help' .or. file_in=='--h') then
        if (params%rank==0) then
            write(*,*) "postprocessing subroutine to read a file from flusi and convert it to wabbit format. command line:"
            write(*,*) "mpi_command -n number_of_processes ./wabbit-post --flusi-to-wabbit source.h5 target.h5 target_blocksize "
        end if
    else
        if (params%rank==0) write(*,*) "ATTENTION: this routine works in parallel only with ghost nodes synchronization --generic-sequence (default)"
        ! get values from command line (filename and desired blocksize)
        call check_file_exists(trim(file_in))
        call get_command_argument(3, file_out)

        ! /todo to be checked (l66-72)

        call get_command_argument(4, dummy)
        read(dummy,*) Bs(1)
        call get_command_argument(5, dummy)
        read(dummy,*) Bs(2)
        call get_command_argument(6, dummy)
        read(dummy,*) Bs(3)
        if (mod(Bs(1),2)==0 .or. mod(Bs(2),2)==0 .or. mod(Bs(3),2)==0) call abort(7844, "ERROR: For WABBIT we need an odd blocksize!")

        ! read attributes such as number of discretisation points, time, domain size
        call get_attributes_flusi(file_in, nxyz, time, domain)
        if (nxyz(1)/=1) then
            params%dim = 3
        else
            params%dim = 2
        end if

        do k=1,3
          level_tmp(k) = log(dble(nxyz(k))/dble(Bs(k)-1)) / log(2.0_rk)
          level(k) = int(level_tmp(k))
        enddo

        ! check the input
        if (nxyz(2)/=nxyz(3)) call abort(8724, "ERROR: nx and ny differ. This is not possible for WABBIT")
        if (mod(nxyz(2),2)/=0) call abort(8324, "ERROR: nx and ny need to be even!")
        if (mod(nxyz(1),(Bs(1)-1))/=0 .or. mod(nxyz(2),(Bs(2)-1))/=0 .or. mod(nxyz(3),(Bs(3)-1))/=0 &
            .or. abs(level(1)-level_tmp(1))>1.e-14 .or. abs(level(2)-level_tmp(2))>1.e-14 .or. abs(level(3)-level_tmp(3))>1.e-14)&
            call abort(2948, "ERROR: I'm afraid your saved blocksize does not match for WABBIT")
        if ((Bs(1)/nxyz(1))/=(Bs(2)/nxyz(2)) .or. (Bs(1)/nxyz(1))/=(Bs(3)/nxyz(3)) .or. (Bs(2)/nxyz(2))/=(Bs(3)/nxyz(3)))&
            call abort(2411921, "ERROR: Blocksizes of flusi and WABBIT are not compatible")

        ! set important parameters
        params%max_treelevel=max(max(level(1),level(2)),max(level(2),level(3)))
        params%domain_size(1) = domain(2)
        params%domain_size(2) = domain(3)
        params%Bs = Bs
        params%n_ghosts = 1_ik
        params%order_predictor = 'multiresolution_4th'
        if (params%dim == 3) then
            lgt_n = 8_ik**params%max_treelevel
            params%domain_size(3) = domain(1)
        else
            lgt_n = 4_ik**params%max_treelevel
        end if
        hvy_n = lgt_n
        params%number_blocks = lgt_n

        call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
            hvy_active, lgt_sortednumlist)

        ! create an equidistant grid (order of light id is important!)
        call create_equidistant_grid( params, lgt_block, hvy_block, hvy_neighbor,&
            lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, &
            params%max_treelevel, .true.)

        ! read the field from flusi file and organize it in WABBITs blocks
        call read_field_flusi_MPI(file_in, hvy_block, lgt_block, hvy_n,&
            hvy_active, params, nxyz)

        ! set refinement status of blocks not lying at the outer edge to 11 (historic fine)
        ! they will therefore send their redundant points to the last blocks in x,y and z-direction
        do k=1,lgt_n
            call get_block_spacing_origin( params, lgt_active(k), lgt_block, x0, dx )
            start_x = nint(x0(1)/dx(1))
            start_y = nint(x0(2)/dx(2))
            if (params%dim == 3) then
                start_z = nint(x0(3)/dx(3))
            else
                start_z = 0_ik
            end if
            lgt_block(lgt_active(k), params%max_treelevel + idx_refine_sts) = refinement_status(start_x, start_y, start_z, nxyz, Bs)
        end do

        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        iteration = 0
        call write_field(file_out, time, iteration, 1, params, lgt_block,&
            hvy_block, lgt_active, lgt_n, hvy_n, hvy_active)
    end if

end subroutine flusi_to_wabbit

function refinement_status(start_x, start_y, start_z, Bs_f, Bs)
  use module_precision
  implicit none
  integer(kind=ik), intent(in) :: start_x, start_y, start_z
  integer(kind=ik), dimension(3), intent(in) :: Bs, Bs_f
  integer(kind=ik) :: refinement_status

  ! if I'm the last block in x,y and/or z-direction,
  ! set my refinement refinement_status to 0, otherwise to 11 (historic fine)
  if (((start_x==Bs_f(1)-Bs(1)+1) .or. (start_y==Bs_f(2)-Bs(2)+1)) .or. (start_z==Bs_f(3)-Bs(3)+1)) then
      refinement_status = 0_ik
  else
      refinement_status = 11_ik
  end if

end function refinement_status
