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
    character(len=5)                  :: dummy
    integer(kind=ik), dimension(3)    :: nxyz
    real(kind=rk), dimension(3)       :: level_tmp
    integer(kind=ik), dimension(3)    :: level
    integer(kind=ik)                  :: status, start_x, start_y, start_z

    !-----------------------------------------------------------------------------------------------------
    params%n_eqn  = 1
    !----------------------------------------------
    call get_command_argument(2, file_in)
    ! does the user need help?
    if (params%rank==0) then
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "Convert a FLUSI hdf5 file to a WABBIT hdf5 file"
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "INPUT: a flusi (or equivalent) generated HDF5 file. It contains an EQUIDISTANT grid with one QTY "
        write(*,*) "       stored (e.g. a velocity component). The dataset in the file must have the name of the file"
        write(*,*) "       (e.g. mask_00000.h5 is expected to hold a dataset named mask)"
        write(*,*) "       The resolution must be ISOTROPIC and EVEN numbers (e.g. 1024x1024x1024)"
        write(*,*) "       You can check using HDF's h5dump -A tool."
        write(*,*) ""
        write(*,*) "OUTPUT: A wabbit compatible file with the same data organized in blocks. No adaptation is performed,"
        write(*,*) "        if you pass a field of zeros thats what you get out. note wabbit uses double precision by "
        write(*,*) "        default, so your file may be exactly twice as large."
        write(*,*) ""
        write(*,*) "USAGE: ./wabbit-post --flusi-to-wabbit [source.h5 [target.h5] [Bsx] [Bsy] [Bsz]"
        write(*,*) "       for 3D data:"
        write(*,*) "       ./wabbit-post --flusi-to-wabbit mask_00000.h5 wabbit_000000000.h5 33 33 33"
        write(*,*) "       for 2D data:"
        write(*,*) "       ./wabbit-post --flusi-to-wabbit mask_00000.h5 wabbit_000000000.h5 33 33 1"
        write(*,*) ""
        write(*,*) "NOTES: * Even though anisotropic resolution would be possible, it is not implemented yet"
        write(*,*) "       * The blocksize BS is isotropic"
        write(*,*) "       * Routine is parallel"
        write(*,*) "       * the resolution of the input data must be an integer multiple of (Bs-1)"
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end if

    ! it should be a much simpler function now, but it has to be updated
    call abort(7666, "This function is to be updated for the new block definition.")

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
    if (mod(Bs(1),2)==0 .or. mod(Bs(2),2)==0 .or. mod(Bs(3),2)==0) then
        write(*,*) Bs
        call abort(7844, "ERROR: For WABBIT we need an odd blocksize!")
    endif

    ! read attributes such as number of discretisation points, time, domain size
    call get_attributes_flusi(file_in, nxyz, time, domain)

    if (params%rank==0) then
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,'(" Input file: ", A)') trim(adjustl(file_in))
        write(*,'(" Output file: ", A)') trim(adjustl(file_out))
        write(*,'(" Resolution in FLUSI input is: ", 3(i5,1x))') nxyz
        write(*,'(" Domain information in FLUSI input is: ", 3(f5.3,1x))') domain
    endif



    if (nxyz(1) /= 1) then
        !-----------------------------------------------------------------------
        ! 3D
        !-----------------------------------------------------------------------
        params%dim = 3

        if (mod(nxyz(1),2)/=0) call abort(8324, "ERROR: nx, ny, nz need to be even!")
        if (mod(nxyz(2),2)/=0) call abort(8324, "ERROR: nx, ny, nz need to be even!")
        if (mod(nxyz(3),2)/=0) call abort(8324, "ERROR: nx, ny, nz need to be even!")

    else
        !-----------------------------------------------------------------------
        ! 2D
        !-----------------------------------------------------------------------
        params%dim = 2

        ! there is a goofy shift. in flusi, 2D runs ALWAYS need to set the x-direction to
        ! 1 point (that is a memory consistency requirement)
        ! Hence, a 2D flusi array has the dimensions 1 x NX x NY (from the wabbit point of view)
        domain(1) = domain(2)
        domain(2) = domain(3)

        nxyz(1) = nxyz(2)
        nxyz(2) = nxyz(3)

        if (mod(nxyz(1),2)/=0) call abort(8324, "ERROR: nx and ny need to be even!")
        if (mod(nxyz(2),2)/=0) call abort(8324, "ERROR: nx and ny need to be even!")
    end if


    level_tmp = 0.0
    level = 0
    do k = 1, params%dim
        if (Bs(k)==1) call abort(11021903, "ERROR: Blocksize cannot be one! (Even though that is an odd number)")
        level_tmp(k) = log(dble(nxyz(k))/dble(Bs(k)-1)) / log(2.0_rk)
        level(k) = int(level_tmp(k))
    enddo


    if (params%rank==0) then
        write(*,'(" Data dimensionality is: ", i5)') params%dim
        write(*,'(" Wabbit level is: (should be integers)", 3(f4.1,1x))') level_tmp
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    endif

    ! this would actually now be possible but needs to be implemented:
    !if (nxyz(1)/=nxyz(2)) call abort(8724, "ERROR: nx and ny differ. This is not possible for WABBIT")

    do k = 1, params%dim
        if (mod(nxyz(k),(Bs(k)-1)) /=0 ) then
            call abort(11021901, "The input data size and your choice of BS are not compatible.")
        endif
        ! non-integer levels indicate the resolution is not a multiple of Bs
        if ( abs(dble(level(k))-level_tmp(k)) > 1.e-10 ) then
            call abort(11021902, "ERROR: I'm afraid your selected blocksize does not match for WABBIT")
        endif
    enddo


    ! set important parameters
    params%max_treelevel=max(max(level(1),level(2)),max(level(2),level(3)))
    params%domain_size = domain
    params%Bs = Bs
    params%n_ghosts = 1_ik
    params%order_predictor = 'multiresolution_4th'
    lgt_n = (2_ik**params%max_treelevel)**params%dim
    hvy_n = lgt_n
    params%number_blocks = lgt_n

    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist)

    ! create an equidistant grid (order of light id is important!)
    call create_equidistant_grid( params, lgt_block, hvy_neighbor,&
    lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, &
    params%max_treelevel, .true., tree_ID=1)

    ! read the field from flusi file and organize it in WABBITs blocks
    call read_field_flusi_MPI(file_in, hvy_block, lgt_block, hvy_n,&
    hvy_active, params, nxyz)
    ! set refinement status of blocks not lying at the outer edge to 11 (historic fine)
    ! they will therefore send their redundant points to the last blocks in x,y and z-direction
    do k = 1, lgt_n
        call get_block_spacing_origin( params, lgt_active(k), lgt_block, x0, dx )
        start_x = nint(x0(1)/dx(1))
        start_y = nint(x0(2)/dx(2))
        if (params%dim == 3) then
            start_z = nint(x0(3)/dx(3))
        else
            start_z = 0_ik
        end if
        lgt_block(lgt_active(k), params%max_treelevel + IDX_REFINE_STS) = refinement_status(start_x, start_y, start_z, nxyz, Bs)
    end do

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    iteration = 0
    call write_field(file_out, time, iteration, 1, params, lgt_block, &
    hvy_block, lgt_active, lgt_n, hvy_n, hvy_active)

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
