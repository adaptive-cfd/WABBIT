!> \brief postprocessing routine that generates a WABBIT-readable .h5 file (a field composed in blocks)
!! from a .h5 file where all data is stored in one block
!-----------------------------------------------------------------------------------------------------
!
subroutine flusi_to_wabbit(params)
    use module_globals
    use module_mesh
    use module_params
    use mpi
    use module_MPI
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout) :: params
    character(len=cshort)                 :: file_in
    character(len=cshort)                 :: file_out
    real(kind=rk)                     :: time
    integer(kind=ik)                  :: iteration

    real(kind=rk), allocatable        :: hvy_block(:, :, :, :, :)!, data_flusi(:,:,:)
    integer(kind=ik)                  :: tree_ID=1, hvy_id, k
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
    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )
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
        write(*,*) "USAGE: ./wabbit-post --flusi-to-wabbit --input=mask_00000.h5 --output=wabbit_000000000.h5 --Bs=33"
        write(*,*) ""
        write(*,*) "NOTES: * Even though anisotropic resolution would be possible, it is not implemented yet"
        write(*,*) "       * The blocksize BS is isotropic"
        write(*,*) "       * Routine is parallel"
        write(*,*) "       * the resolution of the input data must be an integer multiple of (Bs-1)"
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end if

    call get_cmd_arg( "--input", file_in, default="none" )
    call get_cmd_arg( "--output", file_out, default="none" )
    call get_cmd_arg( "--bs", Bs(1), default=32 )

    call check_file_exists(trim(file_in))


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

        bs(1:params%dim) = bs(1)

        if (mod(nxyz(1),2)/=0) call abort(8324, "ERROR: nx, ny, nz need to be even!")
        if (mod(nxyz(2),2)/=0) call abort(8324, "ERROR: nx, ny, nz need to be even!")
        if (mod(nxyz(3),2)/=0) call abort(8324, "ERROR: nx, ny, nz need to be even!")

    else
        !-----------------------------------------------------------------------
        ! 2D
        !-----------------------------------------------------------------------
        params%dim = 2
        bs(1:params%dim) = bs(1)
        bs(3) = 1

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
        level_tmp(k) = log(dble(nxyz(k))/dble(Bs(k))) / log(2.0_rk)
        level(k) = int(level_tmp(k))
    enddo


    if (params%rank==0) then
        write(*,'(" Data dimensionality is: ", i5)') params%dim
        write(*,'(" Wabbit level is: (should be integers)", 3(f4.1,1x))') level_tmp
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    endif

    ! this would actually now be possible but needs to be implemented:
    if (nxyz(1)/=nxyz(2)) call abort(8724, "ERROR: nx and ny differ. This is not possible for WABBIT")

    do k = 1, params%dim
        if ( mod(nxyz(k),Bs(k)) /=0 ) then
            call abort(11021901, "The input data size and your choice of BS are not compatible.")
        endif
        ! non-integer levels indicate the resolution is not a multiple of Bs
        if ( abs(dble(level(k))-level_tmp(k)) > 1.e-10 ) then
            call abort(11021902, "ERROR: I'm afraid your selected blocksize does not match for WABBIT")
        endif
    enddo


    ! set important parameters
    params%Jmax=max(max(level(1),level(2)),max(level(2),level(3)))
    params%domain_size = domain
    params%Bs = Bs
    params%g = 1_ik
    params%order_predictor = 'multiresolution_4th'
    params%wavelet = "CDF40"
    lgt_n(tree_ID) = (2_ik**params%Jmax)**params%dim
    hvy_n(tree_ID) = lgt_n(tree_ID)
    params%number_blocks = lgt_n(tree_ID)

    call allocate_forest(params, hvy_block)

    ! create an equidistant grid (order of light id is important!)
    call createEquidistantGrid_tree( params, hvy_block, params%Jmax, .true., tree_ID)



    ! read the field from flusi file and organize it in WABBITs blocks
    call read_field_flusi_MPI(file_in, hvy_block, params, tree_ID)

    iteration = 0

    call saveHDF5_tree(file_out, time, iteration, 1, params, hvy_block, tree_ID)

end subroutine flusi_to_wabbit
