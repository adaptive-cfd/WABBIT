!> \brief module for all mesh subroutines
! ********************************************************************************************

module module_mesh

    use mpi
    use module_hdf5_wrapper
    use module_forestMetaData
    use module_params               ! global parameters
    use module_timing               ! debug module
    use module_wavelets        ! interpolation routines
    ! use MPI module, since threshold_block needs to synch ghosts
    use module_MPI
    use module_treelib              ! module with evrything related to treecodes (encoding, decoding, neighbors, etc)
    use module_operators
    ! used in executeCoarsening_tree
    use module_helpers
    ! used for multigrid solver
    use module_fft
    use module_poisson
    ! if the threshold_mask option is used, then the mesh module needs to create the mask function here
    ! hence we require the metamodule to be used.
    use module_physics_metamodule

    implicit none

contains

#include "securityZone_tree.f90"
#include "coarseExtensionUpdate_tree.f90"
#include "refineToEquidistant_tree.f90"
#include "InputOutput_Flusi.f90"
#include "InputOutput.f90"
#include "create_active_and_sorted_lists.f90"
#include "createMask_tree.f90"
#include "block_xfer_nonblocking.f90"
#include "updateNeighbors_tree.f90"
#include "find_neighbors.f90"
#include "doesBlockExist_tree.f90"
#include "refine_tree.f90"
#include "respectJmaxJmin_tree.f90"
#include "refinementExecute.f90"
#include "adapt_tree.f90"
#include "coarseningIndicator_tree.f90"
#include "ensureGradedness_tree.f90"
#include "ensure_completeness_block.f90"
#include "executeCoarsening_tree.f90"
#include "merge_blocks.f90"
#include "balanceLoad_tree.f90"
#include "treecode_to_sfc_id_2D.f90"
#include "treecode_to_sfc_id_3D.f90"
#include "treecode_to_hilbertcode_2D.f90"
#include "treecode_to_hilbertcode_3D.f90"
#include "get_block_spacing_origin.f90"
#include "updateFamily_tree.f90"
#include "find_family.f90"
#include "ActiveLevel_tree.f90"
#include "get_free_local_light_id.f90"
#include "quicksort.f90"
#include "updateMetadata_tree.f90"
#include "createEquidistantGrid_tree.f90"
#include "createRandomGrid_tree.f90"
#include "reset_tree.f90"
#include "allocate_forest.f90"
#include "write_block_distribution.f90"
#include "check_lgt_block_synchronization.f90"
#include "remove_nonperiodic_neighbors.f90"
#include "forest.f90"
#include "multigrid_vcycle.f90"

! former module_initialization
#include "setInitialCondition_tree.f90"
#include "setInicondBlocks_tree.f90"
! previously in module_params (but then cannot call setup_wavelet, c'est chiant.)
#include "ini_file_to_params.f90"

! analytical data used for wavelet compression test
subroutine set_block_testing_data(params, u, x0, dx)
    use module_globals

    implicit none

    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    real(kind=rk), intent(inout) :: u(:,:,:)

    integer(kind=ik) :: ix, iy, iz, Bs(1:3), g
    real(kind=rk) :: sigma0, x00, y00, z00, ampli, x, y, z

    Bs = params%Bs
    g  = params%g
    
    sigma0 = 0.3_rk / 15.0_rk
    x00 = params%domain_size(1)/2.0_rk
    y00 = params%domain_size(2)/2.0_rk
    z00 = params%domain_size(3)/2.0_rk
    ampli = 4.0_rk
    
    if (params%dim==2) then
        ! create gauss pulse. 
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! compute x,y coordinates from spacing and origin
                x = real(ix-(g+1), kind=rk) * dx(1) + x0(1) - x00
                y = real(iy-(g+1), kind=rk) * dx(2) + x0(2) - y00
                
                if (x<-params%domain_size(1)/2.0) x = x + params%domain_size(1)
                if (x>params%domain_size(1)/2.0) x = x - params%domain_size(1)
                
                if (y<-params%domain_size(2)/2.0) y = y + params%domain_size(2)
                if (y>params%domain_size(2)/2.0) y = y - params%domain_size(2)
                
                ! set actual inicond gauss blob
                ! here only for the pressure.
                u(ix,iy,:) = 1.0_rk + ampli*exp( -( (x)**2 + (y)**2 ) / (2.0*sigma0**2) )
            end do
        end do
    else
        ! create gauss pulse
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! compute x,y coordinates from spacing and origin
                    x = real(ix-(g+1), kind=rk) * dx(1) + x0(1) - x00
                    y = real(iy-(g+1), kind=rk) * dx(2) + x0(2) - y00
                    z = real(iz-(g+1), kind=rk) * dx(3) + x0(3) - z00
                    
                    if (x<-params%domain_size(1)/2.0) x = x + params%domain_size(1)
                    if (x>params%domain_size(1)/2.0) x = x - params%domain_size(1)
                    
                    if (y<-params%domain_size(2)/2.0) y = y + params%domain_size(2)
                    if (y>params%domain_size(2)/2.0) y = y - params%domain_size(2)
                    
                    if (z<-params%domain_size(3)/2.0) z = z + params%domain_size(3)
                    if (z>params%domain_size(3)/2.0) z = z - params%domain_size(3)
                    
                    ! set actual inicond gauss blob
                    u(ix,iy,iz) = 1.0_rk + ampli*exp( -( (x)**2 + (y)**2 + (z)**2 ) / (2.0*sigma0**2) )
                end do
            end do
        end do
    end if


    sigma0 = 0.4_rk
    x00 = params%domain_size(1)/4.1_rk
    y00 = params%domain_size(2)/4.1_rk
    z00 = params%domain_size(3)/4.1_rk
    ampli = 2.0_rk
    
    if (params%dim==2) then
        ! create gauss pulse. 
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                ! compute x,y coordinates from spacing and origin
                x = real(ix-(g+1), kind=rk) * dx(1) + x0(1) - x00
                y = real(iy-(g+1), kind=rk) * dx(2) + x0(2) - y00
                
                if (x<-params%domain_size(1)/2.0) x = x + params%domain_size(1)
                if (x>params%domain_size(1)/2.0) x = x - params%domain_size(1)
                
                if (y<-params%domain_size(2)/2.0) y = y + params%domain_size(2)
                if (y>params%domain_size(2)/2.0) y = y - params%domain_size(2)
                
                ! set actual inicond gauss blob
                ! here only for the pressure.
                u(ix,iy,:) = u(ix,iy,:) + ampli*exp( -( (x)**2 + (y)**2 ) / (2.0*sigma0**2) )
            end do
        end do
    else
        ! create gauss pulse
        do iz = g+1, Bs(3)+g
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    ! compute x,y coordinates from spacing and origin
                    x = real(ix-(g+1), kind=rk) * dx(1) + x0(1) - x00
                    y = real(iy-(g+1), kind=rk) * dx(2) + x0(2) - y00
                    z = real(iz-(g+1), kind=rk) * dx(3) + x0(3) - z00
                    
                    if (x<-params%domain_size(1)/2.0) x = x + params%domain_size(1)
                    if (x>params%domain_size(1)/2.0) x = x - params%domain_size(1)
                    
                    if (y<-params%domain_size(2)/2.0) y = y + params%domain_size(2)
                    if (y>params%domain_size(2)/2.0) y = y - params%domain_size(2)
                    
                    if (z<-params%domain_size(3)/2.0) z = z + params%domain_size(3)
                    if (z>params%domain_size(3)/2.0) z = z - params%domain_size(3)
                    
                    ! set actual inicond gauss blob
                    u(ix,iy,iz) = u(ix,iy,iz) + ampli*exp( -( (x)**2 + (y)**2 + (z)**2 ) / (2.0*sigma0**2) )
                end do
            end do
        end do
    end if
end subroutine


    logical function runtime_control_stop(  )
        ! reads runtime control command
        use module_ini_files_parser_mpi
        implicit none
        character(len=cshort) :: command
        character(len=cshort) :: file
        type(inifile) :: CTRL_FILE
        logical :: exists
        integer :: mpirank, mpicode

        file ="runtime_control"
        call MPI_Comm_rank(WABBIT_COMM, mpirank, mpicode)

        if (mpirank==0) then
            inquire(file=file, exist=exists)
            if (.not. exists) then
                call Initialize_runtime_control_file()
            endif
        endif

        call MPI_BCAST( exists, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )
        if (.not. exists) then
            runtime_control_stop = .false.
            return
        endif

        ! root reads in the control file
        ! and fetched the command
        call read_ini_file_mpi( CTRL_FILE, file, .false. ) ! false = non-verbose
        call read_param_mpi(CTRL_FILE, "runtime_control","runtime_control", command, "none")
        call clean_ini_file_mpi( CTRL_FILE )

        if (command == "save_stop") then
            runtime_control_stop = .true.
        else
            runtime_control_stop = .false.
        endif

    end function runtime_control_stop


end module module_mesh
