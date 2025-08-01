# Makefile for WABBIT code, adapted from pseudospectators/FLUSI and pseudospectators/UP2D
# Non-module Fortran files to be compiled:
FFILES = treecode_size.f90 array_compare.f90 \
proc_to_lgt_data_start_id.f90 lgt2hvy.f90 hvy2lgt.f90 lgt2proc.f90 init_random_seed.f90 \
init_physics_modules.f90 sparse_to_dense.f90 dense_to_sparse.f90 mult_mask.f90 compute_vorticity_post.f90 compute_scalar_field_post.f90 \
keyvalues.f90 compare_keys.f90 flusi_to_wabbit.f90 post_mean.f90 post_rhs.f90 post_stl2dist.f90 post_add_two_masks.f90 post_prune_tree.f90 \
post_average_snapshots.f90 post_superstl.f90 post_dry_run.f90 performance_test.f90 adaption_test.f90 post_generate_forest.f90 \
post_dump_neighbors.f90 operator_reconstruction.f90 rhs_operator_reconstruction.f90 post_filtertest.f90 post_extract_slice.f90 \
post_evaluate_thresholding.f90 post_unit_test.f90 post_compression_unit_test.f90 post_denoising.f90 post_cvs_invertibility_test.f90 \
post_wavelet_transform.f90 post_denoising_test.f90 post_derivative.f90 proto_GS_multigrid.f90
# Object and module directory:
OBJDIR = OBJ
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = module_forestMetaData.f90 module_globals.f90 module_params.f90 module_timing.f90 module_hdf5_wrapper.f90 \
	module_wavelets.f90 module_mesh.f90 module_time_step.f90 module_mpi.f90 module_unit_test.f90 \
	module_treelib.f90  module_ini_files_parser.f90  module_ini_files_parser_mpi.f90 \
	module_indicators.f90 module_operators.f90 module_navier_stokes.f90 module_ns_penalization.f90 \
	module_physics_metamodule.f90 module_ACM.f90 module_ConvDiff_new.f90 module_bridge_interface.f90 \
	module_bridge.f90 module_navier_stokes_params.f90 module_helpers.f90 module_insects_integration_flusi_wabbit.f90 \
	module_insects.f90 module_funnel.f90 module_navier_stokes_cases.f90\
	module_simple_geometry.f90 module_shock.f90 module_pipe_flow.f90\
	module_MOR.f90 module_sparse_operators.f90 module_stl_file_reader.f90\
	module_t_files.f90 module_saving.f90 module_fft.f90 module_poisson.f90 # module_sync_ghosts_redundantGrid.f90 module_sync_ghosts_uniqueGrid.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = LIB
VPATH += :LIB/MAIN:LIB/MODULE:LIB/INI:LIB/HELPER:LIB/MESH:LIB/IO:LIB/TIME:LIB/EQUATION:LIB/MPI:LIB/TIMING
VPATH += :LIB/PARAMS:LIB/TREE:LIB/INDICATORS:LIB/GEOMETRY:LIB/EQUATION/ACMnew:LIB/TESTING
VPATH += :LIB/OPERATORS:LIB/EQUATION/convection-diffusion:LIB/POSTPROCESSING:LIB/EQUATION/navier_stokes
VPATH += :LIB/EQUATION/navier_stokes:LIB/EQUATION/navier_stokes/case_study:LIB/MPI/BRIDGE
VPATH += :LIB/EQUATION/insects:LIB/BOUNDARYCONDITIONS:LIB/WAVELETS:LIB/FFT:LIB/POISSON

# Set the default compiler if it's not already set
ifndef $(FC)
FC = mpif90
endif

ifndef HDF_SOURCE
HDF_SOURCE = $(HDF_ROOT)
endif

# to print the version number at each run. version number == git hash ID of current commit
GIT_HASH := $(shell git rev-parse HEAD)
BUILD_DATE := $(shell date)

#Place of Sparse BLAS objects
SB_LIB = #-L../../sblas/SOFTWARE -lSparseBLAS_GNU
#Place of Sparse BLAS modules
SB_INCL = #-I../../sblas/SOFTWARE
#-------------------------------------------------------------------------------
# PRAGMAS part.
#-------------------------------------------------------------------------------
# we really should avoid as many PRAGMAS as possible, but some things just cannot
# be done without them. Use them ONLY for machine or compiler dependent settings
# AND NEVER FOR PHYSICS MODULES (e.g. switching the viscosity)
#
# PRAGMAS are added in the compiler specific section, but here is their list:
# DEV : if set, additional tests are done and error messages are issued; these concern
#       errors in development and are thus useful only for people fiddling around with
#       the code.

#-------------------------------------------------------------------------------
# COMPILER-DEPENDEND PART
#-------------------------------------------------------------------------------
# GNU compiler
#-------------------------------------------------------------------------------
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
FFLAGS += -ffree-line-length-none
PPFLAG = -cpp # preprocessor flag
#LDFLAGS = -llapack
FFLAGS += -O3
FFLAGS += -Wuninitialized -fimplicit-none -fbounds-check -pedantic
FFLAGS += -Wall -Wextra -Wconversion -fbacktrace -ffpe-trap=zero,invalid 
FFLAGS += -g3 -g -ggdb  # debugging flags
FFLAGS += -finit-local-zero -finit-real=snan -finit-integer=-99999
FFLAGS += -Wno-unused-variable -Wno-unused-parameter -Wno-unused-dummy-argument # -Wno-unused-function
# HDF_ROOT is set in environment. NOTE: it is an TNT@Tu-berlin oddity that libraries are compiled
# to lib64/ and not lib/ like on all other systems. As a workaround, we use BOTH as linkdirs here.
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_ROOT)/lib -L$(HDF_SOURCE)/fortran/src/.libs -L$(HDF_SOURCE)/fortran/src -lhdf5_fortran -lhdf5
LDFLAGS += -L$(FFT_ROOT)/lib -lfftw3 $(SB_LIB) -lz
FFLAGS += -I$(HDF_ROOT)/include $(SB_INCL) -I$(HDF_SOURCE)/fortran/src $(SB_INCL)
FFLAGS += -I$(FFT_ROOT)/include
# for GNU/gfortran, use -D for example: "PRAGMAS=-DTEST" will turn "#ifdef TEST" to true in the code
# different pragmas are space-separated
PRAGMAS = #-DSBLAS
ifdef DEV
PRAGMAS += -DDEV
endif
# enable / disable FFT depending on the module
ifdef FFT_ROOT
PRAGMAS += -DFFT_ROOT
endif
ifdef MKLROOT
# Use MKL lapack
LDFLAGS += -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl -lm
else
# Use custom LAPACK installation
LDFLAGS += -L${WABBIT_BLAS_ROOT} -L${WABBIT_LAPACK_ROOT} -llapack -lblas -ldl -lm
endif
endif




#-------------------------------------------------------------------------------
# Intel compiler
#-------------------------------------------------------------------------------
mpif90:=$(shell $(FC) --version | head -c 5)
ifeq ($(mpif90),ifort)
PPFLAG = -fpp # preprocessor flag
FFLAGS += -O3 -FR -heap-arrays
##FFLAGS += -xcore-avx2  # use avx2
# timing flags: attention they might disable all optimization!
##FFLAGS += -g -warn all,nounused -traceback -check bounds -check all,noarg_temp_created
FFLAGS += -module $(OBJDIR) # specify directory for modules.
LDFLAGS = -L/usr/X11/lib/ -lX11 #-L/usr/lib64/lapack -llapack
# HDF_ROOT is set in environment. NOTE: it is an TNT@Tu-berlin oddity that libraries are compiled
# to lib64/ and not lib/ like on all other systems. As a workaround, we use BOTH as linkdirs here.
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_ROOT)/lib -L$(HDF_SOURCE)/fortran/src/.libs -L$(HDF_SOURCE)/fortran/src -lhdf5_fortran -lhdf5
LDFLAGS += -L$(FFT_ROOT)/lib -lfftw3 $(SB_LIB) -lz
FFLAGS += -I$(HDF_ROOT)/include $(SB_INCL) -I$(HDF_SOURCE)/fortran/src $(SB_INCL)
FFLAGS += -I$(FFT_ROOT)/include
# for intel, use -D for example: PRAGMAS=-DIFORT will turn #ifdef IFORT to true in the code
# different pragmas are space-separated
PRAGMAS = #
ifdef DEV
PRAGMAS += -DDEV
endif
# enable / disable FFT depending on the module
ifdef FFT_ROOT
PRAGMAS += -DFFT_ROOT
endif
ifdef MKLROOT
# Use MKL lapack
LDFLAGS += -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl -lm
else
# Use custom LAPACK installation
LDFLAGS += -L${WABBIT_BLAS_ROOT} -L${WABBIT_LAPACK_ROOT} -llapack -lblas -ldl -lm
endif
endif

#-------------------------------------------------------------------------------
# IBM compiler
#-------------------------------------------------------------------------------
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
#the xlf compiler is complaining because the arguments to MAX() in the
#code below are of two different types; ct0 is double precision and -0.999
#is single (-0.999d0 would be double). If you want xlf to auto-promote
#single precision constants to doubles, add the compiler option "-qdpc".
FFLAGS += -O4 -qdpc -qmoddir=$(OBJDIR)
FFLAGS += -qnohot
#FFLAGS += -C -qflttrap=overflow:underflow:nanq:zerodivide:qpxstore:enable -qsigtrap -g9 -qnoopt -qfullpath
FFLAGS += -I$(OBJDIR)
PPFLAG=-qsuffix=cpp=f90  #preprocessor flag
# for IBMXLF95 PRAGMAS=-WF,-DIFORT will turn #ifdef IFORT to true in the code
# here different PRAGMAS are comma separated (NO SPACES!!!)
# NOTE first pragma (if any is used) MUST be -WF,
PRAGMAS = #-WF,-DTEST
endif

# add the PRAGMAS to FFLAGS: (for all compilers)
FFLAGS += $(PPFLAG) $(PRAGMAS) -fPIC

.PHONY: version.f90

# Both programs are compiled by default.
all: version.f90 directories wabbit wabbit-post #doc

version.f90: 
	@echo 'character(len=*), parameter :: git_version = "$(GIT_HASH)"' > LIB/MAIN/version.f90
	@echo 'character(len=*), parameter :: build_date = "$(BUILD_DATE)"' >> LIB/MAIN/version.f90

# Compile main programs, with dependencies.
wabbit: main.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

wabbit-post: main_post.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).
# compile module containing all globals
$(OBJDIR)/module_globals.o: module_globals.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_t_files.o: module_t_files.f90 $(OBJDIR)/module_globals.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_stl_file_reader.o: module_stl_file_reader.f90 $(OBJDIR)/module_globals.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_bridge.o: module_bridge.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_insects_integration_flusi_wabbit.o: module_insects_integration_flusi_wabbit.f90  $(OBJDIR)/module_globals.o \
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_insects.o: module_insects.f90 $(OBJDIR)/module_insects_integration_flusi_wabbit.o \
	body_geometry.f90 body_motion.f90 rigid_solid_time_stepper.f90 wings_geometry.f90 \
	wings_motion.f90 stroke_plane.f90 pointcloud.f90 fractal_trees.f90 insect_init_clean.f90 \
	kineloader.f90 active_grid_winglets.f90 $(OBJDIR)/module_t_files.o $(OBJDIR)/module_stl_file_reader.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ini_files_parser.o: module_ini_files_parser.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_bridge.o $(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_params.o: module_params.f90 $(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_globals.o $(OBJDIR)/module_t_files.o \
	ini_file_to_params.f90 $(OBJDIR)/module_helpers.o $(OBJDIR)/module_t_files.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_bridge_interface.o: module_bridge_interface.f90 $(OBJDIR)/module_treelib.o $(OBJDIR)/module_mesh.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_navier_stokes_params.o: module_navier_stokes_params.f90 $(OBJDIR)/module_globals.o\
	$(OBJDIR)/module_helpers.o $(OBJDIR)/module_ini_files_parser_mpi.o \
	inicond_NStokes.f90 filter_block.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ns_penalization.o: module_ns_penalization.f90 $(OBJDIR)/module_navier_stokes_params.o\
	$(OBJDIR)/module_ini_files_parser_mpi.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_funnel.o: module_funnel.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_ns_penalization.o \
	funnel2D.f90 funnel3D.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_pipe_flow.o: module_pipe_flow.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_ns_penalization.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_simple_geometry.o: module_simple_geometry.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_ns_penalization.o
		$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_navier_stokes_cases.o: module_navier_stokes_cases.f90 $(OBJDIR)/module_funnel.o $(OBJDIR)/module_ns_penalization.o\
	$(OBJDIR)/module_shock.o $(OBJDIR)/module_simple_geometry.o $(OBJDIR)/module_pipe_flow.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_navier_stokes.o: module_navier_stokes.f90  $(OBJDIR)/module_ns_penalization.o  $(OBJDIR)/module_navier_stokes_params.o\
	$(OBJDIR)/module_navier_stokes_cases.o $(OBJDIR)/module_sparse_operators.o  $(OBJDIR)/module_operators.o \
	$(OBJDIR)/module_funnel.o  RHS_2D_navier_stokes_periodic.f90 RHS_2D_navier_stokes_bc.f90  $(OBJDIR)/module_t_files.o\
	RHS_2D_cylinder.f90 RHS_3D_navier_stokes.f90 inicond_NStokes.f90 save_data_ns.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_shock.o: module_shock.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_ns_penalization.o \
		$(OBJDIR)/module_navier_stokes_params.o $(OBJDIR)/module_params.o
		$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ACM.o: module_ACM.f90 rhs_ACM.f90 create_mask.f90 sponge.f90 2D_wingsection.f90 save_data_ACM.f90 \
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_globals.o $(OBJDIR)/module_t_files.o \
	$(OBJDIR)/module_helpers.o $(OBJDIR)/module_insects.o statistics_ACM.f90 inicond_ACM.f90 boundcond_ACM.f90 filter_ACM.f90 $(OBJDIR)/module_params.o \
	$(OBJDIR)/module_t_files.o $(OBJDIR)/module_timing.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ConvDiff_new.o: module_ConvDiff_new.f90 rhs_convdiff.f90 statistics_convdiff.f90 \
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_globals.o $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_timing.o: module_timing.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ini_files_parser_mpi.o: module_ini_files_parser_mpi.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_ini_files_parser.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_wavelets.o: module_wavelets.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_globals.o $(OBJDIR)/module_treelib.o \
	conversion_routines.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_fft.o: module_fft.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_globals.o $(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_poisson.o: module_poisson.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_globals.o $(OBJDIR)/module_helpers.o $(OBJDIR)/module_forestMetaData.o \
	$(OBJDIR)/module_treelib.o $(OBJDIR)/module_timing.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_physics_metamodule.o: module_physics_metamodule.f90 $(OBJDIR)/module_globals.o \
	$(OBJDIR)/module_ConvDiff_new.o $(OBJDIR)/module_navier_stokes.o $(OBJDIR)/module_ACM.o $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_hdf5_wrapper.o: module_hdf5_wrapper.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_globals.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_mpi.o: module_mpi.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_wavelets.o \
	$(OBJDIR)/module_treelib.o $(OBJDIR)/module_forestMetaData.o $(OBJDIR)/module_physics_metamodule.o blocks_per_mpirank.f90 reset_ghost_nodes.f90 synchronize_lgt_data.f90 bound_cond_generic.f90 \
	restrict_predict_data.f90 calc_data_bounds.f90 synchronize_ghosts_generic.f90 reconstruction_step.f90 \
	xfer_block_data.f90 block_relations.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_time_step.o: module_time_step.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_mpi.o \
	$(OBJDIR)/module_mesh.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_physics_metamodule.o \
	calculate_time_step.f90 timeStep_tree.f90 RHS_wrapper.f90	 \
	statistics_wrapper.f90 filter_wrapper.f90 krylov.f90  $(OBJDIR)/module_treelib.o \
	runge_kutta_generic.f90 runge_kutta_generic_FSI.f90 runge_kutta_chebychev.f90 runge_kutta_chebychev_FSI.f90 $(OBJDIR)/module_t_files.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_indicators.o: module_indicators.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_operators.o \
	$(OBJDIR)/module_mpi.o $(OBJDIR)/module_wavelets.o refinementIndicator_tree.f90 coarseningIndicator_block.f90 threshold_block.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_helpers.o: module_helpers.f90 $(OBJDIR)/module_globals.o most_common_element.f90 $(OBJDIR)/module_ini_files_parser_mpi.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_mesh.o: module_mesh.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_wavelets.o \
	$(OBJDIR)/module_mpi.o $(OBJDIR)/module_treelib.o $(OBJDIR)/module_physics_metamodule.o $(OBJDIR)/module_indicators.o \
	$(OBJDIR)/module_helpers.o $(OBJDIR)/module_params.o $(OBJDIR)/module_forestMetaData.o $(OBJDIR)/module_poisson.o $(OBJDIR)/module_fft.o \
	refineToEquidistant_tree.f90 \
	InputOutput_Flusi.f90 InputOutput.f90 create_active_and_sorted_lists.f90 createMask_tree.f90 block_xfer_nonblocking.f90 \
	updateNeighbors_tree.f90 find_neighbors.f90 doesBlockExist_tree.f90 refine_tree.f90 respectJmaxJmin_tree.f90 \
	refinementExecute.f90 adapt_tree.f90 coarseningIndicator_tree.f90 ensureGradedness_tree.f90 \
	ensure_completeness_block.f90 executeCoarsening_tree.f90 merge_blocks.f90 balanceLoad_tree.f90 \
	treecode_to_sfc_id_2D.f90 treecode_to_sfc_id_3D.f90 treecode_to_hilbertcode_2D.f90 treecode_to_hilbertcode_3D.f90 get_block_spacing_origin.f90 \
	find_family.f90 ActiveLevel_tree.f90 get_free_local_light_id.f90 quicksort.f90 updateMetadata_tree.f90 createEquidistantGrid_tree.f90 \
	createRandomGrid_tree.f90 reset_tree.f90 allocate_forest.f90 write_block_distribution.f90 check_lgt_block_synchronization.f90 \
	remove_nonperiodic_neighbors.f90 forest.f90 multigrid_vcycle.f90 \
	securityZone_tree.f90 coarseExtensionUpdate_tree.f90 updateFamily_tree.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_unit_test.o: module_unit_test.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_mesh.o $(OBJDIR)/module_time_step.o \
	$(OBJDIR)/module_treelib.o unit_test_treecode.f90 unit_test_Sync.f90 unit_test_ghostSync.f90 unit_test_waveletDecomposition.f90 unit_test_refineCoarsen.f90 \
	unit_test_waveletDecomposition_invertibility.f90 createTestGrids.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_treelib.o: module_treelib.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_globals.o neighborhood.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_operators.o: module_operators.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o \
	$(OBJDIR)/module_forestMetaData.o $(OBJDIR)/module_treelib.o compute_derivative.f90 compute_vorticity.f90 divergence.f90 compute_Qcriterion.f90 \
	componentWiseNorm_tree.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_sparse_operators.o: module_sparse_operators.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o \
	$(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_forestMetaData.o: module_forestMetaData.f90 $(OBJDIR)/module_globals.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_saving.o: module_saving.f90 \
	 $(OBJDIR)/module_physics_metamodule.o $(OBJDIR)/module_mesh.o \
	save_data.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_MOR.o: module_MOR.f90  \
	$(OBJDIR)/module_globals.o $(OBJDIR)/module_ini_files_parser.o $(OBJDIR)/module_mesh.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

clean-fortran:
	rm -rf $(PROGRAMS) $(OBJDIR) a.out wabbit wabbit-post

# when running tests in zero there could be traces - make sure they are cleaned
clean-source:
	rm -f *.key *.h5 *.t *.dat *.htg *.vtm success