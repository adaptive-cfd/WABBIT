##################################################################
##################################################################
#		MAKEFILE - Fortran Routines
##################################################################
#################################################################

# Makefile for WABBIT code, adapted from pseudospectators/FLUSI and pseudospectators/UP2D
# Non-module Fortran files to be compiled:
FFILES = treecode_size.f90 array_compare.f90 \
proc_to_lgt_data_start_id.f90 lgt_id_to_hvy_id.f90 hvy_id_to_lgt_id.f90 lgt_id_to_proc_rank.f90 \
f_xy_2D.f90 f_xyz_3D.f90 init_random_seed.f90 error_msg.f90 \
startup_conditioner.f90 init_physics_modules.f90 sparse_to_dense.f90 mult_mask.f90 \
compute_vorticity_post.f90 keyvalues.f90 compare_keys.f90 flusi_to_wabbit.f90 post_mean.f90 post_rhs.f90
# Object and module directory:
OBJDIR = OBJ
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = module_precision.f90 module_globals.f90 module_params.f90 module_timing.f90 module_hdf5_wrapper.f90 \
	module_interpolation.f90 module_initialization.f90 module_mesh.f90 module_IO.f90 module_time_step.f90 module_mpi.f90 module_unit_test.f90 \
	module_treelib.f90  module_ini_files_parser.f90  module_ini_files_parser_mpi.f90 \
	module_indicators.f90 module_operators.f90 module_navier_stokes.f90 module_ns_penalization.f90 \
	module_physics_metamodule.f90 module_ACM.f90 module_ConvDiff_new.f90 module_bridge_interface.f90 \
	module_bridge.f90 module_navier_stokes_params.f90 module_helpers.f90 module_insects_integration_flusi_wabbit.f90 \
	module_insects.f90 module_boundary_conditions.f90 module_funnel.f90 module_navier_stokes_cases.f90\
	module_simple_geometry.f90 module_shock.f90 module_pipe_flow.f90 module_sparse_operators.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = LIB
VPATH += :LIB/MAIN:LIB/MODULE:LIB/INI:LIB/HELPER:LIB/MESH:LIB/IO:LIB/TIME:LIB/EQUATION:LIB/MPI:LIB/TIMING
VPATH += :LIB/PARAMS:LIB/TREE:LIB/INDICATORS:LIB/GEOMETRY:LIB/EQUATION/ACMnew
VPATH += :LIB/OPERATORS:LIB/EQUATION/convection-diffusion:LIB/POSTPROCESSING:LIB/EQUATION/navier_stokes
VPATH += :LIB/EQUATION/navier_stokes:LIB/EQUATION/navier_stokes/case_study:LIB/MPI/BRIDGE
VPATH += :LIB/EQUATION/insects:LIB/BOUNDARYCONDITIONS

# Set the default compiler if it's not already set
ifndef $(FC)
FC = mpif90
endif


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
# NOHDF5 : if set, all dependency of HDF5 is removed. The code can not save output anymore.
#          useful for developmend on some machines only. use only if absolutely necessary
#          NOT IMPLEMENTED YET
# BLOCKINGSENDRECV

#-------------------------------------------------------------------------------
# COMPILER-DEPENDEND PART
#-------------------------------------------------------------------------------
# GNU compiler
#-------------------------------------------------------------------------------
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
FFLAGS += -O3 -ffree-line-length-none
PPFLAG = -cpp # preprocessor flag
#LDFLAGS = -llapack
# Debug flags for gfortran:
FFLAGS += -Wuninitialized -fimplicit-none -fbounds-check -g -ggdb -pedantic
FFLAGS += -Wall -Wextra -Wconversion -g3 -fbacktrace -ffpe-trap=zero,invalid -finit-real=nan -finit-integer=-99999
FFLAGS += -Wno-unused-variable -Wno-unused-parameter -Wno-unused-dummy-argument # -Wno-unused-function
# HDF_ROOT is set in environment. NOTE: it is an TNT@Tu-berlin oddity that libraries are compiled
# to lib64/ and not lib/ like on all other systems. As a workaround, we use BOTH as linkdirs here.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -L$(HDF_LIB)64 $(SB_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm -lblas -llapack
FFLAGS += -I$(HDF_INC) $(SB_INCL)
# for GNU/gfortran, use -D for example: "PRAGMAS=-DTEST" will turn "#ifdef TEST" to true in the code
# different pragmas are space-separated
PRAGMAS = #-DSBLAS #-DBLOCKINGSENDRECV
endif

#-------------------------------------------------------------------------------
# Intel compiler
#-------------------------------------------------------------------------------
mpif90:=$(shell $(FC) --version | head -c 5)
ifeq ($(mpif90),ifort)
PPFLAG = -fpp # preprocessor flag
FFLAGS += -FR -O3 -heap-arrays
# debug flags: attention they might disable all optimization!
# FFLAGS += -warn all,nounused -traceback -check bounds -debug all -check all,noarg_temp_created
FFLAGS += -module $(OBJDIR) # specify directory for modules.
LDFLAGS = -L/usr/X11/lib/ -lX11 #-L/usr/lib64/lapack -llapack
# HDF_ROOT is set in environment. NOTE: it is an TNT@Tu-berlin oddity that libraries are compiled
# to lib64/ and not lib/ like on all other systems. As a workaround, we use BOTH as linkdirs here.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -L$(HDF_LIB)64 -lhdf5_fortran -lhdf5 -lz -ldl -lm -llapack -lblas
FFLAGS += -I$(HDF_INC)
# for intel, use -D for example: PRAGMAS=-DIFORT will turn #ifdef IFORT to true in the code
# different pragmas are space-separated
PRAGMAS = # -DBLOCKINGSENDRECV
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
PRAGMAS = -WF,-DBLOCKINGSENDRECV
endif

# add the PRAGMAS to FFLAGS: (for all compilers)
FFLAGS += $(PPFLAG) $(PRAGMAS)


# Both programs are compiled by default.
all: directories wabbit wabbit-post #doc

# Compile main programs, with dependencies.
wabbit: main.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).


# compile precision module
$(OBJDIR)/module_precision.o: module_precision.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# compile module containing all globals
$(OBJDIR)/module_globals.o: module_globals.f90 $(OBJDIR)/module_precision.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_bridge.o: module_bridge.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_insects_integration_flusi_wabbit.o: module_insects_integration_flusi_wabbit.f90  $(OBJDIR)/module_precision.o \
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_insects.o: module_insects.f90 $(OBJDIR)/module_insects_integration_flusi_wabbit.o \
	body_geometry.f90 body_motion.f90 rigid_solid_time_stepper.f90 wings_geometry.f90 \
	wings_motion.f90 stroke_plane.f90 pointcloud.f90 fractal_trees.f90 insect_init_clean.f90 \
	kineloader.f90 active_grid_winglets.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ini_files_parser.o: module_ini_files_parser.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_bridge.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_params.o: module_params.f90 $(OBJDIR)/module_ini_files_parser_mpi.o \
	ini_file_to_params.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_sparse_operators.o: module_sparse_operators.f90 $(OBJDIR)/module_params.o \
	$(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_boundary_conditions.o: module_boundary_conditions.f90 \
		$(OBJDIR)/module_params.o $(OBJDIR)/module_globals.o $(OBJDIR)/module_treelib.o
		$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_bridge_interface.o: module_bridge_interface.f90 $(OBJDIR)/module_treelib.o $(OBJDIR)/module_mesh.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_navier_stokes_params.o: module_navier_stokes_params.f90 $(OBJDIR)/module_globals.o\
	$(OBJDIR)/module_helpers.o $(OBJDIR)/module_ini_files_parser_mpi.o \
	inicond_NStokes.f90 filter_block.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ns_penalization.o: module_ns_penalization.f90 $(OBJDIR)/module_navier_stokes_params.o\
	$(OBJDIR)/module_ini_files_parser_mpi.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_funnel.o: module_funnel.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_ns_penalization.o \
	funnel2D.f90 funnel3D.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_pipe_flow.o: module_pipe_flow.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_ns_penalization.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_simple_geometry.o: module_simple_geometry.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_ns_penalization.o
		$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_navier_stokes_cases.o: module_navier_stokes_cases.f90 $(OBJDIR)/module_funnel.o $(OBJDIR)/module_ns_penalization.o\
	$(OBJDIR)/module_shock.o $(OBJDIR)/module_simple_geometry.o $(OBJDIR)/module_pipe_flow.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_navier_stokes.o: module_navier_stokes.f90 $(OBJDIR)/module_ns_penalization.o\
	$(OBJDIR)/module_navier_stokes_params.o $(OBJDIR)/module_sparse_operators.o \
	$(OBJDIR)/module_navier_stokes_cases.o $(OBJDIR)/module_funnel.o RHS_2D_navier_stokes_periodic.f90\
	RHS_2D_navier_stokes_bc.f90 RHS_3D_navier_stokes.f90 RHS_2D_cylinder.f90 inicond_NStokes.f90 save_data_ns.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_shock.o: module_shock.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_ns_penalization.o \
		$(OBJDIR)/module_navier_stokes_params.o
		$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ACM.o: module_ACM.f90 rhs.f90 create_mask.f90 sponge.f90 save_data_ACM.f90 update_grid_qtys_ACM.f90 \
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_globals.o \
	$(OBJDIR)/module_helpers.o $(OBJDIR)/module_insects.o statistics_ACM.f90 inicond_ACM.f90 filter_ACM.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ConvDiff_new.o: module_ConvDiff_new.f90 rhs_convdiff.f90 \
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_globals.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_timing.o: module_timing.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ini_files_parser_mpi.o: module_ini_files_parser_mpi.f90 $(OBJDIR)/module_globals.o $(OBJDIR)/module_ini_files_parser.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_interpolation.o: module_interpolation.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_physics_metamodule.o: module_physics_metamodule.f90 $(OBJDIR)/module_globals.o \
	$(OBJDIR)/module_ConvDiff_new.o $(OBJDIR)/module_navier_stokes.o $(OBJDIR)/module_ACM.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_hdf5_wrapper.o: module_hdf5_wrapper.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_globals.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_initialization.o: module_initialization.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o \
	$(OBJDIR)/module_mesh.o $(OBJDIR)/module_IO.o $(OBJDIR)/module_physics_metamodule.o \
	set_initial_grid.f90 \
	set_inicond_blocks.f90 get_inicond_from_file.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_mpi.o: module_mpi.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_interpolation.o $(OBJDIR)/module_treelib.o\
	synchronize_ghosts.f90 blocks_per_mpirank.f90 reset_ghost_nodes.f90 synchronize_lgt_data.f90 check_redundant_nodes.f90 \
	restrict_predict_data.f90 calc_data_bounds.f90 synchronize_ghosts_generic.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_time_step.o: module_time_step.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_mpi.o\
	$(OBJDIR)/module_mesh.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_physics_metamodule.o \
	calculate_time_step.f90 time_stepper.f90 set_RK_input.f90 RHS_wrapper.f90 final_stage_RK.f90 \
	statistics_wrapper.f90 filter_wrapper.f90 krylov.f90 update_grid_qtys.f90 $(OBJDIR)/module_boundary_conditions.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_indicators.o: module_indicators.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_operators.o \
	refinement_indicator.f90 block_coarsening_indicator.f90 threshold_block.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_helpers.o: module_helpers.f90 $(OBJDIR)/module_globals.o most_common_element.f90 $(OBJDIR)/module_ini_files_parser_mpi.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_mesh.o: module_mesh.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o $(OBJDIR)/module_interpolation.o \
	$(OBJDIR)/module_mpi.o $(OBJDIR)/module_treelib.o $(OBJDIR)/module_indicators.o $(OBJDIR)/module_physics_metamodule.o \
	$(OBJDIR)/module_boundary_conditions.o $(OBJDIR)/module_helpers.o update_neighbors_2D.f90 find_neighbor_edge_2D.f90 does_block_exist.f90 \
	find_neighbor_corner_2D.f90 refine_mesh.f90 respect_min_max_treelevel.f90 refinement_execute_2D.f90 adapt_mesh.f90 threshold_block.f90 \
	ensure_gradedness.f90 ensure_completeness.f90 coarse_mesh.f90 balance_load.f90 set_desired_num_blocks_per_rank.f90 \
	compute_friends_table.f90 compute_affinity.f90 treecode_to_sfc_id_2D.f90 treecode_to_sfc_id_3D.f90 treecode_to_hilbertcode_2D.f90 \
    treecode_to_hilbertcode_3D.f90 update_neighbors_3D.f90 find_neighbor_face_3D.f90 find_neighbor_edge_3D.f90 find_neighbor_corner_3D.f90 \
    refinement_execute_3D.f90 get_block_spacing_origin.f90 update_neighbors.f90 \
	find_sisters.f90 max_active_level.f90 min_active_level.f90 get_free_local_light_id.f90 \
	merge_blocks.f90 create_active_and_sorted_lists.f90 quicksort.f90 grid_coarsening_indicator.f90 \
	create_equidistant_grid.f90 create_random_grid.f90 allocate_grid.f90 reset_grid.f90 block_xfer_nonblocking.f90 block_xfer_blocking.f90 \
		check_lgt_block_synchronization.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_unit_test.o: module_unit_test.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_initialization.o $(OBJDIR)/module_mesh.o $(OBJDIR)/module_time_step.o \
	unit_test_ghost_nodes_synchronization.f90 unit_test_treecode.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_treelib.o: module_treelib.f90 $(OBJDIR)/module_params.o get_neighbor_treecode.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_IO.o: module_IO.f90 $(OBJDIR)/module_mesh.o $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o \
	$(OBJDIR)/module_hdf5_wrapper.o $(OBJDIR)/module_mpi.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_physics_metamodule.o \
	save_data.f90 write_field.f90 read_field.f90 \
	read_mesh.f90 read_attributes.f90 read_file_flusi.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_operators.o: module_operators.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o \
	$(OBJDIR)/module_helpers.o volume_integral.f90 compute_vorticity.f90 divergence.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_sparse_operators.o: module_sparse_operators.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_timing.o \
	$(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

wabbit-post: main_post.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
clean-fortran:
	rm -rf $(PROGRAMS) $(OBJDIR) a.out wabbit wabbit-post
