
# Makefile for WABBIT code, adapted from pseudospectators/FLUSI and pseudospectators/UP2D
# Non-module Fortran files to be compiled:
FFILES = treecode_size.f90 array_compare.f90 \
proc_to_lgt_data_start_id.f90 lgt_id_to_hvy_id.f90 hvy_id_to_lgt_id.f90 lgt_id_to_proc_rank.f90 get_free_light_id.f90 \
f_xy_2D.f90 f_xyz_3D.f90 init_random_seed.f90 error_msg.f90 \
startup_conditioner.f90 init_physics_modules.f90 sparse_to_dense.f90 \
compute_vorticity_post.f90 keyvalues.f90 compare_keys.f90 flusi_to_wabbit.f90 post_mean.f90
# Object and module directory:
OBJDIR = OBJ
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = module_precision.f90 module_params.f90 module_debug.f90 module_hdf5_wrapper.f90 \
	module_interpolation.f90 module_initialization.f90 module_mesh.f90 module_IO.f90 module_time_step.f90 module_mpi.f90 module_unit_test.f90 \
	module_treelib.f90  module_ini_files_parser.f90  module_ini_files_parser_mpi.f90\
	module_indicators.f90 module_operators.f90 module_navier_stokes_new.f90 module_ns_penalization.f90\
	module_physics_metamodule.f90 module_ACM-new.f90 module_ConvDiff_new.f90 module_bridge_interface.f90 module_bridge.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = LIB
VPATH += :LIB/MAIN:LIB/MODULE:LIB/INI:LIB/HELPER:LIB/MESH:LIB/IO:LIB/TIME:LIB/EQUATION:LIB/MPI:LIB/DEBUG
VPATH += :LIB/PARAMS:LIB/TREE:LIB/INDICATORS:LIB/GEOMETRY:LIB/EQUATION/ACMnew
VPATH += :LIB/OPERATORS:LIB/EQUATION/convection-diffusion:LIB/POSTPROCESSING:LIB/EQUATION/navier_stokes
VPATH += :LIB/EQUATION/navier_stokes:LIB/EQUATION/navier_stokes/case_study:LIB/MPI/BRIDGE

# Set the default compiler if it's not already set
ifndef $(FC)
FC = mpif90
endif

#-------------------------------------------------------------------------------
# COMPILER-DEPENDEND PART
#-------------------------------------------------------------------------------
# GNU compiler
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
FFLAGS += -O3 -ffree-line-length-none
PPFLAG= -cpp #preprocessor flag
#LDFLAGS = -llapack
# Debug flags for gfortran:
FFLAGS += -Wuninitialized -fimplicit-none -fbounds-check -g -ggdb #-O
FFLAGS += -Wall -Wextra -Wconversion -g3 -fbacktrace -ffpe-trap=zero,invalid -finit-real=nan
# HDF_ROOT is set in environment. NOTE: it is an TNT@Tu-berlin oddity that libraries are compiled
# to lib64/ and not lib/ like on all other systems. As a workaround, we use BOTH as linkdirs here.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -L$(HDF_LIB)64 -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC)
endif

# Intel compiler
mpif90:=$(shell $(FC) --version | head -c 5)
ifeq ($(mpif90),ifort)
PPFLAG= -fpp
FFLAGS = -FR -O3 -warn all -traceback -check bounds -heap-arrays
FFLAGS += -module $(OBJDIR) # specify directory for modules.
LDFLAGS = -L/usr/X11/lib/ -lX11 #-L/usr/lib64/lapack -llapack
# HDF_ROOT is set in environment. NOTE: it is an TNT@Tu-berlin oddity that libraries are compiled
# to lib64/ and not lib/ like on all other systems. As a workaround, we use BOTH as linkdirs here.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -L$(HDF_LIB)64 -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC)
endif

# Both programs are compiled by default.
all: directories wabbit wabbit-post #doc

# Compile main programs, with dependencies.
wabbit: main.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).

# first compile precision module
$(OBJDIR)/module_precision.o: module_precision.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_bridge.o: module_bridge.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ini_files_parser.o: module_ini_files_parser.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_bridge.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_params.o: module_params.f90 $(OBJDIR)/module_ini_files_parser_mpi.o \
	ini_file_to_params.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_bridge_interface.o: module_bridge_interface.f90 $(OBJDIR)/module_treelib.o $(OBJDIR)/module_mesh.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ns_penalization.o: module_ns_penalization.f90 $(OBJDIR)/module_ini_files_parser_mpi.o\
	funnel.f90 vortex_street.f90 sod_shock_tube.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_navier_stokes_new.o: module_navier_stokes_new.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_ns_penalization.o \
	$(OBJDIR)/module_operators.o RHS_2D_navier_stokes.f90 RHS_3D_navier_stokes.f90 initial_conditions.f90 inicond_shear_layer.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ACM-new.o: module_ACM-new.f90 rhs.f90 create_mask_new.f90 iniconds.f90 sponge_new.f90\
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_precision.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ConvDiff_new.o: module_ConvDiff_new.f90 rhs_convdiff.f90 \
	$(OBJDIR)/module_ini_files_parser_mpi.o $(OBJDIR)/module_precision.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_debug.o: module_debug.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_interpolation.o \
	check_lgt_block_synchronization.f90 write_future_mesh_lvl.f90 write_debug_times.f90 write_block_distribution.f90 write_com_list.f90 \
	write_com_matrix.f90 write_com_matrix_pos.f90 allocate_init_debugging.f90 check_redundant_nodes.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ini_files_parser_mpi.o: module_ini_files_parser_mpi.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_ini_files_parser.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_interpolation.o: module_interpolation.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_physics_metamodule.o: module_physics_metamodule.f90 $(OBJDIR)/module_precision.o \
	$(OBJDIR)/module_ConvDiff_new.o $(OBJDIR)/module_navier_stokes_new.o $(OBJDIR)/module_ACM-new.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_hdf5_wrapper.o: module_hdf5_wrapper.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_initialization.o: module_initialization.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o \
	$(OBJDIR)/module_mesh.o $(OBJDIR)/module_IO.o $(OBJDIR)/module_physics_metamodule.o \
	set_initial_grid.f90 initial_block_distribution_2D.f90 new_block_heavy.f90 \
	inicond_vorticity_filaments.f90 initial_block_distribution_3D.f90 create_equidistant_base_mesh.f90 \
	allocate_grid.f90 reset_grid.f90 \
	set_inicond_blocks.f90 allocate_com_arrays.f90 get_inicond_from_file.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_mpi.o: module_mpi.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_interpolation.o\
	synchronize_ghosts.f90 copy_ghost_nodes_2D.f90 create_send_buffer_2D.f90 write_receive_buffer_2D.f90 create_external_com_list.f90 \
	max_com_num.f90 fill_send_buffer.f90 fill_receive_buffer.f90 RMA_lock_unlock_get_data.f90 RMA_lock_unlock_put_data.f90 \
	isend_irecv_data.f90 copy_ghost_nodes_3D.f90 create_send_buffer_3D.f90 write_receive_buffer_3D.f90 blocks_per_mpirank.f90 \
	reset_ghost_nodes.f90 copy_redundant_nodes_2D.f90 synchronize_lgt_data.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_time_step.o: module_time_step.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_mpi.o \
	$(OBJDIR)/module_mesh.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_physics_metamodule.o \
	filter_block.f90 filter_1D.f90 calculate_time_step.f90 time_stepper.f90 set_RK_input.f90 RHS_wrapper.f90 final_stage_RK.f90 \
	wavelet_filter.f90 get_block_max_velocity_norm.f90 bogey_filter.f90 statistics_wrapper.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_indicators.o: module_indicators.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_mpi.o \
	refinement_indicator.f90 coarsening_indicator.f90 threshold_block.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_mesh.o: module_mesh.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_interpolation.o \
	$(OBJDIR)/module_mpi.o $(OBJDIR)/module_treelib.o $(OBJDIR)/module_indicators.o \
	update_neighbors_2D.f90 find_neighbor_edge_2D.f90 does_block_exist.f90 \
	find_neighbor_corner_2D.f90 refine_mesh.f90 respect_min_max_treelevel.f90 refinement_execute_2D.f90 adapt_mesh.f90 threshold_block.f90 \
	ensure_gradedness.f90 ensure_completeness.f90 coarse_mesh.f90 balance_load.f90 set_desired_num_blocks_per_rank.f90 \
	compute_friends_table.f90 compute_affinity.f90 treecode_to_sfc_id_2D.f90 treecode_to_sfc_id_3D.f90 treecode_to_hilbertcode_2D.f90 \
    treecode_to_hilbertcode_3D.f90 update_neighbors_3D.f90 find_neighbor_face_3D.f90 find_neighbor_edge_3D.f90 find_neighbor_corner_3D.f90 \
    refinement_execute_3D.f90 get_block_spacing_origin.f90 update_neighbors.f90 \
	find_sisters.f90 max_active_level.f90 min_active_level.f90 get_free_local_light_id.f90 gather_blocks_on_proc.f90 \
	merge_blocks.f90 create_active_and_sorted_lists.f90 quicksort.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_unit_test.o: module_unit_test.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_initialization.o $(OBJDIR)/module_mesh.o $(OBJDIR)/module_time_step.o \
	unit_test_ghost_nodes_synchronization.f90 unit_test_wavelet_compression.f90 \
	unit_test_treecode.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_treelib.o: module_treelib.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_IO.o: module_IO.f90 $(OBJDIR)/module_mesh.o $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o \
	$(OBJDIR)/module_hdf5_wrapper.o $(OBJDIR)/module_mpi.o $(OBJDIR)/module_operators.o $(OBJDIR)/module_physics_metamodule.o \
	save_data.f90 write_field.f90 write_vorticity.f90 read_field.f90 \
	read_mesh.f90 check_file_exists.f90 get_attributes.f90 read_file_flusi.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_operators.o: module_operators.f90 $(OBJDIR)/module_mesh.o $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o \
	volume_integral.f90 compute_vorticity.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

wabbit-post: main_post.f90 $(MOBJS) $(OBJS)
		$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR) a.out wabbit wabbit-post
.PHONY: doc test
doc:
	doxygen doc/doc_configuration
	firefox doc/output/html/index.html &
test:
	@cd TESTING/;  ./runtests.sh

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
