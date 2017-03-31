# Makefile for WABBIT code, adapted from pseudospectators/FLUSI and pseudospectators/UP2D
# Non-module Fortran files to be compiled:
FFILES = check_allocation.f90 encoding_2D.f90 int_to_binary.f90 treecode_size.f90 adjacent_block_2D.f90 array_compare.f90 \
proc_to_lgt_data_start_id.f90 lgt_id_to_hvy_id.f90 hvy_id_to_lgt_id.f90 lgt_id_to_proc_rank.f90 get_free_light_id.f90 \
RHS_2D_convection_diffusion.f90 RHS_2D_navier_stokes.f90 encoding_3D.f90 adjacent_block_3D.f90 RHS_3D_convection_diffusion.f90 \
RHS_3D_navier_stokes.f90 f_xy_2D.f90 f_xyz_3D.f90 init_random_seed.f90 error_msg.f90

# Object and module directory:
OBJDIR = OBJ
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = module_precision.f90 module_params.f90 module_debug.f90 module_ini_files_parser.f90 module_hdf5_wrapper.f90 \
	module_interpolation.f90 module_init.f90 module_mesh.f90 module_IO.f90 module_time_step.f90 module_MPI.f90 module_unit_test.f90
# physics modules
MFILED += module_convection_diffusion.f90
MFILED += module_2D_navier_stokes.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = LIB
VPATH += :LIB/MAIN:LIB/MODULE:LIB/INI:LIB/HELPER:LIB/MESH:LIB/IO:LIB/TIME:LIB/EQUATION:LIB/MPI:LIB/DEBUG

# Set the default compiler if it's not already set
FC = mpif90

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
#FFLAGS += -Wuninitialized -O -fimplicit-none -fbounds-check -g -ggdb
FFLAGS += -Wall -Wextra -Wconversion -g3 -fbacktrace -fbounds-check -ffpe-trap=zero -g -ggdb -fimplicit-none
#FFLAGS += -Wconversion -g3 -fbacktrace -fbounds-check -ffpe-trap=zero -g -ggdb -fimplicit-none
# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib64
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC)
endif

# Intel compiler
mpif90:=$(shell $(FC) --version | head -c 5)
ifeq ($(mpif90),ifort)
PPFLAG= -fpp
FFLAGS = -FR -O3 -warn all -traceback -check bounds -heap-arrays

FFLAGS += -module $(OBJDIR) # specify directory for modules.
LDFLAGS = -L/usr/X11/lib/ -lX11 #-L/usr/lib64/lapack -llapack
# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib64
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC)
endif

# Compile main programs, with dependencies.
wabbit: main.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).

# first compile precision module
$(OBJDIR)/module_precision.o: module_precision.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

#compile physics modules
$(OBJDIR)/module_convection_diffusion.o: module_convection_diffusion.f90 $(OBJDIR)/module_precision.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/module_navier_stokes.o: module_navier_stokes.f90 $(OBJDIR)/module_precision.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_params.o: module_params.f90 $(OBJDIR)/module_convection_diffusion.o $(OBJDIR)/module_navier_stokes.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_debug.o: module_debug.f90 $(OBJDIR)/module_params.o \
	check_lgt_block_synchronization.f90 write_future_mesh_lvl.f90 write_debug_times.f90 write_block_distribution.f90 write_com_list.f90 \
	write_com_matrix.f90 write_com_matrix_pos.f90 unit_test_ghost_nodes_synchronization.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_ini_files_parser.o: module_ini_files_parser.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_interpolation.o: module_interpolation.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_hdf5_wrapper.o: module_hdf5_wrapper.f90 $(OBJDIR)/module_params.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_init.o: module_init.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_ini_files_parser.o \
	init_data.f90 allocate_block_list.f90 allocate_block_data.f90 inicond_gauss_blob.f90 initial_block_distribution_2D.f90 new_block_heavy.f90 \
	allocate_work_data.f90 inicond_vorticity_filaments.f90 ini_file_to_params.f90 inicond_zeros.f90 initial_block_distribution_3D.f90 \
	inicond_richtmyer_meshkov.f90 inicond_shear_layer.f90 inicond_sinus_2D.f90 inicond_sinus_3D.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_MPI.o: module_MPI.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_interpolation.o\
	synchronize_ghosts.f90 copy_ghost_nodes_2D.f90 create_send_buffer_2D.f90 write_receive_buffer_2D.f90 synchronize_internal_nodes.f90 \
	max_com_num.f90 fill_send_buffer.f90 fill_receive_buffer.f90 RMA_lock_unlock_get_data.f90 RMA_lock_unlock_put_data.f90 \
	isend_irecv_data.f90 copy_ghost_nodes_3D.f90 create_send_buffer_3D.f90 write_receive_buffer_3D.f90 blocks_per_mpirank.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_time_step.o: module_time_step.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_MPI.o\
	time_step_RK4.f90 filter_block.f90 filter_1D.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_mesh.o: module_mesh.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_interpolation.o \
	$(OBJDIR)/module_time_step.o\
	create_lgt_active_list.f90 create_hvy_active_list.f90 update_neighbors_2D.f90 find_neighbor_edge_2D.f90 does_block_exist.f90 \
	find_neighbor_corner_2D.f90 refine_everywhere.f90 respect_min_max_treelevel.f90 refine_mesh_2D.f90 adapt_mesh.f90 threshold_block.f90 \
	ensure_gradedness.f90 ensure_completeness_2D.f90 coarse_mesh_2D.f90 balance_load_2D.f90 set_desired_num_blocks_per_rank.f90 \
	compute_friends_table.f90 compute_affinity.f90 treecode_to_sfc_id.f90 treecode_to_hilbercode.f90 update_neighbors_3D.f90 \
	find_neighbor_face_3D.f90 find_neighbor_edge_3D.f90 find_neighbor_corner_3D.f90 refine_mesh_3D.f90 ensure_completeness_3D.f90 \
	coarse_mesh_3D.f90 balance_load_3D.f90 treecode_to_3D_z_curve.f90 get_block_spacing_origin.f90 decoding.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_unit_test.o: module_unit_test.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_init.o $(OBJDIR)/module_mesh.o \
	unit_test_ghost_nodes_synchronization.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/module_IO.o: module_IO.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_debug.o $(OBJDIR)/module_hdf5_wrapper.o $(OBJDIR)/module_MPI.o\
	save_data.f90 write_field.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR)/*.o $(OBJDIR)/*.mod a.out

tidy:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod a.out

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
