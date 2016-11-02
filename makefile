# Makefile for WABBIT code, adapted from pseudospectators/FLUSI and pseudospectators/UP2D
# Non-module Fortran files to be compiled:
FFILES = init_data.f90 allocate_block_memory.f90 initial_condition_dense_field.f90 inicond_gauss_blob.f90 matrix_to_block_tree.f90 \
encoding.f90 int_to_binary.f90 initial_block_distribution.f90 new_block_light.f90 treecode_size.f90 new_block_heavy.f90 get_heavy_free_block.f90 \
save_data.f90 write_field.f90 update_neighbors.f90 adjacent_block.f90 does_block_exist.f90 find_block_id.f90 array_compare.f90 interpolate_mesh.f90 \
get_sister_id.f90 delete_block_heavy.f90 delete_block_light.f90 get_light_free_block.f90 find_neighbor_corner.f90 find_neighbor_edge.f90 \
synchronize_ghosts.f90 create_com_list.f90 find_neighborhood.f90 sort_com_list.f90 com_allowed.f90 send_receive_data.f90 calc_dt.f90 \
time_step_RK4.f90 block_count.f90 RHS_2D_convection_diffusion.f90 refine_everywhere.f90 respect_min_max_treelevel.f90  get_light_free_block.f90 \
adapt_mesh.f90 threshold_block.f90 calculate_detail.f90 ensure_gradedness.f90 ensure_completeness.f90 

# Object and module directory:
OBJDIR = OBJ
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = module_params.f90 module_blocks.f90 ini_files_parser.f90 hdf5_wrapper.f90 module_interpolation.f90  
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = LIB
VPATH += :LIB/MAIN:LIB/MODULE:LIB/INI:LIB/HELPER:LIB/MESH:LIB/IO:LIB/TIME:LIB/EQUATION:LIB/MPI

# Set the default compiler if it's not already set
FC = mpif90

#-------------------------------------------------------------------------------
# COMPILER-DEPENDEND PART
#-------------------------------------------------------------------------------
# GNU compiler
#ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
#FFLAGS += -J$(OBJDIR) # specify directory for modules.
#FFLAGS += -O3 -ffree-line-length-none
#PPFLAG= -cpp #preprocessor flag
#LDFLAGS = -llapack
# Debug flags for gfortran:
#FFLAGS += -Wuninitialized -O -fimplicit-none -fbounds-check -g -ggdb
#FFLAGS += -Wuninitialized -Wall -Wextra -Wconversion -g3 -fbacktrace -fbounds-check -ffpe-trap=zero -g -ggdb -fimplicit-none
#FFLAGS += -Wconversion -g3 -fbacktrace -fbounds-check -ffpe-trap=zero -g -ggdb -fimplicit-none
# HDF_ROOT is set in environment.
#HDF_LIB = $(HDF_ROOT)/lib
#HDF_INC = $(HDF_ROOT)/include
#LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
#FFLAGS += -I$(HDF_INC)
#endif

# Intel compiler
#ifort:=$(shell $(FC) --version | head -c 5)
#ifeq ($(ifort),ifort)
#PPFLAG= -fpp #preprocessor flag
#FFLAGS = -FR -O3 -warn all -traceback -check bounds -heap-arrays
#FFLAGS = -FR -O3 -warn all -traceback -check bounds
#FFLAGS += -module $(OBJDIR) # specify directory for modules.
#LDFLAGS = -L/usr/X11/lib/ -lX11 -L/usr/lib64/lapack -llapack
# HDF_ROOT is set in environment.
#HDF_LIB = $(HDF_ROOT)/lib
#HDF_INC = $(HDF_ROOT)/include
#LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
#FFLAGS += -I$(HDF_INC)
#endif

#IBM compiler
#ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
#FFLAGS += -qmoddir=$(OBJDIR)
#FFLAGS += -I$(OBJDIR)
#PPFLAG=-qsuffix=cpp=f90  #preprocessor flag
#endif

# MPI compiler
mpif90:=$(shell $(FC) --version | head -c 5)
ifeq ($(mpif90),ifort)
PPFLAG= -fpp
FFLAGS = -FR -O3 -warn all -traceback -check bounds -heap-arrays

FFLAGS += -module $(OBJDIR) # specify directory for modules.
LDFLAGS = -L/usr/X11/lib/ -lX11 -L/usr/lib64/lapack -llapack
# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC)
endif

PROGRAMS = wabbit

# Both programs are compiled by default.
all: directories $(PROGRAMS)

# Compile main programs, with dependencies.
wabbit: main.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).
$(OBJDIR)/module_blocks.o: module_blocks.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/hdf5_wrapper.o: hdf5_wrapper.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser.o: ini_files_parser.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/module_params.o: module_params.f90 $(OBJDIR)/module_blocks.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/module_interpolation.o: module_interpolation.f90 $(OBJDIR)/module_params.o $(OBJDIR)/module_blocks.o
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
