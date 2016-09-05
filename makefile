# Makefile for WABBIT code, adapted from pseudospectators/FLUSI and pseudospectators/UP2D
# Non-module Fortran files to be compiled:
FFILES = save_data.f90 write_field.f90 giveCertainOrder.f90 \
Dper.f90 D26p.f90 Dnonper.f90 D18j.f90  time_step.f90 calc_dt.f90 local_refinement_status.f90 \
synchronize_ghosts.f90 delete_block.f90 get_sister_id.f90 matrix_to_block_tree.f90 active_blocks_list.f90 \
new_block.f90 does_block_exist.f90 ensure_completeness.f90 adjacent_block.f90 adapt_mesh.f90 blocks_sum.f90 \
encoding.f90 respect_min_max_treelevel.f90 interpolate_mesh.f90 treecode_size.f90 \
find_block_id.f90 ensure_gradedness.f90 block_check.f90 get_free_block.f90 update_neighbors.f90 \
matrix_mult.f90 int_to_binary.f90 factorial.f90 \
print_data.f90 array_compare.f90 fliplr.f90 grad_test.f90 matrix_sum.f90 \
neighbor_search.f90 RHS_2D_block.f90 allocate_block_memory.f90 inicond_dense_field_wrapper.f90 \
inicond_sinus.f90 inicond_gauss_blob.f90 init_empty_file.f90

FFILES += init_data.f90

# Object and module directory:
OBJDIR = OBJ
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = module_params.f90 module_blocks.f90 module_interpolation.f90 hdf5_wrapper.f90 \
ini_files_parser.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = LIB
VPATH += :LIB/DERIVATIVES:LIB/EQUATION:LIB/HELPER:LIB/IO:LIB/MAIN:LIB/MESH:LIB/MODULE:LIB/TIME
VPATH += :LIB/INI

# Set the default compiler if it's not already set
ifndef FC
FC = gfortran
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
LDFLAGS = -llapack
# Debug flags for gfortran:
#FFLAGS += -Wuninitialized -O -fimplicit-none -fbounds-check -g -ggdb
#FFLAGS += -Wuninitialized -Wall -Wextra -Wconversion -g3 -fbacktrace -fbounds-check -ffpe-trap=zero -g -ggdb -fimplicit-none
FFLAGS += -Wconversion -g3 -fbacktrace -fbounds-check -ffpe-trap=zero -g -ggdb -fimplicit-none
# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC)
endif

# Intel compiler
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
PPFLAG= -fpp #preprocessor flag
FFLAGS = -FR -O3 -warn all -traceback -check bounds -check all #-heap-arrays
FFLAGS = -FR -O3 -warn all -traceback -check bounds
FFLAGS += -module $(OBJDIR) # specify directory for modules.
LDFLAGS = -L/usr/X11/lib/ -lX11 -L/usr/lib64/lapack -llapack
# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC)
endif

#IBM compiler
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
FFLAGS += -qmoddir=$(OBJDIR)
FFLAGS += -I$(OBJDIR)
PPFLAG=-qsuffix=cpp=f90  #preprocessor flag
endif

PROGRAMS = main






# Both programs are compiled by default.
all: directories $(PROGRAMS)

# Compile main programs, with dependencies.
main: main.f90 $(MOBJS) $(OBJS)
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
$(OBJDIR)/module_interpolation.o: module_interpolation.f90 $(OBJDIR)/module_params.o prediction_2D.f90 restriction_2D.f90
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
