########################################################################
########################################################################
#	Makefile for Python LIB 
########################################################################
########################################################################
#=======================================================================
#                   define the compiler names
#=======================================================================

CC       = mpicc
F90      = $(FC) 

PYDIR 	 = LIB/PYTHON
moddir   = $(OBJDIR) # Use WABBITs OBJECT Directory
wrap_dir = $(PYDIR)/f90wrap


#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),mpif90)
	FCOMP    = gnu
endif


F90FLAGS = $(FFLAGS)
FPP      = gfortran -E
SOFLAGS  = -fPIC       #     ==> universal for ifort, gfortran, pgi

#=======================================================================
#=======================================================================

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
	LIBTOOL = libtool -static -o
else
	LIBTOOL = ar src
endif

#=======================================================================
#List all source files
#=======================================================================

# Files that create modules:
LIBSRC_FILES = module_precision.f90 module_globals.f90 module_params.f90 module_timing.f90 module_hdf5_wrapper.f90
	#module_interpolation.f90 module_initialization.f90 module_mesh.f90 module_IO.f90 module_time_step.f90 module_mpi.f90 module_unit_test.f90 \
	#module_treelib.f90  module_ini_files_parser.f90  module_ini_files_parser_mpi.f90 \
	#module_indicators.f90 module_operators.f90 module_navier_stokes.f90 module_ns_penalization.f90 \
	#module_physics_metamodule.f90 module_ACM.f90 module_ConvDiff_new.f90 module_bridge_interface.f90 \
	#module_bridge.f90 module_navier_stokes_params.f90 module_helpers.f90 \
	#module_insects.f90 module_boundary_conditions.f90 module_funnel.f90 module_navier_stokes_cases.f90\
	#module_simple_geometry.f90 module_shock.f90 module_pipe_flow.f90 

LIBSRC_OBJECTS := $(LIBSRC_FILES:%.f90=$(OBJDIR)/%.o)


kmap_file        = $(PYDIR)/kind_map

#=======================================================================
#Relevant suffixes
#=======================================================================

.SUFFIXES: .F90 

#=======================================================================
#   Makefile default
#======================================================================

.PHONY: python-clean

python-clean:
	rm -rf OBJ/*.o
	rm -rf  $(wrap_dir)/f90wrap_*.o $(wrap_dir)/f90wrap_*.f90
	rm -rf *wabbit*.so src.linux*
	rm  libsrc.a wabbit.py 
	
.F90.o:
	${F90} ${F90FLAGS} ${SOFLAGS} -fPIC -c $< -o $@
	
.c.o:
	${CC} ${SOFLAGS} ${F90FLAGS} -c $< -o $@
	
.F90.fpp:
	${FPP} ${SOFLAGS} $<  -o 

libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?

#=======================================================================
#shared object (direct compilation, single .py module)
#=======================================================================
python: ${LIBSRC_FILES}
	@echo -e "\n\n \e[1m creating library\e[0m \n\n"
	make libsrc.a
	mkdir -p $(wrap_dir)
	mkdir -p $(moddir)
	@echo -e "\n \n \e[1m create f90wrapper\e[0m \n\n"
	f90wrap -m wabbit $^ -k $(kmap_file) -c -v
	mv f90wrap_*.f90 $(wrap_dir)
	@echo -e "\n\n \e[1m compiling f90 wrapper \e[0m \n\n"
	f2py -c --f90exec=$(F90) -I$(moddir) \
	$(HDF5_FLAGS) -L$(HDF_LIB) -L$(HDF_LIB)64 -lhdf5_fortran -lhdf5 \
	-lz -ldl -lm -lblas -llapack -I$(HDF_INC) $(SB_LIB) $(SB_INCL) --build-dir . \
	-m _wabbit -L. -lsrc $(wrap_dir)/f90wrap*.f90 
########################################################################
########################################################################
