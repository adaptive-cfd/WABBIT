########################################################################
########################################################################
#	Makefile for Python LIB 
########################################################################
########################################################################
#=======================================================================
#                   define the compiler names
#=======================================================================

CC       = mpicc
F90      = mpif90
FPP      = gfortran -E -cpp
PYDIR 	 = LIB/PYTHON
PYTHON   = python
moddir   = $(OBJDIR) # Use WABBITs OBJECT Directory
wrap_dir = $(PYDIR)/f90wrap


#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),mpif90)
	FPP      = gfortran -E
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = -fPIC
    FCOMP    = gfortran
    LIBS     =
endif

ifeq ($(F90),ifort)

	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = -fpscomp logicals -fPIC # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
    LIBS =
endif

CFLAGS = -fPIC #     ==> universal for ifort, gfortran, pgi

#=======================================================================
#=======================================================================

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# ======================================================================
# PROJECT CONFIG, do not put spaced behind the variables
# ======================================================================
# Python module name
PYTHON_MODN = WABBIT
# mapping between Fortran and C types
KIND_MAP = $(PYDIR)/kind_map

#=======================================================================
#
#=======================================================================

#VPATH	=

#=======================================================================
#       List all source files required for the project
#=======================================================================

# names (without suffix), f90 sources
LIBSRC_SOURCES = $(MFILES:%.f90=%) $(FFILES:%.f90=%)

# file names
LIBSRC_FILES = $(addsuffix .f90,${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(OBJS) $(MOBJS)

# only used when cleaning up
LIBSRC_FPP_FILES = $(LIBSRC_SOURCES:%=$(OBJDIR)/%.fpp)

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================

# names for python modules 
PYSRC  = module_precision.f90 module_globals.f90 module_params.f90 module_timing.f90 module_hdf5_wrapper.f90 \
	module_interpolation.f90 module_initialization.f90 module_mesh.f90 module_IO.f90 module_time_step.f90 module_mpi.f90 module_unit_test.f90 \
	module_treelib.f90  module_ini_files_parser.f90  module_ini_files_parser_mpi.f90 \
	module_indicators.f90 module_operators.f90 module_navier_stokes.f90 module_ns_penalization.f90 \
	module_physics_metamodule.f90 module_ACM.f90 module_ConvDiff_new.f90 module_bridge_interface.f90 \
	module_bridge.f90 module_navier_stokes_params.f90 module_helpers.f90 \
	module_insects.f90 module_boundary_conditions.f90 module_funnel.f90 module_navier_stokes_cases.f90\
	$(FFILES)

LIBSRC_WRAP_SOURCES = $(PYSRC:%.f90=%) $(FILES:%.f90=%)

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix .o,${LIBSRC_WRAP_SOURCES})

# fpp files
LIBSRC_WRAP_FPP_FILES =$(LIBSRC_WRAP_SOURCES:%=$(OBJDIR)/%.fpp)

#=======================================================================
#                 Relevant suffixes
#=======================================================================

.SUFFIXES: .f90 .fpp

#=======================================================================
#
#=======================================================================

.PHONY: all clean

.c.o:
	${CC} ${CFLAGS} -c $< -o $@


$(OBJDIR)/%.fpp: %.f90 
	${FPP} ${FPP_F90FLAGS} $<  -o $@


libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?


_${PYTHON_MODN}.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v -P
	mv f90wrap* $(OBJDIR)
	f2py-f90wrap --f90exec=$(F90) --fcompiler=$(FCOMP) --build-dir . -c -m _${PYTHON_MODN} -L. -lsrc $(OBJDIR)/f90wrap*.f90 \
	-I$(OBJDIR) $(HDF5_FLAGS) -L$(HDF_LIB) -L$(HDF_LIB)64 -lhdf5_fortran -lhdf5 \
	-lz -ldl -lm -lblas -llapack -I$(HDF_INC) $(SB_LIB) $(SB_INCL) 


clean-python:
	rm -rf WABBIT/
	rm -rf __pycache__
	rm -rf $(wrap_dir)/* $(wrap_dir)/
	rm -rf *WABBIT*.so src.linux*
	rm -rf _${PYTHON_MODN}.so _${PYTHON_MODN}_pkg.so libsrc.a
	
python: directories _${PYTHON_MODN}.so 
	


########################################################################
# Remarks:
# 1. Libraries must be listed after the objects that use them
# (more precisely, a library will be used only if it contains a 
# symbol that satisfies an undefined reference known at the time
# it is encountered)
# Therefore make sure that -L. -lsrc comes after f90wrap*.f90.
# 2. f90wrap: allocatable subroutine arguments are not currently supported
# To get arround with it read: https://github.com/jameskermode/f90wrap/issues/18
##################################################################
