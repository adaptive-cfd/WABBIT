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
# Modules exluded from python library
excludes_mod = module_insects_integration_flusi_wabbit.f90 \
	   module_sparse_operators.f90
excludes = $(addprefix f90wrap_,$(excludes_mod))
# Files that create modules:
LIBSRC_FILES = $(MFILES)
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
	rm  libwabbit.a wabbit.py 
	
.F90.o:
	${F90} ${F90FLAGS} ${SOFLAGS} -fPIC -c $< -o $@
	
.c.o:
	${CC} ${SOFLAGS} ${F90FLAGS} -c $< -o $@
	
.F90.fpp:
	${FPP} ${SOFLAGS} $<  -o 

libwabbit.a: ${LIBSRC_OBJECTS} ${MOBJS} $(OBJS)
	${LIBTOOL} $@ $?

#=======================================================================
#shared object (direct compilation, single .py module)
#=======================================================================
python: ${LIBSRC_FILES}
	@echo -e "\n\n \e[1m creating library\e[0m \n\n"
	make libwabbit.a
	mkdir -p $(wrap_dir)
	mkdir -p $(moddir)
	@echo -e "\n \n \e[1m create f90wrapper\e[0m \n\n"
	f90wrap -m libwabbit $^ -k $(kmap_file) -c -v
	rm $(excludes)
	mv f90wrap_*.f90 $(wrap_dir)
	@echo -e "\n\n \e[1m compiling f90 wrapper \e[0m \n\n"
	f2py-f90wrap -c --f90exec=$(F90) -I$(moddir) $(wrap_dir)/f90wrap*.f90 \
	$(HDF5_FLAGS) -L$(HDF_LIB) -L$(HDF_LIB)64 -lhdf5_fortran -lhdf5 \
	-lz -ldl -lm -lblas -llapack -I$(HDF_INC) $(SB_LIB) $(SB_INCL) --build-dir . \
	-m _libwabbit -L. -lwabbit 
########################################################################
# Remarks:
# 1. Libraries must be listed after the objects that use them
# (more precisely, a library will be used only if it contains a 
# symbol that satisfies an undefined reference known at the time
# it is encountered)
# Therefore make sure that -L. -lsrc comes after f90wrap*.f90.
##################################################################
