#----------------------------------------

OBJ = $(shell find ./LIB -name "*.f90" -printf %f\ | sed 's/\.f90/.o/g')
SRC = $(shell find ./LIB -name "*.f90")

LDFLAGS = -L/usr/X11/lib/ -lX11 -L/usr/lib64/lapack -llapack
CFLAGS = -FR -warn all -O3 -traceback -check bounds -check all #-heap-arrays

CASE = USER/msr/test004/
FILE = init_data.f90
OBJDIR = OBJ

OBJNEW = $(addprefix $(OBJDIR)/, $OBJ \)

#----------------------------------------
all: main

main: $(OBJ) init_data.o move
	ifort $(LDFLAGS) -o main $(OBJDIR)/*.o

$(OBJ): $(SRC)
	ifort $(CFLAGS) -c $(SRC)

init_data.o: $(CASE)$(FILE) 
	ifort $(CFLAGS) -c $(CASE)$(FILE)
	
$(OBJDIR):
	mkdir $(OBJDIR)
	
move:
	cp -p *.o *.mod *.f90 $(OBJDIR)/
	rm *.o *.mod *.f90

clean:
	rm -r $(OBJDIR) 

$(shell   mkdir -p $(OBJDIR))