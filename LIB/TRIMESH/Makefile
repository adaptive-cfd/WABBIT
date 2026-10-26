
FC = mpifort
FFLAGS = -Wall -Wextra -g -fcheck=all

# Adjust these paths to your project
GLOBALS_SRC = /home/tejada/WABBIT/WABBIT/LIB/MODULE/module_globals.f90
TRIMESH_SRC = module_trimesh.f90
MAIN_SRC = main_trimesh.f90

TARGET = test_trimesh

.PHONY: all run clean

all: $(TARGET)

$(TARGET): module_globals.o module_trimesh.o main_trimesh.o
	$(FC) $(FFLAGS) -o $@ $^

module_globals.o: $(GLOBALS_SRC)
	$(FC) $(FFLAGS) -c $< -o $@

module_stl_file_reader.o: $(STL_SRC) module_globals.o
	$(FC) $(FFLAGS) -c $< -o $@

module_trimesh.o: $(TRIMESH_SRC) module_globals.o \
                  module_stl_file_reader.o
	$(FC) $(FFLAGS) -c $< -o $@

main_trimesh.o: $(MAIN_SRC) module_trimesh.o
	$(FC) $(FFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET) mesh.dat

clean:
	rm -f $(TARGET) *.o *.mod
