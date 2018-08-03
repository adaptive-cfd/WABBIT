# WABBIT
## (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nteractions with (T)urbulence

With WABBIT it is possible to solve partial differential equations on block-based adaptive grids. Calculations in 2D and 3D are possible and is performed fully parallel. As in WABBIT the set of PDE is encapsulated from the rest of the code the PDE implementation is similar to calculatoins with single domain code. Solvable PDEs are of type

<a href="http://www.codecogs.com/eqnedit.php?latex=\partial_t&space;\phi&space;=&space;N(\phi)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\partial_t&space;\phi&space;=&space;N(\phi)" title="\partial_t \phi = N(\phi)" /></a>

and <a href="http://www.codecogs.com/eqnedit.php?latex=N(\phi)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?N(\phi)" title="N(\phi)" /></a> can be defined.

## Getting Started
How to get a copy of WABBIT and compiling the code:

1. Clone from github

```
git clone https://github.com/adaptive-cfd/WABBIT
```

2. Compile the code running make.  
  Note that the mpi library and the HDF5 library have to be installed first and the variables have to be exported (see below).

```
make [FC=[mpif90]]
```

3. Run the testfiles with

```
make test
```

## HDF5 Library

Make sure that the mpi library which is also used for WABBIT is installed (for example by loading mpich3).

This is a short example (/working practice) of how to install hdf5 libary (Tested for version hdf5-1.10.1).

1. download source code from [hdf5](https://www.hdfgroup.org/downloads/hdf5/source-code/ "HDF5 Source Code")

2. open terminal and follow

```
gunzip < hdf5-X.Y.Z.tar.gz | tar xf -
cd hdf5-X.Y.Z
./configure --prefix=path_2_build_dir --enable-fortran --enable-parallel
make
make check                # run test suite.
make install
make check-install        # verify installation.
```

3. export variables:

```
HDF_ROOT=path_2_build_dir
export HDF_ROOT
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF_ROOT}/lib:${HDF_ROOT}/lib64
export LD_RUN_PATH=$LD_LIBRARY_PATH
```

### run WABBIT

Customize the .ini-file and rename file to [your_filename.ini], run WABBIT with option for dimension and .ini-file name

```
wabbit [2D|3D] [your_filename.ini] --memory=2.0GB
```

where the --memory options allows you to approximately control how much memory is globally allocated, i.e., on all ranks. Note that WABBIT does not free memory which is once allocated again during runtime. This is because we intent to use clusters, where the globally available memory is reserved for the exectution.

## Additional Information
For further Information see the documentation. Therefor it is necessary to have [Doxygen](http://www.stack.nl/~dimitri/doxygen/ "Doxygen") installed.

```
make doc
```

## History

* **v0.1** - *2D, periodic boundary, mesh adaption* --- lots of schemes for testing
* **v0.2** - *2D, periodic boundary, mesh adaption* --- change mesh adaption to *paris_meeting*-version, simplify code
* **v0.3** - *2D, periodic boundary, mesh adaption* --- MPI added, split block data in light/heavy data
* **v0.4** - *2D, periodic boundary, mesh adaption* --- reworked data structure, order subroutines and switch to explicit variable/parameter passing
* **v0.5** - *3D, 2D, periodic boundary, mesh adaption* --- now uses 3D
