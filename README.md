# WABBIT
## (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nteractions with (T)urbulence

New in 05/2021: please see this video for an introduction to the code's datastructures: https://www.youtube.com/watch?v=qBBIW2-ktgo


With WABBIT it is possible to solve partial differential equations (PDE) on block-based adaptive grids. Calculations in 2D and 3D are possible and performed fully parallel. The set of PDE is encapsulated from the code handling the adaptive grid, and thus existing monobloc solvers can be adapted easily for this solver. WABBIT can handle PDEs of the following type:

<a href="http://www.codecogs.com/eqnedit.php?latex=\partial_t&space;\phi&space;=&space;N(\phi)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\partial_t&space;\phi&space;=&space;N(\phi)" title="\partial_t \phi = N(\phi)" /></a>

and <a href="http://www.codecogs.com/eqnedit.php?latex=N(\phi)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?N(\phi)" title="N(\phi)" /></a> can be defined. This implementation is handled by the "physics-modules".

## Getting Started
How to get a copy of WABBIT and compiling the code:

1. Clone from github

```
git clone git@github.com:adaptive-cfd/WABBIT.git
```

2. Compile the code running make.

     Note that WABBIT requires:

     + [MPI library](https://www.open-mpi.org/)
     + [HDF5 library](https://www.hdfgroup.org/downloads/hdf5/source-code/ "HDF5 Source Code")
     + [BLAS+LAPACK library](http://ab-initio.mit.edu/wiki/index.php/Template:Installing_BLAS_and_LAPACK)
     + some variables have to be exported (see below).

```
make [FC=[mpif90]]
```

3. Run the unit tests with

```
make test
```

NOTE: since 15 Aug 2023, the unit testing framework has evolved. It now stores full HDF5 files in the TESTING directory, which makes it easier to visualize the reference data and current results, should they be different. We now calculate the L2 error of the field, if the grid is identical. This new framework requires the https://github.com/adaptive-cfd/python-tools repository for comparing two WABBIT HDF5 files. 


## HDF5 Library

Make sure that the mpi library which is also used for WABBIT is installed (for example by loading mpich3).

This is a short example (/working practice) of how to install hdf5 libary (Tested for version hdf5-1.10.1).

1. download source code from [hdf5](https://www.hdfgroup.org/downloads/hdf5/source-code/ "HDF5 Source Code")

2. open terminal and follow  
  (mind that *path_2_build_dir* has to be replaced by the path of the directory of your choice.)

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
   Recommendation: Add the lines to export the variables to your bashrc-file. Otherwise, the export has to be done every time you open a new terminal to compile the code.

### run WABBIT

Customize the template .ini-file and rename file to [your_filename.ini], run WABBIT with option for dimension and .ini-file name

```
mpirun -n 1 ./wabbit [your_filename.ini] --memory=2.0GB
```

where the --memory options allows you to control how much memory is globally allocated, i.e., on all CPUs. Note that WABBIT does not free memory during runtime. This is because we intent to use clusters, where the available memory is reserved for the execution of the code alone. This is quite typical for supercomputing.

## Additional Information
If you are new to WABBIT it is recommended to read the [information for newcomers](https://github.com/adaptive-cfd/WABBIT/issues?q=is%3Aissue+is%3Aopen+label%3A%22for+the+newcomers%22 "newcomer issues")!

In case you have problems with the preparation to use WABBIT, have a look at the informations given in the  [wiki](https://github.com/adaptive-cfd/WABBIT/wiki "additional information for WABBIT on fedora/ubuntu")

For further Information see the documentation. Therefore it is necessary to have [Doxygen](http://www.stack.nl/~dimitri/doxygen/ "Doxygen") installed.

```
make doc
```
## Publications

* ["Wavelet adaptive proper orthogonal decomposition for large-scale flow data"](https://link.springer.com/article/10.1007/s10444-021-09922-2 "Krah2022"); Krah, Engels, Schneider, Reiss; Advances in Comput. Math. Volume 48, Article number: 10, 2022

* ["A wavelet-adaptive method for multiscale simulation of turbulent flows in flying insects"](https://arxiv.org/abs/1912.05371 "Engels2021"); T. Engels, K. Schneider, J. Reiss and M. Farge; Commun. Comput. Phys., 30, 1118-1149, 2021

* ["An Open and Parallel Multiresolution Framework Using Block-Based Adaptive Grids"](https://link.springer.com/chapter/10.1007%2F978-3-319-98177-2_19 "Sroka2018"); Sroka, Engels, Krah, Mutzel, Schneider, Reiss; Active Flow and Combustion Control 2018
