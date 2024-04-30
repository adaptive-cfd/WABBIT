# WABBIT v2.0beta5 (newBiorthogonal)
## (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nteractions with (T)urbulence

> :exclamation: New in 05/2021: please see this video for an introduction to the code's datastructures: https://www.youtube.com/watch?v=qBBIW2-ktgo

With WABBIT it is possible to solve partial differential equations (PDE) on block-based adaptive grids. Simulations in 2D and 3D are possible and performed fully parallel on CPU. The set of PDE is encapsulated from the code handling the adaptive grid, and thus existing monobloc solvers can be adapted easily for this solver. WABBIT can handle PDEs of the following type:

$\partial_t \phi = N\left(\phi\right)$

and $N\left(\phi\right)$ can be defined. This implementation is handled by the "physics-modules". Note the current version of the code does not handle elliptic PDE, such as the Poisson equation that typically arises in incompressible fluid dynamics. Instead, we use a quasi-hyperbolic approximation in that case, the "artificial compressibility method".

## Installation of WABBIT
How to get a copy of WABBIT and compile the code:
```
git clone https://github.com/adaptive-cfd/WABBIT.git
```
Unpack the file and run the compilation and tests with `make`, make sure that all necessary dependencies are loaded:
```
make all
make test
```
> :warning: since 15 Aug 2023, the unit testing framework has evolved. It now stores full HDF5 files in the TESTING directory, which makes it easier to visualize the reference data and current results, should they be different. We now calculate the L2 error of the field, if the grid is identical. This new framework requires the [WABBIT Python Tools](https://github.com/adaptive-cfd/python-tools) repository for comparing two WABBIT HDF5 files.

### Dependencies
WABBIT needs several packages installed in order to compile and run, the main dependencies are:
- [MPI](https://www.open-mpi.org/ "OpenMPI")
- [HDF5](https://github.com/HDFGroup/hdf5/tags "HDF5")
- [BLAS](https://www.netlib.org/blas/ "BLAS") + [LAPACK](https://www.netlib.org/lapack/ "LAPACK")
- Python-scripts found in the "[python-tools](https://github.com/adaptive-cfd/python-tools)" repository - these scripts are used when we perform WABBIT unit tests

Further information on the installation and compilation of all pre-requesites can be found in the wiki under [Install-WABBIT-with-requirements](../../wiki/Install-WABBIT-with-requirements).

A list of all environment variables to be set can be find in the wiki under [Loading-prerequesites](../../Loading-prerequesites). Ensure that the WABBIT-specific variables for HDF5 are set in order for the compilation to finish successfully.

## Running WABBIT
Customize the template `.ini`-file and rename file to `[your_filename.ini]`, run WABBIT and pass it the `.ini`-file as well as the total amount of memory used:
```
mpirun -n 1 ./wabbit [your_filename.ini] --memory=2.0GB
```
alternatively, you can specify the amount of memory per core that the code may use:
```
mpirun -n 1 ./wabbit [your_filename.ini] --mem-per-core=2.0GB
```
where the --memory options allows you to control how much memory is globally allocated, i.e., on all CPUs. Note that WABBIT does not free memory during runtime. This is because the code is intented to run on clusters/supercomputers, where the available memory is reserved for the execution of the code alone. This is quite typical for supercomputing.

## Additional Information
In case you have problems with the preparation to use WABBIT, have a look if you can find anything in the  [wiki](../../ "additional information for WABBIT")

More information can also be found in the documentation build with [`doxygen`](https://www.doxygen.nl/). Therefore it is necessary to have [Doxygen](http://www.stack.nl/~dimitri/doxygen/ "Doxygen") installed. Create the documentation with
``` shell
doxygen doc/doc_configuration
```
and display `doc/output/html/index.html` with your browser. You can also locally display all files with
``` shell
python3 -m http.server --directory doc/output/html
```
and then open in a browser `localhost:8000`

## Publications

- ["Wavelet adaptive proper orthogonal decomposition for large-scale flow data"](https://link.springer.com/article/10.1007/s10444-021-09922-2 "Krah2022"); Krah, Engels, Schneider, Reiss; Advances in Comput. Math. Volume 48, Article number: 10, 2022

- ["A wavelet-adaptive method for multiscale simulation of turbulent flows in flying insects"](https://arxiv.org/abs/1912.05371 "Engels2021"); T. Engels, K. Schneider, J. Reiss and M. Farge; Commun. Comput. Phys., 30, 1118-1149, 2021

- ["An Open and Parallel Multiresolution Framework Using Block-Based Adaptive Grids"](https://link.springer.com/chapter/10.1007%2F978-3-319-98177-2_19 "Sroka2018"); Sroka, Engels, Krah, Mutzel, Schneider, Reiss; Active Flow and Combustion Control 2018
