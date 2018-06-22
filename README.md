# WABBIT
## (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nteractions in (T)urbulence

## Getting Started

how to get a copy of WABBIT and compile the code

### clone from github

```
git clone https://github.com/adaptive-cfd/WABBIT
```

### run make

choose compiler with FC option (from v0.3)

```
make FC=[mpif90]
```

### run WABBIT

customize .ini-file and rename file to [your_filename.ini], run WABBIT with option for dimension and .ini-file name

```
wabbit [2D|3D] [your_filename.ini] --memory=2.0GB
```

where the --memory options allows you to approximately control how much memory is globally allocated, i.e., on all ranks. Note that WABBIT does not free memory which is once allocated again during runtime. This is because we intent to use clusters, where the globally available memory is reserved for the exectution.

## History

* **v0.1** - *2D, periodic boundary, mesh adaption* --- lots of schemes for testing
* **v0.2** - *2D, periodic boundary, mesh adaption* --- change mesh adaption to *paris_meeting*-version, simplify code
* **v0.3** - *2D, periodic boundary, mesh adaption* --- MPI added, split block data in light/heavy data
* **v0.4** - *2D, periodic boundary, mesh adaption* --- reworked data structure, order subroutines and switch to explicit variable/parameter passing
* **v0.5** - *3D, 2D, periodic boundary, mesh adaption* --- now uses 3D
