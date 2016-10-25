# WABBIT
## (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nsects in (T)urbulence

## Getting Started

how to get a copy of WABBIT and compile the code

### clone from github

```
git clone https://github.com/mario-sroka/WABBIT
```

### run make

choose compiler with FC option (to v0.2) 

```
make FC=[gfortran|ifort]
```

choose compiler with FC option (up v0.3) 

```
make FC=[mpif90]
```
### run WABBIT

customize .ini-file and rename file to [your_filename.ini], run WABBIT

```
wabbit [your_filename.ini]
```

## History

* **v0.1** - *2D, periodic boundary, mesh adaption* --- lots of schemes for testing
* **v0.2** - *2D, periodic boundary, mesh adaption* --- change mesh adaption to *paris_meeting*-version, simplify code 
* **v0.3** - *2D, periodic boundary, mesh adaption* --- MPI added