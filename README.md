# WABBIT
## (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nsects in (T)urbulence

## Getting Started

how to get a copy of WABBIT and compile the code

### clone from github

```
git clone https://github.com/mario-sroka/WABBIT
```

### run make

choose the compiler with FC option 

```
make FC=[gfortran|ifort]
```
### run WABBIT

customize .ini-file and rename file to [your_filename.ini], run WABBIT

```
wabbit [your_filename.ini]
```

## History

* **v0.1** - *2D, periodic boundary, mesh adaption* --- lots of schemes for testing
* **v0.2** - *2D, periodic boundary, mesh adaption* --- change mesh adaption to *paris_meeting*-version, simplify code 