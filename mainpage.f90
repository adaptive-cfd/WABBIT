!> \mainpage (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nsects in (T)urbulence
!> \section intro_sec Introduction
!> Wabbit is a Fortran ....
!> \section install_sec Installation
!> \subsection step1 Clone from GitHub
!! \n
!! `git clone https://github.com/mario-sroka/WABBIT`
!> \subsection step2 Run make
!! choose compiler with FC option (to v0.2): \n
!! \n
!! `make FC=[gfortran|ifort]` \n
!! \n
!! choose compiler with FC option (from v0.3): \n
!! \n
!! `make FC=[mpif90]` \n
!> \subsection step3 Run WABBIT
!! customize .ini-file and rename file to [your_filename.ini], run WABBIT with option for dimension and .ini-file name: \n
!! \n
!! `wabbit [2D|3D] [your_filename.ini]`
