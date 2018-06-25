#!/bin/bash

#-------------------------------------------------------------------------------
# WABBIT (FSI) unit test
# This script helps you set up new unit tests.
#       1) You define the parameter file you want
#       2) run it with WABBIT to generate the output files
#       3) with this script, reduce all HDF5 files to four values
#       4) you can delete the *.h5 files
#-------------------------------------------------------------------------------

if  [ "$WABBIT_POST" == "" ]; then
    WABBIT_POST=./wabbit-post
fi

# loop over all HDF5 files an generate keyvalues using WABBIT
for file in *.h5
do
  ${mpi_serial} $WABBIT_POST 2D --keyvalues ${file}
  mv ${file%%.h5}.key ${file%%.h5}.ref
done

echo "delete h5 files? [Y,n]"
read answer
if [ ! "$answer" == "n" ]; then
  rm *.h5
fi
