#!/bin/bash

wavelet="$1"
dir="$2"

happy=0
sad=0
file1="vor_00100.h5"
file2="vor_00200.h5"

echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "This test takes an existing 2D field (which is not equidistant!!) and "
echo "(1) refines it everywhere by one level"
echo "(2) coarsens the result down one level (to the original grid)"
echo "We compare the HDF5 files of both steps with the reference data."
echo "NOTE that on non-equidistant grids, the result of (2) is _NOT_ perfectly identical to the original data."
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""

# run actual test
${mpi_command} ./wabbit-post --sparse-to-dense ./TESTING/wavelets/vor_000020000000.h5 ${file1} --wavelet=${wavelet} --operator=refine-everywhere

${mpi_command} ./wabbit-post --sparse-to-dense ${file1} ${file2} --wavelet=${wavelet} --operator=coarsen-everywhere



wabbit-compare-hdffiles.py ${dir}${file1} ${file1}
result="$?"

if [ $result == "0" ]; then
    echo -e " :) Happy, this looks ok!" $keyfile $reffile
    happy=$((happy+1))
else
    echo -e ":[ Sad, this test failed!" $keyfile $reffile
    sad=$((sad+1))
fi

wabbit-compare-hdffiles.py ${dir}${file2} ${file2}
result="$?"

if [ $result == "0" ]; then
    echo -e " :) Happy, this looks ok!" $keyfile $reffile
    happy=$((happy+1))
else
    echo -e ":[ Sad, this test failed!" $keyfile $reffile
    sad=$((sad+1))
fi


echo -e "\t  happy tests: \t" $happy
echo -e "\t  sad tests: \t" $sad

if [ $sad == 0 ]; then
    exit 0
else
    exit 999
fi
