#!/bin/bash

test=$1
mpi=$2

curdir=$PWD
dir=$(dirname $test)
inifile=$(basename $test)

echo "We will update the unit test "
echo $test "(NOTE THIS SHOULD BE AN INI FILE)"
echo "Dir=" $dir

if [ ! -d "$dir" ]; then
    echo "dir not found.."
    exit 67
fi

echo "inifile=" $inifile

if [ ! -f "$dir"/"$inifile" ]; then
    echo "inifile not found.."
    exit 67
fi

echo "mpi_command=" $mpi
echo "with the current version of the code"
echo "be very very sure what you are doing!! you need to know EXACTLY why the test fails"
echo "and that this is either a bugfix or a new feature."
echo " Do you understand? (yes,no)"
read understood

if [ ! "$understood" == "yes" ];
then
    exit 9
fi


echo $curdir
cd $dir
ln -s $curdir/../wabbit
ln -s $curdir/../wabbit-post
ln -s $curdir/create_ref_files.sh

$mpi ./wabbit 2D $inifile --memory=0.25GB

./create_ref_files.sh

rm -f *times.dat

rm wabbit wabbit-post create_ref_files.sh

