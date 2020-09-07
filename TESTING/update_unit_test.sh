#!/bin/bash

#-------------------------------------------------------------------------------
# Usage:
# ./update_unit_test.sh PATH/TO/UNITEST/test.ini "mpirun -n 4"
#-------------------------------------------------------------------------------

test_inifile=$1
mpi=$2


Color_Off='\e[0m'       # Text Reset
Black='\e[0;30m'        # Black
Red='\e[0;31m'          # Red
Green='\e[0;32m'        # Green
Yellow='\e[0;33m'       # Yellow
Blue='\e[0;34m'         # Blue
Purple='\e[0;35m'       # Purple
Cyan='\e[0;36m'         # Cyan
White='\e[0;37m'        # White

curdir=$PWD
dir=$(dirname $test_inifile)
inifile=$(basename $test_inifile)

echo "We will update the unit test_inifile "
echo -e ${Cyan}${test_inifile}${Color_Off} "(NOTE THIS SHOULD BE AN INI FILE)"
echo -e "Dir=" ${Cyan}${dir}${Color_Off}

if [ ! -d "$dir" ]; then
    echo -e $Red"dir=$dir not found.."$Color_Off
    exit 2
fi

echo -e "inifile=" ${Cyan}${inifile}${Color_Off}

if [ ! -f "$dir"/"$inifile" ]; then
    echo "inifile not found.."
    exit 3
fi

echo -e "mpi_command=" ${Cyan}${mpi}${Color_Off}
echo "with the current version of the code"
echo "be very very sure what you are doing!! you need to know EXACTLY why the test fails"
echo "and that this is either a bugfix or a new feature."
echo " Do you understand? (yes,no)"
read understood

if [ ! "$understood" == "yes" ];
then
    exit 4
fi


echo -e "curren directory" ${Cyan}${curdir}${Color_Off}
cd $dir

if [ ! -f $curdir/../wabbit ]; then
    echo -e $Red"Wabbit executable in $curdir/../wabbit not found.."$Color_Off
    exit 5
fi

if [ ! -f $curdir/create_ref_files.sh ]; then
    echo -e $Red"$curdir/create_ref_files.sh not found.."$Color_Off
    exit 6
fi

ln -s $curdir/../wabbit
ln -s $curdir/../wabbit-post
ln -s $curdir/create_ref_files.sh

if [[ $inifile == "pod.sh" ]]; then
    sh $inifile
else
    $mpi ./wabbit $inifile ${memory}
fi

./create_ref_files.sh

rm -f *times.dat

rm wabbit wabbit-post create_ref_files.sh
