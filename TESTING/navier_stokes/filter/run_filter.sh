#!/bin/bash
#-------------------------------------------------------------------------------
# This routine runs all filters on sods shock tube in the list of available filters
# and produces data (statevector) at time t=0.05 for comparing the filters
# author: P.Krah date: 20.04.2018
#-------------------------------------------------------------------------------

test_dir="./TESTING/"
dir="./TESTING/navier_stokes/filter/"
params=$1

# check if file exists
if [ ! -e $params ]
then
	echo  file $params you passed in argument one does not exists
	echo "call: run_filter params.ini 00000000010000 filter_name"
	exit
fi

# write time
times=$2
# list of prefixes the test generates
prefixes=(Ux p rho)
# list of filters
filter=$3
#filters=(explicit_5pt bogey_shock)


	echo " "
	echo " "
	echo "===================================================================="
	echo -e "running filter:" ${filter}
	echo "===================================================================="
	echo " "
	echo " "
	echo " "

	./${test_dir}/replace_ini_value.sh  ${params} Discretization filter_type ${filter}
	# run actual test
	${mpi_command} ./wabbit 2D ${params} --memory=2GB ${ghosts}

	echo " "
	echo " "
	echo " "

	# loop over all HDF5 files and generate keyvalues using wabbit
	for p in ${prefixes[@]}
	do
	    #current configname
	    configname=${p}"_"${times}
	    #change configname for current filter
	    mv ${configname}".h5" ${filter}_${configname}".h5"
	    #
	    echo "changed: " ${configname}".h5 to" ${filter}_${configname}".h5"
	done
	echo " "
	echo " "
	echo " "
	echo " "
exit 0
