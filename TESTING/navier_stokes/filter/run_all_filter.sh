#!/bin/bash
#-------------------------------------------------------------------------------
# This routine runs all filters on sods shock tube in the list of available filters 
# and produces data (statevector) at time t=0.05 for comparing the filters
# author: P.Krah date: 20.04.2018
#-------------------------------------------------------------------------------

test_dir="./TESTING/"
dir="./TESTING/navier_stokes/filter/"
params=$1
# write time
times=$2
# list of prefixes the test generates
prefixes=(Ux p rho)
# list of filters
filters=(explicit_5pt explicit_7pt explicit_9pt explicit_11pt wavelet bogey_shock)
#filters=(explicit_5pt )


for filter in ${filters[@]}
	do

	echo " "
	echo " "
	echo "===================================================================="
	echo -e "running filter:" ${filter} 
	echo "===================================================================="
	echo " "
	echo " "
	echo " "

	./${test_dir}/replace_ini_value.sh  ${params} filter_type ${filter}
	# run actual test
	${mpi_command} ./wabbit 2D ${params} --memory=2GB

	echo " "
	echo " "
	echo " "

	# loop over all HDF5 files and generate keyvalues using wabbit
	for p in ${prefixes[@]}
	do
	    #current configname
	    configname=${p}"_"${times}
	    #change configname for current filter
	    cp ${configname}".h5" ${filter}_${configname}".h5"
	    #
	    echo "changed: " ${configname}".h5 to" ${filter}_${configname}".h5"
	done
	echo " "
	echo " "
	echo " "
	echo " "
done

exit 0