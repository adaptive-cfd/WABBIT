#!/bin/bash
#-------------------------------------------------------------------------------
# WABBIT unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# what parameter file
test_dir="./TESTING/"
dir="./TESTING/navier_stokes/filter/"
params=${dir}"filter_test.ini"
happy=0
sad=0
echo "testing filters"

# list of prefixes the test generates
prefixes=(Ux p rho)
# list of filters
filters=(no_filter explicit_5pt explicit_7pt explicit_9pt explicit_11pt wavelet bogey_shock)
## change this list in run_all_filter.sh to!!
# list of possible times (no need to actually have them)
times=(000000050000)




${dir}/run_all_filter.sh ${params}


	echo "============================"
	echo "run done, analyzing data now"
	echo "============================"

for filter in ${filters[@]}
do
	# loop over all HDF5 files and generate keyvalues using wabbit
	for p in ${prefixes[@]}
	do
	  for t in ${times[@]}
	  do
	    echo "--------------------------------------------------------------------"
	    #configname
	    configname=${filter}_${p}"_"${t}
	    # *.h5 file coming out of the code
	   	file=${configname}".h5"
	    # will be transformed into this *.key file
	    keyfile=${configname}".key"
	    # which we will compare to this *.ref file
	    reffile=${dir}/${configname}".ref"
	    if [ -f $file ]; then
	        # get four characteristic values describing the field
	        ${mpi_command} ./wabbit-post 2D --keyvalues ${file}
	        # and compare them to the ones stored
	        if [ -f $reffile ]; then
	        	${mpi_command} ./wabbit-post 2D --compare-keys $keyfile $reffile
	            result=$(cat return); rm return
	            if [ $result == "0" ]; then
	              echo -e ":) Happy, this looks okay!" $keyfile $reffile
	              happy=$((happy+1))
	            else
	              echo -e ":[ Sad, this is failed!" $keyfile $reffile
	              sad=$((sad+1))
	            fi
	        else
	            sad=$((sad+1))
	            echo -e ":[ Sad: Reference file not found"
	        fi
	    else
	        sad=$((sad+1))
	        echo -e ":[ Sad: output file not found"
	    fi
	    echo " "
	    echo " "

	  done
	done
done

echo -e "\t happy tests:\t" $happy
echo -e "\t sad tests: \t" $sad

#-------------------------------------------------------------------------------
#                               RETURN
#-------------------------------------------------------------------------------
if [ $sad == 0 ]
then
  exit 0
else
  exit 999
fi
