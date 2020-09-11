#!/bin/bash
#-------------------------------------------------------------------------------
# WABBIT unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# what parameter file
dir="./TESTING/wabbit_post/pod/"
happy=0
sad=0
echo "testing wabbit_post -POD"

# list of prefixes the test generates
prefixes=(u mode1 reconst1-005)
# list of possible times (no need to actually have them)
times=( 000000000001 000000000002 000000000010 000000000012)

# run actual test
sh ${dir}/pod.sh $dir

echo "============================"
echo "run done, analyzing data now"
echo "============================"

# loop over all HDF5 files and generate TESTING/wabbit_post/pod/pod_test.shkeyvalues using wabbit
for p in ${prefixes[@]}
do
  for t in ${times[@]}
  do
    echo "--------------------------------------------------------------------"
    # *.h5 file coming out of the code
    file=${p}"_"${t}".h5"
    # will be transformed into this *.key file
    keyfile=${p}"_"${t}".key"
    # which we will compare to this *.ref file
    reffile=${dir}${p}"_"${t}".ref"

    if [ -f $file ]; then
        # get four characteristic values describing the field
        ${mpi_serial} ./wabbit-post --keyvalues ${file}
        # and compare them to the ones stored
        if [ -f $reffile ]; then
            ${mpi_serial} ./wabbit-post --compare-keys $keyfile $reffile
            result=$(cat return); rm return
            if [ $result == "0" ]; then
              echo -e " :) Happy, this looks ok!" $keyfile $reffile
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


echo -e "\t  happy tests: \t" $happy
echo -e "\t  sad tests: \t" $sad

#-------------------------------------------------------------------------------
#                               RETURN
#-------------------------------------------------------------------------------
if [ $sad == 0 ]
then
  exit 0
else
  exit 999
fi
