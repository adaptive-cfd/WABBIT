#!/bin/bash

#--------------------------------
# WABBIT run unit tests
#--------------------------------
rootdir="../"
cd $rootdir
echo "   "
echo -e "\t \033[4m WABBIT: run all existing unit tests \033[0m"
echo "   "

if [ -z "$nprocs" ]; then
    echo "unset"
    nprocs=4
fi

if [ -z "$mpi_command" ]; then
    export mpi_command="nice mpiexec -n ${nprocs}"
fi

fail_color=$'\033[31;1m'
pass_color=$'\033[92;1m'
end_color=$'\033[0m'

# list all tests here. For readability, we include header lines like ---acm---
# they structure teh output on the screen. Note the three dashes mark those headers
tests=("---convection---"
       "TESTING/conv/blob_3D_adaptive/blob3d-adaptive.sh"
       "TESTING/conv/blob_3D_equi/blob3d.sh"
       "TESTING/conv/blob_3D_nonequi/blob3d.sh"
       "TESTING/conv/blob_convection_equi/blob-convection-equi.sh"
       "TESTING/conv/blob_convection_2nd_serial/blob-conv-adaptive-serial.sh"
       "TESTING/conv/blob_convection/blob-convection-adaptive.sh"
       "TESTING/conv/blob_convection_2nd/blob-convection-adaptive-2nd.sh"
       "TESTING/conv/blob_convection_2nd_zcurve/blob-conv-adaptive-2nd-zcurve.sh"
       "---navier-stokes---"
       "TESTING/navier_stokes/pressure_blob/pressure_blob.sh"
       "TESTING/navier_stokes/funnel/ion_funnel_light.sh"
       "TESTING/navier_stokes/filter/bogey_shock/bogey_shock.sh"
       "TESTING/navier_stokes/filter/explicit_5pt/explicit_5pt.sh"
       "TESTING/navier_stokes/filter/explicit_7pt/explicit_7pt.sh"
       "TESTING/navier_stokes/filter/explicit_9pt/explicit_9pt.sh"
       "TESTING/navier_stokes/filter/explicit_11pt/explicit_11pt.sh"
       "---acm---"
       "TESTING/acm/acm_cyl_equi/acm_cylinder_equi.sh"
       "TESTING/acm/acm_cyl_nonequi/acm_cylinder_nonequi.sh"
       "TESTING/acm/acm_cyl_adaptive/acm_cylinder_adaptive.sh"
       "TESTING/acm/acm_cyl_adaptive/acm_cylinder_adaptive_zcurve.sh"
       )

happy_sum=0
sad_sum=0
numtests=0

echo "employed command for parallel exec: " $mpi_command
echo "to modify the command, set \$nprocs, \$mpi_command in shell"
echo "   "

if [ $nprocs != 4 ]; then
    echo "$fail_color WARNING $end_color"
    echo "your tests might fail because the keyvalues for load balancing may differ if you don't use nprocs=4 for testing"
fi

T="$(date +%s)"

for ts in ${tests[@]}
do
    if [[ $ts == "---"* ]]; then
        # this is a header line - we just need to print it on the screen
        echo $ts
    else
        # output usually sent to screen will be written to log file.
        logfile=${ts%%.sh}.log
        rm -f $logfile
        touch $logfile

        echo "Running test: " ${ts}
        echo "Writing output to: " ${logfile}

        # run the actual test
        ./${ts} > $logfile

        if [ $? == 0 ]; then
            printf "%s \n" "${pass_color}pass${end_color}"
            happy_sum=$((happy_sum+1))
            summary[$numtests]=0
        else
            printf "%s \n" "${fail_color}fail${end_color}"
            sad_sum=$((sad_sum+1))
            summary[$numtests]=1
        fi
        numtests=$((numtests+1))
        rm -f *.key *.h5 *.t *.dat
        printf "\n"
    fi
done

echo
T="$(($(date +%s)-T))"
echo "Time used in tests: ${T} seconds"

echo " "
echo "All in all we have: "
echo " "

numtests=0
for ts in ${tests[@]}
do
    if [[ $ts != "---"* ]]; then
        if [ ${summary[$numtests]} == 0 ]; then
            printf "%-80s %s \n" ${ts} "$pass_color ok $end_color"
        else
            printf "%-80s %s \n" ${ts} "$fail_color X $end_color"
        fi
        numtests=$((numtests+1))
    fi
done
echo " "

echo -e "\t $pass_color sum happy tests: $end_color \t" $happy_sum
echo -e "\t $fail_color sum sad tests: $end_color \t" $sad_sum
