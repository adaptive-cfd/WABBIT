#!/bin/bash

#--------------------------------
# WABBIT run unit tests
#--------------------------------

echo "   "
echo -e "\t \033[4m WABBIT: run all existing unit tests \033[0m"
echo "   "

fail_color=$'\033[31;1m'
pass_color=$'\033[92;1m'
end_color=$'\033[0m'

if [ -z "$nprocs" ]; then
    echo "unset"
    nprocs=4
fi

if [ -z "$mpi_command" ]; then
    export mpi_command="nice mpirun -n ${nprocs} "
fi

if [ "${memory}" == "" ]; then
    export memory="--memory=5.0GB"
fi

if [ "$(which wabbit-compare-hdffiles.py)" == "" ]; then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "Since 15 Aug 2023, the unit testing framework has evolved. It now stores "
    echo "full HDF5 files in the TESTING directory, which makes it easier to visualize the reference data"
    echo "and current results, should they be different. We now calculate the L2 error of the field, if"
    echo "the grid is identical. This new framework requires the "
    echo "https://github.com/adaptive-cfd/python-tools repository for comparing two WABBIT HDF5 files. "
    echo ""
    echo " You do not seem to have the command ${fail_color}wabbit-compare-hdffiles.py${end_color} available! Either you do not have"
    echo " the repository, or its directory is not in you \$PYTHONPATH"
    echo ""
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "${fail_color}Cannot run unit tests !!${end_color}"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    exit 881
fi

#export mpi_serial="nice mpirun -n 1 "



# list all tests here. For readability, we include header lines like ---acm---
# they structure teh output on the screen. Note the three dashes mark those headers
tests=( "---post---"
        # "TESTING/wabbit_post/pod/pod_test.sh"
        "---wavelets---"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_2D_CDF20.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_2D_CDF22.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_2D_CDF40.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_2D_CDF42.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_2D_CDF44.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_2D_CDF60.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_2D_CDF62.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_3D_CDF20.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_3D_CDF22.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_3D_CDF40.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_3D_CDF42.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_3D_CDF44.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_3D_CDF60.sh"
        "TESTING/wavelets/equi_refineCoarsen_FWT-IWT_3D_CDF62.sh"
        "TESTING/wavelets/ghost_nodes_2D_2nd.sh"
        "TESTING/wavelets/ghost_nodes_2D_4th.sh"
        "TESTING/wavelets/ghost_nodes_2D_6th.sh"
        "TESTING/wavelets/ghost_nodes_3D_2nd.sh"
        "TESTING/wavelets/ghost_nodes_3D_4th.sh"
        "TESTING/wavelets/ghost_nodes_3D_6th.sh"
        "TESTING/wavelets/adaptive_CDF20/run.sh"
        "TESTING/wavelets/adaptive_CDF22/run.sh"
        "TESTING/wavelets/adaptive_CDF40/run.sh"
        "TESTING/wavelets/adaptive_CDF42/run.sh"
        "TESTING/wavelets/adaptive_CDF44/run.sh"
        "TESTING/wavelets/adaptive_CDF62/run.sh"
        "---convection---"
        "TESTING/conv/blob_2D_equi_4th/blob-convection-adaptive.sh"
        "TESTING/conv/blob_2D_adaptive_CDF20/blob-convection-adaptive.sh"
        "TESTING/conv/blob_2D_adaptive_CDF22/blob-convection-adaptive.sh"
        "TESTING/conv/blob_2D_adaptive_CDF40/blob-convection-adaptive.sh"
        "TESTING/conv/blob_2D_adaptive_CDF42/blob-convection-adaptive.sh"
        "TESTING/conv/blob_2D_adaptive_CDF44/blob-convection-adaptive.sh"
        "TESTING/conv/blob_3D_equi_2nd/blob3d.sh"
        "TESTING/conv/blob_3D_equi_4th/blob3d.sh"
        "TESTING/conv/blob_3D_adaptive_CDF40/blob3d-adaptive.sh"
        "TESTING/conv/blob_3D_adaptive_CDF22/blob3d-adaptive.sh"
        "TESTING/conv/blob_3D_adaptive_CDF44/blob3d-adaptive.sh"

        # "---navier-stokes---"
        # "TESTING/navier_stokes/simple_geometry/rhombus/rhombus_moving_shock.sh"
        # "TESTING/navier_stokes/simple_geometry/triangle/triangle_adapt.sh"
        # "TESTING/navier_stokes/pressure_blob/2D/pressure_blob.sh"
        # "TESTING/navier_stokes/pressure_blob/3D/pressure_blob3D.sh"
        # "TESTING/navier_stokes/funnel/2D/ion_funnel_light.sh"
        # "TESTING/navier_stokes/funnel/3D/ion_funnel_3D.sh"
        # "TESTING/navier_stokes/filter/bogey_shock/2D/bogey_shock.sh"
        # "TESTING/navier_stokes/filter/bogey_shock/3D/bogey_shock3D.sh"
        # "TESTING/navier_stokes/filter/explicit_5pt/explicit_5pt.sh"
        # "TESTING/navier_stokes/filter/explicit_7pt/explicit_7pt.sh"
        # "TESTING/navier_stokes/filter/explicit_9pt/explicit_9pt.sh"
        # "TESTING/navier_stokes/filter/explicit_11pt/explicit_11pt.sh"
        "---acm---"
        "TESTING/acm/acm_cyl_adaptive_CDF44/run.sh"
        )

happy_sum=0
sad_sum=0
numtests=0

echo "employed command for parallel exec: " $mpi_command
echo "memory flag for wabbit is: " $memory
echo "to modify the command, export \$memory, \$mpi_command in shell "
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
        T2="$(date +%s)"
        ./${ts} > $logfile
        err="$?"
        T2="$(($(date +%s)-T2))"


        if [ $err == 0 ]; then
            printf "%s \t" "${pass_color} pass ${end_color}"
            happy_sum=$((happy_sum+1))
            summary[$numtests]=0
        else
            printf "%s \t" "${fail_color} fail ${end_color}"
            sad_sum=$((sad_sum+1))
            summary[$numtests]=1
        fi
        echo "Time used in test: ${T2} seconds"
        numtests=$((numtests+1))
        rm -f *.key *.h5 *.t *.dat
        printf "\n"
    fi
done

echo ""
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
