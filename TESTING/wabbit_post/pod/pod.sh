#!/bin/bash
#-------------------------------------------------------------------------------
# THIS IS NOT AN INI FILE it is used as a workarround for update_unit_test
#-------------------------------------------------------------------------------
# what parameter file
dir=$1

${mpi_command} ./wabbit-post --generate_forest --Jmax=5
${mpi_command} ./wabbit-post --POD --components=1 --adapt=0.01 --nmodes=15 --list=${dir}u_list.txt ${memory}
${mpi_command} ./wabbit-post --POD-reconstruct --time_coefficients=a_coefs.txt --nmodes=5 --save_all --adapt=0.1 --mode-list=${dir}mode_list.txt ${memory}
