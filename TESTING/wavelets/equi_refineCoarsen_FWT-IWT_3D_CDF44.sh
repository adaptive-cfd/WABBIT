#!/bin/bash

##################
wavelet="CDF44"
##################

${mpi_command} ./wabbit-post --refine-coarsen-test --wavelet=${wavelet} ${memory} --dim=3

if [ ! "$?" == 0 ]; then
    exit 999
fi


${mpi_command} ./wabbit-post --wavelet-decomposition-unit-test --wavelet=${wavelet} ${memory} --dim=3

if [ "$?" == 0 ]; then
    exit 0
else
    exit 999
fi
