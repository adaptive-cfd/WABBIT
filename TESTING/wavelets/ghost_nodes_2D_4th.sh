#!/bin/bash

# run actual test
${mpi_command} ./wabbit-post --ghost-nodes-test --wavelet=CDF40 ${memory} --dim=2

if [ "$?" == 0 ]; then
    exit 0
else
    exit 999
fi
