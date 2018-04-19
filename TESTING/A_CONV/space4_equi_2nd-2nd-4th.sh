#!/bin/bash

ini=blob-convection.ini
dx=(2 3 4 5)
test_number=4
test_name="equi_2nd-2nd-4th"

# delete all data:
rm -rf dx${test_number}_*
# run tests
for ddx in ${dx[@]}
do
	dir=dx${test_number}_${test_name}_${ddx}
	mkdir $dir
	cd $dir
	cp ../$ini .
	ln -s ../../../wabbit

	# time
	../replace_ini_value.sh $ini Time time_max 1.0
	../replace_ini_value.sh $ini Time dt_fixed 0.0
	../replace_ini_value.sh $ini Time dt_max 2.0e-3
	../replace_ini_value.sh $ini Time CFL 1.2

	# order
	../replace_ini_value.sh $ini Discretization order_discretization FD_2nd_central
	../replace_ini_value.sh $ini Discretization order_predictor multiresolution_2nd

	# blocks
	../replace_ini_value.sh $ini Blocks adapt_mesh 0
	../replace_ini_value.sh $ini Blocks adapt_inicond 0
	../replace_ini_value.sh $ini Blocks inicond_refinements 0
	../replace_ini_value.sh $ini Blocks number_block_nodes 17
	../replace_ini_value.sh $ini Blocks number_ghost_nodes 4
	../replace_ini_value.sh $ini Blocks eps 1.0e-3
	../replace_ini_value.sh $ini Blocks max_treelevel $ddx
	../replace_ini_value.sh $ini Blocks min_treelevel $ddx

	# other
	../replace_ini_value.sh $ini ConvectionDiffusion nu 0.0
	../replace_ini_value.sh $ini ConvectionDiffusion blob_width 0.01

	$mpi ./wabbit 2D $ini --memory=3.0GB
	cd ..
done
