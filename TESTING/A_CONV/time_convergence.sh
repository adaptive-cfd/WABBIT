#!/bin/bash


#---------------------------------------------------------------------------------------------------

ini=adv-adaptive-RK2.ini
dt=(5.0e-2 2.0e-2 1.5e-2 1.2e-2 1.1e-2 1.0e-2 0.9e-2 5.0e-3 4.0e-3 3.0e-3 )
test_number=1
test_name="adaptive-4th-4th-2nd"

# delete all data:
rm -rf dt${test_number}_*

# run tests
for ddt in ${dt[@]}
do
	dir=dt${test_number}_${test_name}_${ddt}
	mkdir $dir
	cd $dir
	cp ../$ini .
	ln -s ../../../wabbit

	# time
	../replace_ini_value.sh $ini Time time_max 1.0
	../replace_ini_value.sh $ini Time dt_fixed $ddt
	../replace_ini_value.sh $ini Time CFL 9.9

	# order
	../replace_ini_value.sh $ini Discretization order_discretization FD_4th_central_optimized
	../replace_ini_value.sh $ini Discretization order_predictor multiresolution_4th

	# blocks
	../replace_ini_value.sh $ini Blocks adapt_mesh 0
	../replace_ini_value.sh $ini Blocks adapt_inicond 0
	../replace_ini_value.sh $ini Blocks inicond_refinements 0
	../replace_ini_value.sh $ini Blocks number_block_nodes 33
	../replace_ini_value.sh $ini Blocks number_ghost_nodes 4
	../replace_ini_value.sh $ini Blocks eps 1.0e-4
	../replace_ini_value.sh $ini Blocks max_treelevel 12
	../replace_ini_value.sh $ini Blocks min_treelevel 2

	# other
	../replace_ini_value.sh $ini ConvectionDiffusion nu 0.0
	../replace_ini_value.sh $ini ConvectionDiffusion blob_width 0.01

	$mpi ./wabbit 2D $ini --memory=3.0GB
	cd ..
done
