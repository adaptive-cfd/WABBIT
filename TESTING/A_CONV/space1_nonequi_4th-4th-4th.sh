#!/bin/bash


ini=blob-convection.ini
dx=(0 1 2 3 4)
test_number=1
test_name="nonequi_4th-4th-4th"
# delete all data:
rm -rf dx${test_number}_*
# run tests
for ddx in ${dx[@]}
do
	dir=dx${test_number}_${test_name}_${ddx}
	mkdir $dir
	cd $dir
	cp ../$ini .
	ln -s ../../wabbit

	# time
	replace_ini_value.sh $ini time_max 1.0
	replace_ini_value.sh $ini dt_fixed 0.0
	replace_ini_value.sh $ini dt_max 2.0e-3
	replace_ini_value.sh $ini CFL 1.2

	# order
	replace_ini_value.sh $ini order_discretization FD_4th_central_optimized
	replace_ini_value.sh $ini order_predictor multiresolution_4th

	# blocks
	replace_ini_value.sh $ini adapt_mesh 0 
	replace_ini_value.sh $ini adapt_inicond 1
	replace_ini_value.sh $ini inicond_refinements $ddx
	replace_ini_value.sh $ini number_block_nodes 17
	replace_ini_value.sh $ini number_ghost_nodes 4
	replace_ini_value.sh $ini eps 1.0e-4
	replace_ini_value.sh $ini max_treelevel 13
	replace_ini_value.sh $ini min_treelevel 2

	# other
	replace_ini_value.sh $ini nu 0.0
	replace_ini_value.sh $ini blob_width 0.01
	
	$mpi ./wabbit 2D $ini --memory=0.5GB
	cd ..
done

