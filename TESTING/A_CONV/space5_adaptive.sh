#!/bin/bash

#-------------
test1=1
test2=1
#-------------

if [ "$test1" == "1" ]; then
	ini=adv-adaptive.ini
	eps=(1.00000000e-07   1.34339933e-07   1.80472177e-07   2.42446202e-07
	   3.25702066e-07   4.37547938e-07   5.87801607e-07   7.89652287e-07
	   1.06081836e-06   1.42510267e-06   1.91448198e-06   2.57191381e-06
	   3.45510729e-06   4.64158883e-06   6.23550734e-06   8.37677640e-06
	   1.12533558e-05   1.51177507e-05   2.03091762e-05   2.72833338e-05
	   3.66524124e-05   4.92388263e-05   6.61474064e-05   8.88623816e-05
	   1.19377664e-04   1.60371874e-04   2.15443469e-04   2.89426612e-04
	   3.88815518e-04   5.22334507e-04   7.01703829e-04   9.42668455e-04
	   1.26638017e-03   1.70125428e-03   2.28546386e-03   3.07029063e-03
	   4.12462638e-03   5.54102033e-03   7.44380301e-03   1.00000000e-02)
	number=1
	pre=adapt
	# delete all data:
	rm -r ${pre}${number}_*

	i=0
	for j in ${eps[@]}
	do
		dir=${pre}${number}_${i}_${j}
		mkdir $dir
		echo $dir
		cd $dir
		cp ../$ini .

		ln -s ../../../wabbit

		../replace_ini_value.sh $ini Saving N_fields_saved 1
		../replace_ini_value.sh $ini Saving field_names phi


		../replace_ini_value.sh $ini Discretization order_discretization FD_4th_central_optimized
		../replace_ini_value.sh $ini Discretization order_predictor multiresolution_4th

		../replace_ini_value.sh $ini Blocks adapt_mesh 1
		../replace_ini_value.sh $ini Blocks adapt_inicond 1
		../replace_ini_value.sh $ini Blocks eps $j
		../replace_ini_value.sh $ini Blocks number_block_nodes 17
		../replace_ini_value.sh $ini Blocks number_ghost_nodes 4
		../replace_ini_value.sh $ini Blocks max_treelevel 13
		../replace_ini_value.sh $ini Blocks min_treelevel 1

		../replace_ini_value.sh $ini ConvectionDiffusion nu 0.0
		../replace_ini_value.sh $ini ConvectionDiffusion blob_width 0.01

		../replace_ini_value.sh $ini Time time_max 1.0
		../replace_ini_value.sh $ini Time CFL 1.0

		$mpi ./wabbit 2D $ini --memory=3.0GB
		i=$((i+1))
		cd ..
	done
fi

#--------------------------------------------------------------------------------

if [ "$test2" == "1" ]; then
	ini=adv-adaptive.ini
	eps=( 2.06913808e-05   2.97635144e-05
	   4.28133240e-05   6.15848211e-05   8.85866790e-05   1.27427499e-04
	   1.83298071e-04   2.63665090e-04   3.79269019e-04   5.45559478e-04
	   7.84759970e-04   1.12883789e-03   1.62377674e-03   2.33572147e-03
	   3.35981829e-03   4.83293024e-03   6.95192796e-03   1.00000000e-02)
	number=2
	pre=adapt
	# delete all data:
	rm -r ${pre}${number}_*

	i=0
	for j in ${eps[@]}
	do
		dir=${pre}${number}_${i}_${j}
		mkdir $dir
		echo $dir
		cd $dir
		cp ../$ini .

		ln -s ../../../wabbit

		../replace_ini_value.sh $ini Saving N_fields_saved 1
		../replace_ini_value.sh $ini Saving field_names phi

		../replace_ini_value.sh $ini Discretization order_discretization FD_2nd_central
		../replace_ini_value.sh $ini Discretization order_predictor multiresolution_2nd

		../replace_ini_value.sh $ini Blocks adapt_mesh 1
		../replace_ini_value.sh $ini Blocks adapt_inicond 1
		../replace_ini_value.sh $ini Blocks eps $j
		../replace_ini_value.sh $ini Blocks number_block_nodes 17
		../replace_ini_value.sh $ini Blocks number_ghost_nodes 4
		../replace_ini_value.sh $ini Blocks max_treelevel 13
		../replace_ini_value.sh $ini Blocks min_treelevel 1

		../replace_ini_value.sh $ini ConvectionDiffusion nu 0.0
		../replace_ini_value.sh $ini ConvectionDiffusion blob_width 0.01

		../replace_ini_value.sh $ini Time time_max 1.0
		../replace_ini_value.sh $ini Time CFL 1.0

		$mpi ./wabbit 2D $ini --memory=3.0GB
		i=$((i+1))
		cd ..
	done
fi

#--------------------------------------------------------------------------------

