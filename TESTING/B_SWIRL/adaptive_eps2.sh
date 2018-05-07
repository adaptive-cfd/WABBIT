#!/bin/bash


#-------------
test1=1
test2=1
#-------------

if [ "$test1" == "1" ]; then
	ini=swirl-tough.ini
	eps=( 1.00000000e-02   6.15848211e-03   3.79269019e-03
         2.33572147e-03   1.43844989e-03   8.85866790e-04
         5.45559478e-04   3.35981829e-04   2.06913808e-04
         1.27427499e-04   7.84759970e-05   4.83293024e-05
         2.97635144e-05   1.83298071e-05   1.12883789e-05
         6.95192796e-06   4.28133240e-06   2.63665090e-06
         1.62377674e-06   1.00000000e-06)
	number=2
	pre=adaptive
	# delete all data:
	rm -r ${pre}${number}_*

	i=0
	for j in ${eps[@]}
	do
		dir=${pre}${number}_swirl-easy_${i}_${j}
		mkdir $dir
		echo $dir
		cd $dir
		cp ../$ini .

		ln -s ../../../wabbit

		../replace_ini_value.sh $ini Saving N_fields_saved 1
		../replace_ini_value.sh $ini Saving field_names phi
		../replace_ini_value.sh $ini Time write_time 1.25

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
		../replace_ini_value.sh $ini ConvectionDiffusion  blob_width 0.01
# T=4 "easy" test
		../replace_ini_value.sh $ini Time time_max 4.0
		../replace_ini_value.sh $ini Time CFL 1.0

		$mpi ./wabbit 2D $ini --memory=3.0GB
		i=$((i+1))
		cd ..
	done
fi

# -------------------------------------------------------------------------------
# adaptive, but with fixed Jmax

if [ "$test2" == "1" ]; then
	ini="swirl-tough.ini"


	JMAX=(4 5 6 7)
	EPS=( 1.00000000e-02   6.15848211e-03   3.79269019e-03
         2.33572147e-03   1.43844989e-03   8.85866790e-04
         5.45559478e-04   3.35981829e-04   2.06913808e-04
         1.27427499e-04   7.84759970e-05   4.83293024e-05
         2.97635144e-05   1.83298071e-05   1.12883789e-05
         6.95192796e-06   4.28133240e-06   2.63665090e-06
         1.62377674e-06   1.00000000e-06)

	for jmax in ${JMAX[@]}
	do
		number=$jmax
		pre=adapt_jmax
		# delete all data:
		rm -r ${pre}${number}_*

		i=0
		for eps in ${EPS[@]}
		do
			dir=${pre}${number}_swirl-easy_${i}_${eps}
			if [ ! -f ${dir}/phi_000004000000.h5 ]; then
				mkdir $dir
				echo $dir
				cd $dir
				cp ../$ini .

				ln -s ../../../wabbit

				../replace_ini_value.sh $ini Saving N_fields_saved 1
				../replace_ini_value.sh $ini Saving field_names phi
				../replace_ini_value.sh $ini Time write_time 1.25

				../replace_ini_value.sh $ini Discretization order_discretization FD_4th_central_optimized
				../replace_ini_value.sh $ini Discretization order_predictor multiresolution_4th

				../replace_ini_value.sh $ini Blocks adapt_mesh 1
				../replace_ini_value.sh $ini Blocks adapt_inicond 1
				../replace_ini_value.sh $ini Blocks eps "$eps"
				../replace_ini_value.sh $ini Blocks number_block_nodes 17
				../replace_ini_value.sh $ini Blocks number_ghost_nodes 4
				../replace_ini_value.sh $ini Blocks max_treelevel "$jmax"
				../replace_ini_value.sh $ini Blocks min_treelevel 1

				../replace_ini_value.sh $ini ConvectionDiffusion nu 0.0
				../replace_ini_value.sh $ini ConvectionDiffusion blob_width 0.01
# T=4 "easy" test
				../replace_ini_value.sh $ini Time time_max 4.0
				../replace_ini_value.sh $ini Time CFL 1.0

				$mpi ./wabbit 2D $ini --memory=3.0GB
				cd ..
			else
				echo "we skip!"
			fi
			i=$((i+1))
		done
	done
fi
