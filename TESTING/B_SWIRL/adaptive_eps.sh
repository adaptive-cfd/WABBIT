#!/bin/bash

export mpi="nice mpirun -n 3"

#-------------
test1=1
#-------------

if [ "$test1" == "1" ]; then
	ini=swirl-tough.ini
	eps=(5.0e-5 1.0e-4 5.0e-4 1.0e-3)
	number=1
	pre=adaptive
	# delete all data:
	rm -r ${pre}${number}_*

	i=0
	for j in ${eps[@]}
	do
		dir=${pre}${number}_swirl-tough_${i}_${j}
		mkdir $dir
		echo $dir
		cd $dir
		cp ../$ini .

		ln -s ../../wabbit

		./replace_ini_value.sh $ini N_fields_saved 1
		./replace_ini_value.sh $ini field_names phi

		./replace_ini_value.sh $ini order_discretization FD_4th_central_optimized
		./replace_ini_value.sh $ini order_predictor multiresolution_4th

		./replace_ini_value.sh $ini adapt_mesh 1
		./replace_ini_value.sh $ini adapt_inicond 1
		./replace_ini_value.sh $ini eps $j
	 
		./replace_ini_value.sh $ini number_block_nodes 17
		./replace_ini_value.sh $ini number_ghost_nodes 4
		./replace_ini_value.sh $ini max_treelevel 13
		./replace_ini_value.sh $ini min_treelevel 1
		./replace_ini_value.sh $ini nu 0.0
		./replace_ini_value.sh $ini time_max 10.0
		./replace_ini_value.sh $ini CFL 1.0

		./replace_ini_value.sh $ini blob_width 0.01

		cleanhere -f
	
		$mpi ./wabbit 2D $ini --memory=0.75GB
		i=$((i+1))
		cd ..
	done
fi

