#!/bin/bash


#-------------
test1=1
#-------------


# -------------------------------------------------------------------------------
# adaptive, but with fixed Jmax

if [ "$test1" == "1" ]; then
	ini="swirl-tough.ini"


	JMAX=(3 4 5 6 7 8 9)
	eps=1.0e-7

	for jmax in ${JMAX[@]}
	do
		number=$jmax
		pre=fixeps_adapt_jmax
		# delete all data:
		rm -r ${pre}${number}_*

		dir=${pre}${number}_swirl-easy_${jmax}_${eps}

		if [ ! -f ${dir}/phi_000004000000.h5 ]; then
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
			../replace_ini_value.sh $ini Blocks eps "$eps"
			../replace_ini_value.sh $ini Blocks number_block_nodes 17
			../replace_ini_value.sh $ini Blocks number_ghost_nodes 4
			../replace_ini_value.sh $ini Blocks max_treelevel "$jmax"
			../replace_ini_value.sh $ini Blocks min_treelevel 1

			# T=4 "easy" test
			../replace_ini_value.sh $ini Time time_max 4.0
			../replace_ini_value.sh $ini Time CFL 1.0
			../replace_ini_value.sh $ini Time write_time 1.00

			$mpi ./wabbit 2D $ini --memory=3.5GB
			cd ..
		else
			echo "we skip!"
		fi
	done
fi
