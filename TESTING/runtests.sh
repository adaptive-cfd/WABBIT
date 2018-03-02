#!/bin/bash
#--------------------------------
# WABBIT run unit tests
#--------------------------------
rootdir="../"
cd $rootdir
echo -e "\t \033[4m WABBIT: run all existing unit tests \033[0m" 

fail_color=$'\033[31;1m'
pass_color=$'\033[92;1m'
end_color=$'\033[0m'

tests=("TESTING/navier_stokes/pressure_blob/pressure_blob.sh" "TESTING/acm/acm_cyl/acm_cylinder.sh")
happy=0
sad=0
i=0
for ts in ${tests[@]}
do
    sh ${ts}
    if [ $? == 0 ]; then
	happy=$((happy+1))
	summary[$i]=0
    else
	sad=$((sad+1))
	summary[$i]=1
    fi
    i=$((i+1))
done

echo " "
echo "All in all we have: "
echo " "
i=0
for ts in ${summary[@]}
do
    if [ ${summary[$i]} == 0 ]; then
	printf "%-80s %s \n" ${tests[$i]} "$pass_color ok $end_color"
    else
        printf "%-80s %s \n" ${tests[$i]} "$fail_color X $end_color"
    fi
    i=$((i+1))
done
echo " "

echo -e "\t $pass_color sum happy tests: $end_color \t" $happy
echo -e "\t $fail_color sum sad tests: $end_color \t" $sad
