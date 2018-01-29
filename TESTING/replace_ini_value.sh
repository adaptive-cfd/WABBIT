#!/bin/bash


if [ "$1" == "" ]; then
	echo "NO ARGUMENT GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini CFL 2.0 [section]"
	exit 6
fi
file=$1

if [ "$2" == "" ]; then
	echo "NO PARAMETER GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini CFL 2.0 [section]"
	exit 6
fi
key=$2

if [ "$3" == "" ]; then
	echo "NO NEW VALUE GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini CFL 2.0 [section]"
	exit 6
fi
value=$3



if [ "$section" == "" ]; then
	# no section argument given, just replace first occurance

	current=$(grep "$key=" $file | head -n 1)
	new="$key=$value;"

	sed  -i "s+$current+$new+g" $file

	new=$(grep "$key=" $file | head -n 1)
	echo Changed: $current to $new

else
	# with section argument
	section=$4
	found=0
	
	while read -r line
	do		
		search="\[${section}\]"

		# found section begin
		if [ ! "$(echo $line | grep $search)" == "" ]; then
			found=1
		fi

		# if section found, look for string to replace
		if [ "$found" == 1 ]; then
			if [ ! "$(echo $line | grep $key\=)" == "" ]; then
				current=$line
				if [ ! "$current" == "" ]; then
					new="$key=$value;"

					sed  -i "s+$current+$new+g" $file

					new=$(grep "$key=" $file | head -n 1)
					echo Changed: $current to $new
				fi
			fi
		fi
	done < $file

fi
