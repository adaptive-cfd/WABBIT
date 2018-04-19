#!/bin/bash


if [ "$1" == "" ]; then
	echo "NO ARGUMENT GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini SECTION KEYWORD VALUE"
	echo "replace_ini_value.sh PARAMS.ini Time time_max 3.0"
	exit 6
fi
file=$1

if [ "$2" == "" ]; then
	echo "NO SECTION GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini SECTION KEYWORD VALUE"
	echo "replace_ini_value.sh PARAMS.ini Time time_max 3.0"
	exit 6
fi
section=$2

if [ "$3" == "" ]; then
	echo "NO KEYWORD GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini SECTION KEYWORD VALUE"
	echo "replace_ini_value.sh PARAMS.ini Time time_max 3.0"
	exit 6
fi
key=$3

if [ "$3" == "" ]; then
	echo "NO VALUE GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini SECTION KEYWORD VALUE"
	echo "replace_ini_value.sh PARAMS.ini Time time_max 3.0"
	exit 6
fi
value=$4


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
				exit 0
			fi
		fi
	fi
done < $file
