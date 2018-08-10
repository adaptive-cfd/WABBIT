#!/bin/bash
##########################################################################################
# This Script performs some changes in the input ini-file:
# 	- it removes the section [DomainSize]
# 	- it adds the relevant DomainSize variables to the new section [Domain]
##########################################################################################

if [ "$1" == "" ]; then
	echo "NO ARGUMENT GIVEN!"
	echo "usage: replace_ini_value.sh PARAMS.ini"
	echo "replace_ini_value.sh PARAMS.ini Time time_max 3.0"
	exit 6
fi
file=$1

sed  -i "s+\[Dimensionality\]+[Domain]+g" $file

Lx=$(grep 'Lx=' $file | egrep -o '[[:digit:]]+\.*[[:digit:]]*')
Ly=$(grep 'Ly=' $file | egrep -o '[[:digit:]]+\.*[[:digit:]]*')
Lz=$(grep 'Lz=' $file | egrep -o '[[:digit:]]+\.*[[:digit:]]*')

# add new lines after dim with domain_size
sed -i '/dim=[2-3]/a domain_size='$Lx' '$Ly' '$Lz';' $file
# remove Domain_size section
sed -i '/\[DomainSize\]/,/Lz=.*/d' $file

echo "created new section [domain] and removed [Domain], [DomainSize]"
