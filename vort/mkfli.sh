#!/bin/bash
#
# check for existence of arguments
if ([ -z $1 ] || [ -z $2 ]) ; then {
	echo "usage: sh mkfli.sh nx ny"
	exit
} ; fi
#
# convert to something sensible
for fil in `ls im????`
do
	r8toras $1 $2 /home/gammie/bin/john.pal < $fil > $fil.ras
	convert $fil.ras $fil.ppm
	rm $fil.ras
done
#
# make a list file
ls im????.ppm > tmp.lis
#
# make the fli
ppm2fli -g$1x$2 -N tmp.lis rho.fli
#
# tidy up
rm im????.ppm
#
