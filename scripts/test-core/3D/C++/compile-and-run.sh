#!/bin/bash
file=$1
mode=$2
rm -r *.plot
if [[ $mode == 'start' ]]; then
	./$file.bin
else
	mpic++ $file -O2 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -ftemplate-depth-30 -o $file.bin
	./$file.bin
fi
#./gnuplot-all-blocks_test.sh $file
