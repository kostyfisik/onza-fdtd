#!/bin/bash
file=$1
rm -r *.png
for i in *.plot; do
if [[ $file == '3D_test_blitz.cc' ]]; then
 	sed "s/name/`echo $i`/g" time-blitz.plt | gnuplot
else
	sed "s/name/`echo $i`/g" time-stdc.plt | gnuplot
fi
done

