#!/bin/bash
for i in *.dat; do
 sed "s/name/`echo $i`/g" mean_time.plt | gnuplot
 sed "s/name/`echo $i`/g" efficiency.plt | gnuplot
done
for i in *.fdat; do
 sed "s/name/`echo $i`/g" mean_time_fixed.plt | gnuplot
done
