#!/bin/bash
for i in *.dat; do
 sed "s/name/`echo $i`/g" mean-time.plt | gnuplot
 sed "s/name/`echo $i`/g" efficiency.plt | gnuplot
done
for i in *.fdat; do
 sed "s/name/`echo $i`/g" mean-time-fixed-mpirun-parameters.plt | gnuplot
done
