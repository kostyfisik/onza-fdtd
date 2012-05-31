#!/bin/bash
for i in *.onza ; do
 sed "s/calc1/`echo $i`/g" onza.plt | gnuplot
# sed "s/calc1/`echo $i`/g" data2.plt | gnuplot
done
