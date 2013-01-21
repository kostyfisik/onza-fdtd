#!/bin/bash
for i in *.f95; do
	gfortran -O2 -o $i.bin $i
	./$i.bin
done
