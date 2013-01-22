#!/bin/bash
mode=$1 # valgrind - b
path=$PWD; path_fortran=$path/Fortran; path_valgrind=$path/valgrind
# Standard launch: mode not specified
##############################
if [ $# -eq "0" ]; then
echo ">> Standard launch of 3D Fortran tests"
cd $path_fortran
for i in *.f95; do
	gfortran -O2 -o $i.bin $i
	./$i.bin
done
echo ">> Done running standard 3D Fortran tests"
##############################
elif [ $# -eq "1" ]; then
	##############################
	if [[ $mode == "b" ]]; then
	echo ">> Valgrind mode selected"
	echo ">> Launches 3D Fortran tests and valgrind profiling for 0.1x timesteps"
	cd $path_fortran
	for i in *.f95; do
		gfortran -O2 -o $i.bin $i
		./$i.bin
	done
	cd $path_valgrind
	for i in *.f95; do
		gfortran -O2 -o $i.bin $i
		./valgrind.sh $i.bin
	done
	echo ">> Done running 3D Fortran tests and valgrind profiling for 0.1x timesteps"
	##############################	
	else
	echo ">> Unknown mode!"
	fi
else
echo ">> Please specify the mode correctly!"
fi

