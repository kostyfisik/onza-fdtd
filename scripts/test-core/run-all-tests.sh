#!/bin/bash
mode=$1 # blocks - a, valgrind - b, both - c
Usage=">> This script runs all FDTD tests located in ~/test-core.\n
	>> Usage: $./run-all-tests.sh mode\n
	mode:\n
	a - enables block calculation,\n
	b - enables valgrind profiling,\n
	c - enables both options."
if [[ $# > 1 ]]; then
	echo -e $Usage
else

path_FDTD_test=$PWD;   path_2D=$path_FDTD_test/2D;   path_3D=$path_FDTD_test/3D
path_2D_cpp=$path_2D/C,C++; path_2D_fortran=$path_2D/Fortran
path_3D_cpp=$path_3D/C,C++; path_3D_fortran=$path_3D/Fortran

# Running 2D tests
cd $path_2D_cpp
chmod +x *.sh
./compile-and-run-2D-C,C++-tests.sh $mode
cd $path_2D_fortran
chmod +x *.sh
tmp_1=$mode
##############################
if [[ $tmp_1 == "a" ]]; then
echo -e ">> Block calculation option is not available for Fortran.
	Will use standard launch instead."
./compile-and-run-2D-Fortran-tests.sh
##############################
elif [[ $tmp_1 == "c" ]]; then
	echo ">> Block calculation option is not available for Fortran. 
	Will use only valgrind launch."
	tmp_1="b"
	##############################
fi
./compile-and-run-2D-Fortran-tests.sh $tmp_1
# End Running 2D tests

# Running 3D tests
cd $path_3D_cpp
chmod +x *.sh
./compile-and-run-3D-C,C++-tests.sh $mode
cd $path_3D_fortran
chmod +x *.sh
tmp_2=$mode
##############################
if [[ $tmp_2 == "a" ]]; then
echo -e ">> Block calculation option is not available for Fortran.
	Will use standard launch instead."
./compile-and-run-3D-Fortran-tests.sh
##############################
elif [[ $tmp_2 == "c" ]]; then
	echo ">> Block calculation option is not available for Fortran. 
	Will use only valgrind launch."
	tmp_2="b"
	##############################
fi
./compile-and-run-3D-Fortran-tests.sh $tmp_2
# End Running 3D tests
fi


