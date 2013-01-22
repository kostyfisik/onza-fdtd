#!/bin/bash
mode=$1 # blocks - a, valgrind - b, both - c
blitz=2D_test_blitz.cc
path=$PWD; path_blitz=$path/blitz; path_c=$path/c; path_valgrind=$path/valgrind
# Standard launch: mode not specified
##############################
if [ $# -eq "0" ]; then
echo ">> Standard launch of 2D C/C++ tests"
cd $path_blitz
mpic++ $blitz -O2 -Dblk=0 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -ftemplate-depth-30 -o $blitz.bin
./$blitz.bin
cd $path_c
for i in *.c; do
	gcc -O3 -fstrict-aliasing -std=c99 -lm $i -o $i.bin
	./$i.bin
done
echo ">> Done running standard 2D C/C++ tests"
##############################
elif [ $# -eq "1" ]; then
	##############################
	if [[ $mode == "a" ]]; then
	echo ">> Blocks mode selected"
	echo ">> Launches 2D C/C++ tests with blocks option"
	cd $path_blitz
	mpic++ $blitz -O2 -Dblk=1 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -ftemplate-depth-30 -o $blitz.bin
	./$blitz.bin
	cd $path_c
	for i in *.c; do
	gcc -O3 -fstrict-aliasing -std=c99 -lm $i -o $i.bin
	./$i.bin
	done
	echo ">> Done running 2D C/C++ tests with blocks option"
	##############################
	elif [[ $mode == "b" ]]; then
	echo ">> Valgrind mode selected"
	echo ">> Launches standard 2D C/C++ tests and valgrind profiling for 0.1x timesteps"
	cd $path_blitz
	mpic++ $blitz -O2 -Dblk=0 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -ftemplate-depth-30 -o $blitz.bin
	./$blitz.bin
	cd $path_c
	for i in *.c; do
		gcc -O3 -fstrict-aliasing -std=c99 -lm $i -o $i.bin
	./$i.bin
	done
	cd $path_valgrind
	mpic++ $blitz -O2 -Dblk=0 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -ftemplate-depth-30 -o $blitz.bin
	./valgrind.sh $blitz.bin
	for i in *.c; do
		gcc -O3 -fstrict-aliasing -std=c99 -lm $i -o $i.bin
		./valgrind.sh $i.bin
	done
	echo ">> Done running standard 2D C/C++ tests and valgrind profiling for 0.1x timesteps"
	##############################
	elif [[ $mode == "c" ]]; then
	echo ">> Blocks + Valgrind mode selected"
	echo ">> Launches 2D C/C++ tests with blocks option and valgrind profiling for 0.1x timesteps"
	cd $path_blitz
	mpic++ $blitz -O2 -Dblk=1 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -ftemplate-depth-30 -o $blitz.bin
	./$blitz.bin
	cd $path_c
	for i in *.c; do
	gcc -O3 -fstrict-aliasing -std=c99 -lm $i -o $i.bin
	./$i.bin
	done
	cd $path_valgrind
	mpic++ $blitz -O2 -Dblk=1 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -ftemplate-depth-30 -o $blitz.bin
	./valgrind.sh $blitz.bin
	for i in *.c; do
	gcc -O3 -fstrict-aliasing -std=c99 -lm $i -o $i.bin
	./valgrind.sh $i.bin
	done
	echo ">> Done running 2D C/C++ tests with blocks option and valgrind profiling for 0.1x timesteps"
	##############################
	else
	echo ">> Unknown mode!"
	fi
else
echo ">> Please specify the mode correctly!"
fi
