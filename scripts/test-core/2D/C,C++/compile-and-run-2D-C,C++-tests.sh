#!/bin/bash
mode=$1 # Control sum - a, valgrind - b
##############################
max_size=2048; max_steps=100;
let "C = max_size**2 * max_steps"
##############################
path=$PWD; path_blitz=$path/blitz; path_c=$path/c; path_valgrind=$path/valgrind
stencil=2D_test_stencil.cc
range=2D_test_range.cc
stencil_blocks=2D_test_stencil_blocks.cc
range_blocks=2D_test_range_blocks.cc
c=2D_test_c.c
# Standard launch: mode not specified
##############################
if [ $# -eq "0" ]; then
echo ">> Standard launch of 2D C/C++ tests"
size=32
while [ $size -le $max_size ]
do
    let "steps = C / size**2"
    cd $path_blitz
    mpic++ $stencil -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $stencil.bin
    ./$stencil.bin
    mpic++ $range -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $range.bin
    ./$range.bin
    cd $path_c
    mpicc -O3 -Dsz=$size -Dst=$steps -Dcs=0 -fstrict-aliasing -std=c99 -lm $c -o $c.bin
    ./$c.bin
    rm -r *.bin
    let "size *= 2"
done
echo ">> Done running standard 2D C/C++ tests"
##############################
elif [ $# -eq "1" ]; then
	##############################
    if [[ $mode == "a" ]]; then
	echo ">> Control sum mode selected"
	echo ">> Launches 2D C/C++ tests with control sum option"
	size=32
        while [ $size -le $max_size ]
        do
            let "steps = C / size**2"
	    cd $path_blitz
            mpic++ $stencil -O2 -Dsz=$size -Dst=$steps -Dcs=1 -ftemplate-depth-30 -o $stencil.bin
            ./$stencil.bin
            mpic++ $range -O2 -Dsz=$size -Dst=$steps -Dcs=1 -ftemplate-depth-30 -o $range.bin
            ./$range.bin
	    cd $path_c
            mpicc -O3 -Dsz=$size -Dst=$steps -Dcs=1 -fstrict-aliasing -std=c99 -lm $c -o $c.bin
            ./$c.bin
            rm -r *.bin
            let "size *= 2"
        done
	echo ">> Done running 2D C/C++ tests with Control sum option"
	##############################
	elif [[ $mode == "b" ]]; then
	echo ">> Valgrind mode selected"
	echo ">> Launches standard 2D C/C++ tests and valgrind profiling for 0.1x timesteps"
	size=32
        while [ $size -le $max_size ]
        do
            let "steps = C / size**2"
            cd $path_blitz
            mpic++ $stencil -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $stencil.bin
            ./$stencil.bin
            mpic++ $range -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $range.bin
            ./$range.bin
            cd $path_c
            mpicc -O3 -Dsz=$size -Dst=$steps -Dcs=0 -fstrict-aliasing -std=c99 -lm $c -o $c.bin
            ./$c.bin
            rm -r *.bin
            let "size *= 2"
        done
	cd $path_valgrind
        size=32
        while [ $size -le $max_size ]
        do
            let "steps = C / (size**2 * 10)"
            mpic++ $stencil -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $stencil.bin
            ./valgrind.sh $stencil.bin
            mpic++ $range -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $range.bin
            ./valgrind.sh $range.bin
	    mpicc -O3 -Dsz=$size -Dst=$steps -Dcs=0 -fstrict-aliasing -std=c99 -lm $c -o $c.bin
            ./valgrind $c.bin
            rm -r *.bin
            let "size *= 2"
        done
	echo ">> Done running standard 2D C/C++ tests and valgrind profiling for 0.1x timesteps"
	##############################
	else
	echo ">> Unknown mode!"
	fi
else
echo ">> Please specify the mode correctly!"
fi
