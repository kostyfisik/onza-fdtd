#!/bin/bash
mode=$1 # control sum - a, valgrind - b
##############################
max_size=128; max_steps=100;
let "C = max_size**3 * max_steps"
##############################
path=$PWD; path_blitz=$path/blitz; path_c=$path/c; path_valgrind=$path/valgrind
stencil=3D_test_stencil.cc
range=3D_test_range.cc
stencil_blocks=3D_test_stencil_blocks.cc
range_blocks=3D_test_range_blocks.cc
c=3D_test_c.c
# Standard launch: mode not specified
##############################
if [ $# -eq "0" ]; then
echo ">> Standard launch of 3D C/C++ tests"
size=8
while [ $size -le $max_size ]
do
    let "steps = C / size**3"
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
echo ">> Done running standard 3D C/C++ tests"
##############################
elif [ $# -eq "1" ]; then
	##############################
    if [[ $mode == "a" ]]; then
	echo ">> Control sum mode selected"
	echo ">> Launches 3D C/C++ tests with Control sum option"
	size=8
        while [ $size -le $max_size ]
        do
            let "steps = C / size**3"
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
	echo ">> Done running 3D C/C++ tests with Control sum option"
	##############################
	elif [[ $mode == "b" ]]; then
	echo ">> Valgrind mode selected"
	echo ">> Launches standard 3D C/C++ tests and valgrind profiling for 0.1x timesteps"
	size=32
        while [ $size -le $max_size ]
        do
            let "steps = C / size**3"
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
        size=8
        while [ $size -le $max_size ]
        do
            let "steps = C / (size**3 * 10)"
            mpic++ $stencil -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $stencil.bin
            ./valgrind.sh $stencil.bin
            mpic++ $range -O2 -Dsz=$size -Dst=$steps -Dcs=0 -ftemplate-depth-30 -o $range.bin
            ./valgrind.sh $range.bin
	    mpicc -O3 -Dsz=$size -Dst=$steps -Dcs=0 -fstrict-aliasing -std=c99 -lm $c -o $c.bin
            ./valgrind $c.bin
            rm -r *.bin
            let "size *= 2"
        done
	echo ">> Done running standard 3D C/C++ tests and valgrind profiling for 0.1x timesteps"
	##############################
	else
	echo ">> Unknown mode!"
	fi
else
echo ">> Please specify the mode correctly!"
fi
