#!/bin/bash
mode=$1 # control sum - a
##############################
path=$PWD; path_fortran=$path/Fortran; path_valgrind=$path/valgrind
file=2D_fortran.f95
max_size=2048; max_steps=100;
let "C = max_size**2 * max_steps"
##############################
# Standard launch: mode not specified
##############################
if [ $# -eq "0" ]; then
echo ">> Standard launch of 2D Fortran tests"
size=32
cd $path_fortran
while [ $size -le $max_size ]
do
    let "steps = C / size**2"
    gfortran $file -cpp -Dsz=$size -Dst=$steps -O2 -Dcs=0 -o $file.bin
    ./$file.bin
    rm -r *.bin
    let "size *= 2"
done
echo ">> Done running standard 2D Fortran tests"
##############################
elif [ $# -eq "1" ]; then
	##############################
    if [[ $mode == "a" ]]; then
	echo ">> Control sum mode selected"
	echo ">> Launches 2D Fortran tests with Control sum option"
	size=32
        cd $path_fortran
        while [ $size -le $max_size ]
        do
            let "steps = C / size**2"
	    gfortran $file -cpp -Dsz=$size -Dst=$steps -O2 -Dcs=1 -o $file.bin
            ./$file.bin
            rm -r *.bin
            let "size *= 2"
        done
	echo ">> Done running 2D Fortran tests with Control sum option"
        fi
	##############################
else
echo ">> Please specify the mode correctly!"
fi

