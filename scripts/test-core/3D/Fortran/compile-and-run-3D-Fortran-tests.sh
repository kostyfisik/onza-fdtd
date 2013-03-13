#!/bin/bash
mode=$1 # control sum - a
##############################
path=$PWD; path_fortran=$path/Fortran; path_valgrind=$path/valgrind
file=3D_fortran.f95
max_size=128; max_steps=100;
let "C = max_size**3 * max_steps"
##############################
# Standard launch: mode not specified
##############################
if [ $# -eq "0" ]; then
echo ">> Standard launch of 3D Fortran tests"
size=8
cd $path_fortran
while [ $size -le $max_size ]
do
    let "steps = C / size**3"
    gfortran $file -cpp -Dsz=$size -Dst=$steps -O2 -Dcs=0 -o $file.bin
    ./$file.bin
    rm -r *.bin
    let "size *= 2"
done
echo ">> Done running standard 3D Fortran tests"
##############################
elif [ $# -eq "1" ]; then
	##############################
    if [[ $mode == "a" ]]; then
	echo ">> Control sum mode selected"
	echo ">> Launches 3D Fortran tests with Control sum option"
	size=8
        cd $path_fortran
        while [ $size -le $max_size ]
        do
            let "steps = C / size**3"
	    gfortran $file -cpp -Dsz=$size -Dst=$steps -O2 -Dcs=1 -o $file.bin
            ./$file.bin
            rm -r *.bin
            let "size *= 2"
        done
	echo ">> Done running 3D Fortran tests with Control sum option"
        fi
	##############################
else
echo ">> Please specify the mode correctly!"
fi
