#!/bin/bash
mode=$1 # Control sum - a, valgrind - b
Usage=">> This script runs all FDTD tests located in ~/test-core.\n
	>> Usage: $./run-all-tests.sh mode\n
	mode:\n
	a - enables control sum calculation,\n
	b - enables valgrind profiling."
if [[ $# > 1 ]]; then
	echo -e $Usage
else
  if [[ $# == 1 ]]; then
      echo -e ">> Control sum mode selected.\n"
      dummy=1
  else
      echo -e ">> Time-only mode selected.\n"
      dummy=0
  fi
path_FDTD_test=$PWD;   path_2D=$path_FDTD_test/2D;   path_3D=$path_FDTD_test/3D
path_2D_cpp=$path_2D/C,C++; path_2D_fortran=$path_2D/Fortran
path_3D_cpp=$path_3D/C,C++; path_3D_fortran=$path_3D/Fortran

# Running 2D tests
echo -e "\n>> Running 2D tests\n"
max_size=2048; max_steps=100;
path=$path_2D; path_blitz=$path/C,C++/blitz; path_c=$path/C,C++/c; path_fortran=$path/Fortran/Fortran;
stencil=2D_test_stencil.cc; range=2D_test_range.cc; c=2D_test_c.c; fortran=2D_fortran.f95
let "C = max_size**2 * max_steps"
size=32
while [ $size -le $max_size ]
do
        let "steps = C / size**2"
        #########################
        cd $path_fortran
        gfortran $fortran -cpp -Dsz=$size -Dst=$steps -O2 -Dcs=$dummy -o $fortran.bin
        ./$fortran.bin
        rm -r *.bin
        #########################
        cd $path_blitz
        mpic++ $stencil -O2 -Dsz=$size -Dst=$steps -Dcs=$dummy -ftemplate-depth-30 -o $stencil.bin
        ./$stencil.bin
        mpic++ $range -O2 -Dsz=$size -Dst=$steps -Dcs=$dummy -ftemplate-depth-30 -o $range.bin
        ./$range.bin
        rm -r *.bin
        cd $path_c
        mpicc -O3 -Dsz=$size -Dst=$steps -Dcs=$dummy -fstrict-aliasing -std=c99 -lm $c -o $c.bin
        ./$c.bin
        rm -r *.bin
        #########################
        let "size *= 2"
done
# End Running 2D tests

# Running 3D tests
echo -e "\n>> Running 3D tests\n"
max_size=128; max_steps=100;
path=$path_3D; path_blitz=$path/C,C++/blitz; path_c=$path/C,C++/c; path_fortran=$path/Fortran/Fortran;
stencil=3D_test_stencil.cc; range=3D_test_range.cc; c=3D_test_c.c; fortran=3D_fortran.f95
let "C = max_size**3 * max_steps"
size=8
while [ $size -le $max_size ]
do
        let "steps = C / size**3"
        #########################
        cd $path_fortran
        gfortran $fortran -cpp -Dsz=$size -Dst=$steps -O2 -Dcs=$dummy -o $fortran.bin
        ./$fortran.bin
        rm -r *.bin
        #########################
        cd $path_blitz
        mpic++ $stencil -O2 -Dsz=$size -Dst=$steps -Dcs=$dummy -ftemplate-depth-30 -o $stencil.bin
        ./$stencil.bin
        mpic++ $range -O2 -Dsz=$size -Dst=$steps -Dcs=$dummy -ftemplate-depth-30 -o $range.bin
        ./$range.bin
        rm -r *.bin
        cd $path_c
        mpicc -O3 -Dsz=$size -Dst=$steps -Dcs=$dummy -fstrict-aliasing -std=c99 -lm $c -o $c.bin
        ./$c.bin
        rm -r *.bin
        #########################
        let "size *= 2"
done
# End Running 3D tests
fi


