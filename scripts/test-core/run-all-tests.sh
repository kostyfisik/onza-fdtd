#!/bin/bash

path_FDTD_test=$PWD;   path_2D=$path_FDTD_test/2D;   path_3D=$path_FDTD_test/3D
path_2D_cpp=$path_2D/C++; path_2D_fortran=$path_2D/Fortran
path_2D_fortran_fast=$path_2D_fortran/fast_version
path_3D_cpp=$path_3D/C++; path_3D_fortran=$path_3D/Fortran
path_3D_fortran_fast=$path_3D_fortran/fast_version

# Running 2D tests
cd $path_2D_cpp
echo ">> Running 2D C++ tests"
chmod +x *.sh
./compile-and-run.sh 2D_test_blitz.cc
./compile-and-run.sh 2D_test_std.cc
echo ">> Done running 2D C++ tests"
cd $path_2D_fortran
echo ">> Running 2D Fortran tests"
cd $path_2D_fortran_fast
chmod +x *.sh
./compile-and-run-tests.sh
echo ">> Done Running 2D Fortran tests"
# End Running 2D tests

# Running 3D tests
cd $path_3D_cpp
echo ">> Running 3D C++ tests"
chmod +x *.sh
./compile-and-run.sh 3D_test_blitz.cc
./compile-and-run.sh 3D_test_std.cc
echo ">> Done running 3D C++ tests"
cd $path_3D_fortran
echo ">> Running 3D Fortran tests"
cd $path_3D_fortran_fast
chmod +x *.sh
./compile-and-run-tests.sh
echo ">> Done Running 3D Fortran tests"
# End Running 3D tests


