#!/bin/bash
isNEW=$1
corepath=$PWD
echo $corepath
rm -r $corepath/bin/* >/dev/null 2>&1
if [[ ! $isNEW ]]; then
    isNEW="new"
fi
cd $corepath/build/clang
if [ $isNEW = "new" ]; then
    rm -r $corepath/build/clang/*
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/clang"
else
    cmake $corepath
fi
export OMPI_CC=clang
export OMPI_CXX=clang++
make -j4 
make install
cd $corepath/build/gcc
if [ $isNEW = "new" ]; then
    rm -r $corepath/build/gcc/*
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/gcc"
else
    cmake $corepath
fi
export OMPI_CC=gcc
export OMPI_CXX=g++
make -j4 
make install
cd $corepath/bin
cp clang/run-onza-fdtd run-onza-fdtd-clang
cp gcc/run-onza-fdtd run-onza-fdtd-gcc
echo
echo "Executing clang version"
echo
mpirun -np 2 ./run-onza-fdtd-clang
echo
echo "Executing gcc version"
echo
mpirun -np 2 ./run-onza-fdtd-gcc
