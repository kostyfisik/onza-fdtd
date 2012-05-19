#!/bin/bash
isNEW=$1
corepath=$PWD
MPIsize=16
# MPIoptions=--bind-to-core
echo $corepath
find -L src  -name "*.cc" -o  -name "*.h" | xargs etags
rm -r $corepath/bin/* >/dev/null 2>&1
if [[ ! $isNEW ]]; then
    isNEW="new"
fi
#export OMPI_CXXFLAGS=-O0 # Debug mode
export OMPI_CXXFLAGS=-O2 # Benchmark mode
#export OMPI_CXXFLAGS=-O3 # Benchmark mode
cd $corepath/build/clang
export OMPI_CC=clang
export OMPI_CXX=clang++
if [ $isNEW = "new" ]; then
    rm -r $corepath/build/clang/*
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/clang"
else
    cmake $corepath
fi
make -j4 
make install
cd $corepath/build/gcc
export OMPI_CC=gcc
export OMPI_CXX=g++
if [ $isNEW = "new" ]; then
    rm -r $corepath/build/gcc/*
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/gcc"
else
    cmake $corepath
fi
make -j4 
make install
cd $corepath/bin
cp clang/run-onza-fdtd run-onza-fdtd-clang
cp gcc/run-onza-fdtd run-onza-fdtd-gcc
HOST=`cat /etc/hostname`
if [ $HOST == "head.phoif.ifmo.ru" ]
then
    echo
    echo "Executing clang version"
    echo
    time salloc -n $MPIsize -p max1hour mpirun $MPIoptions ./run-onza-fdtd-clang
    echo
    echo "Executing gcc version"
    echo
    time salloc -n $MPIsize -p max1hour mpirun $MPIoptions ./run-onza-fdtd-gcc
    echo
else
    echo
    echo "Executing clang version"
    echo
    time mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd-clang
    echo
    echo "Executing gcc version"
    echo
    time mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd-gcc
    echo
fi
