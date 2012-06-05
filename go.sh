#!/bin/bash
isNEW=$1
corepath=$PWD
MPIsize=2
#MPInodes="-N 16"
# MPIoptions=--bind-to-core
echo $corepath
find -L src  -name "*.cc" -o  -name "*.h" | xargs etags
rm -r $corepath/bin/* >/dev/null 2>&1
if [[ ! $isNEW ]]; then
    isNEW="new"
fi
#export OMPI_CXXFLAGS="-ftemplate-depth-30 -DBZ_DEBUG" # Debug mode
#export OMPI_CXXFLAGS="-O2 -ftemplate-depth-30" # Benchmark mode
#export OMPI_CXXFLAGS="-O3 -ftemplate-depth-30" # Benchmark mode
export OMPI_CXXFLAGS="-O3 -ffast-math -ftemplate-depth-30 -march=native -mtune=native -mfpmath=both -malign-double -mstackrealign" # Benchmark mode. May bring no seed up, only errors.
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
HOST=`cat /etc/hostname`
if [ $HOST == "head.phoif.ifmo.ru" ]
then
    echo
    echo "Executing clang version"
    echo
    cd $corepath/bin/clang
    time salloc $MPInodes -n $MPIsize -p max1hour mpirun $MPIoptions ./run-onza-fdtd
    echo
    echo "Executing gcc version"
    echo
    cd $corepath/bin/gcc
    time salloc $MPInodes -n $MPIsize -p max1hour mpirun $MPIoptions ./run-onza-fdtd
    echo
else
    echo
    echo "Executing clang version"
    echo
    cd $corepath/bin/clang
    time mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd
    echo
    echo "Executing gcc version"
    echo
    cd $corepath/bin/gcc
    time mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd
    echo
fi
echo "Prepare *.png from gnuplot for clang ..."
cd $corepath/bin/clang
cp $corepath/data/gnuplot/* ./
time ./gnuplot-all.sh