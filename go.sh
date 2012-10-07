#!/bin/bash
isNEW=$1
corepath=$PWD
MPIsize=1
MPInodes="-N 2"
#MPInodes="-N 1 --ntasks-per-socket=2"
MPIoptions=--bind-to-core
echo $corepath
find -L src  -name "*.cc" -o  -name "*.h" | xargs etags
rm -r $corepath/bin/* >/dev/null 2>&1
if [[ ! $isNEW ]]; then
    isNEW="new"
fi
#export OMPI_CXXFLAGS="-ftemplate-depth-30 -DBZ_DEBUG" # Debug mode
#export OMPI_CXXFLAGS="" # Debug mode
export OMPI_CXXFLAGS="-O2 -ftemplate-depth-30" # Benchmark mode
#export OMPI_CXXFLAGS="-O3 -ftemplate-depth-30" # Benchmark mode
#export OMPI_CXXFLAGS="-O3 -ffast-math -ftemplate-depth-30 -march=native -mtune=native -mfpmath=both -malign-double -mstackrealign -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5" # Benchmark mode. May bring no seed up, only errors.
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
# cd $corepath/build/gcc
# export OMPI_CC=gcc
# export OMPI_CXX=g++
# if [ $isNEW = "new" ]; then
#     rm -r $corepath/build/gcc/*
#     CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/gcc"
# else
#     cmake $corepath
# fi
# make -j4 
# make install
HOST=`cat /etc/hostname`
echo
echo "Executing clang version"
echo
cd $corepath/bin/clang
cp $corepath/data/default-onza.config $corepath/bin/clang/onza.config
if [ $HOST == "head.phoif.ifmo.ru" ]
then
    echo "Waiting for shared file system to distibute files"
    sleep 7
    # salloc $MPInodes -n $MPIsize -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # echo "(2) Nodes 16   procs 128"
    # salloc -N 16 -n 128 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # echo "(2) Nodes 16   procs 64"
    # salloc -N 16 -n 64 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # echo "(2) Nodes 16   procs 32"
    # salloc -N 16 -n 32 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # echo "(1) Nodes 16   procs 16"
    # salloc -N 16 -n 16 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # # echo "(2) Nodes 8   procs 16"
    # # salloc -N 8 -n 16 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # # echo
    # echo "(1) Nodes 8   procs 8"
    # salloc -N 8 -n 8 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # # echo "(2) Nodes 4   procs 8"
    # # salloc -N 4 -n 8 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # # echo
    # echo "(1) Nodes 4   procs 4"
    # salloc -N 4 -n 4 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # # echo "(2) Nodes 2   procs 4"
    # # salloc -N 2 -n 4 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # # echo
    # echo "(1) Nodes 2   procs 2"
    # salloc -N 2 -n 2 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # # echo "(8) Nodes 1   procs 8"
    # # salloc -N 1 -n 8 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # # echo
    # # echo "(4) Nodes 1   procs 4"
    # # salloc -N 1 -n 4 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # # echo
    # # echo "(2) Nodes 1   procs 2"
    # # salloc -N 1 -n 2 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # # echo
    echo "(1) Nodes 1   procs 1"
    salloc -N 1 -n 1 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    # mpirun -H n05 -n 1  $MPIoptions ./run-onza-fdtd onza.config
    # echo
    # echo "Executing gcc version"
    # echo
    # cd $corepath/bin/gcc
    # time salloc $MPInodes -n $MPIsize -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    echo
else
    mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd onza.config
    echo
    # echo "Executing gcc version"
    # echo
    # cd $corepath/bin/gcc
    # time mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd onza.config
    # echo
fi
echo "Prepare *.png from gnuplot for clang ..."
cd $corepath/bin/clang
cp $corepath/data/gnuplot/* ./
./gnuplot-all.sh >/dev/null  2>&1