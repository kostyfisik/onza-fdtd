#!/bin/bash
NEW=$1
rm -r bin/* >/dev/null 2>&1
corepath=$PWD
echo $corepath
cd build/clang
if [ $NEW = "new" ]; then
    rm -r *
    CC=clang CXX=clang++ VERBOSE=1 cmake ../.. -DCMAKE_INSTALL_PREFIX="$corepath/bin/clang"
else
    cmake ../..
fi
make -j4 
make install
cd ../gcc
if [ $NEW = "new" ]; then
    rm -r *
    CC=gcc CXX=g++ VERBOSE=1 cmake ../.. -DCMAKE_INSTALL_PREFIX="$corepath/bin/gcc"
else
    cmake ../..
fi
make -j4 
make install
cd ../../bin
cp clang/run-onza-fdtd run-onza-fdtd-clang
cp gcc/run-onza-fdtd run-onza-fdtd-gcc
echo "Executing clang version"
./run-onza-fdtd-clang
echo "Executing gcc version"
./run-onza-fdtd-gcc
