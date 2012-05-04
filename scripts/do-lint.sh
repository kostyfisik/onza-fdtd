#!/bin/bash
FILES="../src/onza-fdtd.cc
../src/common.h
../src/mpi-decomposition/halo-exchange-process.h
../src/mpi-decomposition/halo-exchange-process.cc
../src/simulation-core/basic-fdtd.h
../src/simulation-core/basic-fdtd.cc"
for file in $FILES
do
 echo 
 ./cpplint.py $file
done