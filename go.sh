#!/bin/bash
# @file   go.sh
# @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
# @copyright 2012 Ladutenko Konstantin
# @section LICENSE
# This file is part of Onza FDTD.
#
# Onza FDTD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Onza FDTD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Onza FDTD.  If not, see <http://www.gnu.org/licenses/>.
#
# @brief  Compile Onza FDTD and run simulation.
#
# Usage:
#
isNew=$1
isTest="no"
corepath=$PWD
bin_path=corepath
MPIsize=1
MPInodes="-N 2"
#MPInodes="-N 1 --ntasks-per-socket=2"
MPIoptions=--bind-to-core
echo $corepath
find -L src  -name "*.cc" -o  -name "*.h" | xargs etags
rm -r $corepath/bin/* >/dev/null 2>&1
if [[ ! $isNew ]]; then
    isNew="new"
fi
if [ $isNew = "test" ]; then
    isNew="new2"
    isTest="yes"
fi
#export OMPI_CXXFLAGS="-ftemplate-depth-30 -DBZ_DEBUG" # Debug mode
#export OMPI_CXXFLAGS="" # Debug mode
export OMPI_CXXFLAGS="-O2 -ftemplate-depth-30" # Benchmark mode
#export OMPI_CXXFLAGS="-O3 -ftemplate-depth-30" # Benchmark mode
#export OMPI_CXXFLAGS="-O3 -ffast-math -ftemplate-depth-30 -march=native -mtune=native -mfpmath=both -malign-double -mstackrealign -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5" # Benchmark mode. May bring no speed up, only errors.
#export OMPI_LDFLAG="-pg"
cd $corepath/build/clang
export OMPI_CC=clang
export OMPI_CXX=clang++
if [ $isNew = "new" ]; then
    rm -r $corepath/build/clang/*
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/clang"
    make -j4 
    make install
    bin_path=$corepath/bin/clang
fi
if [ $isNew = "old" ]; then
    cmake $corepath
    make -j4 
    make install
    bin_path=$corepath/bin/clang
fi
cd $corepath/build/gcc
export OMPI_CC=gcc
export OMPI_CXX=g++
if [ $isNew = "new2" ]; then
    rm -r $corepath/build/gcc/*
    # export CXXFLAGS='-pg'
    # export CXXLDFLAGS='-pg'
    #CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/gcc"
    echo "Compile with gprof"
    ## TODO See
    ## http://www.open-mpi.org/community/lists/users/2009/04/9039.php
    ## to use gprof with mpi
    ## setenv GMON_OUT_PREFIX gout
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $corepath -DCMAKE_INSTALL_PREFIX="$corepath/bin/gcc" -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg
    make -j4 
    make install
    bin_path=$corepath/bin/gcc
fi
if [ $isNew = "old2" ]; then
    cmake $corepath
    make -j4 
    make install
    bin_path=$corepath/bin/gcc
fi

HOST=`cat /etc/hostname`
if [ $isTest = "no" ]; then
    cd $bin_path
    cp $corepath/data/default-onza.config $bin_path/onza.config
    if [ $HOST == "head.phoif.ifmo.ru" ]
    then
        echo "Waiting for shared file system to distibute files"
        sleep 2
        # Grpof output for each process in separate file
        export GMON_OUT_PREFIX='gmon.out'
        echo "(1) Nodes XX  procs XX"
        salloc -N 8 -n 16 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
        gprof --no-flat-profile run-onza-fdtd gmon.out* | grep '^index'
        for file in gmon.out*;
        do
            echo $file
            gprof --no-flat-profile run-onza-fdtd $file | grep 'DoBorderStep' | grep '^\['
            gprof --no-flat-profile run-onza-fdtd $file | grep 'DoStep' | grep '^\['
        done
        echo Average
        gprof --no-flat-profile run-onza-fdtd gmon.out* | grep 'DoBorderStep' | grep '^\['
        gprof --no-flat-profile run-onza-fdtd gmon.out* | grep 'DoStep' | grep '^\['
        rm gmon.out*
        echo
    else
        mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd onza.config
        echo
    fi
    rm *.onza
    # echo "Prepare *.png from gnuplot ..."
    # cd $bin_path
    # cp $corepath/data/gnuplot/* ./
    # ./gnuplot-all.sh >/dev/null  2>&1
fi