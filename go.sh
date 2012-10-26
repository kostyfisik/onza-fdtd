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
USAGE="\nUsage: Onza can be (re-)compiled and run with call \n
$./go.sh [mode] [config]\n
e.g. from top-level Onza FDTD dir $./go.sh new\n
\n
Possible modes are:\n
\n
help        \t\t- show usage message.\n
<empty>, new \t- start totaly new build of Onza with Clang\n
             \t\t\x20 compiler. Use it for first time compilation.\n
new2         \t\t- same as new, but using gcc compiler.\n
new3         \t\t- same as new2, but with -O3 and -ffast-math flags.
old, old2    \t- recomile updated source files using Clang or gcc.\n
test         \t\t- compile with gcc benchmark flags, run a number of\n
             \t\t\x20 self tests, report about passing or failing them.\n
prof         \t\t- compile with gcc with gprof flag.\n
old2prof     \t- recompile with gprof. flag\n
custom       \t- use new mode and config file should be defined\n
debug        \t- use debug flag fro Blitz++
\n
Possible predefined config files:\n
\n
testX1Dzero      \t\t- 1D FDTD (used in test to check border exchange)\n
testTMz2Dspeedup \t- 2D FDTD (used in test to check speed up in 2D case)\n
test3Dsimple     \t\t- 3D FDTD\n"
#############################################################################
if [[ ! $MPIsize ]]; then MPIsize=2; fi
MPInodes="-N 2"
#MPInodes="-N 1 --ntasks-per-socket=2"
MPIoptions="--bind-to-core"
#############################################################################
#   Parse input parameters   
#############################################################################
mode=$1
configFile=$2
wrong_params=$3
if [[ $wrong_params ]]; then
    echo ================ !ERROR! =================
    echo Should be not more than two input parameters: mode and config file
    echo ================ !ERROR! =================
    echo -e  $USAGE
    exit 1
fi
# Mode names
mode_new="new"    # clang
mode_new2="new2"  # gcc
mode_new3="new3"  # gcc
mode_old="old"    
mode_old2="old2"
mode_test="test"
mode_prof="prof"
mode_old2prof="old2prof"
mode_custom="custom"
mode_debug="debug"
# Default values
compiler_gcc="gcc"
compiler_clang="clang"
yes="yes"
no="no"
isNew=$yes
isTest=$no
isProfile=$no
usedCompiler=$compiler_gcc # or clang
path_onza=$PWD
path_bin=$path_onza/bin
path_build=$path_onza/build
if [[ ! $mode ]]; then
    mode="new"
fi
if [[ $mode = "help" || \
    $mode = "--help" || \
    $mode = "-h" ]]; then
 echo -e $USAGE
 exit 0
fi 
# Check mode
if [[ $mode != $mode_new && \
    $mode != $mode_new2 && \
    $mode != $mode_new3 && \
    $mode != $mode_old && \
    $mode != $mode_old2 && \
    $mode != $mode_test && \
    $mode != $mode_prof && \
    $mode != $mode_old2prof && \
    $mode != $mode_custom && \
    $mode != $mode_debug ]]; then
    echo Using default mode: Full build with clang compiler.
    if [[ $configFile ]]; then
        echo ================ !ERROR! =================
        echo Undefined mode
        echo ================ !ERROR! =================
        echo -e  $USAGE
        exit 1
    fi  
    configFile=$mode
    mode=$mode_new
fi 
if [[ $mode = $mode_test && $configFile ]]; then
    echo ================ !ERROR! =================
    echo Test mode do not support external config file
    echo ================ !ERROR! =================
    exit 1
fi
# Check config file(s)
path_testX1Dzero=$path_onza/data/testX1Dzero.config
path_testTMz2Dspeedup=$path_onza/data/testTMz2Dspeedup.config
path_test3Dsimple=$path_onza/data/test3Dsimple.config
path_default_config=$path_onza/data/default-onza.config
if [[ $mode = $mode_test ]]; then
    echo Check for tests config files...
    if [[ ! -r $path_testX1Dzero || \
        ! -r $path_testTMz2Dspeedup || \
        ! -r $path_test3Dsimple ]];
    then
        echo ================ !ERROR! =================
        echo Test mode was not able to access some of config files for tests
        echo $path_testX1Dzero
        echo $path_testTMz2Dspeedup
        echo $path_test3Dsimple
        echo ================ !ERROR! =================
        exit 1
    fi
fi
if [ ! $configFile ]; then
    echo Setting default config file $path_default_config
    configFile=$path_default_config
fi
if  [[ ! -a $configFile ]]; then
    echo ================ !ERROR! =================
    echo Was not able to found config file using path: $configFile
    echo ================ !ERROR! =================
    exit 1
fi
# Convert relative path for custom config file to absolute.
firstChar=${configFile:0:1}
if [[ ! $firstChar = "/" ]]; then
    configFile=$path_onza/$configFile
    echo Change config file path to absolute path: $configFile
fi
echo ============ Current script settings =============
echo mode: $mode
echo config file: $configFile
echo base dir path: $path_onza
# Check directory structure
path_src=$path_onza/src
if [[ ! -d $path_src ]]; then
    echo ================ !ERROR! =================
    echo No source folder $path_src
    echo ================ !ERROR! =================
    exit
fi
if [[ -a $path_build ]]; then
    echo Found build folder.
else
    echo Creating build folder...
    mkdir build
fi
if [[ -a $path_bin ]]; then
    echo Found bin folder.
else
    echo Creating build folder...
    mkdir bin
fi

#############################################################################
#   Compile settings
#############################################################################
echo ============ Compile settings =============
if [[ $mode = $mode_old || $mode = $mode_old2 || $mode = $mode_old2prof ]]; then
    isNew=$no
    echo Recompile mode is on.
else
    echo Cleaning build path.
    rm -r $path_build/* >/dev/null 2>&1
fi 
echo Cleaning bin
rm -r $path_bin/* >/dev/null 2>&1
# Profiling mode should use gcc compiler.
flag_cmake_profile=
if [[ $mode = $mode_prof || $mode = $mode_old2prof ]]; then
    echo Using gprof.
    isProfile=$yes
    flag_profile="-pg";
    flag_cmake_profile="-DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg"
fi 
# Set compiler
if [[ $mode = $mode_new || $mode = $mode_old ]]; then
    echo Using \'clang\' compiler.
    usedCompiler=$compiler_clang
    export OMPI_CC=clang
    export OMPI_CXX=clang++
else
    echo Using \'gcc\' compiler.
    export OMPI_CC=gcc
    export OMPI_CXX=g++
fi 
# Select OMPI_CXXFLAGS
flags_O2="-O2 -ftemplate-depth-30"
flags_debug="-ftemplate-depth-30 -DBZ_DEBUG"
flags_O3="-O3 -ffast-math -ftemplate-depth-30 -march=native -mtune=native -mfpmath=both -malign-double -mstackrealign -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5"
export OMPI_CXXFLAGS=$flags_O2
if [[ $mode = $mode_debug ]]; then
    echo Using debug mode for Biltz++.
    export OMPI_CXXFLAGS=$flags_debug
fi
if [[ $mode = $mode_new3 ]]; then
    export OMPI_CXXFLAGS=$flags_O3
fi
#############################################################################
#   Build
#############################################################################
cd $path_build
if [[ $isNew = $yes ]]; then
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $path_onza -DCMAKE_INSTALL_PREFIX="$path_bin" $flag_cmake_profile
fi
# if [[ $isNew = $no ]]; then
#     cmake $path_onza
# fi
# if [[ $isNew = "new2" ]]; then
#     echo "Compile with gprof"
#     ## TODO See
#     ## http://www.open-mpi.org/community/lists/users/2009/04/9039.php
#     ## to use gprof with mpi
#     ## setenv GMON_OUT_PREFIX gout
#     CC=mpicc CXX=mpic++ VERBOSE=1 cmake $path_onza -DCMAKE_INSTALL_PREFIX="$path_bin" -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg
# fi
make -j4 
make install
#############################################################################
#   Run
#############################################################################
HOST=`cat /etc/hostname`
if [[ $isProfile = $yes ]]; then 
    # Grpof output for each process in separate file
    export GMON_OUT_PREFIX='gmon.out'
fi
if [[ $isTest = $no ]]; then
    cd $path_bin
    cp $configFile $path_bin/onza.config
    if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
        echo "Waiting for shared file system to distibute files"
        sleep 6
        salloc -N 8 -n 16 -p max1hour mpirun $MPIoptions ./run-onza-fdtd onza.config
    elif  [[ $HOST == "deb00" || $HOST == "dmmrkovich-birzha" ]]; then
        echo "(1) Nodes XX  procs XX"
        mpirun -np 4 $MPIoptions ./run-onza-fdtd onza.config
    else
        mpirun -np $MPIsize $MPIoptions ./run-onza-fdtd onza.config
    fi
fi  # end of if [[ $isTest = $no ]]

if [[ $isProfile = $yes ]]; then 
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
    gprof --no-flat-profile run-onza-fdtd gmon.out* > average-profile
    rm gmon.out*
fi  # end of if isProfile
rm *.onza
# echo "Prepare *.png from gnuplot ..."
# cd $path_bin
# cp $path_onza/data/gnuplot/* ./
# ./gnuplot-all.sh >/dev/null  2>&1
    