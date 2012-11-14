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
<empty>, new \t- start totaly new build of Onza with Clang.\n
             \t\t\x20 compiler. Use it for first time compilation.\n
new2         \t\t- same as new, but using gcc compiler.\n
new3         \t\t- same as new2, but with -O3 and -ffast-math flags.\n
old, old2    \t- recomile updated source files using Clang or gcc.\n
test         \t\t- compile with gcc benchmark flags, run a number of\n
             \t\t\x20 self tests, report about passing or failing them.\n
prof         \t\t- compile with gcc with gprof flag.\n
old2prof     \t- recompile with gprof flag.\n
custom       \t- use 'new' mode and config file should be defined.\n
build        \t\t- just build with gcc compiler.\n
debug        \t\t- use debug flag fro Blitz++.
\n\n
Possible predefined config files:\n
\n
testX1Dzero      \t\t- 1D FDTD (used in test to check border exchange)\n
testTMz2Dspeedup \t- 2D FDTD (used in test to check speed up in 2D case)\n
test3Dsimple     \t\t- 3D FDTD\n
\n
Possible enviroment parameters:\n
\n
Onza_MPI_size   \t- total number of MPI processes
Onza_MPI_nodes  \t- total number of MPI nodes for cluster enviroment
\n
In the case of enviroment parameters are not (or set with value 'unset')\n
set some default vale (depending on executing host) will be used.\n"
#############################################################################
if [[ ! $Onza_MPI_size ]]; then Onza_MPI_size="unset"; fi
if [[ ! $Onza_MPI_nodes ]]; then Onza_MPI_nodes="unset"; fi
MPIoptions="--bind-to-core"
#############################################################################
#   Parse input parameters   
#############################################################################
mode=$1; configFile=$2; wrong_params=$3
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
mode_new3="new3"  # gcc with -O3
mode_old="old";       mode_old2="old2";         mode_test="test"; 
mode_prof="prof";     mode_old2prof="old2prof"
mode_custom="custom"; mode_debug="debug";       mode_build="build"
# Default values
compiler_gcc="gcc"; compiler_clang="clang"
usedCompiler=$compiler_gcc # or clang
yes="yes";        no="no"
isNew=$yes;       isTest=$no ;               isProfile=$no
path_onza=$PWD;   path_bin=$path_onza/bin;   path_build=$path_onza/build
path_src=$path_onza/src
# Should be same as in cmake config in $path_src
onza_bin="run-onza-fdtd"
if [[ ! $mode ]]; then  mode="new"; fi
if [[ $mode = "help" || $mode = "--help" || $mode = "-h" ]]; then
 echo -e $USAGE
 exit 0
fi 
# Check mode
if [[ $mode != $mode_new && $mode != $mode_new2 && $mode != $mode_new3 && \
    $mode != $mode_old && $mode != $mode_old2 && \
    $mode != $mode_test && \
    $mode != $mode_prof && $mode != $mode_old2prof && \
    $mode != $mode_custom && \
    $mode != $mode_build && \
    $mode != $mode_debug ]]; then
    # So mode may be miss spelled or contains config file path.
    if [[ $configFile ]]; then
        echo ================ !ERROR! =================
        echo Undefined mode
        echo ================ !ERROR! =================
        echo -e  $USAGE
        exit 1
    fi  
    echo Using default mode: Full build with clang compiler.
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
tests="testX1Dzero testTMz2Dspeedup test3Dsimple"
for test in $tests; do
    path_test_config="path_${test}_config"
    eval path_${test}_config=$path_onza/data/$test.config  
    echo ${!path_test_config}
done
# path_testX1Dzero_config=$path_onza/data/testX1Dzero.config
# path_testTMz2Dspeedup_config=$path_onza/data/testTMz2Dspeedup.config
# path_test3Dsimple_config=$path_onza/data/test3Dsimple.config
path_default_config=$path_onza/data/default-onza.config
if [[ $mode = $mode_test ]]; then
    echo Check for tests config files...
    for test in $tests; do
        path_test_config="path_${test}_config"
        eval path_${test}_config=$path_onza/data/$test.config  
        echo ${!path_test_config}
        if [[ ! -r ${!path_test_config} ]];
        then
            echo ================ !ERROR! =================
            echo Test mode was not able to access some of config files for tests
            echo ${!path_test_config}
            echo ================ !ERROR! =================
            exit 1
        fi
    done
fi
# Should be tested after checking $test_mode 
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
if [[ ! -d $path_src ]]; then
    echo ================ !ERROR! =================
    echo No source folder $path_src
    echo ================ !ERROR! =================
    exit
fi
force_new_build=$no
if [[ -a $path_build ]]; then
    echo Found build folder.
else
    echo Creating build folder...
    mkdir build
    force_new_build=$yes
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
if [[ ( $mode = $mode_old || $mode = $mode_old2 || $mode = $mode_old2prof || \
    $mode = $mode_test ) && $force_new_build = $no ]]; then
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
    CC=mpicc CXX=mpic++ VERBOSE=1 cmake $path_onza \
        -DCMAKE_INSTALL_PREFIX="$path_bin" $flag_cmake_profile
fi
make -j4 
make install
#############################################################################
#   Run
#############################################################################
if [[ $mode = $mode_build ]]; then
    cd $path_bin
    mv run-onza-fdtd onza-fdtd.bin
    exit 0; 
fi
HOST=`cat /etc/hostname`
echo "Executing on host -> $HOST <-"
if [[ $mode = $mode_test ]]; then isTest=$yes; fi
if [[ $isProfile = $yes ]]; then 
    # Grpof output for each process in separate file
    export GMON_OUT_PREFIX='gmon.out'
fi
if [[ $isTest = $no ]]; then
    cd $path_bin
    cp $configFile $path_bin/onza.config
    if [[ $Onza_MPI_size = "unset" && $configFile = $path_testX1Dzero_config ]]; then
        Onza_MPI_size=2
        Onza_MPI_nodes=1
    fi
    if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
        echo "Waiting for shared file system to distibute files..."
        sleep 2
        if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=16; fi
        if [[ $Onza_MPI_nodes = "unset" ]]; then Onza_MPI_nodes=8; fi
        echo "(1) Nodes $Onza_MPI_nodes procs $Onza_MPI_size"
        salloc -N $Onza_MPI_nodes -n $Onza_MPI_size -p max1hour \
            mpirun $MPIoptions ./$onza_bin onza.config
    elif  [[ $HOST == "deb00" || $HOST == "dmmrkovich-birzha" ]]; then
        if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=4; fi            
        if [[ $Onza_MPI_size = "1-4" ]]; then
            Onza_MPI_size=4
            echo "(1) Nodes 1  procs $Onza_MPI_size"
            mpirun -np $Onza_MPI_size $MPIoptions ./$onza_bin onza.config
            Onza_MPI_size=2
            echo "(1) Nodes 1  procs $Onza_MPI_size"
            mpirun -np $Onza_MPI_size $MPIoptions ./$onza_bin onza.config
            Onza_MPI_size=1
            echo "(1) Nodes 1  procs $Onza_MPI_size"
            mpirun -np $Onza_MPI_size $MPIoptions ./$onza_bin onza.config
        else
            echo "(1) Nodes 1  procs $Onza_MPI_size"
            mpirun -np $Onza_MPI_size $MPIoptions ./$onza_bin onza.config
        fi
    else
        if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=2; fi
        echo "(1) Nodes 1  procs $Onza_MPI_size"
        mpirun -np $Onza_MPI_size $MPIoptions ./$onza_bin onza.config
    fi
fi  # end of if [[ $isTest = $no ]]
if [[ $isProfile = $yes ]]; then 
    gprof  $onza_bin gmon.out* | grep '^index'
    for file in gmon.out*;
    do
        echo $file
        gprof  $onza_bin $file | grep 'DoBorderStep' | grep '^\['
        gprof  $onza_bin $file | grep 'DoStep' | grep '^\['
    done
    echo Average
    gprof  $onza_bin gmon.out* | grep 'DoBorderStep' | grep '^\['
    gprof  $onza_bin gmon.out* | grep 'DoStep' | grep '^\['
    gprof --no-flat-profile $onza_bin gmon.out* > average-profile
    gprof $onza_bin gmon.out* > flat-average-profile
    rm gmon.out*
fi  # end of if isProfile
if [[ $isTest = $yes ]]; then
    echo "Prepare files for tests ..."
    cd $path_bin
    for test in $tests; do
        # For each test make dir, copy config and binary to it, run.
        cd $path_bin
        mkdir $test
        path_test="$path_bin/$test"  
        path_test_config="path_${test}_config"
        cp ${!path_test_config} ${path_test}/onza.config
        cp $onza_bin $path_test
        cd $path_test
        backup_Onza_MPI_size=$Onza_MPI_size
        backup_Onza_MPI_nodes=$Onza_MPI_nodes
        if [[ $test = "testX1Dzero" ]]; then
            Onza_MPI_size=2
            Onza_MPI_nodes=1
        fi
        echo "============ Running test $test ============="
        if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
            echo "Waiting for shared file system to distibute files..."
            sleep 2
            if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=16; fi
            if [[ $Onza_MPI_nodes = "unset" ]]; then Onza_MPI_nodes=8; fi
            echo "(1) Nodes $Onza_MPI_nodes procs $Onza_MPI_size"
            salloc -N $Onza_MPI_nodes -n $Onza_MPI_size -p max1hour \
                mpirun $MPIoptions ./$onza_bin onza.config
            if [[ $test = "testTMz2Dspeedup" ]]; then
                echo "(*******) Nodes 16 procs 128 (1024 x 1024, step 1000, ~15.4s)"
                salloc -N 16 -n 128 -p max1hour \
                    mpirun $MPIoptions ./$onza_bin onza.config
                echo "(*******) Nodes 16 procs 16 (1024 x 1024, step 1000, ~7s)"
                salloc -N 16 -n 16 -p max1hour \
                    mpirun $MPIoptions ./$onza_bin onza.config
                echo "(*******) Nodes 1 procs 8 (1024 x 1024, step 1000, ~30.2s)"
                salloc -N 1 -n 8 -p max1hour \
                    mpirun $MPIoptions ./$onza_bin onza.config
                echo "(*******) Nodes 1 procs 1 (1024 x 1024, step 1000, ~47.6s)"
                salloc -N 1 -n 1 -p max1hour \
                    mpirun $MPIoptions ./$onza_bin onza.config
            fi

        elif  [[ $HOST == "deb00" || $HOST == "dmmrkovich-birzha" ]]; then
            echo "(*******) Procs 1"
            mpirun -np 1 $MPIoptions ./$onza_bin onza.config
            echo "(*******) Procs 2"
            mpirun -np 2 $MPIoptions ./$onza_bin onza.config
            echo "(*******) Procs 4"
            mpirun -np 4 $MPIoptions ./$onza_bin onza.config
        else
            if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=2; fi
            echo "(1) Nodes 1  procs $Onza_MPI_size"
            mpirun -np $Onza_MPI_size $MPIoptions ./$onza_bin onza.config
        fi
        if [[ $test = "testX1Dzero" ]]; then
            echo "Prepare *.png from gnuplot ..."
            cp $path_onza/data/gnuplot/* ./
            ./gnuplot-all.sh >/dev/null  2>&1
            cp *.png $path_bin/
            # cp *02??-*.png $path_bin/
            rm *.png
        fi
        Onza_MPI_size=$backup_Onza_MPI_size
        Onza_MPI_nodes=$backup_Onza_MPI_nodes
        rm *.onza  >/dev/null  2>&1
    done  # end of for test in $tests; do
fi  # end of if [[ $isTest = $yes ]]
if [[ $configFile = $path_testX1Dzero_config ]]; then
    echo "Prepare *.png from gnuplot ..."
    ls *
    cp $path_onza/data/gnuplot/* ./
    ./gnuplot-all.sh >/dev/null  2>&1
    # mkdir tmpdir
    # cp *0241-* $path_bin/tmpdir
    # rm $path_bin/*    
    # cp $path_bin/tmpdir/* ./
    # rm -r $path_bin/tmpdir
fi
#rm *.onza  >/dev/null  2>&1
    