#!/bin/bash
# @file   go.sh
# @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
# @copyright 2012 Konstantin Ladutenko
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
usage_msg="\nUsage: Onza can be (re-)compiled and run with call \n
$./go.sh [mode] [config]\n
e.g. from top-level Onza FDTD dir $./go.sh new1\n
\n
Possible modes (build + run) are:\n
\n
help        \t\t- show usage message.\n
<empty>, new1 \t- start totaly new build of Onza with Clang.\n
             \t\t\x20 compiler. Use it for first time compilation.\n
new2         \t\t- same as new1, but using gcc compiler.\n
new3         \t\t- same as new2, but with -O3 -ffast-math -flto\n
             \t\t  flags (needs gcc > 4.5). You can run build-gcc.sh in\n
             \t\t  onza-fdtd/scripts/build-additional-soft to build gcc compiler\n
old1, old2    \t- recomile updated source files using Clang or gcc.\n
test         \t\t- compile with gcc benchmark flags, run a number of\n
             \t\t\x20 self tests, report about passing or failing them.\n
prof         \t\t- compile with gcc with gprof flag.\n
old2prof     \t- recompile with gprof flag.\n
pgo          \t\t- compile with gcc with in pgo mode.\n
custom       \t- use 'new1' mode and config file should be defined.\n
build        \t\t- just build with gcc compiler (no run).\n
debug        \t\t- use debug flag fro Blitz++.
\n\n
Possible config files:\n
\n
Any valid file. If it can be treated as config will be checked with onza\n
or\n
build - only to build with predefined mode (no run).
\n
Possible enviroment parameters:\n
\n
Onza_MPI_size   \t- total number of MPI processes.\n
Onza_MPI_nodes  \t- total number of MPI nodes for cluster enviroment.\n
\n
In the case of enviroment parameters are not (or set with value 'unset')\n
set some default vale (depending on executing host) will be used.\n"
#############################################################################
if [[ ! $Onza_MPI_size ]]; then Onza_MPI_size="unset"; fi
if [[ ! $Onza_MPI_nodes ]]; then Onza_MPI_nodes="unset"; fi
MPI_options="--bind-to-core"
HOST=`hostname`
#############################################################################
#   Parse input parameters   
#############################################################################
mode=$1; config_file=$2; wrong_params=$3
if [[ $wrong_params ]]; then
    echo ================ !ERROR! =================
    echo Should be not more than two input parameters: mode and config file
    echo ================ !ERROR! =================
    echo -e  $usage_msg
    exit 1
fi
# Mode names
mode_new1="new1"    # clang
mode_new2="new2"  # gcc
mode_new3="new3"  # gcc with -O3
mode_old1="old1";     mode_old2="old2";         mode_test="test"; 
mode_prof="prof";     mode_old2prof="old2prof"; mode_pgo="pgo"
mode_custom="custom"; mode_debug="debug";       mode_build="build"
# Default values
yes="yes";        no="no"
compiler_gcc="gcc"; compiler_clang="clang"
usedCompiler=$compiler_gcc # or clang
useGCC47=$yes  # use gcc 4.7 if it is available in build area of scripts folder
isNew=$yes;       isTest=$no ;  isProfile=$no ; isPGO=$no
isBuildOnly=$no;
if [[ $mode = $mode_build || $config_file = $mode_build ]]; then
    isBuildOnly=$yes    
fi
path_onza=$PWD;   path_bin=$path_onza/bin;   path_build=$path_onza/build
path_src=$path_onza/src
# Should be same as in cmake config in $path_src
onza_bin="run-onza-fdtd"
if [[ ! $mode ]]; then  mode="new1"; fi
if [[ $mode = "help" || $mode = "--help" || $mode = "-h" ]]; then
 echo -e $usage_msg
 exit 0
fi 
# Check mode
if [[ $mode != $mode_new1 && $mode != $mode_new2 && $mode != $mode_new3 && \
    $mode != $mode_old1 && $mode != $mode_old2 && \
    $mode != $mode_test && \
    $mode != $mode_prof && $mode != $mode_old2prof && \
    $mode != $mode_pgo && \
    $mode != $mode_custom && \
    $mode != $mode_build && \
    $mode != $mode_debug ]]; then
    # So mode may be miss spelled or contains config file path.
    if [[ $config_file ]]; then
        echo ================ !ERROR! =================
        echo Undefined mode
        echo ================ !ERROR! =================
        echo -e  $usage_msg
        exit 1
    fi  
    echo Using default mode: Full build with clang compiler.
    config_file=$mode 
    mode=$mode_new1
fi 
if [[ ( $mode = $mode_test || $mode = $mode_build ) && $config_file ]]; then
    echo ================ !ERROR! =================
    echo Test and build modes do not support external config file
    echo ================ !ERROR! =================
    exit 1
fi
# Check config file(s)
path_default_config=$path_onza/data/default-onza.config
# Test containing "X1D-zero" in its name is a special optin in TuneOnzaOptionsMPI
tests="self-test-X1D-zero self-test-TMz2D-speedup self-test-3D-simple"
if [[ $mode = $mode_test ]]; then
    echo Check for tests config files...
    for test in $tests; do
        test_name=${test//-/_}
        path_test_config="path_${test_name}_config"
        eval path_${test_name}_config=$path_onza/data/$test.config  
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
if [[ ! $config_file && $mode != $mode_build ]]; then
    echo Setting default config file $path_default_config
    config_file=$path_default_config
fi
if  [[ ! -a $config_file && $isBuildOnly != $yes ]]; then
    echo ================ !ERROR! =================
    echo Was not able to found config file using path: $config_file
    echo ================ !ERROR! =================
    exit 1
fi
# Convert relative path for custom config file to absolute.
firstChar=${config_file:0:1}
if [[ $firstChar != "/" && $isBuildOnly != $yes ]]; then
    config_file=$path_onza/$config_file
    echo Change config file path to absolute path: $config_file
fi
echo ============ Current script settings =============
echo mode: $mode
if [[ $mode != $mode_test && $isBuildOnly != $yes ]]; then
    echo config file: $config_file
fi
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
if [[ ( $mode = $mode_old1 || $mode = $mode_old2 || $mode = $mode_old2prof || \
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
if [[ $mode = $mode_prof || $mode = $mode_old2prof ]]; then
    echo Using gprof.
    isProfile=$yes
    flag_cmake_profile="-DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg"
fi 
# PGO mode should use gcc compiler.
if [[ $mode = $mode_pgo ]]; then
    echo Using pgo.
    isPGO=$yes
fi 
# Set compiler
if [[ $mode = $mode_new1 || $mode = $mode_old1 ]]; then
    echo Using \'clang\' compiler.
    usedCompiler=$compiler_clang
    export OMPI_CC=clang
    export OMPI_CXX=clang++
    export CC=clang
    export CXX=clang++	
else
    path_gcc47=$path_onza/scripts/build-additional-soft/gcc-4.7/output/bin/
    if [[ -a $path_gcc47/gcc-4.7 && $useGCC47 = $yes ]]; then
        echo Using \'gcc-4.7.2\' compiler.
        export OMPI_CC=$path_gcc47/gcc-4.7
        export OMPI_CXX=$path_gcc47/g++-4.7
    else
        echo Using gcc compiler.
        path_gcc47=
        export OMPI_CC=gcc
        export OMPI_CXX=g++
    fi
fi 
if  [[ $HOST == "rh-lum.metalab.ifmo.ru" ]]; then
    echo Setting MPI path on rh-lum.metalab.ifmo.ru !
    ompi_path_bin=/usr/lib64/openmpi/bin/
    ompi_path_lib=/usr/lib64/openmpi/lib/
    if [ -d "$ompi_path_bin" ] && [[ ":$PATH:" != *":$ompi_path_bin:"* ]]; then
        PATH="${PATH:+"$PATH:"}$ompi_path_bin"
    fi
    export LD_LIBRARY_PATH=$ompi_path_lib
fi

# Select OMPI_CXXFLAGS
flags_O2="-O2 -Wall -std=c++11 -DNDEBUG"
flags_debug="-DBZ_DEBUG -Wall -std=c++11 "
flags_O3="-O3 -ffast-math -march=native -mtune=native -mfpmath=both -malign-double -mstackrealign -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5  -Wall -std=c++11  -DNDEBUG" 
# TODO option -flto   -- Do we need it?
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
function BuildOnza {
    local path_tmp=`pwd`
    cd $path_build
    if [[ $flag_cmake_profile && $flag_cmake_pgo ]]; then
        echo ================ !ERROR! =================
        echo Config with cmake flags $flag_cmake_profile and $flag_cmake_pgo
        echo ================ !ERROR! =================
        exit        
    fi
    flag_cmake=
    if [[ $flag_cmake_profile ]]; then
        flag_cmake=$flag_cmake_profile
    fi
    if [[ $flag_cmake_pgo ]]; then
        flag_cmake=$flag_cmake_pgo
    fi
    if [[ $isNew = $yes ]]; then
        CC=mpicc CXX=mpic++ VERBOSE=1 cmake $path_onza \
            -DCMAKE_INSTALL_PREFIX="$path_bin" $flag_cmake
    fi
    make -j4 
    make install
    cd $path_tmp
}  # end of function BuildOnza
BuildOnza
#############################################################################
#   Run
#############################################################################
if [[ $isBuildOnly = $yes ]]; then
    cd $path_bin
    mv run-onza-fdtd onza-fdtd.bin
    exit 0; 
fi
echo "Executing on host -> $HOST <-"
if [[ $mode = $mode_test ]]; then isTest=$yes; fi
if [[ $isProfile = $yes ]]; then 
    # Grpof output for each process in separate file
    export GMON_OUT_PREFIX='gmon.out'
fi
# Task and host dependant MPI parameters tuning.
function TuneOnzaOptionsMPI {
    if [[ $Onza_MPI_size = "unset" && \
         "$config_file" = *"X1D-zero"* ]]; then
        Onza_MPI_size=2
        Onza_MPI_nodes=1
    fi
    if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
        if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=16; fi
        if [[ $Onza_MPI_nodes = "unset" ]]; then Onza_MPI_nodes=8; fi
    elif  [[ $HOST == "deb00" || $HOST == "dmmrkovich-birzha" ]]; then
        if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=4; fi            
    elif  [[ $HOST == "rh-lum.metalab.ifmo.ru" ]]; then
        if [[ $JADE_MPI_size = "unset" ]]; then JADE_MPI_size=12; fi            
    else
        if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=2; fi
    fi
}  # end of TuneOnzaOptionsMPI
# Host dependant MPI run of Onza
function RunOnza {
    local path_tmp=`pwd`
    cd $path_bin
    cp $config_file $path_bin/onza.config
    if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
        echo "Waiting for shared file system to distibute files..."
        sleep 2
        echo "(1) Nodes $Onza_MPI_nodes procs $Onza_MPI_size"
        salloc -N $Onza_MPI_nodes -n $Onza_MPI_size -p max1day \
            mpirun $MPI_options ./$onza_bin onza.config
    else
        echo "(1) Nodes 1  procs $Onza_MPI_size"
        mpirun -np $Onza_MPI_size $MPI_options ./$onza_bin onza.config
    fi    
    cd $path_tmp
}  # end of RunOnza
if [[ $isTest = $no ]]; then
    TuneOnzaOptionsMPI
    RunOnza
    if [[ $isPGO = $yes ]]; then
        # export OMPI_CXXFLAGS="$flags_O2 -fprofile-generate"
        # export OMPI_LDFLAGS=$OMPI_CXXFLAGS
        flag_cmake_pgo="-DCMAKE_CXX_FLAGS=-fprofile-generate -DCMAKE_EXE_LINKER_FLAGS=-fprofile-generate"
        rm -r $path_bin/*
        BuildOnza
        echo ============================================
        echo Run PGO profiling
        sleep 2
        RunOnza

        cp $path_build/src/mpi-decomposition/CMakeFiles/mpi-decomposition.dir/__/simulation-core/basic-fdtd.cc.gcda $path_build/src/simulation-core/CMakeFiles/simulation-core.dir/basic-fdtd.cc.gcda
        cp $path_build/src/mpi-decomposition/CMakeFiles/mpi-decomposition.dir/__/profiling/timer.cc.gcda $path_build/src/profiling/CMakeFiles/profiling.dir/timer.cc.gcda
        # export OMPI_CXXFLAGS="$flags_O2 -fprofile-use"
        # export OMPI_LDFLAGS=$OMPI_CXXFLAGS
        flag_cmake_pgo="-DCMAKE_CXX_FLAGS=-fprofile-use -DCMAKE_EXE_LINKER_FLAGS=-fprofile-use"
        rm -r $path_bin/*
        BuildOnza
        echo ============================================
        echo Run after PGO compiled
        sleep 2
        RunOnza 
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
        test_name=${test//-/_}
        path_test_config="path_${test_name}_config"
        path_test="$path_bin/$test"  
        #path_test_config="path_${test}_config"
        cp ${!path_test_config} ${path_test}/onza.config
        cp $onza_bin $path_test
        cd $path_test
        backup_Onza_MPI_size=$Onza_MPI_size
        backup_Onza_MPI_nodes=$Onza_MPI_nodes
        if [[ $test = "self-test-X1D-zero" ]]; then
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
                mpirun $MPI_options ./$onza_bin onza.config
            if [[ $test = "self-test-TMz2D-speedup" ]]; then
                echo "(*******) Nodes 16 procs 128 (1024 x 1024, step 1000, ~15.4s)"
                salloc -N 16 -n 128 -p max1hour \
                    mpirun $MPI_options ./$onza_bin onza.config
                echo "(*******) Nodes 16 procs 16 (1024 x 1024, step 1000, ~7s)"
                salloc -N 16 -n 16 -p max1hour \
                    mpirun $MPI_options ./$onza_bin onza.config
                echo "(*******) Nodes 1 procs 8 (1024 x 1024, step 1000, ~30.2s)"
                salloc -N 1 -n 8 -p max1hour \
                    mpirun $MPI_options ./$onza_bin onza.config
                echo "(*******) Nodes 1 procs 1 (1024 x 1024, step 1000, ~47.6s)"
                salloc -N 1 -n 1 -p max1hour \
                    mpirun $MPI_options ./$onza_bin onza.config
            fi

        elif  [[ $HOST == "deb00" || $HOST == "dmmrkovich-birzha" ]]; then
            echo "(*******) Procs 1"
            mpirun -np 1 $MPI_options ./$onza_bin onza.config
            echo "(*******) Procs 2"
            mpirun -np 2 $MPI_options ./$onza_bin onza.config
            echo "(*******) Procs 4"
            mpirun -np 4 $MPI_options ./$onza_bin onza.config
        else
            if [[ $Onza_MPI_size = "unset" ]]; then Onza_MPI_size=2; fi
            echo "(1) Nodes 1  procs $Onza_MPI_size"
            mpirun -np $Onza_MPI_size $MPI_options ./$onza_bin onza.config
        fi
        if [[ $test = "self-test-X1D-zero" ]]; then
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
if [[ $config_file = $path_self_test_X1D_zero_config ]]; then
    echo "Prepare *.png from gnuplot ..."
    cp $path_onza/data/gnuplot/* ./
    ./gnuplot-all.sh >/dev/null  2>&1
    # mkdir tmpdir
    # cp *0241-* $path_bin/tmpdir
    # rm $path_bin/*    
    # cp $path_bin/tmpdir/* ./
    # rm -r $path_bin/tmpdir
fi
#rm *.onza  >/dev/null  2>&1
    
