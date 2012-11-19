#!/bin/bash
USAGE=">> Compiles and runs parallel programs\n"
binary_file=$1; config_file=$2; fixed_N=$3; fixed_n=$4;
#echo $1
#echo $2
#if [[ $extra_params ]]; then
#    echo -e ">> ================ !ERROR! =================\n"
#    echo -e ">> Script takes two input parameters: file name and config file\n"
#    echo -e ">> ================ !ERROR! =================\n"
#    echo -e  $USAGE
#    exit 1
#fi

path_onza=$PWD/../..;   path_bin=$path_onza/bin;   path_build=$path_onza/build
path_src=$path_onza/src; path_test_parallel=$path_onza/scripts/test-parallel
path_scripts=$path_onza/scripts; path_data=$path_onza/data

if [[ $binary_file == "onza-fdtd" ]]; then
#    echo ">> Processing $1 with config $2\n"
    cd $path_onza
    ./go.sh build
    cd $path_test_parallel
    cp ./exec_script.py $path_bin/exec_script.py
    cp ./efficiency.plt $path_bin/efficiency.plt
    cp ./mean_time.plt $path_bin/mean_time.plt
    cp ./mean_time_fixed.plt $path_bin/mean_time_fixed.plt
    cp ./gnuplot-all.sh $path_bin/gnuplot-all.sh
    cd $path_data
    cp ./$config_file $path_bin/$config_file
    cd $path_bin
    
    ./exec_script.py $binary_file.bin $config_file $fixed_N $fixed_n
    ./gnuplot-all.sh
else
    echo "Not onza"
   # mpic++ $binary_file -o $file_name.bin
    
fi