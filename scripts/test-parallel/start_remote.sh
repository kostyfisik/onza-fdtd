#!/bin/bash
binary_file=$1
config_file=$2
fixed_N=$3
fixed_n=$4
./rsy.sh
ssh -t phoif "cd ~/onza-fdtd/scripts/test-parallel && ./start_local.sh $binary_file $config_file $fixed_N $fixed_n"
echo "Rsync from phoif to local..."
rm -r ../bin/*
rsync --exclude run-onza-fdtd -ave ssh phoif:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null
