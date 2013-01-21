#!/bin/bash
binary_file=$1
config_file=$2
fixed_N=$3
fixed_n=$4
./rsync-to-phoif.sh
ssh -t phoif "cd ~/onza-fdtd/scripts/test-core && ./run-all-tests.sh"
echo "Rsync from phoif to local..."
rm -r ../bin/*
rsync --exclude run-onza-fdtd -ave ssh phoif:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null
