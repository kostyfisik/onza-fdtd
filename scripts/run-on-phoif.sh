#!/bin/bash
mode=$1
configFile=$2
wrong_params=$3
./rsync-to-phoif.sh
ssh phoif "cd ~/onza-fdtd && ./go.sh $mode $configFile $wrong_params"
echo "Rsync from phoif to local..."
rm -r ../bin/*
rsync --exclude run-onza-fdtd -ave ssh phoif:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null