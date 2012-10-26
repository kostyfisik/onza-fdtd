#!/bin/bash
mode=$1
configFile=$2
wrong_params=$3
./rsync-to-deb00.sh
ssh deb00 "cd ~/onza-fdtd && ./go.sh $mode $configFile $wrong_params"
echo "Rsync from deb00 to local..."
rm -r ../bin/*
rsync --exclude run-onza-fdtd -ave ssh deb00:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null