#!/bin/bash
mode=$1
configFile=$2
wrong_params=$3
if [[ ! $MPIsize ]]; then MPIsize="unset"; fi
if [[ ! $MPInodes ]]; then MPInodes="unset"; fi
./rsync-to-deb00.sh
ssh -t deb00 "cd ~/onza-fdtd && Onza_MPI_size=$Onza_MPI_size Onza_MPI_nodes=$Onza_MPI_nodes ./go.sh $mode $configFile $wrong_params"
echo "Rsync from deb00 to local..."
rm -r ../bin/*
rsync --exclude run-onza-fdtd -ave ssh deb00:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null