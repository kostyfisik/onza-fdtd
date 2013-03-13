#!/bin/bash
./rsync-to-phoif.sh
mode=$1 # blocks - a, valgrind - b, both - c
Usage=">> This script runs all FDTD tests located in ~/test-core.\n
	>> Usage: $./run-all-tests.sh mode\n
	mode:\n
	a - enables block calculation,\n
	b - enables valgrind profiling,\n
	c - enables both options."
if [[ $# > 1 ]]; then
	echo -e $Usage
else
ssh -t phoif "cd ~/onza-fdtd/scripts/test-core && salloc --nodelist=n07 ./run-all-tests.sh $mode"
echo "Rsync from phoif to local..."
rm -r ../bin/*
rsync --exclude run-onza-fdtd -ave ssh phoif:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null
fi
