#!/bin/bash
./rsync-to-phoif.sh
ssh -t phoif "cd ~/onza-fdtd/scripts/test-core && salloc --nodelist=n07 ./run-all-tests.sh"
echo "Rsync from phoif to local..."
rm -r ../bin/*
rsync --exclude run-onza-fdtd -ave ssh phoif:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null
