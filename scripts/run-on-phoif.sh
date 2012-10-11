#!/bin/bash
./rsync-to-phoif.sh
ssh phoif "cd ~/onza-fdtd && ./go.sh new2"
echo "Rsync from phoif to local..."
rm -r ../bin/*
rsync -ave ssh phoif:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null