#!/bin/bash
./rsync-to-deb00.sh
ssh deb00 "cd ~/onza-fdtd && ./go.sh"
echo "Rsync from deb00 to local..."
rm -r ../bin/*
rsync -ave ssh deb00:~/onza-fdtd/bin/ ~/onza-fdtd/bin/ >/dev/null