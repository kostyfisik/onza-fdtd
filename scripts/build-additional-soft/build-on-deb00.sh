#!/bin/bash
rsync --exclude run-onza-fdtd -ave ssh ~/onza-fdtd/scripts/build-additional-soft/ deb00:~/onza-fdtd/scripts/build-additional-soft/
ssh -t deb00 "cd ~/onza-fdtd/scripts/build-additional-soft && ./build-gcc.sh"
