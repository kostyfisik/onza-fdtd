#!/bin/bash
rsync --exclude run-onza-fdtd -ave ssh ~/onza-fdtd/scripts/build-additional-soft/ phoif:~/onza-fdtd/scripts/build-additional-soft/
ssh -t phoif "cd ~/onza-fdtd/scripts/build-additional-soft && ./build-gcc.sh"