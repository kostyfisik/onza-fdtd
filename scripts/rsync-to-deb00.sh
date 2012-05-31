#!/bin/bash
cd ..
rsync -ave ssh --exclude=bin --exclude=build --exclude=doc ./ deb00:~/onza-fdtd/