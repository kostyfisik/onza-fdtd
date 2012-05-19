#!/bin/bash
cd ..
rsync -ave ssh --exclude=bin --exclude=build --exclude=doc ./ phoif:~/onza-fdtd/