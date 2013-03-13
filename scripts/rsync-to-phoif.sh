#!/bin/bash
cd ..
rsync --delete -ave ssh --exclude=bin --exclude=build* --exclude=doc ~/onza-fdtd/ phoif:~/onza-fdtd/